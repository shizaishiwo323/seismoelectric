import importlib.util
import inspect
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

MODULE_PATH = ROOT / "schakel2011_sommerfeld.py"
spec = importlib.util.spec_from_file_location("schakel2011_sommerfeld", MODULE_PATH)
model = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = model
spec.loader.exec_module(model)


class Schakel2011SommerfeldSourceTest(unittest.TestCase):
    def test_config_parameters_are_declared_locally_not_inherited_from_base_config(self):
        source = inspect.getsource(model.ZeroOffsetSchakelConfig)
        cfg = model.ZeroOffsetSchakelConfig()

        self.assertFalse(issubclass(model.ZeroOffsetSchakelConfig, model.base.SEConfig))
        self.assertNotIn("base.SEConfig", source)
        for name in ("z_s", "receiver_z_min", "receiver_z_max", "receiver_spacing", "waveform_nt"):
            self.assertIn(f"{name}:", source)
            self.assertEqual(getattr(cfg, name), model.ZeroOffsetSchakelConfig.__dataclass_fields__[name].default)

    def test_default_source_is_documented_as_causal_ricker(self):
        cfg = model.ZeroOffsetSchakelConfig()
        doc = model.schakel_source_A_spectrum.__doc__

        self.assertEqual(cfg.schakel_source_mode, "causal_ricker")
        self.assertIn("causal Ricker", doc)
        self.assertNotIn("visual Fig. 4 pressure-pulse approximation", doc)

    def test_default_source_spectrum_matches_base_causal_ricker(self):
        cfg = model.ZeroOffsetSchakelConfig()
        frequencies = np.array([0.0, cfg.f0, 2.0 * cfg.f0], dtype=float)

        actual = model.schakel_source_A_spectrum(frequencies, cfg)
        expected = cfg.source_pressure_amp * model.base.causal_ricker_source_spectrum(
            2.0 * np.pi * frequencies,
            cfg,
        )

        np.testing.assert_allclose(actual, expected)
        self.assertAlmostEqual(abs(actual[0]), abs(expected[0]))

    def test_parameter_table_marks_source_mode_as_causal_ricker_assumption(self):
        cfg = model.ZeroOffsetSchakelConfig()

        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            model.save_parameter_table(cfg, outdir)
            params = pd.read_csv(outdir / "parameters_used.csv")

        row = params[params["parameter"] == "schakel_source_mode"].iloc[0]
        self.assertEqual(row["value"], "causal_ricker")
        self.assertIn("causal Ricker", row["meaning"])

    def test_waveform_time_axis_starts_at_zero_and_still_contains_T0(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.receiver_z_min = -0.001
        cfg.receiver_z_max = 0.001
        cfg.receiver_spacing = 0.001
        cfg.waveform_nt = 40
        cfg.schakel_gamma_max = 1.0
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        _, t, _ = model.synthesize_waveforms_schakel2011(row, cfg, n_frequencies=3, n_theta=6)
        T0 = cfg.z_s / np.sqrt(cfg.K_fl / cfg.rho_fl)

        self.assertAlmostEqual(float(t[0]), 0.0)
        self.assertLess(float(t[0]), T0)
        self.assertGreater(float(t[-1]), T0)

    def test_interpolated_interface_receiver_is_excluded_from_peak_metrics(self):
        cfg = model.ZeroOffsetSchakelConfig()
        ts = pd.DataFrame({"valid_poroelastic": [True]})
        row = pd.DataFrame(
            {
                "Porosity": [0.24],
                "Permeability_mD": [100.0],
                "Tortuosity": [2.0],
                "OutletHConc": [1.0e-10],
            }
        )
        z = np.array([-0.001, 0.0, 0.001])
        t = np.array([0.0, 1.0e-6])
        u = np.array(
            [
                [1.0, 2.0],
                [1000.0, 2000.0],
                [-3.0, -4.0],
            ]
        )

        with patch.object(model, "synthesize_waveforms_schakel2011", return_value=(z, t, u)):
            out = model.compute_peak_amplitude_schakel2011(ts, row, cfg, n_frequencies=3, n_theta=6)

        self.assertEqual(float(out["Amax_waveform_schakel2011"].iloc[0]), 4.0)
        self.assertEqual(float(out["Amax_waveform_schakel2011_RE"].iloc[0]), 2.0)
        self.assertEqual(float(out["Amax_waveform_schakel2011_TE"].iloc[0]), 4.0)

        diag = model.waveform_spatial_peak_diagnostics_schakel2011(z, t, u, cfg)
        interface = diag[diag["side"] == "interface"].iloc[0]
        self.assertFalse(bool(interface["include_in_quantitative_summary"]))
        self.assertEqual(interface["interface_value_policy"], "plot_only_linear_interpolation")

    def test_convergence_diagnostics_record_peak_time_amplitude_and_polarity(self):
        cfg = model.ZeroOffsetSchakelConfig()
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        def fake_synthesize(_row, _cfg, n_frequencies=None, n_theta=48, integration_method="fixed"):
            scale = float(n_frequencies + n_theta)
            z = np.array([-0.001, 0.0, 0.001])
            t = np.array([0.0, 1.0e-6, 2.0e-6])
            u = np.array(
                [
                    [0.0, scale, 0.25 * scale],
                    [0.0, 999.0 * scale, 999.0 * scale],
                    [0.0, -0.5 * scale, -1.5 * scale],
                ]
            )
            return z, t, u

        diag = model.waveform_convergence_diagnostics(
            row,
            cfg,
            levels=(4, 8),
            integration_method="fixed",
            synthesize_func=fake_synthesize,
        )

        self.assertEqual(set(diag["side"]), {"all", "R_E", "T_E"})
        self.assertTrue(np.all(diag["quantitative_excludes_interface"]))
        all_rows = diag[diag["side"] == "all"].sort_values("n_frequencies")
        self.assertEqual(list(all_rows["n_theta"]), [4, 8])
        self.assertEqual(float(all_rows.iloc[0]["peak_time_us"]), 2.0)
        self.assertEqual(float(all_rows.iloc[0]["peak_signed_polarity"]), -1.0)
        self.assertTrue(np.isfinite(float(all_rows.iloc[1]["peak_abs_relative_change_from_previous"])))


if __name__ == "__main__":
    unittest.main()
