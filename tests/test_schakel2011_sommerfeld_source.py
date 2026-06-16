import importlib.util
import inspect
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace

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
        for name in (
            "z_s",
            "receiver_z_min",
            "receiver_z_max",
            "receiver_spacing",
            "waveform_nt",
            "plot_time_min_us",
            "plot_time_max_us",
        ):
            self.assertIn(f"{name}:", source)
            self.assertEqual(getattr(cfg, name), model.ZeroOffsetSchakelConfig.__dataclass_fields__[name].default)

    def test_plot_time_axis_cli_overrides_are_applied_to_config(self):
        cfg = model.ZeroOffsetSchakelConfig()
        args = SimpleNamespace(
            z_s=None,
            z_s_mm=None,
            receiver_z_min_mm=None,
            receiver_z_max_mm=None,
            receiver_spacing_mm=None,
            plot_time_min_us=0.0,
            plot_time_max_us=160.0,
            f0=None,
            upper_fluid_conductivity_mode=None,
            transducer_radius_mm=None,
            schakel_bandpass_low_hz=None,
            schakel_bandpass_high_hz=None,
            schakel_gamma_max=None,
            source_mode=None,
            convergence_levels=None,
        )

        model._apply_common_overrides(args, cfg)

        self.assertEqual(cfg.plot_time_min_us, 0.0)
        self.assertEqual(cfg.plot_time_max_us, 160.0)

    def test_default_source_is_documented_as_causal_ricker(self):
        cfg = model.ZeroOffsetSchakelConfig()
        doc = model.schakel_source_A_spectrum.__doc__

        self.assertEqual(cfg.schakel_source_mode, "causal_ricker")
        self.assertIn("causal Ricker", doc)
        self.assertNotIn("visual Fig. 4 pressure-pulse approximation", doc)

    def test_default_snapshot_theta_quadrature_is_192(self):
        run_sig = inspect.signature(model.run_simulation)
        main_source = inspect.getsource(model.main)

        self.assertEqual(run_sig.parameters["n_theta"].default, 192)
        self.assertIn('parser.add_argument("--n-theta", type=int, default=192)', main_source)

    def test_causal_ricker_pressure_spectrum_is_unscaled_pressure_time_integral(self):
        cfg = model.ZeroOffsetSchakelConfig()
        frequencies = np.array([0.5, 1.0, 1.5], dtype=float) * cfg.f0
        omega = 2.0 * np.pi * frequencies
        n_time = 1024

        actual = model.causal_ricker_pressure_spectrum(frequencies, cfg, pressure_peak_pa=1.0, n_time=n_time)

        duration = max(cfg.source_duration_cycles, cfg.source_peak_cycles + 2.0) / cfg.f0
        tau = np.linspace(0.0, duration, n_time)
        src = model.base.ricker(tau - cfg.source_peak_cycles / cfg.f0, cfg.f0)
        ramp_len = max(4, min(len(tau), int(0.25 / cfg.f0 / max(tau[1] - tau[0], 1e-15))))
        ramp = np.ones_like(src)
        ramp[:ramp_len] = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, ramp_len)))
        src *= ramp
        integrate = getattr(np, "trapezoid", np.trapz)
        expected = integrate(src[None, :] * np.exp(-1j * omega[:, None] * tau[None, :]), tau, axis=1)

        np.testing.assert_allclose(actual, expected)
        self.assertLess(float(np.max(np.abs(actual))), 2.0e-6)

    def test_default_source_A_spectrum_uses_one_pa_ricker_pressure_and_reference_distance(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.z_s = 0.1
        cfg.source_pressure_amp = 1.0
        frequencies = np.array([0.0, cfg.f0, 2.0 * cfg.f0], dtype=float)

        actual = model.schakel_source_A_spectrum(frequencies, cfg)
        expected = (
            model.causal_ricker_pressure_spectrum(frequencies, cfg, pressure_peak_pa=1.0)
            * abs(cfg.z_s)
            * model._bandpass_taper(frequencies, cfg)
        )

        np.testing.assert_allclose(actual, expected)
        self.assertLess(float(np.max(np.abs(actual))), 2.0e-7)

    def test_source_pressure_amplitude_scales_linearly_in_pa(self):
        cfg_1pa = model.ZeroOffsetSchakelConfig()
        cfg_1pa.source_pressure_amp = 1.0
        cfg_50kpa = model.ZeroOffsetSchakelConfig()
        cfg_50kpa.source_pressure_amp = 50_000.0
        frequencies = np.array([cfg_1pa.f0], dtype=float)

        one_pa = model.schakel_source_A_spectrum(frequencies, cfg_1pa)
        fifty_kpa = model.schakel_source_A_spectrum(frequencies, cfg_50kpa)

        np.testing.assert_allclose(fifty_kpa, one_pa * 50_000.0)

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

    def test_waveform_time_axis_extends_to_configured_plot_max_us(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.receiver_z_min = -0.001
        cfg.receiver_z_max = 0.001
        cfg.receiver_spacing = 0.001
        cfg.waveform_nt = 40
        cfg.schakel_gamma_max = 1.0
        cfg.plot_time_max_us = 160.0
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        _, t, _ = model.synthesize_waveforms_schakel2011(row, cfg, n_frequencies=3, n_theta=6)

        self.assertAlmostEqual(float(t[-1]), 160.0e-6)

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
