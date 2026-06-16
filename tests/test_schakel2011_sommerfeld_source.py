import importlib.util
import inspect
import sys
import tempfile
import unittest
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


if __name__ == "__main__":
    unittest.main()
