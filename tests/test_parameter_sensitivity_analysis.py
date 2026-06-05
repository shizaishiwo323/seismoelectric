import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "parameter_sensitivity_analysis.py"
spec = importlib.util.spec_from_file_location("parameter_sensitivity_analysis", MODULE_PATH)
psa = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = psa
spec.loader.exec_module(psa)


class ParameterSensitivityAnalysisTest(unittest.TestCase):
    def test_pore_volume_to_surface_ratio_uses_porosity_grain_volume_and_surface_area(self):
        row = pd.Series({
            "Porosity": 0.25,
            "GrainVolume_cm3": 0.75,
            "SurfaceArea_cm2": 2.0,
        })

        ratio_m = psa.pore_volume_to_surface_ratio_m(row)

        self.assertTrue(np.isclose(ratio_m, 0.125e-2))

    def test_one_at_a_time_row_changes_only_requested_parameter_group(self):
        base = pd.Series({
            "Porosity": 0.2,
            "Permeability_mD": 10.0,
            "Tortuosity": 2.0,
            "OutletHConc": 1.0e-10,
        })
        target = pd.Series({
            "Porosity": 0.4,
            "Permeability_mD": 30.0,
            "Tortuosity": 3.0,
            "OutletHConc": 5.0e-10,
        })

        row = psa.make_one_at_a_time_row(base, target, "permeability")

        self.assertEqual(row["Porosity"], base["Porosity"])
        self.assertEqual(row["Permeability_mD"], target["Permeability_mD"])
        self.assertEqual(row["Tortuosity"], base["Tortuosity"])
        self.assertEqual(row["OutletHConc"], base["OutletHConc"])

    def test_log_contribution_table_keeps_signed_effects_and_interaction_residual(self):
        metrics = {
            "baseline": 10.0,
            "full_target": 100.0,
            "porosity": 20.0,
            "permeability": 5.0,
        }

        table = psa.build_log_contribution_table(metrics, ["porosity", "permeability"])

        by_name = table.set_index("component")
        self.assertTrue(np.isclose(by_name.loc["porosity", "delta_log10_metric"], np.log10(2.0)))
        self.assertTrue(np.isclose(by_name.loc["permeability", "delta_log10_metric"], np.log10(0.5)))
        expected_residual = np.log10(10.0) - np.log10(2.0) - np.log10(0.5)
        self.assertTrue(np.isclose(by_name.loc["nonlinear_interaction_residual", "delta_log10_metric"], expected_residual))
        self.assertEqual(by_name.loc["porosity", "sign"], "positive")
        self.assertEqual(by_name.loc["permeability", "sign"], "negative")

    def test_schakel_table_config_matches_reference_sensitivity_settings(self):
        cfg, baseline = psa.schakel_reference_config()

        self.assertTrue(np.isclose(cfg.K_s, 40.0e9))
        self.assertTrue(np.isclose(cfg.eps_s, 4.0))
        self.assertTrue(np.isclose(baseline["phi"], 0.24))
        self.assertTrue(np.isclose(baseline["k0_m2"], 0.390e-12))
        self.assertTrue(np.isclose(baseline["alpha_inf"], 2.3))
        self.assertTrue(np.isclose(baseline["omega"], 1.0e6))
        self.assertTrue(np.isclose(baseline["theta_deg"], 45.0))


if __name__ == "__main__":
    unittest.main()
