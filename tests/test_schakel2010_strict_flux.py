import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "schakel2010_strict_sensitivity.py"
spec = importlib.util.spec_from_file_location("schakel2010_strict_sensitivity", MODULE_PATH)
strict = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = strict
spec.loader.exec_module(strict)


class Schakel2010StrictFluxTest(unittest.TestCase):
    def test_table_iii_energy_flux_coefficients_at_omega_1e6_theta_45(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=45.0)

        self.assertTrue(np.isclose(result["RE_EE"], 5.4312e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["TE_TM_TM"], -1.6581e-6, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["RE_Pr_Pr"], -5.6982e-2, rtol=3.0e-3))

    def test_table_iii_energy_flux_coefficients_at_omega_1e6_theta_30(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=30.0)

        self.assertTrue(np.isclose(result["RE_EE"], 1.3503e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["TE_TM_TM"], -5.0880e-7, rtol=3.0e-3))
        self.assertTrue(np.isclose(result["RE_Pr_Pr"], -3.2877e-1, rtol=3.0e-3))

    def test_strict_coefficients_are_not_squared_potential_proxies(self):
        result = strict.reference_energy_flux_coefficients(omega=1.0e6, theta_deg=45.0)

        self.assertGreater(abs(result["R_E"]) ** 2, 1.0e16)
        self.assertLess(result["RE_EE"], 1.0e-5)
        self.assertLess(abs(result["TE_TM_TM"]), 1.0e-4)


if __name__ == "__main__":
    unittest.main()
