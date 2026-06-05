import importlib.util
import sys
import unittest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "liu2018_fig2b_dipole_comparison.py"
spec = importlib.util.spec_from_file_location("liu2018_fig2b_dipole_comparison", MODULE_PATH)
comparison = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = comparison
spec.loader.exec_module(comparison)


class Liu2018Fig2bComparisonTest(unittest.TestCase):
    def test_default_quadrature_matches_reference_spectral_model(self):
        self.assertEqual(comparison.DEFAULT_SPECTRAL_N_OMEGA, 48)
        self.assertEqual(comparison.DEFAULT_SPECTRAL_N_K, 401)
        self.assertEqual(comparison.PRIDE_COLUMN_PREFIX, "pride")

    def test_receiver_line_dipole_geometry_follows_liu_eq4(self):
        D = 0.045
        z_peak = D / np.sqrt(2.0)
        z = np.array([-z_peak, 0.0, z_peak])

        geometry = comparison.liu_receiver_line_dipole_geometry(z, D)

        expected_theta = np.array([
            180.0 - np.degrees(np.arctan(np.sqrt(2.0))),
            90.0,
            np.degrees(np.arctan(np.sqrt(2.0))),
        ])
        self.assertTrue(np.allclose(geometry["theta_deg"], expected_theta))
        self.assertEqual(geometry["dipole_signed"][1], 0.0)
        self.assertTrue(np.isclose(abs(geometry["dipole_signed"][0]), abs(geometry["dipole_signed"][2])))

    def test_liu_semicircle_axis_has_bottom_diameter(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        try:
            comparison.configure_liu_semicircle_axis(ax)

            self.assertEqual(ax.get_theta_direction(), 1.0)
            self.assertTrue(np.isclose(ax.get_theta_offset(), 0.0))
            self.assertTrue(np.isclose(ax.get_thetamin(), 0.0))
            self.assertTrue(np.isclose(ax.get_thetamax(), 180.0))
        finally:
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
