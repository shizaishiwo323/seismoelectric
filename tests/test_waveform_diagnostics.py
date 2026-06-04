import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "seismoelectric_offset_liu2018_spectral.py"
spec = importlib.util.spec_from_file_location("seismoelectric_offset_liu2018_spectral", MODULE_PATH)
se = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = se
spec.loader.exec_module(se)


class WaveformDiagnosticsTest(unittest.TestCase):
    def test_liu_interface_coefficient_uses_incidence_wavenumber_magnitude(self):
        self.assertEqual(se.liu_interface_coefficient_kx(-12.5), 12.5)
        self.assertEqual(se.liu_interface_coefficient_kx(0.0), 0.0)
        self.assertEqual(se.liu_interface_coefficient_kx(7.0), 7.0)

    def test_liu_electrical_potential_coefficients_convert_schakel_potentials(self):
        coeff = {
            "R_E": 2.0 + 3.0j,
            "T_TM": 4.0 - 1.0j,
            "alpha_TM": 0.5 + 0.25j,
            "k3_E": 6.0 - 2.0j,
            "k3_TM": 8.0 + 1.0j,
        }

        re_potential, te_potential = se.liu_electrical_potential_coefficients(coeff, -2.0)

        self.assertTrue(np.allclose(re_potential, -(coeff["k3_E"] / 2.0) * coeff["R_E"]))
        self.assertTrue(
            np.allclose(
                te_potential,
                (coeff["k3_TM"] / 2.0) * coeff["alpha_TM"] * coeff["T_TM"],
            )
        )

    def test_liu_electrical_potential_coefficients_require_nonzero_limit_sample(self):
        coeff = {
            "R_E": 0.0j,
            "T_TM": 0.0j,
            "alpha_TM": 1.0 + 0.0j,
            "k3_E": 1.0 + 0.0j,
            "k3_TM": 1.0 + 0.0j,
        }

        with self.assertRaises(ValueError):
            se.liu_electrical_potential_coefficients(coeff, 0.0)

    def test_spectral_synthesis_passes_nonnegative_kx_to_interface_solver(self):
        cfg = se.SEConfig()
        cfg.receiver_z_min = -0.001
        cfg.receiver_z_max = 0.001
        cfg.receiver_spacing = 0.001
        cfg.waveform_nt = 4
        cfg.spectral_f_min_factor = 1.0
        cfg.spectral_f_max_factor = 1.1
        cfg.source_kb_m_inv = 1.0e9
        row = pd.Series({
            "Porosity": 0.24,
            "Permeability_mD": 100.0,
            "Tortuosity": 2.0,
            "OutletHConc": 1.0e-10,
        })
        seen_kx = []
        original = se.se_coefficients

        def fake_se_coefficients(*args, **kwargs):
            seen_kx.append(kwargs["kx_override"])
            return {
                "R_E": 1.0 + 0.0j,
                "T_TM": 1.0 + 0.0j,
                "alpha_TM": 1.0 + 0.0j,
                "k3_fl": 1.0 + 0.0j,
                "k3_E": 1.0 + 0.0j,
                "k3_TM": 1.0 + 0.0j,
            }

        try:
            se.se_coefficients = fake_se_coefficients
            se.synthesize_waveforms_spectral(row, cfg, n_omega=2, n_k=3)
        finally:
            se.se_coefficients = original

        self.assertTrue(seen_kx)
        self.assertTrue(all(kx >= 0.0 for kx in seen_kx))

    def test_spectral_synthesis_evaluates_normal_incidence_potential_limit(self):
        cfg = se.SEConfig()
        cfg.receiver_z_min = -0.001
        cfg.receiver_z_max = 0.001
        cfg.receiver_spacing = 0.001
        cfg.waveform_nt = 4
        cfg.spectral_f_min_factor = 1.0
        cfg.spectral_f_max_factor = 1.1
        cfg.source_kb_m_inv = 1.0e9
        row = pd.Series({
            "Porosity": 0.24,
            "Permeability_mD": 100.0,
            "Tortuosity": 2.0,
            "OutletHConc": 1.0e-10,
        })
        seen_kx = []
        original = se.se_coefficients

        def fake_se_coefficients(*args, **kwargs):
            seen_kx.append(kwargs["kx_override"])
            return {
                "R_E": kwargs["kx_override"] + 0.0j,
                "T_TM": kwargs["kx_override"] + 0.0j,
                "alpha_TM": 1.0 + 0.0j,
                "k3_fl": 1.0 + 0.0j,
                "k3_E": 1.0 + 0.0j,
                "k3_TM": 1.0 + 0.0j,
            }

        try:
            se.se_coefficients = fake_se_coefficients
            se.synthesize_waveforms_spectral(row, cfg, n_omega=2, n_k=3)
        finally:
            se.se_coefficients = original

        self.assertTrue(seen_kx)
        self.assertTrue(all(kx > 0.0 for kx in seen_kx))

    def test_finite_offset_waveform_display_indices_exclude_interface_trace(self):
        z = np.array([-0.01, -0.005, 0.0, 0.005, 0.01])
        cfg = se.SEConfig()
        cfg.receiver_spacing = 0.005
        cfg.offset_D = 0.045

        finite_offset_idx = se.waveform_display_indices(z, cfg)
        cfg.offset_D = 0.0
        zero_offset_idx = se.waveform_display_indices(z, cfg)

        self.assertNotIn(2, finite_offset_idx)
        self.assertIn(2, zero_offset_idx)

    def test_waveform_spatial_peak_diagnostics_separate_re_and_ttm_sides(self):
        z = np.array([-0.02, -0.01, 0.0, 0.01, 0.02])
        t = np.array([0.0, 1.0, 2.0])
        U = np.array(
            [
                [0.0, 2.0, 0.0],
                [0.0, 5.0, 0.0],
                [0.0, 7.0, 0.0],
                [0.0, 3.0, 0.0],
                [0.0, 1.0, 0.0],
            ]
        )

        diag = se.waveform_spatial_peak_diagnostics(z, t, U, offset_D=0.045)

        self.assertEqual(list(diag["side"]), ["R_E", "R_E", "interface", "T_E", "T_E"])
        self.assertTrue(np.allclose(diag["peak_abs"], [2.0, 5.0, 7.0, 3.0, 1.0]))
        self.assertTrue(np.allclose(diag["distance_from_interface_m"], [0.02, 0.01, 0.0, 0.01, 0.02]))
        self.assertTrue(np.isclose(diag["liu_dipole_peak_distance_m"].iloc[0], 0.045 / np.sqrt(2.0)))


if __name__ == "__main__":
    unittest.main()
