import importlib.util
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "seismoelectric_offset_liu2018_spectral.py"
spec = importlib.util.spec_from_file_location("seismoelectric_offset_liu2018_spectral", MODULE_PATH)
se = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = se
spec.loader.exec_module(se)


class WaveformDiagnosticsTest(unittest.TestCase):
    def test_complex_sqrt_branch_attenuates_with_schakel_spatial_convention(self):
        k = se.complex_sqrt_branch(1.0 - 1.0j)

        self.assertGreater(k.real, 0.0)
        self.assertLess(k.imag, 0.0)
        self.assertLess(abs(np.exp(-1j * k)), 1.0)

    def test_liu_phase_helpers_follow_exp_iomega_t_convention(self):
        spatial_phase = 0.7
        omega = 2.3
        t = np.array([0.0, 0.4])

        self.assertTrue(np.allclose(se.liu_spatial_phase(spatial_phase), np.exp(-1j * spatial_phase)))
        self.assertTrue(np.allclose(se.liu_time_phase(omega, t), np.exp(1j * omega * t)))

    def test_causal_ricker_source_spectrum_uses_exp_minus_iomega_tau(self):
        cfg = se.SEConfig()
        omega = np.array([0.7, 1.3]) * 2.0 * np.pi * cfg.f0
        n_time = 512

        actual = se.causal_ricker_source_spectrum(omega, cfg, n_time=n_time)

        duration = max(cfg.source_duration_cycles, cfg.source_peak_cycles + 2.0) / cfg.f0
        tau = np.linspace(0.0, duration, n_time)
        src = se.ricker(tau - cfg.source_peak_cycles / cfg.f0, cfg.f0)
        ramp_len = max(4, min(len(tau), int(0.25 / cfg.f0 / max(tau[1] - tau[0], 1e-15))))
        ramp = np.ones_like(src)
        ramp[:ramp_len] = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, ramp_len)))
        src *= ramp
        integrate = getattr(np, "trapezoid", np.trapz)
        expected = integrate(src[None, :] * np.exp(-1j * omega[:, None] * tau[None, :]), tau, axis=1)
        expected /= np.max(np.abs(expected))

        self.assertTrue(np.allclose(actual, expected))

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

    def test_zero_offset_waveform_display_indices_follow_receiver_spacing(self):
        z = np.arange(-20.0e-3, 20.0e-3 + 0.5e-3, 1.0e-3)
        cfg = se.SEConfig()
        cfg.receiver_spacing = 1.0e-3
        cfg.offset_D = 0.0

        display_idx = se.waveform_display_indices(z, cfg)

        np.testing.assert_array_equal(display_idx, np.arange(len(z)))

    def test_waveform_gather_uses_configured_time_axis_limits_us(self):
        cfg = se.SEConfig()
        cfg.plot_time_min_us = 0.0
        cfg.plot_time_max_us = 160.0
        z = np.array([-0.001, 0.001])
        t = np.linspace(0.0, 200.0e-6, 8)
        U = np.vstack([np.sin(np.linspace(0.0, np.pi, len(t))), -np.sin(np.linspace(0.0, np.pi, len(t)))])
        row = pd.Series({"Time_s": 1.0, "Porosity": 0.25})

        with tempfile.TemporaryDirectory() as tmp, patch("matplotlib.axes.Axes.set_xlim") as set_xlim:
            se.plot_waveform_gather(z, t, U, row, cfg, Path(tmp))

        set_xlim.assert_any_call(0.0, 160.0)

    def test_spectral_waveform_time_axis_extends_to_configured_plot_max_us(self):
        cfg = se.SEConfig()
        cfg.receiver_z_min = -0.001
        cfg.receiver_z_max = 0.001
        cfg.receiver_spacing = 0.001
        cfg.waveform_nt = 8
        cfg.spectral_f_min_factor = 1.0
        cfg.spectral_f_max_factor = 1.1
        cfg.source_kb_m_inv = 1.0e9
        cfg.plot_time_max_us = 160.0
        row = pd.Series({
            "Porosity": 0.24,
            "Permeability_mD": 100.0,
            "Tortuosity": 2.0,
            "OutletHConc": 1.0e-10,
        })

        _, t, _ = se.synthesize_waveforms_spectral(row, cfg, n_omega=2, n_k=3)

        self.assertAlmostEqual(float(t[-1]), 160.0e-6)

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

    def test_peak_amplitude_spectral_records_all_reflected_and_transmitted_peaks(self):
        cfg = se.SEConfig()
        ts = pd.DataFrame({"valid_poroelastic": [True, True]})
        df = pd.DataFrame({"Time_s": [1.0, 2.0]})
        z = np.array([-0.01, 0.0, 0.01])
        t = np.array([0.0, 1.0])
        original = se.synthesize_waveforms_spectral

        def fake_synthesize(row, *_args, **_kwargs):
            if row["Time_s"] == 1.0:
                U = np.array([[0.0, 2.0], [0.0, 9.0], [0.0, 3.0]])
            else:
                U = np.array([[0.0, 4.0], [0.0, 5.0], [0.0, 12.0]])
            return z, t, U

        try:
            se.synthesize_waveforms_spectral = fake_synthesize
            out = se.compute_peak_amplitude_spectral(ts, df, cfg)
        finally:
            se.synthesize_waveforms_spectral = original

        self.assertTrue(np.allclose(out["Amax_waveform_spectral"], [9.0, 12.0]))
        self.assertTrue(np.allclose(out["Amax_waveform_spectral_RE"], [2.0, 4.0]))
        self.assertTrue(np.allclose(out["Amax_waveform_spectral_TE"], [3.0, 12.0]))
        self.assertTrue(np.allclose(out["Amax_waveform_spectral_RE_norm"], [1.0, 2.0]))
        self.assertTrue(np.allclose(out["Amax_waveform_spectral_TE_norm"], [1.0, 4.0]))


if __name__ == "__main__":
    unittest.main()
