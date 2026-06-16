import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

BASE_PATH = ROOT / "seismoelectric_offset_liu2018_spectral.py"
base_spec = importlib.util.spec_from_file_location("seismoelectric_offset_liu2018_spectral", BASE_PATH)
base_model = importlib.util.module_from_spec(base_spec)
sys.modules[base_spec.name] = base_model
base_spec.loader.exec_module(base_model)

MODULE_PATH = ROOT / "schakel2011_sommerfeld.py"
spec = importlib.util.spec_from_file_location("schakel2011_sommerfeld", MODULE_PATH)
model = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = model
spec.loader.exec_module(model)


class ZeroOffsetSchakel2011RTTest(unittest.TestCase):
    def test_default_config_forces_zero_offset_and_preserves_receiver_geometry(self):
        cfg = model.ZeroOffsetSchakelConfig()

        self.assertEqual(cfg.offset_D, 0.0)
        self.assertEqual(cfg.receiver_z_min, -0.02)
        self.assertEqual(cfg.receiver_z_max, 0.02)
        self.assertEqual(cfg.receiver_spacing, 0.001)

    def test_sommerfeld_pressure_coefficient_uses_schakel_pressure_normalization(self):
        coeff = {"R_E": 6.0 + 8.0j}

        actual = model.pressure_normalized_re_from_coeff(coeff, rho_f=1000.0, omega=2.0)

        self.assertEqual(actual, (6.0 + 8.0j) / (1000.0 * 4.0))

    def test_schakel_source_spectrum_preserves_causal_source_phase(self):
        cfg = model.ZeroOffsetSchakelConfig()
        frequencies = np.linspace(cfg.schakel_bandpass_low_hz, cfg.schakel_bandpass_high_hz, 256)
        source = model.schakel_source_A_spectrum(frequencies, cfg) * 0.5 / cfg.z_s
        tau = np.linspace(-5.0e-6, 8.0e-6, 1200)

        reconstructed = 2.0 * np.real(
            np.trapezoid(
                source[None, :] * np.exp(1j * 2.0 * np.pi * frequencies[None, :] * tau[:, None]),
                frequencies,
                axis=1,
            )
        )

        peak_tau = float(tau[np.argmax(np.abs(reconstructed))])
        self.assertGreater(peak_tau, 0.0)
        self.assertGreater(np.nanmax(np.abs(reconstructed[tau >= 0.0])), np.nanmax(np.abs(reconstructed[tau < 0.0])))

    def test_waveform_synthesis_returns_existing_output_shape(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.receiver_z_min = -0.002
        cfg.receiver_z_max = 0.002
        cfg.receiver_spacing = 0.002
        cfg.waveform_nt = 24
        cfg.schakel_gamma_max = 1.5
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        z, t, u = model.synthesize_waveforms_schakel2011(
            row,
            cfg,
            n_frequencies=3,
            n_theta=6,
            integration_method="fixed",
        )

        self.assertEqual(u.shape, (len(z), len(t)))
        self.assertEqual(len(z), 3)
        self.assertEqual(len(t), cfg.waveform_nt)
        self.assertTrue(np.all(np.isfinite(u)))

    def test_waveform_has_porous_side_response_and_polarity_reversal(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.receiver_z_min = -0.002
        cfg.receiver_z_max = 0.002
        cfg.receiver_spacing = 0.002
        cfg.waveform_nt = 160
        cfg.schakel_gamma_max = 1.5
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        z, t, u = model.synthesize_waveforms_schakel2011(
            row,
            cfg,
            n_frequencies=24,
            n_theta=10,
            integration_method="fixed",
        )

        fluid = u[np.where(z < 0.0)[0][0]]
        porous = u[np.where(z > 0.0)[0][0]]
        self.assertGreater(float(np.nanmax(np.abs(porous))), 0.0)
        diag = model.interface_em_polarity_diagnostics(z, t, u, cfg, distances_m=(0.002,))
        self.assertEqual(len(diag), 1)
        self.assertTrue(bool(diag["polarity_reversed"].iloc[0]))
        self.assertGreater(float(diag["time_after_T0_us"].iloc[0]), 0.0)
        self.assertLess(float(diag["time_after_T0_us"].iloc[0]), 8.0)

    def test_zero_offset_peak_amplitude_increases_toward_interface_after_T0(self):
        cfg = model.ZeroOffsetSchakelConfig()
        cfg.receiver_z_min = -0.004
        cfg.receiver_z_max = 0.004
        cfg.receiver_spacing = 0.002
        cfg.waveform_nt = 180
        cfg.schakel_gamma_max = 1.5
        row = pd.Series(
            {
                "Porosity": 0.24,
                "Permeability_mD": 100.0,
                "Tortuosity": 2.0,
                "OutletHConc": 1.0e-10,
            }
        )

        z, t, u = model.synthesize_waveforms_schakel2011(
            row,
            cfg,
            n_frequencies=24,
            n_theta=10,
            integration_method="fixed",
        )
        T0 = cfg.z_s / np.sqrt(cfg.K_fl / cfg.rho_fl)
        post = u[:, t >= T0]
        peaks = np.nanmax(np.abs(post), axis=1)

        self.assertGreater(peaks[np.where(np.isclose(z, -0.002))[0][0]], peaks[np.where(np.isclose(z, -0.004))[0][0]])
        self.assertGreater(peaks[np.where(np.isclose(z, 0.002))[0][0]], peaks[np.where(np.isclose(z, 0.004))[0][0]])
        if np.any(t < T0):
            self.assertGreater(float(np.nanmax(np.abs(post))), float(np.nanmax(np.abs(u[:, t < T0]))))
        else:
            self.assertGreaterEqual(float(t[0]), T0)

    def test_fixed_receiver_peak_uses_nearest_noninterface_electrodes(self):
        z = np.array([-0.004, -0.002, 0.0, 0.002, 0.004])
        t = np.array([0.0, 1.0])
        u = np.array(
            [
                [9.0, 1.0],
                [2.0, -3.0],
                [0.0, 0.0],
                [4.0, -5.0],
                [8.0, 1.0],
            ]
        )

        summary = model.fixed_nearest_receiver_peak_summary_schakel2011(z, t, u)

        self.assertEqual(summary["R_E"]["peak_abs"], 3.0)
        self.assertEqual(summary["R_E"]["receiver_z_m"], -0.002)
        self.assertEqual(summary["T_E"]["peak_abs"], 5.0)
        self.assertEqual(summary["T_E"]["receiver_z_m"], 0.002)

    def test_waveform_plot_values_are_millivolts(self):
        values = np.array([0.001, -0.0025, 0.0])

        values_mv = model.electric_potential_to_millivolts(values)

        np.testing.assert_allclose(values_mv, np.array([1.0, -2.5, 0.0]))

    def test_pipeline_writes_same_named_result_types(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            input_csv = tmpdir / "rt.csv"
            input_csv.write_text(
                "Time_s,Porosity,Permeability_mD,Tortuosity,OutletHConc\n"
                "0,0.24,100,2,1e-10\n"
                "60,0.25,120,2.1,1e-10\n",
                encoding="utf-8",
            )
            outdir = tmpdir / "out"

            model.run_simulation(
                input_csv,
                outdir,
                n_frequencies=3,
                n_theta=6,
                peak_n_frequencies=3,
                peak_n_theta=6,
                integration_method="fixed",
                snapshot_target_phi=0.24,
            )

            expected = [
                "seismoelectric_timeseries_results.csv",
                "parameters_used.csv",
                "run_summary.csv",
                "coefficients_vs_dissolution_time.png",
                "coefficients_vs_dissolution_time_logy.png",
                "dynamic_coefficients_vs_dissolution_time.png",
                "dynamic_coefficients_vs_dissolution_time_logy.png",
                "peak_amplitude_RE_TE_vs_dissolution_time.png",
                "peak_amplitude_RE_TE_vs_dissolution_time_logy.png",
                "waveform_snapshot_schakel2011.csv",
                "waveform_snapshot_schakel2011.npz",
                "waveform_snapshot_schakel2011.png",
                "waveform_interface_em_polarity_diagnostics.csv",
                "waveform_t0_causality_diagnostics.csv",
                "waveform_spatial_peak_diagnostics.csv",
                "waveform_spatial_peak_diagnostics.png",
                "formula_audit.md",
            ]
            for name in expected:
                self.assertTrue((outdir / name).exists(), name)

            removed = [
                "peak_amplitude_vs_dissolution_time.png",
                "transmitted_peak_amplitude_vs_dissolution_time.png",
            ]
            for name in removed:
                self.assertFalse((outdir / name).exists(), name)

            summary = pd.read_csv(outdir / "run_summary.csv", index_col=0).iloc[:, 0]
            self.assertEqual(summary["waveform_mode"], "schakel2011_sommerfeld_zerooffset")
            self.assertAlmostEqual(float(summary["offset_D_m"]), 0.0)

            ts = pd.read_csv(outdir / "seismoelectric_timeseries_results.csv")
            self.assertIn("Amax_waveform_schakel2011_RE_fixed_receiver_z_m", ts.columns)
            self.assertIn("Amax_waveform_schakel2011_TE_fixed_receiver_z_m", ts.columns)


if __name__ == "__main__":
    unittest.main()
