import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = Path(__file__).resolve().parents[1] / "schakel2011_geophysics_reproduction.py"
spec = importlib.util.spec_from_file_location("schakel2011_geophysics_reproduction", MODULE_PATH)
schakel2011 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = schakel2011
spec.loader.exec_module(schakel2011)


class Schakel2011GeophysicsReproductionTest(unittest.TestCase):
    def test_table_parameters_match_geophysics_paper(self):
        params = schakel2011.geophysics_parameters()

        self.assertAlmostEqual(params.K_s, 49.9e9)
        self.assertAlmostEqual(params.K_f, 2.2e9)
        self.assertAlmostEqual(params.eta, 1.0e-3)
        self.assertAlmostEqual(params.rho_f, 998.0)
        self.assertAlmostEqual(params.Lambda, 1.229e-5)
        self.assertAlmostEqual(params.eps_f, 80.1)
        self.assertAlmostEqual(params.eps_s, 4.0)
        self.assertAlmostEqual(params.K_b, 6.6e9)
        self.assertAlmostEqual(params.G, 5.5e9)
        self.assertAlmostEqual(params.phi, 0.345)
        self.assertAlmostEqual(params.rho_s, 2.212e3)
        self.assertAlmostEqual(params.k0_m2, 3.1e-12)
        self.assertAlmostEqual(params.alpha_inf, 2.1)
        self.assertAlmostEqual(params.source_z_m, -0.15)
        self.assertAlmostEqual(params.transducer_radius_m, 0.0142875)

        salinity = schakel2011.salinity_case("1e-3")
        self.assertAlmostEqual(salinity.concentration_mol_l, 1.0e-3)
        self.assertAlmostEqual(salinity.conductivity_s_m, 1.20e-2)
        self.assertAlmostEqual(salinity.zeta_v, -61.5e-3)
        self.assertAlmostEqual(salinity.amplitude_scale, 0.19)

    def test_piston_directivity_zero_angle_limit_and_symmetry(self):
        params = schakel2011.geophysics_parameters()
        frequency = 500_000.0

        zero = schakel2011.piston_directivity(frequency, 0.0, params)
        pos = schakel2011.piston_directivity(frequency, 0.25, params)
        neg = schakel2011.piston_directivity(frequency, -0.25, params)

        self.assertAlmostEqual(zero, 0.5, places=12)
        self.assertAlmostEqual(pos, neg, places=12)
        self.assertGreater(pos, 0.0)
        self.assertLess(pos, 0.5)

    def test_digitized_source_spectrum_peaks_near_500_khz(self):
        source = schakel2011.digitized_source_fig4()

        self.assertIn("time_us", source.columns)
        self.assertIn("pressure_kpa", source.columns)
        self.assertIn("frequency_mhz", source.columns)
        self.assertIn("spectrum_mpa", source.columns)
        peak_freq = source.loc[source["spectrum_mpa"].idxmax(), "frequency_mhz"]
        self.assertGreater(peak_freq, 0.45)
        self.assertLess(peak_freq, 0.60)
        self.assertGreater(source["pressure_kpa"].max(), 45.0)
        self.assertLess(source["pressure_kpa"].min(), -45.0)

    def test_fig11_digitization_includes_repeated_experiment_series(self):
        peaks = schakel2011.digitized_peak_amplitudes_fig11()

        self.assertGreaterEqual(peaks["series"].nunique(), 3)
        self.assertIn("0-1 hours", set(peaks["series"]))
        self.assertIn("18-24 hours", set(peaks["series"]))
        self.assertIn("40-41 hours", set(peaks["series"]))
        self.assertGreaterEqual(len(peaks[peaks["panel"] == "z"]), 20)
        self.assertGreaterEqual(len(peaks[peaks["panel"] == "r"]), 18)

    def test_pressure_normalized_re_is_finite_for_geophysics_case(self):
        coeff = schakel2011.pressure_normalized_reflection_coefficient(500_000.0, 20.0, "1e-3")

        self.assertTrue(np.isfinite(coeff.real))
        self.assertTrue(np.isfinite(coeff.imag))
        self.assertGreater(abs(coeff), 0.0)

    def test_model_peak_amplitude_decays_away_from_interface(self):
        result = schakel2011.model_waveforms_for_positions(
            positions_m=[-0.003, -0.013, -0.023, -0.043],
            radial_m=0.0,
            salinity_key="1e-3",
            n_frequencies=41,
            n_theta=28,
        )

        peaks = result.groupby("z_m")["model_scaled_mv"].apply(lambda x: float(x.max() - x.min()))
        self.assertGreater(peaks.loc[-0.003], peaks.loc[-0.013])
        self.assertGreater(peaks.loc[-0.013], peaks.loc[-0.023])
        self.assertGreater(peaks.loc[-0.023], peaks.loc[-0.043])
        self.assertLess(peaks.loc[-0.003], 5.0)

    def test_frequency_sampling_is_stable_after_absolute_time_phase(self):
        coarse = schakel2011.model_waveforms_for_positions(
            positions_m=[-0.013],
            radial_m=0.0,
            salinity_key="1e-3",
            n_frequencies=61,
            n_theta=42,
        )
        fine = schakel2011.model_waveforms_for_positions(
            positions_m=[-0.013],
            radial_m=0.0,
            salinity_key="1e-3",
            n_frequencies=121,
            n_theta=42,
        )

        coarse_peak = float(coarse["model_scaled_mv"].max() - coarse["model_scaled_mv"].min())
        fine_peak = float(fine["model_scaled_mv"].max() - fine["model_scaled_mv"].min())
        self.assertAlmostEqual(coarse_peak, fine_peak, delta=0.01)

    def test_fig6_model_uses_unscaled_nearest_axis_receiver(self):
        waveform, spectrum = schakel2011.model_fig6_waveform_and_spectrum(
            n_frequencies=41,
            n_theta=28,
        )

        self.assertIn("model_unscaled_mv", waveform.columns)
        self.assertIn("model_display_mv", waveform.columns)
        self.assertIn("spectrum_unscaled_mv", spectrum.columns)
        self.assertIn("spectrum_display_mv", spectrum.columns)
        self.assertAlmostEqual(float(waveform["z_m"].iloc[0]), -0.003)
        self.assertAlmostEqual(float(waveform["r_m"].iloc[0]), 0.0)
        self.assertEqual(str(waveform["salinity_key"].iloc[0]), "1e-3")
        self.assertGreater(float(waveform["model_unscaled_mv"].max()), 0.3)
        self.assertLess(float(waveform["model_unscaled_mv"].min()), -0.3)
        self.assertGreater(float(waveform["model_display_mv"].max()), 0.8)
        self.assertLess(float(waveform["model_display_mv"].min()), -0.5)
        self.assertGreater(float(spectrum["spectrum_unscaled_mv"].max()), 0.01)
        self.assertGreater(float(spectrum["spectrum_display_mv"].max()), 25.0)
        self.assertLess(float(spectrum["spectrum_display_mv"].max()), 30.0)
        peak_frequency = float(spectrum.loc[spectrum["spectrum_unscaled_mv"].idxmax(), "frequency_mhz"])
        self.assertGreater(peak_frequency, 0.35)
        self.assertLess(peak_frequency, 0.65)
        self.assertLess(float(spectrum["frequency_mhz"].iloc[0]), 0.2)
        self.assertGreater(float(spectrum["frequency_mhz"].iloc[-1]), 0.8)

    def test_fig6_paper_digitization_captures_published_axes(self):
        waveform, spectrum = schakel2011.digitized_paper_fig6()

        self.assertIn("time_ms", waveform.columns)
        self.assertIn("paper_electric_potential_mv", waveform.columns)
        self.assertIn("frequency_mhz", spectrum.columns)
        self.assertIn("paper_electric_potential_spectrum_mv", spectrum.columns)
        self.assertLessEqual(float(waveform["time_ms"].min()), 0.055)
        self.assertGreaterEqual(float(waveform["time_ms"].max()), 0.145)
        self.assertGreater(float(waveform["paper_electric_potential_mv"].max()), 0.75)
        self.assertLess(float(waveform["paper_electric_potential_mv"].min()), -1.0)
        self.assertGreater(float(spectrum["paper_electric_potential_spectrum_mv"].max()), 25.0)
        peak_frequency = float(
            spectrum.loc[spectrum["paper_electric_potential_spectrum_mv"].idxmax(), "frequency_mhz"]
        )
        self.assertGreater(peak_frequency, 0.45)
        self.assertLess(peak_frequency, 0.55)

    def test_source_spectrum_preserves_complex_phase_without_extra_scale(self):
        params = schakel2011.geophysics_parameters()
        source = schakel2011._source_A_spectrum(np.array([500_000.0, 600_000.0]), params)
        module_source = MODULE_PATH.read_text(encoding="utf-8")

        self.assertTrue(np.iscomplexobj(source))
        self.assertNotIn("source_spectrum_scale", module_source)
        self.assertNotIn("0.92 * norm", module_source)

        frequencies, spectrum = schakel2011._source_fft_table()
        mask = (frequencies >= params.bandpass_low_hz) & (frequencies <= params.bandpass_high_hz)
        time_s = np.linspace(0.0, 3.6e-6, 401)
        reconstructed = 2.0 * np.real(
            np.trapezoid(
                spectrum[mask, None] * np.exp(1j * 2.0 * np.pi * frequencies[mask, None] * time_s[None, :]),
                frequencies[mask],
                axis=0,
            )
        )
        self.assertGreater(reconstructed.max() - reconstructed.min(), 70_000.0)
        self.assertLess(reconstructed.max() - reconstructed.min(), 100_000.0)

    def test_forward_model_source_uses_ricker_wavelet_requested_for_unknown_experiment(self):
        params = schakel2011.geophysics_parameters()

        self.assertTrue(hasattr(params, "ricker_peak_frequency_hz"))
        frequencies = np.array([0.0, params.ricker_peak_frequency_hz, 2.0 * params.ricker_peak_frequency_hz])

        spectrum = schakel2011._source_A_spectrum(frequencies, params)

        self.assertEqual(params.forward_source_kind, "ricker")
        self.assertAlmostEqual(abs(spectrum[0]), 0.0, places=15)
        self.assertGreater(abs(spectrum[1]), abs(spectrum[2]))
        self.assertIn("ricker", schakel2011._source_A_spectrum.__doc__.lower())

    def test_ricker_parameters_are_marked_as_model_assumptions_not_paper_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            schakel2011._write_parameters(outdir)
            params = pd.read_csv(outdir / "parameters_used.csv")

        ricker_rows = params[params["name"].str.contains("ricker|forward_source_kind", regex=True)]

        self.assertGreaterEqual(len(ricker_rows), 5)
        self.assertTrue(ricker_rows["source"].str.contains("modeling assumption", case=False).all())
        self.assertFalse(ricker_rows["source"].str.contains("Geophysics 2011 Table 1", case=False).any())

    def test_pipeline_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            schakel2011.run_reproduction(outdir, n_frequencies=31, n_theta=24)

            expected = [
                "fig2_reproduction.png",
                "fig4_reproduction.png",
                "fig5_reproduction.png",
                "fig6_reproduction.png",
                "fig7_reproduction.png",
                "fig9_reproduction.png",
                "fig11_reproduction.png",
                "fig11_source_normalization_diagnostics.png",
                "parameters_used.csv",
                "digitized_source_fig4.csv",
                "digitized_source_directivity_fig5.csv",
                "digitized_paper_fig6_waveform.csv",
                "digitized_paper_fig6_spectrum.csv",
                "digitized_peak_amplitudes_fig11.csv",
                "digitized_paper_present_model_fig11.csv",
                "model_waveform_fig6.csv",
                "model_spectrum_fig6.csv",
                "reproduction_residuals_fig6.csv",
                "model_waveforms_fig7.csv",
                "model_waveforms_fig9.csv",
                "model_peak_amplitudes_fig11.csv",
                "formula_audit.md",
                "source_phase_convention_audit.md",
                "sommerfeld_integral_audit.md",
                "digitization_metadata.csv",
                "model_dipole_approximation_fig11.csv",
                "sommerfeld_convergence.csv",
                "frequency_sampling_sensitivity.csv",
                "reproduction_residuals_fig11.csv",
                "reproduction_residuals_fig11_paper_model.csv",
                "reproduction_residual_summary.csv",
                "source_normalization_diagnostics.csv",
            ]
            for name in expected:
                self.assertTrue((outdir / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
