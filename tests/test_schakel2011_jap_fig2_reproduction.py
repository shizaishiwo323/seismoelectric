import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "schakel2011_jap_fig2_reproduction.py"
spec = importlib.util.spec_from_file_location("schakel2011_jap_fig2_reproduction", MODULE_PATH)
jap = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = jap
spec.loader.exec_module(jap)


class Schakel2011JAPFig2ReproductionTest(unittest.TestCase):
    def test_table_ii_parameters_and_geometry(self):
        params = jap.jap_parameters()
        self.assertAlmostEqual(params.K_b, 0.93e9)
        self.assertAlmostEqual(params.G, 0.88e9)
        self.assertAlmostEqual(params.phi, 0.52)
        self.assertAlmostEqual(params.rho_s, 2570.0)
        self.assertAlmostEqual(params.k0_m2, 3.4e-12)
        self.assertAlmostEqual(params.alpha_inf, 1.7)
        self.assertAlmostEqual(params.sigma_f_s_m, 4.8e-2)
        self.assertAlmostEqual(params.zeta_v, -4.0e-2)
        self.assertAlmostEqual(params.sample_width_m, 3.21e-2)

        self.assertEqual(list(jap.JAP_POSITIONS["position"]), [1, 2, 3, 4, 5, 6])
        self.assertTrue(np.allclose(jap.JAP_POSITIONS.loc[:4, "z_m"], [-0.023, -0.018, -0.013, -0.008, -0.003]))
        self.assertAlmostEqual(float(jap.JAP_POSITIONS.loc[5, "z_m"]), 0.010)

    def test_jap_config_reproduces_table_lambda_and_zeta_setup(self):
        params = jap.jap_parameters()
        cfg = jap._jap_config(params)
        state = jap.strict.dynamic_coefficients(
            params.phi,
            params.k0_m2,
            params.alpha_inf,
            jap._cH(cfg),
            2.0 * np.pi * 500_000.0,
            cfg,
            C_override_molL=cfg.C_molL,
        )
        self.assertAlmostEqual(float(state["Lambda"]), params.Lambda, delta=1e-10)
        self.assertAlmostEqual(cfg.pH, 6.0)
        self.assertGreater(cfg.C_molL, 0.0)

    def test_source_spectrum_uses_eq1_on_axis_piston_calibration(self):
        params = jap.jap_parameters()
        frequency = np.array([500_000.0])
        freqs, fft_amp = jap._source_fft_table()
        raw = np.interp(frequency, freqs, np.real(fft_amp)) + 1j * np.interp(frequency, freqs, np.imag(fft_amp))

        calibrated = jap._source_A_spectrum(frequency, params)

        self.assertTrue(np.allclose(calibrated, raw * params.source_z_to_interface_m / 0.5 * jap._bandpass_taper(frequency, params)))

    def test_sommerfeld_coefficients_are_finite(self):
        coeff = jap._fluid_incident_coefficients(500_000.0, 1000.0)
        pf = jap._pf_interface_coefficients(500_000.0, 1000.0)
        self.assertTrue(np.isfinite(coeff["R_E"].real))
        self.assertTrue(np.isfinite(coeff["T_TM"].imag))
        self.assertTrue(np.isfinite(pf["R_Pf"].real))
        self.assertTrue(np.isfinite(pf["R_TM_potential"].imag))

    def test_eq7_front_return_tm_term_is_modulated_by_back_pf_reflection(self):
        params = jap.jap_parameters()
        omega = 2.0 * np.pi * 500_000.0
        kx = omega / params.cP_fluid_m_s * 0.35
        z_receiver = 0.010
        coeff = jap._fluid_incident_coefficients(500_000.0, kx)
        pf = jap._pf_interface_coefficients(500_000.0, kx)

        terms = jap._porous_eq8_terms(500_000.0, kx, z_receiver)

        expected = (
            pf["R_Pf"]
            * np.exp(-1j * coeff["k3_Pf"] * params.sample_width_m)
            * pf["R_TM_potential"]
            * np.exp(-1j * coeff["k3_TM"] * z_receiver)
        )
        self.assertAlmostEqual(terms["se_front_return_tm"], expected)

    def test_fig2_waveforms_match_key_paper_features(self):
        result = jap.synthesize_jap_fig2(n_frequencies=161, n_theta=64, integration_method="fixed")
        summary = jap._arrival_summary(result).set_index("position")

        fluid = summary.loc[[1, 2, 3, 4, 5], "max_abs_mv"].to_numpy()
        self.assertTrue(np.all(np.diff(fluid) > 0.0))
        self.assertGreater(summary.loc[5, "max_abs_mv"], 0.12)
        self.assertLess(summary.loc[5, "max_abs_mv"], 0.25)
        self.assertLess(summary.loc[6, "max_abs_mv"], 0.5)
        self.assertAlmostEqual(summary.loc[5, "global_peak_time_ms"], 0.10, delta=0.003)
        self.assertAlmostEqual(summary.loc[6, "second_0p105ms_peak_time_ms"], 0.1057, delta=0.003)
        self.assertLess(summary.loc[5, "early_0p10ms_peak_mv"], 0.0)
        self.assertGreater(summary.loc[6, "early_0p10ms_peak_mv"], 0.0)
        self.assertGreater(summary.loc[6, "second_minus_early_us"], 2.0)
        self.assertLess(summary.loc[6, "second_minus_early_us"], 8.0)
        pos6 = result[result["position"] == 6]
        early_peak = float(pos6[(pos6["time_ms"] >= 0.096) & (pos6["time_ms"] <= 0.108)]["model_mv"].abs().max())
        late_tail = float(pos6[pos6["time_ms"] >= 0.112]["model_mv"].abs().max())
        self.assertLess(late_tail / early_peak, 0.35)

    def test_fig2_display_waveforms_suppress_out_of_arrival_sidelobes(self):
        result = jap.synthesize_jap_fig2(n_frequencies=161, n_theta=64, integration_method="fixed")
        for position in range(1, 7):
            part = result[result["position"] == position]
            peak = float(part["model_mv"].abs().max())
            pre_arrival = float(part[part["time_ms"] < 0.094]["model_mv"].abs().max())
            self.assertLess(pre_arrival / peak, 0.05, position)

    def test_model_amplitudes_are_close_to_published_predicted_panels_without_fitting(self):
        result = jap.synthesize_jap_fig2(n_frequencies=161, n_theta=64, integration_method="fixed")
        digitized = jap.digitized_published_fig2_traces(n_time=361)
        published = digitized[digitized["column"] == "predicted"]
        for position in range(1, 7):
            model_part = result[result["position"] == position]["model_mv"]
            paper_part = published[published["position"] == position]["digitized_mv"]
            model_p2p = float(model_part.max() - model_part.min())
            paper_p2p = float(paper_part.max() - paper_part.min())
            self.assertGreater(model_p2p / paper_p2p, 0.65, position)
            self.assertLess(model_p2p / paper_p2p, 1.5, position)

    def test_digitized_published_fig2_traces_include_predicted_and_observed_columns(self):
        traces = jap.digitized_published_fig2_traces(n_time=181)

        self.assertEqual(set(traces["column"]), {"predicted", "observed"})
        self.assertEqual(set(traces["position"]), {1, 2, 3, 4, 5, 6})
        self.assertEqual(len(traces), 12 * 181)
        self.assertTrue(np.allclose(traces.groupby(["column", "position"])["time_ms"].min(), 0.055))
        self.assertTrue(np.allclose(traces.groupby(["column", "position"])["time_ms"].max(), 0.145))

        predicted = traces[traces["column"] == "predicted"]
        observed = traces[traces["column"] == "observed"]
        self.assertLessEqual(predicted["digitized_mv"].abs().max(), 0.55)
        self.assertLessEqual(observed["digitized_mv"].abs().max(), 0.25)
        self.assertGreater(predicted.groupby("position")["digitized_mv"].apply(lambda x: x.max() - x.min()).max(), 0.15)
        self.assertGreater(observed.groupby("position")["digitized_mv"].apply(lambda x: x.max() - x.min()).max(), 0.05)

    def test_pipeline_writes_expected_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            jap.run_reproduction(outdir, n_frequencies=31, n_theta=24, integration_method="fixed")
            for name in [
                "parameters_used.csv",
                "jap_fig2_positions.csv",
                "jap_fig2_model_waveforms.csv",
                "jap_fig2_arrival_summary.csv",
                "digitized_published_fig2_waveforms.csv",
                "digitization_metadata_fig2.csv",
                "jap_fig2_source_crop.png",
                "jap_fig2_reproduction.png",
                "jap_fig2_model_vs_digitized.png",
                "formula_audit_jap_fig2.md",
            ]:
                self.assertTrue((outdir / name).exists(), name)


if __name__ == "__main__":
    unittest.main()
