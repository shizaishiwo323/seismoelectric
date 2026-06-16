import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

MODULE_PATH = ROOT / "zerooffset_waveform_display_search.py"
spec = importlib.util.spec_from_file_location("zerooffset_waveform_display_search", MODULE_PATH)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


class ZeroOffsetWaveformDisplaySearchTest(unittest.TestCase):
    def test_candidate_grid_uses_symmetric_zero_offset_close_views(self):
        candidates = module.default_candidate_grid()

        self.assertGreaterEqual(len(candidates), 6)
        for candidate in candidates:
            self.assertLess(candidate.receiver_z_min_m, 0.0)
            self.assertGreater(candidate.receiver_z_max_m, 0.0)
            self.assertAlmostEqual(abs(candidate.receiver_z_min_m), candidate.receiver_z_max_m)
            self.assertGreater(candidate.receiver_spacing_m, 0.0)
            self.assertEqual(candidate.waveform_t_before_s, 0.0)
            self.assertGreater(candidate.waveform_t_after_s, 0.0)

    def test_visibility_score_prefers_near_interface_contrast_and_valid_polarity(self):
        z = np.array([-0.02, -0.01, -0.002, 0.002, 0.01, 0.02])
        t = np.linspace(5.0e-5, 6.0e-5, 100)
        pulse = np.sin(np.linspace(0.0, 2.0 * np.pi, len(t)))
        clear = np.array([0.05, 0.12, 1.0, -0.6, -0.10, -0.04])[:, None] * pulse[None, :]
        flat = np.ones((len(z), len(t))) * 0.05

        clear_score = module.score_display_candidate(z, t, clear)
        flat_score = module.score_display_candidate(z, t, flat)

        self.assertGreater(clear_score["visibility_score"], flat_score["visibility_score"])
        self.assertTrue(clear_score["polarity_reversal_near_interface"])
        self.assertGreater(clear_score["near_to_far_peak_ratio"], 5.0)

    def test_enhanced_plot_writes_zoom_gather_and_heatmap(self):
        z = np.array([-0.004, -0.002, 0.002, 0.004])
        t = np.linspace(5.0e-5, 5.8e-5, 80)
        pulse = np.sin(np.linspace(0.0, 4.0 * np.pi, len(t)))
        u = np.array([0.2, 1.0, -0.7, -0.15])[:, None] * pulse[None, :]

        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            module.plot_enhanced_waveform_gather(
                z,
                t,
                u,
                outdir,
                name="demo",
                title="Demo",
                wiggle_scale_fraction=0.75,
            )

            self.assertTrue((outdir / "demo.png").exists())
            self.assertTrue((outdir / "demo_heatmap.png").exists())
            self.assertTrue((outdir / "demo_global_heatmap.png").exists())

    def test_reference_style_plots_write_pixpin_and_schakel_fig6_outputs(self):
        z = np.linspace(-0.045, 0.0, 10)
        t = np.linspace(1.01e-4, 1.06e-4, 220)
        center = 1.032e-4
        width = 4.0e-7
        pulse = np.exp(-0.5 * ((t - center) / width) ** 2) - 0.55 * np.exp(-0.5 * ((t - center - 7e-7) / width) ** 2)
        amplitudes = np.linspace(0.2, 1.0, len(z))
        u = amplitudes[:, None] * pulse[None, :]

        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            module.plot_pixpin_reference_style_gather(
                z,
                t,
                u,
                outdir,
                name="pixpin",
                title="PixPin style",
            )
            module.plot_schakel_fig6_style_panels(
                t,
                u[-1],
                outdir,
                name="fig6",
            )

            self.assertTrue((outdir / "pixpin.png").exists())
            self.assertTrue((outdir / "fig6.png").exists())
            self.assertTrue((outdir / "fig6_time_trace.csv").exists())
            self.assertTrue((outdir / "fig6_frequency_spectrum.csv").exists())

    def test_reference_style_diagnostics_flag_visible_traces_and_electric_potential_spectrum(self):
        z = np.linspace(-0.045, 0.0, 10)
        t = np.linspace(5.36e-5, 6.16e-5, 300)
        center = 5.79e-5
        width = 2.2e-7
        pulse = np.exp(-0.5 * ((t - center) / width) ** 2)
        u = np.linspace(0.1, 1.0, len(z))[:, None] * pulse[None, :]

        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            module.plot_pixpin_reference_style_gather(z, t, u, outdir, name="pixpin", title="")
            module.plot_schakel_fig6_style_panels(
                t,
                u[-1],
                outdir,
                name="fig6",
                spectrum_frequency_hz=np.linspace(1.0, 1.0e6, 200),
                spectrum_amplitude=np.exp(-0.5 * ((np.linspace(1.0, 1.0e6, 200) - 5.0e5) / 1.7e5) ** 2),
            )

            peaks = pd.read_csv(outdir / "pixpin_trace_peaks.csv")
            self.assertGreaterEqual(int((peaks["display_peak_abs"] > 0.7).sum()), 8)
            trace = pd.read_csv(outdir / "fig6_time_trace.csv")
            self.assertLess(float(trace["active_width_us"].iloc[0]), 4.0)
            spec = pd.read_csv(outdir / "fig6_frequency_spectrum.csv")
            self.assertIn("electric_potential_spectrum_display_mV", spec.columns)
            self.assertEqual(spec["spectrum_quantity"].iloc[0], "reflected_electric_potential")

    def test_strict_fig6_frequency_count_pushes_periodic_replica_outside_window(self):
        cfg = module.zero_model.ZeroOffsetSchakelConfig()

        n_freq = module.strict_fig6_frequency_count(
            cfg,
            display_window_s=100.0e-6,
            requested_n_frequencies=96,
        )
        frequencies = module.zero_model._frequency_grid(cfg, n_freq)
        replica_period_s = 1.0 / float(frequencies[1] - frequencies[0])

        self.assertGreater(n_freq, 96)
        self.assertGreaterEqual(replica_period_s, 2.0 * 100.0e-6)

    def test_schakel_fig6_display_diagnostics_reject_background_highpass_and_replicas(self):
        t = np.linspace(5.36e-5, 1.536e-4, 800)
        center = 5.80e-5
        width = 2.5e-7
        pulse = np.exp(-0.5 * ((t - center) / width) ** 2)
        trace = pulse + 0.015 * np.sin(2.0 * np.pi * 0.2e6 * t)

        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            module.plot_schakel_fig6_style_panels(
                t,
                trace,
                outdir,
                name="fig6",
                time_xlim_s=(float(t[0]), float(t[-1])),
            )

            diag = pd.read_csv(outdir / "fig6_time_trace.csv")
            self.assertEqual(diag["display_processing"].iloc[0], "median_removed_only")
            self.assertLessEqual(int(diag["pulse_groups_above_20pct"].iloc[0]), 1)

    def test_cleanup_removes_stale_candidate_images_without_touching_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp)
            stale = outdir / "candidate_rank1_old_case.png"
            stale_heatmap = outdir / "candidate_rank2_old_case_heatmap.png"
            score = outdir / "display_candidate_scores.csv"
            stale.write_text("old", encoding="utf-8")
            stale_heatmap.write_text("old", encoding="utf-8")
            score.write_text("keep", encoding="utf-8")

            module.clear_previous_display_outputs(outdir)

            self.assertFalse(stale.exists())
            self.assertFalse(stale_heatmap.exists())
            self.assertTrue(score.exists())


if __name__ == "__main__":
    unittest.main()
