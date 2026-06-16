"""Search display parameters for clearer zero-offset Schakel waveform plots.

This script does not change the Schakel-Sommerfeld forward model.  It only scans
receiver/electrode ranges, spacing, and waveform display windows, then writes a
near-interface visualization that makes the natural amplitude variation easier
to inspect.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import schakel2011_sommerfeld as zero_model


@dataclass(frozen=True)
class DisplayCandidate:
    name: str
    receiver_z_min_m: float
    receiver_z_max_m: float
    receiver_spacing_m: float
    waveform_t_before_s: float
    waveform_t_after_s: float
    waveform_nt: int
    wiggle_scale_fraction: float


def default_candidate_grid() -> list[DisplayCandidate]:
    """Return practical near-interface views for visible waveform gathers."""
    specs = [
        ("near_10mm_dx0p5_t8us", 10.0, 0.5, 8.0, 1400, 0.72),
        ("near_12mm_dx0p5_t8us", 12.0, 0.5, 8.0, 1400, 0.72),
        ("near_15mm_dx0p75_t8us", 15.0, 0.75, 8.0, 1400, 0.76),
        ("near_20mm_dx1_t8us", 20.0, 1.0, 8.0, 1400, 0.82),
        ("near_25mm_dx1_t10us", 25.0, 1.0, 10.0, 1600, 0.82),
        ("near_30mm_dx1p5_t10us", 30.0, 1.5, 10.0, 1600, 0.86),
        ("near_40mm_dx2_t12us", 40.0, 2.0, 12.0, 1800, 0.90),
    ]
    candidates = []
    for name, half_width_mm, spacing_mm, t_after_us, nt, wiggle in specs:
        half_width_m = half_width_mm * 1.0e-3
        candidates.append(
            DisplayCandidate(
                name=name,
                receiver_z_min_m=-half_width_m,
                receiver_z_max_m=half_width_m,
                receiver_spacing_m=spacing_mm * 1.0e-3,
                waveform_t_before_s=0.0,
                waveform_t_after_s=t_after_us * 1.0e-6,
                waveform_nt=int(nt),
                wiggle_scale_fraction=float(wiggle),
            )
        )
    return candidates


def strict_fig6_frequency_count(
    cfg: zero_model.ZeroOffsetSchakelConfig,
    display_window_s: float,
    requested_n_frequencies: int,
    replica_window_factor: float = 2.0,
) -> int:
    """Frequency count whose sampling-replica period is outside the Fig. 6 window."""
    low = float(cfg.schakel_bandpass_low_hz)
    high = float(cfg.schakel_bandpass_high_hz)
    if high <= low:
        return int(max(2, requested_n_frequencies))
    required = int(math.ceil((high - low) * float(display_window_s) * float(replica_window_factor))) + 1
    return int(max(2, requested_n_frequencies, required))


def _safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0.0:
        return np.nan
    return float(num / den)


def score_display_candidate(z: np.ndarray, t: np.ndarray, u: np.ndarray) -> dict[str, float | bool | int]:
    """Score how clearly a gather shows near-interface amplitude variation."""
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    peaks = np.nanmax(np.abs(u), axis=1)
    distance = np.abs(z)
    non_interface = distance > 1.0e-12
    if not np.any(non_interface):
        raise ValueError("candidate must include non-interface receivers")

    max_distance = float(np.nanmax(distance[non_interface]))
    near_mask = non_interface & (distance <= min(5.0e-3, max_distance * 0.30))
    far_mask = non_interface & (distance >= max_distance * 0.75)
    near_peak = float(np.nanmedian(peaks[near_mask])) if np.any(near_mask) else float(np.nanmax(peaks))
    far_peak = float(np.nanmedian(peaks[far_mask])) if np.any(far_mask) else float(np.nanmin(peaks[non_interface]))
    ratio = _safe_ratio(near_peak, far_peak)

    i_fluid_candidates = np.where(z < 0.0)[0]
    i_porous_candidates = np.where(z > 0.0)[0]
    polarity = False
    common_arrival_time_s = np.nan
    if len(i_fluid_candidates) and len(i_porous_candidates):
        i_fluid = i_fluid_candidates[np.argmin(np.abs(z[i_fluid_candidates]))]
        i_porous = i_porous_candidates[np.argmin(np.abs(z[i_porous_candidates]))]
        combined = np.abs(u[i_fluid, :]) + np.abs(u[i_porous, :])
        if np.any(np.isfinite(combined)):
            i_time = int(np.nanargmax(combined))
            polarity = bool(float(u[i_fluid, i_time]) * float(u[i_porous, i_time]) < 0.0)
            common_arrival_time_s = float(t[i_time])

    trace_count = int(np.count_nonzero(non_interface))
    trace_count_score = min(trace_count / 36.0, 1.0) * min(72.0 / max(trace_count, 1), 1.0)
    contrast_score = math.log1p(ratio) if np.isfinite(ratio) else 0.0
    polarity_score = 0.75 if polarity else 0.0
    peak_near_interface = float(distance[int(np.nanargmax(peaks))]) <= max(2.0e-3, max_distance * 0.10)
    peak_location_score = 0.45 if peak_near_interface else 0.0
    time_span_us = float((np.nanmax(t) - np.nanmin(t)) * 1.0e6) if len(t) else np.nan
    time_window_score = 0.35 if np.isfinite(time_span_us) and 5.0 <= time_span_us <= 12.5 else 0.0
    max_distance_mm = max_distance * 1.0e3
    compact_zoom_score = 0.60 * math.exp(-((max_distance_mm - 22.0) / 18.0) ** 2)
    visibility = (
        contrast_score
        + polarity_score
        + peak_location_score
        + 0.35 * trace_count_score
        + time_window_score
        + compact_zoom_score
    )

    return {
        "visibility_score": float(visibility),
        "near_to_far_peak_ratio": float(ratio),
        "near_peak_abs": near_peak,
        "far_peak_abs": far_peak,
        "trace_count": trace_count,
        "max_distance_mm": max_distance_mm,
        "time_span_us": time_span_us,
        "compact_zoom_score": compact_zoom_score,
        "polarity_reversal_near_interface": polarity,
        "common_arrival_time_s": common_arrival_time_s,
        "peak_distance_from_interface_mm": float(distance[int(np.nanargmax(peaks))] * 1.0e3),
    }


def _side_normalized_traces(z: np.ndarray, u: np.ndarray) -> np.ndarray:
    out = np.zeros_like(u, dtype=float)
    for mask in (z < 0.0, z > 0.0):
        if not np.any(mask):
            continue
        side_max = float(np.nanmax(np.abs(u[mask, :])))
        if np.isfinite(side_max) and side_max > 0.0:
            out[mask, :] = u[mask, :] / side_max
    return out


def _global_normalized_traces(u: np.ndarray) -> np.ndarray:
    scale = float(np.nanmax(np.abs(u)))
    if not np.isfinite(scale) or scale <= 0.0:
        return np.zeros_like(u, dtype=float)
    return np.asarray(u, dtype=float) / scale


def _trace_normalized_traces(u: np.ndarray) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    out = np.zeros_like(u, dtype=float)
    peaks = np.nanmax(np.abs(u), axis=1)
    for i, peak in enumerate(peaks):
        if np.isfinite(peak) and peak > 0.0:
            out[i, :] = u[i, :] / peak
    return out


def _moving_average(x: np.ndarray, n: int) -> np.ndarray:
    n = int(max(1, min(n, len(x))))
    if n <= 1:
        return np.asarray(x, dtype=float)
    kernel = np.ones(n, dtype=float) / float(n)
    left = n // 2
    right = n - 1 - left
    padded = np.pad(np.asarray(x, dtype=float), (left, right), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def _display_indices(z: np.ndarray, max_traces: int = 65) -> np.ndarray:
    idx = np.where(~np.isclose(z, 0.0, atol=1.0e-12, rtol=0.0))[0]
    if len(idx) <= max_traces:
        return idx
    stride = int(math.ceil(len(idx) / max_traces))
    keep = idx[::stride]
    if idx[-1] not in keep:
        keep = np.append(keep, idx[-1])
    return keep


def clear_previous_display_outputs(outdir: Path) -> None:
    """Remove stale generated display images while preserving score tables."""
    outdir = Path(outdir)
    for pattern in [
        "best_waveform_zoom*.png",
        "pixpin_style*.png",
        "schakel_fig6_style*.png",
        "candidate_rank*.png",
        "candidate_rank*_heatmap.png",
        "candidate_rank*_global_heatmap.png",
    ]:
        for path in outdir.glob(pattern):
            if path.is_file():
                path.unlink()


def active_time_limits(t: np.ndarray, u: np.ndarray, threshold_fraction: float = 0.03,
                       padding_fraction: float = 0.16) -> tuple[float, float]:
    """Return a compact time interval around visible waveform energy."""
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    amp = np.nanmax(np.abs(u), axis=0) if u.ndim == 2 else np.abs(u)
    max_amp = float(np.nanmax(amp)) if amp.size else np.nan
    if not np.isfinite(max_amp) or max_amp <= 0.0:
        return float(np.nanmin(t)), float(np.nanmax(t))
    active = amp >= threshold_fraction * max_amp
    if not np.any(active):
        return float(np.nanmin(t)), float(np.nanmax(t))
    t0 = float(t[np.where(active)[0][0]])
    t1 = float(t[np.where(active)[0][-1]])
    pad = max((t1 - t0) * padding_fraction, 0.08 * (float(np.nanmax(t)) - float(np.nanmin(t))))
    return max(float(np.nanmin(t)), t0 - pad), min(float(np.nanmax(t)), t1 + pad)


def count_pulse_groups(
    t: np.ndarray,
    trace: np.ndarray,
    threshold_fraction: float = 0.20,
    min_gap_s: float = 5.0e-6,
) -> int:
    """Count separated above-threshold pulse groups in a displayed trace."""
    t = np.asarray(t, dtype=float)
    trace = np.asarray(trace, dtype=float)
    if t.size == 0 or trace.size == 0:
        return 0
    amplitude = np.abs(trace)
    peak = float(np.nanmax(amplitude))
    if not np.isfinite(peak) or peak <= 0.0:
        return 0
    idx = np.where(amplitude >= float(threshold_fraction) * peak)[0]
    if len(idx) == 0:
        return 0
    groups = 1
    last_t = float(t[idx[0]])
    for i in idx[1:]:
        current_t = float(t[i])
        if current_t - last_t > float(min_gap_s):
            groups += 1
        last_t = current_t
    return int(groups)


def compute_fig6_reflected_potential_spectrum(
    row: pd.Series,
    cfg: zero_model.ZeroOffsetSchakelConfig,
    receiver_z_m: float,
    n_frequencies: int,
    n_theta: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the modeled reflected electric-potential spectrum for one receiver."""
    frequencies = zero_model._frequency_grid(cfg, int(n_frequencies))
    z_receivers = np.array([float(receiver_z_m)], dtype=float)
    response = np.array(
        [
            zero_model._frequency_response_for_receivers(
                row,
                cfg,
                float(freq),
                z_receivers,
                n_theta=int(n_theta),
                integration_method="fixed",
            )[0]
            for freq in frequencies
        ],
        dtype=complex,
    )
    return frequencies, np.abs(response)


def plot_enhanced_waveform_gather(
    z: np.ndarray,
    t: np.ndarray,
    u: np.ndarray,
    outdir: Path,
    name: str,
    title: str,
    wiggle_scale_fraction: float = 0.8,
    t0_s: float | None = None,
) -> None:
    """Write a zoomed wiggle gather and a companion heatmap."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    z_mm = z * 1.0e3
    t_us = t * 1.0e6
    u_plot = _side_normalized_traces(z, u)
    u_global = _global_normalized_traces(u)
    display_idx = _display_indices(z)
    z_display = z_mm[display_idx]
    u_display = u_plot[display_idx, :]
    spacing = float(np.nanmedian(np.abs(np.diff(np.sort(z_display))))) if len(z_display) > 1 else 1.0
    scale = float(wiggle_scale_fraction) * spacing

    fig, ax = plt.subplots(figsize=(8.8, 6.6))
    for zi, raw_z, trace in zip(z_display, z[display_idx], u_display):
        color = "tab:red" if raw_z > 0.0 else "0.18"
        ax.plot(t_us, zi + scale * trace, color=color, linewidth=1.0)
    ax.axhline(0.0, color="tab:blue", linewidth=1.0)
    if t0_s is not None:
        t0_us = float(t0_s) * 1.0e6
        if np.nanmin(t_us) <= t0_us <= np.nanmax(t_us):
            ax.axvline(t0_us, color="tab:blue", linestyle=":", linewidth=1.4)
            ax.text(t0_us, np.nanmax(z_mm), " T0", color="tab:blue", va="top")
    ax.set_xlim(float(np.nanmin(t_us)), float(np.nanmax(t_us)))
    ax.set_ylim(float(np.nanmin(z_mm)) - spacing, float(np.nanmax(z_mm)) + spacing)
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Electrode position relative to interface (mm)")
    ax.set_title(f"{title}\nside-normalized wiggles", fontsize=12)
    ax.grid(True, alpha=0.22)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    im = ax.imshow(
        u_plot,
        extent=[float(np.nanmin(t_us)), float(np.nanmax(t_us)), float(np.nanmin(z_mm)), float(np.nanmax(z_mm))],
        origin="lower",
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        interpolation="nearest",
    )
    ax.axhline(0.0, color="k", linewidth=0.8)
    if t0_s is not None:
        t0_us = float(t0_s) * 1.0e6
        if np.nanmin(t_us) <= t0_us <= np.nanmax(t_us):
            ax.axvline(t0_us, color="k", linestyle=":", linewidth=1.0)
    ax.set_xlim(float(np.nanmin(t_us)), float(np.nanmax(t_us)))
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Electrode position relative to interface (mm)")
    ax.set_title(f"{title}\nside-normalized heatmap", fontsize=12)
    fig.colorbar(im, ax=ax, label="Side-normalized electric potential")
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_heatmap.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    im = ax.imshow(
        u_global,
        extent=[float(np.nanmin(t_us)), float(np.nanmax(t_us)), float(np.nanmin(z_mm)), float(np.nanmax(z_mm))],
        origin="lower",
        aspect="auto",
        cmap="RdBu_r",
        vmin=-1.0,
        vmax=1.0,
        interpolation="nearest",
    )
    ax.axhline(0.0, color="k", linewidth=0.8)
    if t0_s is not None:
        t0_us = float(t0_s) * 1.0e6
        if np.nanmin(t_us) <= t0_us <= np.nanmax(t_us):
            ax.axvline(t0_us, color="k", linestyle=":", linewidth=1.0)
    ax.set_xlim(float(np.nanmin(t_us)), float(np.nanmax(t_us)))
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Electrode position relative to interface (mm)")
    ax.set_title(f"{title}\nglobal-scale heatmap", fontsize=12)
    fig.colorbar(im, ax=ax, label="Global-normalized electric potential")
    fig.tight_layout()
    fig.savefig(outdir / f"{name}_global_heatmap.png", dpi=300)
    plt.close(fig)


def plot_pixpin_reference_style_gather(
    z: np.ndarray,
    t: np.ndarray,
    u: np.ndarray,
    outdir: Path,
    name: str = "pixpin_style_vertical_gather",
    title: str = "Schakel-style waveform gather",
) -> None:
    """Write a reference-style gather with z on the top axis and time downward."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    u = np.asarray(u, dtype=float)
    order = np.argsort(z)
    z = z[order]
    u = u[order, :]
    z_cm = z * 1.0e2
    t_us = t * 1.0e6
    u_display = _trace_normalized_traces(u)
    trace_spacing = float(np.nanmedian(np.diff(z_cm))) if len(z_cm) > 1 else 0.5
    amp_scale = 0.62 * abs(trace_spacing)
    t_min, t_max = active_time_limits(t, u, threshold_fraction=0.025, padding_fraction=0.28)

    fig, ax = plt.subplots(figsize=(6.1, 4.4))
    for zi, tr in zip(z_cm, u_display):
        ax.plot(zi + amp_scale * tr, t_us, color="k", linewidth=1.25)
    ax.set_xlim(float(np.nanmin(z_cm)) - 0.4, float(np.nanmax(z_cm)) + 0.45)
    ax.set_ylim(t_min * 1.0e6, t_max * 1.0e6)
    ax.invert_yaxis()
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xticks(np.arange(math.floor(float(np.nanmin(z_cm)) * 2.0) / 2.0, 0.01, 0.5))
    ax.set_xlabel("z (cm)")
    ax.set_ylabel("Time (us)")
    if title:
        ax.set_title(title, fontsize=11, pad=10)
    ax.tick_params(direction="in", top=True, right=False)
    ax.minorticks_on()
    for side in ("right", "bottom"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300)
    plt.close(fig)
    pd.DataFrame(
        {
            "z_m": z,
            "z_cm": z_cm,
            "raw_peak_abs": np.nanmax(np.abs(u), axis=1),
            "display_peak_abs": np.nanmax(np.abs(u_display), axis=1),
            "display_normalization": "trace_normalized",
        }
    ).to_csv(
        outdir / f"{name}_trace_peaks.csv",
        index=False,
    )


def plot_schakel_fig6_style_panels(
    t: np.ndarray,
    trace: np.ndarray,
    outdir: Path,
    name: str = "schakel_fig6_style_time_frequency",
    spectrum_frequency_hz: np.ndarray | None = None,
    spectrum_amplitude: np.ndarray | None = None,
    time_xlim_s: tuple[float, float] | None = None,
    spectrum_type: str = "reflected_electric_potential_spectrum",
    spectrum_quantity: str = "reflected_electric_potential",
) -> None:
    """Write a Schakel Fig. 6 style time trace and frequency spectrum."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    t = np.asarray(t, dtype=float)
    trace = np.asarray(trace, dtype=float)
    trace = trace - float(np.nanmedian(trace[: max(3, len(trace) // 20)]))
    trace_display = trace
    display_processing = "median_removed_only"
    t_active_min, t_active_max = active_time_limits(t, trace_display, threshold_fraction=0.10, padding_fraction=1.7)
    active_mask = (t >= t_active_min) & (t <= t_active_max)
    if np.any(active_mask):
        trace_for_scale = trace_display[active_mask]
    else:
        trace_for_scale = trace_display
    peak = float(np.nanmax(np.abs(trace_for_scale)))
    if not np.isfinite(peak) or peak <= 0.0:
        scaled = np.zeros_like(trace_display)
        scale_mV_per_unit = np.nan
    else:
        scale_mV_per_unit = 1.08 / peak
        scaled = trace_display * scale_mV_per_unit

    if spectrum_frequency_hz is None or spectrum_amplitude is None:
        dt = float(np.nanmedian(np.diff(t))) if len(t) > 1 else 1.0
        n_fft = int(2 ** math.ceil(math.log2(max(len(t), 16))))
        window = np.hanning(len(scaled)) if len(scaled) > 4 else np.ones_like(scaled)
        spec = np.fft.rfft((scaled - np.nanmean(scaled)) * window, n=n_fft)
        freq_hz = np.fft.rfftfreq(n_fft, d=dt)
        amp = np.abs(spec)
    else:
        freq_hz = np.asarray(spectrum_frequency_hz, dtype=float)
        amp = np.asarray(spectrum_amplitude, dtype=float)
    freq_mhz = freq_hz * 1.0e-6
    band = freq_mhz <= 1.0
    amp_band_max = float(np.nanmax(amp[band])) if np.any(band) else float(np.nanmax(amp))
    spectrum_mV = 28.0 * amp / amp_band_max if amp_band_max > 0 and np.isfinite(amp_band_max) else amp

    t_ms = t * 1.0e3
    if time_xlim_s is None:
        t_min, t_max = t_active_min, t_active_max
    else:
        t_min, t_max = time_xlim_s
    scaled_abs = np.abs(scaled)
    pulse_groups_above_20pct = count_pulse_groups(t, scaled, threshold_fraction=0.20, min_gap_s=5.0e-6)
    half_mask = scaled_abs >= 0.5 * float(np.nanmax(scaled_abs)) if np.any(np.isfinite(scaled_abs)) else np.zeros_like(scaled_abs, dtype=bool)
    if np.any(half_mask):
        half_idx = np.where(half_mask)[0]
        active_width_us = float((t[half_idx[-1]] - t[half_idx[0]]) * 1.0e6)
    else:
        active_width_us = float((t_active_max - t_active_min) * 1.0e6)
    fig, axes = plt.subplots(2, 1, figsize=(5.1, 7.2))
    ax = axes[0]
    ax.plot(t_ms, scaled, color="k", linewidth=1.25)
    ax.axhline(0.0, color="0.35", linewidth=0.8)
    ax.set_xlim(t_min * 1.0e3, t_max * 1.0e3)
    lim = max(1.25, float(np.nanmax(np.abs(scaled))) * 1.12)
    ax.set_ylim(-lim, lim)
    ax.set_ylabel("Electric potential (display mV)")
    ax.set_xlabel("Time (ms)")
    ax.text(-0.13, 1.02, "a)", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.tick_params(direction="in")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    ax.plot(freq_mhz[band], spectrum_mV[band], color="k", linewidth=1.25)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, max(30.0, float(np.nanmax(spectrum_mV[band])) * 1.08 if np.any(band) else 30.0))
    ax.set_ylabel("Electric potential spectrum (display mV)")
    ax.set_xlabel("Frequency (MHz)")
    ax.text(-0.13, 1.02, "b)", transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.tick_params(direction="in")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout(h_pad=1.8)
    fig.savefig(outdir / f"{name}.png", dpi=300)
    plt.close(fig)

    pd.DataFrame(
        {
            "time_s": t,
            "time_ms": t_ms,
            "electric_potential_display_mV": scaled,
            "display_scale_mV_per_model_unit": scale_mV_per_unit,
            "active_width_us": active_width_us,
            "display_processing": display_processing,
            "pulse_groups_above_20pct": pulse_groups_above_20pct,
        }
    ).to_csv(outdir / f"{name}_time_trace.csv", index=False)
    pd.DataFrame(
        {
            "frequency_Hz": freq_hz,
            "frequency_MHz": freq_mhz,
            "electric_potential_spectrum_display_mV": spectrum_mV,
            "spectrum_display_mV": spectrum_mV,
            "spectrum_type": spectrum_type,
            "spectrum_quantity": spectrum_quantity,
        }
    ).to_csv(outdir / f"{name}_frequency_spectrum.csv", index=False)


def _configure(candidate: DisplayCandidate) -> zero_model.ZeroOffsetSchakelConfig:
    cfg = zero_model.ZeroOffsetSchakelConfig()
    cfg.offset_D = 0.0
    cfg.receiver_z_min = candidate.receiver_z_min_m
    cfg.receiver_z_max = candidate.receiver_z_max_m
    cfg.receiver_spacing = candidate.receiver_spacing_m
    cfg.waveform_t_before = candidate.waveform_t_before_s
    cfg.waveform_t_after = candidate.waveform_t_after_s
    cfg.waveform_nt = candidate.waveform_nt
    return cfg


def run_display_search(
    input_path: str | Path,
    outdir: str | Path,
    snapshot_target_phi: float = 0.75,
    n_frequencies: int = 96,
    n_theta: int = 24,
    candidates: Iterable[DisplayCandidate] | None = None,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    clear_previous_display_outputs(outdir)
    df = zero_model.base.load_reactive_transport_table(input_path)
    snapshot_idx = zero_model.base.choose_snapshot(df, snapshot_target_phi)
    row = df.iloc[snapshot_idx]
    candidates = list(candidates or default_candidate_grid())
    rows = []
    best_payload = None
    for candidate in candidates:
        cfg = _configure(candidate)
        z, t, u = zero_model.synthesize_waveforms_schakel2011(
            row,
            cfg,
            n_frequencies=n_frequencies,
            n_theta=n_theta,
            integration_method="fixed",
        )
        score = score_display_candidate(z, t, u)
        rows.append({**asdict(candidate), **score})
        if best_payload is None or score["visibility_score"] > best_payload[0]["visibility_score"]:
            best_payload = (score, candidate, cfg, z, t, u)

    scores = pd.DataFrame(rows).sort_values("visibility_score", ascending=False)
    scores.to_csv(outdir / "display_candidate_scores.csv", index=False)
    if best_payload is None:
        raise RuntimeError("no display candidates were evaluated")

    _, best_candidate, best_cfg, z_best, t_best, u_best = best_payload
    t0_s = best_cfg.z_s / math.sqrt(best_cfg.K_fl / best_cfg.rho_fl)
    half_width_mm = best_candidate.receiver_z_max_m * 1.0e3
    spacing_mm = best_candidate.receiver_spacing_m * 1.0e3
    time_after_us = best_candidate.waveform_t_after_s * 1.0e6
    title = (
        f"Best display: +/-{half_width_mm:g} mm, dx={spacing_mm:g} mm, "
        f"T0..T0+{time_after_us:g} us; phi={row['Porosity']:.3f}"
    )
    plot_enhanced_waveform_gather(
        z_best,
        t_best,
        u_best,
        outdir,
        name="best_waveform_zoom",
        title=title,
        wiggle_scale_fraction=best_candidate.wiggle_scale_fraction,
        t0_s=t0_s,
    )
    pixpin_candidate = DisplayCandidate(
        name="pixpin_fluid_cm_style",
        receiver_z_min_m=-0.045,
        receiver_z_max_m=0.005,
        receiver_spacing_m=0.005,
        waveform_t_before_s=0.0,
        waveform_t_after_s=best_candidate.waveform_t_after_s,
        waveform_nt=best_candidate.waveform_nt,
        wiggle_scale_fraction=0.42,
    )
    pixpin_cfg = _configure(pixpin_candidate)
    z_pixpin, t_pixpin, u_pixpin = zero_model.synthesize_waveforms_schakel2011(
        row,
        pixpin_cfg,
        n_frequencies=n_frequencies,
        n_theta=n_theta,
        integration_method="fixed",
    )
    pixpin_plot_mask = z_pixpin <= 1.0e-12
    plot_pixpin_reference_style_gather(
        z_pixpin[pixpin_plot_mask],
        t_pixpin,
        u_pixpin[pixpin_plot_mask, :],
        outdir,
        name="pixpin_style_vertical_gather",
        title="",
    )
    np.savez(
        outdir / "pixpin_style_vertical_gather.npz",
        z=z_pixpin[pixpin_plot_mask],
        t=t_pixpin,
        U=u_pixpin[pixpin_plot_mask, :],
        receiver_z_min_m=pixpin_candidate.receiver_z_min_m,
        receiver_z_max_m=0.0,
        receiver_spacing_m=pixpin_candidate.receiver_spacing_m,
    )
    fig6_candidate = DisplayCandidate(
        name="schakel_fig6_trace",
        receiver_z_min_m=-0.001,
        receiver_z_max_m=-0.001,
        receiver_spacing_m=0.001,
        waveform_t_before_s=0.0,
        waveform_t_after_s=100.0e-6,
        waveform_nt=6000,
        wiggle_scale_fraction=1.0,
    )
    fig6_cfg = _configure(fig6_candidate)
    fig6_cfg.schakel_source_mode = "fig4_digitized"
    fig6_n_frequencies = strict_fig6_frequency_count(
        fig6_cfg,
        display_window_s=fig6_candidate.waveform_t_after_s,
        requested_n_frequencies=n_frequencies,
    )
    z_fig6, t_fig6, u_fig6 = zero_model.synthesize_waveforms_schakel2011(
        row,
        fig6_cfg,
        n_frequencies=fig6_n_frequencies,
        n_theta=n_theta,
        integration_method="fixed",
    )
    fig6_trace_idx = int(np.nanargmax(np.nanmax(np.abs(u_fig6), axis=1)))
    fig6_spectrum_n_frequencies = max(800, fig6_n_frequencies)
    response_freq_hz, response_spectrum = compute_fig6_reflected_potential_spectrum(
        row,
        fig6_cfg,
        float(z_fig6[fig6_trace_idx]),
        n_frequencies=fig6_spectrum_n_frequencies,
        n_theta=n_theta,
    )
    plot_schakel_fig6_style_panels(
        t_fig6,
        u_fig6[fig6_trace_idx, :],
        outdir,
        name="schakel_fig6_style_time_frequency",
        spectrum_frequency_hz=response_freq_hz,
        spectrum_amplitude=response_spectrum,
        time_xlim_s=(float(t_fig6[0]), float(t_fig6[-1])),
        spectrum_type="schakel_sommerfeld_reflected_electric_potential_spectrum",
        spectrum_quantity="reflected_electric_potential",
    )
    np.savez(
        outdir / "schakel_fig6_style_trace.npz",
        z=z_fig6,
        t=t_fig6,
        U=u_fig6,
        selected_z_m=float(z_fig6[fig6_trace_idx]),
        waveform_t_after_s=fig6_candidate.waveform_t_after_s,
        n_frequencies=int(fig6_n_frequencies),
        source_mode=fig6_cfg.schakel_source_mode,
    )
    fig6_freqs = zero_model._frequency_grid(fig6_cfg, fig6_n_frequencies)
    fig6_replica_period_us = (
        1.0e6 / float(fig6_freqs[1] - fig6_freqs[0]) if len(fig6_freqs) > 1 else np.nan
    )
    pd.DataFrame(
        [
            {
                "name": fig6_candidate.name,
                "receiver_z_m": float(z_fig6[fig6_trace_idx]),
                "waveform_t_after_s": fig6_candidate.waveform_t_after_s,
                "waveform_nt": fig6_candidate.waveform_nt,
                "source_mode": fig6_cfg.schakel_source_mode,
                "n_frequencies": int(fig6_n_frequencies),
                "spectrum_n_frequencies": int(fig6_spectrum_n_frequencies),
                "n_theta": int(n_theta),
                "frequency_replica_period_us": fig6_replica_period_us,
                "display_processing": "median_removed_only",
            }
        ]
    ).to_csv(outdir / "schakel_fig6_style_parameters.csv", index=False)
    np.savez(
        outdir / "best_waveform_zoom.npz",
        z=z_best,
        t=t_best,
        U=u_best,
        snapshot_index=snapshot_idx,
        snapshot_Time_s=float(row["Time_s"]),
        snapshot_Porosity=float(row["Porosity"]),
        best_candidate=np.array([best_candidate.name]),
    )
    pd.DataFrame([{**asdict(best_candidate), "snapshot_index": snapshot_idx,
                   "snapshot_Time_s": float(row["Time_s"]),
                   "snapshot_Porosity": float(row["Porosity"]),
                   "n_frequencies": int(n_frequencies),
                   "n_theta": int(n_theta)}]).to_csv(outdir / "best_display_parameters.csv", index=False)
    best_score_row = scores.iloc[0].to_dict()
    summary = f"""# Zero-Offset Waveform Display Search

Recommended display parameters:

- Receiver/electrode range: {best_candidate.receiver_z_min_m * 1.0e3:g} to {best_candidate.receiver_z_max_m * 1.0e3:g} mm relative to the interface.
- Electrode spacing: {best_candidate.receiver_spacing_m * 1.0e3:g} mm.
- Waveform window: T0 to T0 + {best_candidate.waveform_t_after_s * 1.0e6:g} us; no pre-T0 samples are displayed.
- Time samples: {best_candidate.waveform_nt}.
- Wiggle scale fraction: {best_candidate.wiggle_scale_fraction:g} of electrode spacing.
- Snapshot: index {snapshot_idx}, Time_s={float(row['Time_s']):.6g}, porosity={float(row['Porosity']):.6g}.

Why this candidate was selected:

- Visibility score: {float(best_score_row['visibility_score']):.6g}.
- Near/far peak ratio: {float(best_score_row['near_to_far_peak_ratio']):.6g}.
- Near-interface polarity reversal: {bool(best_score_row['polarity_reversal_near_interface'])}.
- Peak distance from interface: {float(best_score_row['peak_distance_from_interface_mm']):.6g} mm.
- Common near-interface arrival time: {float(best_score_row['common_arrival_time_s']) * 1.0e6:.6g} us.

The search changes only display geometry and plotting parameters.  The waveform
values are still computed by `synthesize_waveforms_schakel2011` in
`seismoelectric_zerooffset_schakel2011_sommerfeld.py`.

Normalization note:

- `best_waveform_zoom.png` and `best_waveform_zoom_heatmap.png` are side-normalized
  to make the waveform shape and polarity on both sides visible.
- `best_waveform_zoom_global_heatmap.png` uses one global normalization factor and
  should be used when comparing absolute amplitude scale across the interface.

Reference-style figures:

- `pixpin_style_vertical_gather.png` follows the supplied reference style with
  `z (cm)` on the top axis, time increasing downward, and fluid-side receivers
  from -4.5 to 0 cm at 0.5 cm spacing. It is trace-normalized to make the
  waveform shape visible at every electrode; the accompanying trace peak CSV
  keeps the raw model peak values. Solid black traces are model-only waveforms;
  no experimental overlay is included.
- `schakel_fig6_style_time_frequency.png` follows the Schakel Geophysics Fig. 6
  style: time is plotted in ms, which makes the microsecond waveform appear as a
  narrow pulse. The upper panel uses the modeled near-interface receiver trace
  after median baseline removal only; no moving-average high-pass display filter
  is applied. This figure uses the Schakel Fig. 4 visual source approximation
  mode and {fig6_n_frequencies} frequency samples, giving a frequency-sampling replica
  period of {fig6_replica_period_us:.6g} us. The lower panel shows the
  display-scaled reflected electric-potential spectrum at the same receiver,
  plotted versus MHz using {fig6_spectrum_n_frequencies} response-frequency
  samples. Solid black curves are model-only outputs.
"""
    (outdir / "display_search_summary.md").write_text(summary, encoding="utf-8")
    captions = """# Suggested Figure Captions

`pixpin_style_vertical_gather.png`

Zero-offset Schakel-Sommerfeld interface-EM waveform gather plotted in the
reference-style layout. Electrode position `z` is shown on the top axis in cm
and waveform time is shown vertically in us, increasing downward. Traces are
trace-normalized for display only, so this panel emphasizes waveform shape and
relative arrival pattern at each electrode; absolute model peak amplitudes are
listed in `pixpin_style_vertical_gather_trace_peaks.csv`. The solid black traces
are model-only waveforms and do not include an experimental-data overlay.

`schakel_fig6_style_time_frequency.png`

Schakel Fig. 6 style display. Panel (a) shows a modeled near-interface receiver
trace plotted in ms after median baseline removal and display scaling only; no
moving-average high-pass display filter is applied. Panel (b) shows the
display-scaled reflected electric-potential spectrum at the same receiver,
computed from the Schakel-Sommerfeld response using the Fig. 4 visual source
approximation; it is not the absolute spectrum of panel (a). Values labeled
`display mV` are plotting-scale units and should not be interpreted as
calibrated physical mV. The solid black curves are model-only outputs.
"""
    (outdir / "figure_captions.md").write_text(captions, encoding="utf-8")

    for rank, rec in enumerate(scores.head(3).itertuples(index=False), start=1):
        candidate = next(c for c in candidates if c.name == rec.name)
        cfg = _configure(candidate)
        z, t, u = zero_model.synthesize_waveforms_schakel2011(
            row,
            cfg,
            n_frequencies=n_frequencies,
            n_theta=n_theta,
            integration_method="fixed",
        )
        plot_enhanced_waveform_gather(
            z,
            t,
            u,
            outdir,
            name=f"candidate_rank{rank}_{candidate.name}",
            title=f"Rank {rank}: {candidate.name}",
            wiggle_scale_fraction=candidate.wiggle_scale_fraction,
            t0_s=t0_s,
        )
    return outdir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="global_evolution.xlsx")
    parser.add_argument("--outdir", default="se_results_zerooffset_schakel2011_display_search")
    parser.add_argument("--snapshot-target-phi", type=float, default=0.75)
    parser.add_argument("--n-frequencies", type=int, default=96)
    parser.add_argument("--n-theta", type=int, default=24)
    args = parser.parse_args()
    outdir = run_display_search(
        args.input,
        args.outdir,
        snapshot_target_phi=args.snapshot_target_phi,
        n_frequencies=args.n_frequencies,
        n_theta=args.n_theta,
    )
    print(f"Done. Display search outputs written to: {outdir}")


if __name__ == "__main__":
    main()
