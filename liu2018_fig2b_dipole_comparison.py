#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liu 2018 Fig. 2(b)-style comparison for the current seismoelectric model.

This script compares the two quantities used in Liu et al. (2018) Fig. 2(b):

1) "Pride theory": the full electrical-potential waveform produced by
   seismoelectric_offset_liu2018_spectral.py using Liu Eq. (1)-(2) spectral
   synthesis and Schakel/Pride seismoelectric coupling coefficients.
2) "Electric dipole": the independent Liu Eq. (4) explanation model,
   electric potential proportional to cos(theta)/r^2.

The dipole factor is not multiplied into the Pride-theory waveform. Liu uses it
to explain the amplitude pattern of the full modeling result, not to generate
that result.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_SPECTRAL_N_OMEGA = 48
DEFAULT_SPECTRAL_N_K = 401
PRIDE_COLUMN_PREFIX = "pride"


def load_model_module(script_path: Path):
    spec = importlib.util.spec_from_file_location("se_offset_model", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import model script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["se_offset_model"] = module
    spec.loader.exec_module(module)
    return module


def signed_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    m = np.nanmax(np.abs(x))
    if np.isfinite(m) and m > 0:
        return x / m
    return np.full_like(x, np.nan, dtype=float)


def side_abs_normalize(values: np.ndarray, z: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    for mask in (z < 0, z > 0):
        if not np.any(mask):
            continue
        m = np.nanmax(np.abs(values[mask]))
        if np.isfinite(m) and m > 0:
            out[mask] = np.abs(values[mask]) / m
    out[np.abs(z) < 1e-15] = 0.0
    return out


def side_signed_normalize(values: np.ndarray, z: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    for mask in (z < 0, z > 0):
        if not np.any(mask):
            continue
        m = np.nanmax(np.abs(values[mask]))
        if np.isfinite(m) and m > 0:
            out[mask] = values[mask] / m
    out[np.abs(z) < 1e-15] = 0.0
    return out


def liu_receiver_line_dipole_geometry(z: np.ndarray, offset_D: float) -> dict[str, np.ndarray]:
    """Return Liu Eq. (4) geometry for a finite-offset vertical receiver line."""
    z = np.asarray(z, dtype=float)
    D = abs(float(offset_D))
    r = np.sqrt(D**2 + z**2)
    theta_deg = np.degrees(np.arccos(np.divide(z, r, out=np.zeros_like(z), where=r > 0)))
    dipole_signed = np.divide(z, r**3, out=np.zeros_like(z), where=r > 0)
    return {
        "r": r,
        "theta_deg": theta_deg,
        "dipole_signed": dipole_signed,
    }


def configure_liu_semicircle_axis(ax) -> None:
    """Configure a polar axis like Liu Fig. 2(b): upper semicircle, diameter below."""
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_rlim(0, 1.05)


def build_comparison_table(z: np.ndarray, U: np.ndarray, cfg, prefix: str) -> pd.DataFrame:
    peak_idx = np.nanargmax(np.abs(U), axis=1)
    signed_peak = U[np.arange(len(z)), peak_idx]
    peak_abs = np.abs(signed_peak)

    geometry = liu_receiver_line_dipole_geometry(z, cfg.offset_D)
    r = geometry["r"]
    theta_deg = geometry["theta_deg"]
    dipole_signed = geometry["dipole_signed"]

    modeled_signed_norm = signed_normalize(signed_peak)
    dipole_signed_norm = signed_normalize(dipole_signed)

    modeled_abs_norm = peak_abs / np.nanmax(peak_abs) if np.nanmax(peak_abs) > 0 else np.nan
    dipole_abs = np.abs(dipole_signed)
    dipole_abs_norm = dipole_abs / np.nanmax(dipole_abs) if np.nanmax(dipole_abs) > 0 else np.nan
    modeled_abs_side_norm = side_abs_normalize(signed_peak, z)
    modeled_signed_side_norm = side_signed_normalize(signed_peak, z)
    dipole_abs_side_norm = side_abs_normalize(dipole_signed, z)
    dipole_signed_side_norm = side_signed_normalize(dipole_signed, z)

    return pd.DataFrame({
        "z_m": z,
        "z_mm": z * 1e3,
        "offset_D_mm": cfg.offset_D * 1e3,
        "r_mm": r * 1e3,
        "theta_deg": theta_deg,
        f"{prefix}_signed_peak": signed_peak,
        f"{prefix}_peak_abs": peak_abs,
        f"{prefix}_signed_norm": modeled_signed_norm,
        f"{prefix}_abs_norm": modeled_abs_norm,
        f"{prefix}_signed_side_norm": modeled_signed_side_norm,
        f"{prefix}_abs_side_norm": modeled_abs_side_norm,
        "dipole_signed": dipole_signed,
        "dipole_abs": dipole_abs,
        "dipole_signed_norm": dipole_signed_norm,
        "dipole_abs_norm": dipole_abs_norm,
        "dipole_signed_side_norm": dipole_signed_side_norm,
        "dipole_abs_side_norm": dipole_abs_side_norm,
        "side": np.where(z < 0, "fluid/reflected", np.where(z > 0, "porous/transmitted", "interface")),
    }).sort_values("theta_deg")


def plot_comparison(tab: pd.DataFrame, row: pd.Series, outdir: Path) -> None:
    theta_full = np.linspace(0.0, 180.0, 721)
    # Liu Eq. (4) sampled over the full radiation-angle range for the
    # explanation model.  Along a finite-offset receiver line, r = D/sin(theta),
    # so |cos(theta)|/r^2 is proportional to |cos(theta)| sin^2(theta).
    theta_rad_full = np.radians(theta_full)
    dipole_full_abs = np.abs(np.cos(theta_rad_full)) * np.sin(theta_rad_full) ** 2
    dipole_full_signed = np.cos(theta_rad_full) * np.sin(theta_rad_full) ** 2
    dipole_full_abs = dipole_full_abs / np.nanmax(dipole_full_abs)
    dipole_full_signed = dipole_full_signed / np.nanmax(np.abs(dipole_full_signed))

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.plot(theta_full, dipole_full_abs, color="tab:red", linewidth=2.0,
            label=r"Electric dipole $|\cos\theta|/r^2$")
    ax.plot(tab["theta_deg"], tab["pride_abs_side_norm"], "o-", color="0.15", linewidth=1.4,
            markersize=3.6, label="Pride theory")
    ax.axvline(90.0, color="tab:blue", linestyle=":", linewidth=1.2)
    ax.set_xlabel(r"Radiation angle $\theta$ (deg)")
    ax.set_ylabel("Side-normalized peak amplitude")
    ax.set_title(f"Liu Fig. 2b-style amplitude comparison, phi={float(row['Porosity']):.3f}")
    ax.set_xlim(0, 180)
    ax.set_ylim(-0.03, 1.08)
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "liu2018_fig2b_amplitude_vs_angle.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.plot(theta_full, dipole_full_signed, color="tab:red", linewidth=2.0,
            label=r"Electric dipole $\cos\theta/r^2$")
    ax.plot(tab["theta_deg"], tab["pride_signed_side_norm"], "o-", color="0.15", linewidth=1.4,
            markersize=3.6, label="Pride theory")
    ax.axhline(0.0, color="0.4", linewidth=0.9)
    ax.axvline(90.0, color="tab:blue", linestyle=":", linewidth=1.2)
    ax.set_xlabel(r"Radiation angle $\theta$ (deg)")
    ax.set_ylabel("Side-normalized signed peak")
    ax.set_title(f"Signed polarity comparison, phi={float(row['Porosity']):.3f}")
    ax.set_xlim(0, 180)
    ax.set_ylim(-1.08, 1.08)
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "liu2018_fig2b_signed_vs_angle.png", dpi=300)
    plt.close(fig)

    polar = tab[np.isfinite(tab["theta_deg"])].copy()
    fig = plt.figure(figsize=(7.6, 4.8))
    ax = fig.add_subplot(111, projection="polar")
    theta_rad = np.radians(polar["theta_deg"].to_numpy(float))
    ax.plot(theta_rad_full, dipole_full_abs, color="tab:red", linewidth=2.0,
            label="Electric dipole")
    ax.plot(theta_rad, polar["pride_abs_side_norm"], color="0.15", linewidth=1.7,
            label="Pride theory")
    configure_liu_semicircle_axis(ax)
    ax.set_title("Liu Fig. 2(b): normalized electrical-potential amplitude", pad=18)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18), ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "liu2018_fig2b_polar_directivity.png", dpi=300)
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-script", type=str, default=str(here / "seismoelectric_offset_liu2018_spectral.py"))
    parser.add_argument("--input", type=str, default=str(here / "global_evolution.xlsx"))
    parser.add_argument("--outdir", type=str, default=str(here / "liu2018_fig2b_comparison"))
    parser.add_argument("--snapshot-target-phi", type=float, default=0.75)
    parser.add_argument("--offset-D-mm", type=float, default=45.0)
    parser.add_argument("--receiver-z-min-mm", type=float, default=-100.0)
    parser.add_argument("--receiver-z-max-mm", type=float, default=100.0)
    parser.add_argument("--receiver-spacing-mm", type=float, default=1.0)
    parser.add_argument("--spectral-n-omega", type=int, default=DEFAULT_SPECTRAL_N_OMEGA)
    parser.add_argument("--spectral-n-k", type=int, default=DEFAULT_SPECTRAL_N_K)
    args = parser.parse_args()

    model_script = Path(args.model_script)
    se = load_model_module(model_script)
    cfg = se.SEConfig()
    cfg.offset_D = args.offset_D_mm * 1e-3
    cfg.receiver_z_min = args.receiver_z_min_mm * 1e-3
    cfg.receiver_z_max = args.receiver_z_max_mm * 1e-3
    cfg.receiver_spacing = args.receiver_spacing_mm * 1e-3
    cfg.waveform_nt = 1200

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = se.load_reactive_transport_table(args.input)
    idx = se.choose_snapshot(df, target_phi=args.snapshot_target_phi)
    row = df.iloc[idx]

    z_pride, t_pride, U_pride = se.synthesize_waveforms_spectral(
        row, cfg, n_omega=args.spectral_n_omega, n_k=args.spectral_n_k
    )

    tab = build_comparison_table(z_pride, U_pride, cfg, prefix=PRIDE_COLUMN_PREFIX)
    tab["snapshot_index"] = idx
    tab["snapshot_Time_s"] = float(row["Time_s"])
    tab["snapshot_Porosity"] = float(row["Porosity"])
    tab.to_csv(outdir / "liu2018_fig2b_comparison.csv", index=False)
    np.savez_compressed(outdir / "liu2018_fig2b_waveforms_pride_theory.npz", z_m=z_pride, t_s=t_pride, U=U_pride)
    plot_comparison(tab, row, outdir)

    summary = {
        "model_script": str(model_script),
        "input": str(args.input),
        "outdir": str(outdir),
        "snapshot_index": int(idx),
        "snapshot_Time_s": float(row["Time_s"]),
        "snapshot_Porosity": float(row["Porosity"]),
        "offset_D_mm": args.offset_D_mm,
        "receiver_z_min_mm": args.receiver_z_min_mm,
        "receiver_z_max_mm": args.receiver_z_max_mm,
        "receiver_spacing_mm": args.receiver_spacing_mm,
        "spectral_n_omega": args.spectral_n_omega,
        "spectral_n_k": args.spectral_n_k,
        "note": "Pride theory is the reference script's full Liu Eq. (1)-(2) spectral model; Liu Eq. (4) is an independent comparison curve.",
    }
    pd.Series(summary).to_csv(outdir / "run_summary.csv")

    print("Done. Outputs written to:", outdir)
    print("Snapshot index:", idx, "Time_s:", float(row["Time_s"]), "Porosity:", float(row["Porosity"]))
    print("Peak Pride-theory normalized at theta:")
    imax = int(np.nanargmax(tab["pride_abs_norm"].to_numpy(float)))
    print(tab.iloc[imax][["theta_deg", "z_mm", "pride_abs_side_norm", "dipole_abs_side_norm"]].to_string())


if __name__ == "__main__":
    main()
