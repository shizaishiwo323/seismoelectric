#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liu 2018 Fig. 2(b)-style comparison for the current seismoelectric model.

This script compares two quantities along a finite-offset receiver line:

1) "Current full formula model": the waveform produced by the current
   seismoelectric_offset_liu2018_spectral.py model.  In that model the Schakel
   R_E/T_TM spectral synthesis is combined with Liu Eq. (4) radiation geometry
   for finite-offset receiver-line amplitudes.
2) "Dipole explanation model": Liu Eq. (4), u proportional to cos(theta)/r^2.

For audit, the script also saves a diagnostic waveform synthesized without the
final dipole-geometry weighting.  That diagnostic is not used as the Fig. 2(b)
main comparison because it is missing the radiation-directivity layer that Liu
uses to explain the finite-offset amplitude pattern.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


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


def synthesize_waveforms_spectral_without_dipole(se, row: pd.Series, cfg,
                                                 n_omega: int,
                                                 n_k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Copy of the Liu spectral synthesis without final dipole-geometry weighting."""
    z_receivers = np.arange(
        cfg.receiver_z_min,
        cfg.receiver_z_max + 0.5 * cfg.receiver_spacing,
        cfg.receiver_spacing,
    )
    vf = math.sqrt(cfg.K_fl / cfg.rho_fl)
    t0 = cfg.z_s / vf
    t = np.linspace(max(0.0, t0 - cfg.waveform_t_before), t0 + cfg.waveform_t_after, cfg.waveform_nt)
    x = float(cfg.offset_D)

    phi = float(np.clip(row["Porosity"], cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    tau = max(float(row["Tortuosity"]), 1.0 + 1e-6)
    cH = float(row["OutletHConc"])
    c_override = se.optional_float(row, "ElectrolyteConcentration_molL")
    sigma_f_override = se.optional_float(row, "FluidConductivity_S_m")

    f_min = max(1.0, cfg.spectral_f_min_factor * cfg.f0)
    f_max = max(f_min * 1.01, cfg.spectral_f_max_factor * cfg.f0)
    omega_grid = 2.0 * math.pi * np.linspace(f_min, f_max, int(n_omega))
    domega = float(omega_grid[1] - omega_grid[0])
    omega_weights = se._trapz_weights(int(n_omega), domega)

    theta0 = math.radians(cfg.source_beam_theta_deg)
    kb = max(float(cfg.source_kb_m_inv), 1e-12)
    source_spectrum = se.causal_ricker_source_spectrum(omega_grid, cfg)
    inv_2pi2 = 1.0 / (2.0 * math.pi) ** 2
    U = np.zeros((len(z_receivers), len(t)), dtype=complex)

    for iw, omega in enumerate(omega_grid):
        k_ac = omega / vf
        k_center = k_ac * math.sin(theta0)
        k_lim = max(1e-12, cfg.spectral_k_limit_factor * k_ac)
        k_grid = np.linspace(-k_lim, k_lim, int(n_k))
        dk = float(k_grid[1] - k_grid[0])
        k_weights = se._trapz_weights(int(n_k), dk)
        phase_t = np.exp(-1j * omega * t)

        for ik, kx in enumerate(k_grid):
            Akw = source_spectrum[iw] * np.exp(-((kx - k_center) / kb) ** 2)
            if not np.isfinite(Akw) or abs(Akw) < 1e-14:
                continue
            try:
                coeff = se.se_coefficients(
                    phi,
                    k0_m2,
                    tau,
                    cH,
                    omega,
                    None,
                    cfg,
                    kx_override=float(kx),
                    C_override_molL=c_override,
                    sigma_f_override=sigma_f_override,
                )
            except Exception:
                continue

            RE = coeff["R_E"]
            TTM = coeff["T_TM"]
            if not (
                np.isfinite(RE.real)
                and np.isfinite(RE.imag)
                and np.isfinite(TTM.real)
                and np.isfinite(TTM.imag)
            ):
                continue

            kzi = coeff["k3_fl"]
            kEr = coeff["k3_E"]
            kEt = coeff["k3_TM"]
            weight = 2.0 * Akw * omega_weights[iw] * k_weights[ik] * inv_2pi2
            common_x = np.exp(1j * kx * x)

            z_neg = z_receivers < 0
            z_pos = ~z_neg
            amp_z = np.empty(len(z_receivers), dtype=complex)
            amp_z[z_neg] = RE * np.exp(1j * (kzi * cfg.z_s - kEr * z_receivers[z_neg])) * common_x
            amp_z[z_pos] = TTM * np.exp(1j * (kzi * cfg.z_s + kEt * z_receivers[z_pos])) * common_x
            U += weight * amp_z[:, None] * phase_t[None, :]

    return z_receivers, t, np.real(U)


def build_comparison_table(z: np.ndarray, U: np.ndarray, cfg, prefix: str) -> pd.DataFrame:
    peak_idx = np.nanargmax(np.abs(U), axis=1)
    signed_peak = U[np.arange(len(z)), peak_idx]
    peak_abs = np.abs(signed_peak)

    D = abs(float(cfg.offset_D))
    r = np.sqrt(D**2 + z**2)
    theta_deg = np.degrees(np.arccos(np.divide(z, r, out=np.zeros_like(z), where=r > 0)))
    dipole_signed = np.divide(z, r**3, out=np.zeros_like(z), where=r > 0)

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
            label=r"Dipole model $|\cos\theta|/r^2$")
    ax.plot(tab["theta_deg"], tab["current_abs_side_norm"], "o-", color="0.15", linewidth=1.4,
            markersize=3.6, label="Current full formula model")
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
            label=r"Dipole model $\cos\theta/r^2$")
    ax.plot(tab["theta_deg"], tab["current_signed_side_norm"], "o-", color="0.15", linewidth=1.4,
            markersize=3.6, label="Current full formula model")
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
    fig = plt.figure(figsize=(6.4, 6.0))
    ax = fig.add_subplot(111, projection="polar")
    theta_rad = np.radians(polar["theta_deg"].to_numpy(float))
    ax.plot(theta_rad_full, dipole_full_abs, color="tab:red", linewidth=2.0,
            label="Dipole")
    ax.plot(theta_rad, polar["current_abs_side_norm"], "o-", color="0.15", linewidth=1.2,
            markersize=3.2, label="Current full model")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    ax.set_rlim(0, 1.05)
    ax.set_title("Normalized amplitude directivity")
    ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.05))
    fig.tight_layout()
    fig.savefig(outdir / "liu2018_fig2b_polar_directivity.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.plot(theta_full, dipole_full_abs, color="tab:red", linewidth=2.0,
            label=r"Dipole model $|\cos\theta|/r^2$")
    ax.plot(tab["theta_deg"], tab["without_dipole_abs_side_norm"], "o-", color="tab:purple",
            linewidth=1.2, markersize=3.0, label="Spectral kernel without dipole layer")
    ax.axvline(90.0, color="tab:blue", linestyle=":", linewidth=1.2)
    ax.set_xlabel(r"Radiation angle $\theta$ (deg)")
    ax.set_ylabel("Normalized peak amplitude")
    ax.set_title("Diagnostic: spectral kernel before Liu dipole radiation layer")
    ax.set_xlim(0, 180)
    ax.set_ylim(-0.03, 1.08)
    ax.grid(True, alpha=0.28)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "liu2018_fig2b_without_dipole_diagnostic.png", dpi=300)
    plt.close(fig)


def main() -> None:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-script", type=str, default=str(here / "seismoelectric_offset_liu2018_spectral.py"))
    parser.add_argument("--input", type=str, default=str(here / "global_evolution.xlsx"))
    parser.add_argument("--outdir", type=str, default=str(here / "liu2018_fig2b_comparison"))
    parser.add_argument("--snapshot-target-phi", type=float, default=0.75)
    parser.add_argument("--offset-D-mm", type=float, default=45.0)
    parser.add_argument("--receiver-z-min-mm", type=float, default=-80.0)
    parser.add_argument("--receiver-z-max-mm", type=float, default=80.0)
    parser.add_argument("--receiver-spacing-mm", type=float, default=2.5)
    parser.add_argument("--spectral-n-omega", type=int, default=80)
    parser.add_argument("--spectral-n-k", type=int, default=81)
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

    z, t, U_without = synthesize_waveforms_spectral_without_dipole(
        se, row, cfg, n_omega=args.spectral_n_omega, n_k=args.spectral_n_k
    )
    z_current, t_current, U_current = se.synthesize_waveforms_spectral(
        row, cfg, n_omega=args.spectral_n_omega, n_k=args.spectral_n_k
    )
    if not np.allclose(z, z_current):
        raise RuntimeError("Current-model and diagnostic receiver grids do not match.")

    tab_current = build_comparison_table(z_current, U_current, cfg, prefix="current")
    tab_without = build_comparison_table(z, U_without, cfg, prefix="without_dipole")
    keep_without = [
        "theta_deg",
        "without_dipole_signed_peak",
        "without_dipole_peak_abs",
        "without_dipole_signed_norm",
        "without_dipole_abs_norm",
        "without_dipole_signed_side_norm",
        "without_dipole_abs_side_norm",
    ]
    tab = tab_current.merge(tab_without[keep_without], on="theta_deg", how="left")
    tab["snapshot_index"] = idx
    tab["snapshot_Time_s"] = float(row["Time_s"])
    tab["snapshot_Porosity"] = float(row["Porosity"])
    tab.to_csv(outdir / "liu2018_fig2b_comparison.csv", index=False)
    np.savez_compressed(outdir / "liu2018_fig2b_waveforms_current_model.npz", z_m=z_current, t_s=t_current, U=U_current)
    np.savez_compressed(outdir / "liu2018_fig2b_waveforms_without_dipole.npz", z_m=z, t_s=t, U=U_without)
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
        "note": "Main comparison uses the current full model; without-dipole spectral kernel is saved as a diagnostic.",
    }
    pd.Series(summary).to_csv(outdir / "run_summary.csv")

    print("Done. Outputs written to:", outdir)
    print("Snapshot index:", idx, "Time_s:", float(row["Time_s"]), "Porosity:", float(row["Porosity"]))
    print("Peak current-model normalized at theta:")
    imax = int(np.nanargmax(tab["current_abs_norm"].to_numpy(float)))
    print(tab.iloc[imax][["theta_deg", "z_mm", "current_abs_side_norm", "dipole_abs_side_norm"]].to_string())


if __name__ == "__main__":
    main()
