#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone parameter sensitivity analysis for the seismoelectric forward model.

The script keeps two analysis branches separate:

1. Schakel-style reproduction scans using the Table I reference settings and
   the current Schakel boundary-value solver.
2. Reactive-transport research-data sensitivity, including one-at-a-time
   parameter substitutions and a signed waveform-peak contribution bar chart.

The Schakel paper sensitivity figures are vertical energy-flux coefficients.
The current main model exposes the conversion potentials R_E and T_TM, so this
script reports potential magnitudes and squared-magnitude energy proxies rather
than claiming to reproduce the full orthodox/interference flux decomposition.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import replace
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seismoelectric_offset_liu2018_spectral as se


PARAMETER_GROUPS: Mapping[str, Tuple[str, ...]] = {
    "porosity": ("Porosity",),
    "permeability": ("Permeability_mD",),
    "tortuosity": ("Tortuosity",),
    "fluid_chemistry": (
        "OutletHConc",
        "ElectrolyteConcentration_molL",
        "FluidConductivity_S_m",
    ),
}


def h_concentration_for_ph(pH: float, cfg: se.SEConfig) -> float:
    """Return an OutletHConc value that produces the requested pH."""
    c_mol_l = 10.0 ** (-float(pH))
    if cfg.outlet_h_unit == "mol_cm3":
        return c_mol_l / 1000.0
    return c_mol_l


def schakel_reference_config() -> Tuple[se.SEConfig, Dict[str, float]]:
    """Reference settings following Schakel and Smeulders Table I and Fig. 4-7."""
    cfg = se.SEConfig()
    cfg.K_b = 5.8e9
    cfg.G = 3.4e9
    cfg.K_s = 40.0e9
    cfg.K_f = 2.22e9
    cfg.K_fl = 2.22e9
    cfg.eta = 1.0e-3
    cfg.rho_f = 1000.0
    cfg.rho_fl = 1000.0
    cfg.rho_s = 2760.0
    cfg.eps_f = 80.0
    cfg.eps_s = 4.0
    cfg.eps_fl = 80.0
    cfg.C_background_molL = 1.0e-3
    cfg.temperature = 295.0
    cfg.M_similarity = 1.0
    cfg.coeff_theta_deg = 45.0
    cfg.upper_fluid_conductivity_mode = "constant"

    baseline = {
        "phi": 0.24,
        "k0_m2": 0.390e-12,
        "alpha_inf": 2.3,
        "pH": 7.0,
        "cH": h_concentration_for_ph(7.0, cfg),
        "C_molL": 1.0e-3,
        "omega": 1.0e6,
        "theta_deg": 45.0,
    }
    return cfg, baseline


def coefficient_record(
    phi: float,
    k0_m2: float,
    alpha_inf: float,
    cH: float,
    omega: float,
    theta_deg: float,
    cfg: se.SEConfig,
    C_override_molL: float | None = None,
    sigma_f_override: float | None = None,
) -> Dict[str, float]:
    coeff = se.se_coefficients(
        phi,
        k0_m2,
        alpha_inf,
        cH,
        omega,
        theta_deg,
        cfg,
        C_override_molL=C_override_molL,
        sigma_f_override=sigma_f_override,
    )
    re_abs = abs(coeff["R_E"])
    ttm_abs = abs(coeff["T_TM"])
    return {
        "Porosity_used": float(phi),
        "k0_m2": float(k0_m2),
        "Permeability_mD": float(k0_m2 / 9.869233e-16),
        "Tortuosity": float(alpha_inf),
        "OutletHConc": float(cH),
        "C_molL": float(coeff["C_molL"]),
        "pH": float(coeff["pH"]),
        "eta_Pa_s": float(cfg.eta),
        "K_f_Pa": float(cfg.K_f),
        "omega_rad_s": float(omega),
        "theta_deg": float(theta_deg),
        "Lambda_m": float(coeff["Lambda"]),
        "debye_d_m": float(coeff["debye_d"]),
        "L_abs": abs(coeff["L"]),
        "sigma_abs": abs(coeff["sigma"]),
        "RE_real": float(np.real(coeff["R_E"])),
        "RE_imag": float(np.imag(coeff["R_E"])),
        "RE_abs": float(re_abs),
        "TTM_real": float(np.real(coeff["T_TM"])),
        "TTM_imag": float(np.imag(coeff["T_TM"])),
        "TTM_abs": float(ttm_abs),
        "RE_energy_proxy": float(re_abs**2),
        "TTM_energy_proxy_abs": float(ttm_abs**2),
        "matrix_cond": float(coeff["matrix_cond"]),
    }


def schakel_scan_definitions() -> Dict[str, np.ndarray]:
    return {
        "electrolyte_concentration_molL": np.logspace(-6, 0, 80),
        "viscosity_Pa_s": np.logspace(-6, 0, 80),
        "permeability_m2": np.logspace(-16, -10, 100),
        "tortuosity": np.linspace(1.0, 5.0, 80),
        "pH": np.linspace(2.0, 10.0, 80),
        "pore_fluid_bulk_modulus_Pa": np.logspace(7, 11, 80),
        "porosity": np.linspace(0.10, 0.85, 80),
    }


def run_single_schakel_scan(
    scan_name: str,
    values: Sequence[float],
    cfg: se.SEConfig,
    baseline: Mapping[str, float],
) -> pd.DataFrame:
    rows: List[Dict[str, float]] = []
    for value in values:
        local_cfg = replace(cfg)
        phi = baseline["phi"]
        k0_m2 = baseline["k0_m2"]
        alpha_inf = baseline["alpha_inf"]
        cH = baseline["cH"]
        C_override = baseline["C_molL"]

        if scan_name == "electrolyte_concentration_molL":
            C_override = float(value)
        elif scan_name == "viscosity_Pa_s":
            local_cfg.eta = float(value)
        elif scan_name == "permeability_m2":
            k0_m2 = float(value)
        elif scan_name == "tortuosity":
            alpha_inf = float(value)
        elif scan_name == "pH":
            cH = h_concentration_for_ph(float(value), local_cfg)
        elif scan_name == "pore_fluid_bulk_modulus_Pa":
            local_cfg.K_f = float(value)
        elif scan_name == "porosity":
            phi = float(value)
        else:
            raise ValueError(f"Unknown Schakel scan: {scan_name}")

        record = coefficient_record(
            phi,
            k0_m2,
            alpha_inf,
            cH,
            baseline["omega"],
            baseline["theta_deg"],
            local_cfg,
            C_override_molL=C_override,
        )
        record["scan"] = scan_name
        record["scan_value"] = float(value)
        rows.append(record)
    return pd.DataFrame(rows)


def normalize_scan_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for scan_name, group in out.groupby("scan"):
        for col in ["RE_abs", "TTM_abs", "RE_energy_proxy", "TTM_energy_proxy_abs", "L_abs"]:
            ref = group[col].iloc[0]
            norm_col = col + "_norm_to_scan_start"
            idx = out["scan"] == scan_name
            out.loc[idx, norm_col] = group[col] / ref if np.isfinite(ref) and ref != 0 else np.nan
    return out


def plot_schakel_two_by_two(
    df: pd.DataFrame,
    scan_names: Sequence[str],
    outpath: Path,
    title: str,
    xlabels: Mapping[str, str],
    xscale: Mapping[str, str],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.5), constrained_layout=True)
    metrics = [("RE_energy_proxy", "R_E energy proxy"), ("TTM_energy_proxy_abs", "T_TM energy proxy")]
    for col, scan_name in enumerate(scan_names):
        part = df[df["scan"] == scan_name]
        for row, (metric, ylabel) in enumerate(metrics):
            ax = axes[row, col]
            ax.plot(part["scan_value"], part[metric], color="#1f77b4", linewidth=1.8)
            ax.set_xlabel(xlabels[scan_name])
            ax.set_ylabel(ylabel)
            if xscale.get(scan_name) == "log":
                ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="both", alpha=0.25)
    fig.suptitle(title)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_schakel_porosity(df: pd.DataFrame, outpath: Path) -> None:
    part = df[df["scan"] == "porosity"]
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.7), constrained_layout=True)
    axes[0].plot(part["scan_value"], part["RE_energy_proxy"], linewidth=1.8)
    axes[1].plot(part["scan_value"], part["TTM_energy_proxy_abs"], linewidth=1.8)
    axes[0].set_ylabel("R_E energy proxy")
    axes[1].set_ylabel("T_TM energy proxy")
    for ax in axes:
        ax.set_xlabel("porosity")
        ax.set_yscale("log")
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle("Schakel-style porosity sensitivity")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_permeability_competition(df: pd.DataFrame, outpath: Path) -> None:
    part = df[df["scan"] == "permeability_m2"]
    fig, ax1 = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    ax2 = ax1.twinx()
    ax1.plot(part["scan_value"], part["L_abs"], color="#1f77b4", linewidth=1.8, label="|L|")
    ax1.plot(part["scan_value"], part["RE_energy_proxy"], color="#2ca02c", linewidth=1.5, label="R_E proxy")
    ax2.plot(part["scan_value"], part["Lambda_m"], color="#d62728", linewidth=1.5, label="Lambda")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax1.set_xlabel("permeability k0 (m^2)")
    ax1.set_ylabel("|L| and R_E proxy")
    ax2.set_ylabel("pore volume-to-surface ratio Lambda (m)")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="best")
    ax1.grid(True, which="both", alpha=0.25)
    fig.suptitle("Permeability and Lambda competition")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_schakel_reproduction(outdir: Path) -> pd.DataFrame:
    cfg, baseline = schakel_reference_config()
    scans = schakel_scan_definitions()
    frames = [run_single_schakel_scan(name, values, cfg, baseline) for name, values in scans.items()]
    out = normalize_scan_columns(pd.concat(frames, ignore_index=True))
    out.to_csv(outdir / "schakel_sensitivity_coefficients.csv", index=False)

    summary_rows = []
    for scan_name, part in out.groupby("scan"):
        for metric in ["RE_energy_proxy", "TTM_energy_proxy_abs", "L_abs"]:
            idx = part[metric].idxmax()
            summary_rows.append({
                "scan": scan_name,
                "metric": metric,
                "min": part[metric].min(),
                "max": part[metric].max(),
                "max_over_min": part[metric].max() / part[metric].min() if part[metric].min() > 0 else np.nan,
                "scan_value_at_max": out.loc[idx, "scan_value"],
            })
    pd.DataFrame(summary_rows).to_csv(outdir / "schakel_sensitivity_summary.csv", index=False)

    xlabels = {
        "electrolyte_concentration_molL": "electrolyte concentration C (mol/L)",
        "viscosity_Pa_s": "viscosity eta (Pa s)",
        "permeability_m2": "permeability k0 (m^2)",
        "tortuosity": "tortuosity alpha_inf",
        "pH": "pH",
        "pore_fluid_bulk_modulus_Pa": "pore-fluid bulk modulus K_f (Pa)",
    }
    xscale = {
        "electrolyte_concentration_molL": "log",
        "viscosity_Pa_s": "log",
        "permeability_m2": "log",
        "pore_fluid_bulk_modulus_Pa": "log",
    }
    plot_schakel_two_by_two(
        out,
        ["electrolyte_concentration_molL", "viscosity_Pa_s"],
        outdir / "fig4_like_concentration_viscosity.png",
        "Schakel-style sensitivity: concentration and viscosity",
        xlabels,
        xscale,
    )
    plot_schakel_two_by_two(
        out,
        ["permeability_m2", "tortuosity"],
        outdir / "fig5_like_permeability_tortuosity.png",
        "Schakel-style sensitivity: permeability and tortuosity",
        xlabels,
        xscale,
    )
    plot_schakel_two_by_two(
        out,
        ["pH", "pore_fluid_bulk_modulus_Pa"],
        outdir / "fig6_like_ph_bulk_modulus.png",
        "Schakel-style sensitivity: pH and pore-fluid bulk modulus",
        xlabels,
        xscale,
    )
    plot_schakel_porosity(out, outdir / "fig7_like_porosity.png")
    plot_permeability_competition(out, outdir / "permeability_lambda_competition.png")
    return out


def pore_volume_to_surface_ratio_m(row: pd.Series) -> float:
    """Estimate pore-volume/surface-area ratio from RTM grain volume and surface area."""
    try:
        phi = float(row["Porosity"])
        grain_volume_cm3 = float(row["GrainVolume_cm3"])
        surface_area_cm2 = float(row["SurfaceArea_cm2"])
    except Exception:
        return np.nan
    if not np.isfinite(phi) or not (0.0 < phi < 1.0):
        return np.nan
    if not np.isfinite(grain_volume_cm3) or grain_volume_cm3 <= 0.0:
        return np.nan
    if not np.isfinite(surface_area_cm2) or surface_area_cm2 <= 0.0:
        return np.nan
    total_volume_cm3 = grain_volume_cm3 / (1.0 - phi)
    pore_volume_cm3 = phi * total_volume_cm3
    return (pore_volume_cm3 / surface_area_cm2) * 1.0e-2


def make_one_at_a_time_row(base: pd.Series, target: pd.Series, parameter_group: str) -> pd.Series:
    if parameter_group not in PARAMETER_GROUPS:
        raise ValueError(f"Unknown parameter group: {parameter_group}")
    row = base.copy()
    for col in PARAMETER_GROUPS[parameter_group]:
        if col in target.index:
            row[col] = target[col]
    return row


def row_coefficient_metrics(row: pd.Series, cfg: se.SEConfig) -> Dict[str, float]:
    phi = float(np.clip(float(row["Porosity"]), cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    alpha_inf = max(float(row["Tortuosity"]), 1.0 + 1e-6)
    cH = float(row["OutletHConc"])
    C_override = se.optional_float(row, "ElectrolyteConcentration_molL")
    sigma_f_override = se.optional_float(row, "FluidConductivity_S_m")
    rec = coefficient_record(
        phi,
        k0_m2,
        alpha_inf,
        cH,
        2.0 * math.pi * cfg.f0,
        cfg.coeff_theta_deg,
        cfg,
        C_override_molL=C_override,
        sigma_f_override=sigma_f_override,
    )
    rec["coefficient_combined_log_metric"] = math.log10(max(rec["RE_abs"], 1.0e-300))
    return rec


def waveform_peak_metric(
    row: pd.Series,
    cfg: se.SEConfig,
    n_omega: int,
    n_k: int,
) -> float:
    _, _, U = se.synthesize_waveforms_spectral(row, cfg, n_omega=n_omega, n_k=n_k)
    val = float(np.nanmax(np.abs(U)))
    return val if np.isfinite(val) else np.nan


def build_log_contribution_table(metrics: Mapping[str, float], parameter_groups: Sequence[str]) -> pd.DataFrame:
    baseline = max(float(metrics["baseline"]), 1.0e-300)
    full_target = max(float(metrics["full_target"]), 1.0e-300)
    full_delta = math.log10(full_target / baseline)
    rows = []
    running = 0.0
    for group in parameter_groups:
        value = max(float(metrics[group]), 1.0e-300)
        delta = math.log10(value / baseline)
        running += delta
        rows.append({
            "component": group,
            "metric": value,
            "delta_log10_metric": delta,
            "sign": "positive" if delta > 0 else "negative" if delta < 0 else "zero",
        })
    residual = full_delta - running
    rows.append({
        "component": "nonlinear_interaction_residual",
        "metric": np.nan,
        "delta_log10_metric": residual,
        "sign": "positive" if residual > 0 else "negative" if residual < 0 else "zero",
    })
    rows.append({
        "component": "full_observed_change",
        "metric": full_target,
        "delta_log10_metric": full_delta,
        "sign": "positive" if full_delta > 0 else "negative" if full_delta < 0 else "zero",
    })
    return pd.DataFrame(rows)


def plot_research_parameter_trajectories(df: pd.DataFrame, outpath: Path) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(10.0, 8.0), constrained_layout=True)
    x = df["Time_s"]
    plots = [
        ("Porosity", "porosity"),
        ("Permeability_mD", "permeability (mD)"),
        ("Tortuosity", "tortuosity"),
        ("C_molL", "electrolyte C (mol/L)"),
        ("Lambda_m", "Lambda from Schakel (m)"),
        ("PoreVolumeToSurface_m", "RTM pore volume/surface (m)"),
    ]
    for ax, (col, label) in zip(axes.ravel(), plots):
        if col not in df.columns:
            ax.axis("off")
            continue
        ax.plot(x, df[col], marker="o", markersize=2.5, linewidth=1.2)
        if col in {"Permeability_mD", "C_molL", "Lambda_m", "PoreVolumeToSurface_m"}:
            ax.set_yscale("log")
        ax.set_xlabel("dissolution time (s)")
        ax.set_ylabel(label)
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle("Reactive-transport parameter trajectories")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_oat_coefficient_sensitivity(df: pd.DataFrame, outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.5), constrained_layout=True)
    for component, part in df.groupby("component"):
        ax.plot(part["Time_s"], part["RE_abs_norm"], linewidth=1.5, label=component)
    ax.set_xlabel("dissolution time (s)")
    ax.set_ylabel("R_E amplitude normalized to baseline")
    ax.set_yscale("log")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.suptitle("One-at-a-time research-data coefficient sensitivity")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_contribution_bar(table: pd.DataFrame, outpath: Path) -> None:
    part = table[~table["component"].eq("full_observed_change")].copy()
    colors = ["#2ca02c" if v > 0 else "#d62728" if v < 0 else "#7f7f7f" for v in part["delta_log10_metric"]]
    fig, ax = plt.subplots(figsize=(8.0, 4.2), constrained_layout=True)
    ax.bar(part["component"], part["delta_log10_metric"], color=colors)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_ylabel("contribution to log10 waveform peak")
    ax.set_xlabel("parameter component")
    ax.tick_params(axis="x", rotation=25)
    ax.grid(True, axis="y", alpha=0.25)
    fig.suptitle("Signed contribution to interface EM waveform peak")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_research_data_sensitivity(
    input_path: Path,
    outdir: Path,
    cfg: se.SEConfig,
    contribution_n_omega: int,
    contribution_n_k: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = se.load_reactive_transport_table(input_path)
    ts = se.compute_time_series(df, cfg)
    diagnostics = df.copy()
    diagnostics["PoreVolumeToSurface_m"] = diagnostics.apply(pore_volume_to_surface_ratio_m, axis=1)
    diagnostics = diagnostics.merge(
        ts[["Time_s", "C_molL", "pH", "Lambda_m", "L_abs", "RE_abs", "TTM_abs", "valid_poroelastic"]],
        on="Time_s",
        how="left",
    )
    diagnostics.to_csv(outdir / "research_parameter_diagnostics.csv", index=False)
    plot_research_parameter_trajectories(diagnostics, outdir / "research_parameter_trajectories.png")

    valid = df[(df["Porosity"] > cfg.phi_min) & (df["Porosity"] < cfg.phi_max_valid)]
    if valid.empty:
        raise ValueError("No valid poroelastic rows are available for research-data sensitivity.")
    base = valid.iloc[0]
    target = valid.iloc[-1]

    oat_rows = []
    baseline_metric = row_coefficient_metrics(base, cfg)["RE_abs"]
    for _, target_row in valid.iterrows():
        full_rec = row_coefficient_metrics(target_row, cfg)
        full_rec.update({"Time_s": target_row["Time_s"], "component": "full_path"})
        oat_rows.append(full_rec)
        for group in PARAMETER_GROUPS:
            hybrid = make_one_at_a_time_row(base, target_row, group)
            rec = row_coefficient_metrics(hybrid, cfg)
            rec.update({"Time_s": target_row["Time_s"], "component": group})
            oat_rows.append(rec)
    oat = pd.DataFrame(oat_rows)
    oat["RE_abs_norm"] = oat["RE_abs"] / baseline_metric if baseline_metric > 0 else np.nan
    oat.to_csv(outdir / "research_one_at_a_time_coefficient_sensitivity.csv", index=False)
    plot_oat_coefficient_sensitivity(oat, outdir / "research_one_at_a_time_RE_sensitivity.png")

    waveform_metrics = {
        "baseline": waveform_peak_metric(base, cfg, contribution_n_omega, contribution_n_k),
        "full_target": waveform_peak_metric(target, cfg, contribution_n_omega, contribution_n_k),
    }
    for group in PARAMETER_GROUPS:
        hybrid = make_one_at_a_time_row(base, target, group)
        waveform_metrics[group] = waveform_peak_metric(hybrid, cfg, contribution_n_omega, contribution_n_k)
    contribution = build_log_contribution_table(waveform_metrics, list(PARAMETER_GROUPS))
    contribution["baseline_Time_s"] = float(base["Time_s"])
    contribution["target_Time_s"] = float(target["Time_s"])
    contribution["contribution_n_omega"] = int(contribution_n_omega)
    contribution["contribution_n_k"] = int(contribution_n_k)
    contribution.to_csv(outdir / "waveform_peak_contribution_table.csv", index=False)
    pd.DataFrame([waveform_metrics]).to_csv(outdir / "waveform_peak_metrics.csv", index=False)
    plot_contribution_bar(contribution, outdir / "waveform_peak_contribution_bar.png")
    return diagnostics, oat, contribution


def write_readme(outdir: Path) -> None:
    text = """# Parameter Sensitivity Outputs

This directory separates two tasks.

- `schakel_reproduction/`: Schakel-style parameter sweeps using Table I reference settings. Figures use the current boundary-value solver's `R_E` and `T_TM` potentials plus squared-magnitude energy proxies, not the full orthodox/interference vertical-flux decomposition.
- `research_data/`: Sensitivity analysis for `global_evolution.xlsx`. One-at-a-time curves isolate model-input groups along the dissolution path. The waveform contribution bar uses signed changes in `log10(Amax)` between the first valid and last valid poroelastic snapshots; the residual is the nonlinear interaction not assigned to one parameter.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="global_evolution.xlsx")
    parser.add_argument("--outdir", default="sensitivity_results")
    parser.add_argument("--skip-schakel", action="store_true")
    parser.add_argument("--skip-research", action="store_true")
    parser.add_argument("--contribution-n-omega", type=int, default=12)
    parser.add_argument("--contribution-n-k", type=int, default=61)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    schakel_dir = outdir / "schakel_reproduction"
    research_dir = outdir / "research_data"
    schakel_dir.mkdir(parents=True, exist_ok=True)
    research_dir.mkdir(parents=True, exist_ok=True)
    write_readme(outdir)

    if not args.skip_schakel:
        run_schakel_reproduction(schakel_dir)
    if not args.skip_research:
        run_research_data_sensitivity(
            Path(args.input),
            research_dir,
            se.SEConfig(),
            contribution_n_omega=args.contribution_n_omega,
            contribution_n_k=args.contribution_n_k,
        )


if __name__ == "__main__":
    main()
