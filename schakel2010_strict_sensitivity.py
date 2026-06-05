#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Strict Schakel & Smeulders (2010) Fig. 4-7 sensitivity reproduction.

This script computes the two vertical energy-flux coefficients plotted in
Schakel & Smeulders (2010):

    R_E^{E,E}      reflected EM orthodox flux / incident acoustic flux
    T_E^{TM,TM}   transmitted TM orthodox flux / incident acoustic flux

The boundary-value solution itself is reused from
`seismoelectric_offset_liu2018_spectral.py`, but the plotted quantities are
not squared potential proxies.  They are obtained from Schakel Eq. (47)-(49)
for the individual reflected EM and transmitted TM waves.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import replace
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import parameter_sensitivity_analysis as proxy
import seismoelectric_offset_liu2018_spectral as se


def incident_acoustic_flux(coeff: Mapping[str, complex], omega: float, cfg: se.SEConfig) -> float:
    """Incident fluid P-wave vertical flux for unit incident scalar potential."""
    return 0.5 * cfg.rho_fl * omega**2 * float(np.real(coeff["k3_fl"]))


def reflected_em_flux_coefficient(coeff: Mapping[str, complex], omega: float, cfg: se.SEConfig) -> float:
    """Schakel reflected EM orthodox flux coefficient R_E^{E,E}."""
    pin = incident_acoustic_flux(coeff, omega, cfg)
    r_e = coeff["R_E"]
    e1 = -coeff["k3_E"] * r_e
    h2 = -(coeff["s2_E"] / cfg.mu0) * r_e
    flux = 0.5 * np.real(e1 * np.conj(h2))
    return float(flux / pin)


def reflected_acoustic_flux_coefficient(coeff: Mapping[str, complex]) -> float:
    """Schakel reflected acoustic orthodox flux coefficient R_E^{Pr,Pr}."""
    return float(-abs(coeff["R_M"]) ** 2)


def transmitted_tm_flux_coefficient(
    coeff: Mapping[str, complex],
    omega: float,
    cfg: se.SEConfig,
) -> float:
    """Schakel transmitted TM orthodox flux coefficient T_E^{TM,TM}.

    The field amplitudes follow the scaling used in Appendix B.  Substitution
    into Eq. (49) gives the orthodox TM contribution plotted as -T_E^{TM,TM}
    in Schakel Fig. 4-7.
    """
    pin = incident_acoustic_flux(coeff, omega, cfg)
    t_tm = coeff["T_TM"]
    k1 = coeff["k1"]
    k3_tm = coeff["k3_TM"]
    s2_tm = coeff["s2_TM"]
    alpha_tm = coeff["alpha_TM"]
    beta_tm = coeff["beta_TM"]

    e1 = -k3_tm * alpha_tm * t_tm
    h2 = alpha_tm * s2_tm / cfg.mu0 * t_tm
    u1 = k3_tm * t_tm
    u3 = k1 * t_tm
    # U3 is included for completeness; the TM pressure term is zero.
    _u3_fluid = beta_tm * k1 * t_tm
    tau31 = -2.0 * cfg.G * (k1**2 - 0.5 * omega**2 * s2_tm) * t_tm
    tau33 = -2.0 * cfg.G * k1 * k3_tm * t_tm

    flux = 0.5 * np.real(
        e1 * np.conj(h2)
        + tau31 * np.conj(1j * omega * u1)
        + tau33 * np.conj(1j * omega * u3)
    )
    return float(flux / pin)


def strict_energy_flux_coefficients(
    phi: float,
    k0_m2: float,
    alpha_inf: float,
    cH: float,
    omega: float,
    theta_deg: float,
    cfg: se.SEConfig,
    C_override_molL: float | None = None,
) -> Dict[str, float | complex]:
    coeff = se.se_coefficients(
        phi,
        k0_m2,
        alpha_inf,
        cH,
        omega,
        theta_deg,
        cfg,
        C_override_molL=C_override_molL,
    )
    return {
        **coeff,
        "RE_EE": reflected_em_flux_coefficient(coeff, omega, cfg),
        "TE_TM_TM": transmitted_tm_flux_coefficient(coeff, omega, cfg),
        "minus_TE_TM_TM": -transmitted_tm_flux_coefficient(coeff, omega, cfg),
        "RE_Pr_Pr": reflected_acoustic_flux_coefficient(coeff),
        "Pin_flux": incident_acoustic_flux(coeff, omega, cfg),
    }


def reference_energy_flux_coefficients(omega: float = 1.0e6, theta_deg: float = 45.0) -> Dict[str, float | complex]:
    cfg, baseline = proxy.schakel_reference_config()
    return strict_energy_flux_coefficients(
        baseline["phi"],
        baseline["k0_m2"],
        baseline["alpha_inf"],
        baseline["cH"],
        omega,
        theta_deg,
        cfg,
        C_override_molL=baseline["C_molL"],
    )


def scan_definitions() -> Dict[str, np.ndarray]:
    return {
        "electrolyte_concentration_molL": np.logspace(-6, 0, 100),
        "viscosity_Pa_s": np.logspace(-6, 0, 100),
        "permeability_m2": np.logspace(-16, -10, 120),
        "tortuosity": np.linspace(1.0, 5.0, 100),
        "pH": np.linspace(2.0, 10.0, 100),
        "pore_fluid_bulk_modulus_Pa": np.logspace(7, 11, 100),
        "porosity": np.linspace(0.10, 0.85, 100),
    }


def run_single_scan(scan_name: str, values: Sequence[float]) -> pd.DataFrame:
    cfg, baseline = proxy.schakel_reference_config()
    rows = []
    for value in values:
        local_cfg = replace(cfg)
        phi = baseline["phi"]
        k0_m2 = baseline["k0_m2"]
        alpha_inf = baseline["alpha_inf"]
        cH = baseline["cH"]
        c_override = baseline["C_molL"]

        if scan_name == "electrolyte_concentration_molL":
            c_override = float(value)
        elif scan_name == "viscosity_Pa_s":
            local_cfg.eta = float(value)
        elif scan_name == "permeability_m2":
            k0_m2 = float(value)
        elif scan_name == "tortuosity":
            alpha_inf = float(value)
        elif scan_name == "pH":
            # Schakel Fig. 6 varies pH directly above neutral.  The main
            # reactive-transport model protects field data with a neutral-water
            # floor, but the reference reproduction must not clip this scan at 7.
            local_cfg.H_min_molL = 1.0e-12
            cH = proxy.h_concentration_for_ph(float(value), local_cfg)
        elif scan_name == "pore_fluid_bulk_modulus_Pa":
            local_cfg.K_f = float(value)
        elif scan_name == "porosity":
            phi = float(value)
        else:
            raise ValueError(f"Unknown scan: {scan_name}")

        coeff = strict_energy_flux_coefficients(
            phi,
            k0_m2,
            alpha_inf,
            cH,
            baseline["omega"],
            baseline["theta_deg"],
            local_cfg,
            C_override_molL=c_override,
        )
        rows.append({
            "scan": scan_name,
            "scan_value": float(value),
            "Porosity_used": phi,
            "k0_m2": k0_m2,
            "Permeability_mD": k0_m2 / 9.869233e-16,
            "Tortuosity": alpha_inf,
            "C_molL": float(coeff["C_molL"]),
            "pH": float(coeff["pH"]),
            "eta_Pa_s": local_cfg.eta,
            "K_f_Pa": local_cfg.K_f,
            "Lambda_m": float(coeff["Lambda"]),
            "debye_d_m": float(coeff["debye_d"]),
            "L_abs": abs(coeff["L"]),
            "R_E_abs": abs(coeff["R_E"]),
            "T_TM_abs": abs(coeff["T_TM"]),
            "RE_EE": coeff["RE_EE"],
            "TE_TM_TM": coeff["TE_TM_TM"],
            "minus_TE_TM_TM": coeff["minus_TE_TM_TM"],
            "RE_Pr_Pr": coeff["RE_Pr_Pr"],
            "matrix_cond": float(coeff["matrix_cond"]),
        })
    return pd.DataFrame(rows)


def plot_two_by_two(
    df: pd.DataFrame,
    scan_names: Sequence[str],
    outpath: Path,
    title: str,
    xlabels: Mapping[str, str],
    xscale: Mapping[str, str],
    yscale: str = "linear",
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.3), constrained_layout=True)
    metrics = [("RE_EE", r"$R_E^{E,E}$"), ("minus_TE_TM_TM", r"$-T_E^{TM,TM}$")]
    for col, scan_name in enumerate(scan_names):
        part = df[df["scan"] == scan_name]
        for row, (metric, ylabel) in enumerate(metrics):
            ax = axes[row, col]
            ax.plot(part["scan_value"], part[metric], color="black", linewidth=1.3)
            ax.set_xlabel(xlabels[scan_name])
            ax.set_ylabel(ylabel)
            if xscale.get(scan_name) == "log":
                ax.set_xscale("log")
            if yscale == "log":
                ax.set_yscale("log")
            ax.grid(False)
    fig.suptitle(title)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_porosity(df: pd.DataFrame, outpath: Path) -> None:
    part = df[df["scan"] == "porosity"]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), constrained_layout=True)
    axes[0].plot(part["scan_value"], part["RE_EE"], color="black", linewidth=1.3)
    axes[1].plot(part["scan_value"], part["minus_TE_TM_TM"], color="black", linewidth=1.3)
    axes[0].set_ylabel(r"$R_E^{E,E}$")
    axes[1].set_ylabel(r"$-T_E^{TM,TM}$")
    for ax in axes:
        ax.set_xlabel("porosity")
    fig.suptitle("Strict Schakel Fig. 7 reproduction")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_strict_sensitivity(outdir: Path) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    frames = [run_single_scan(name, values) for name, values in scan_definitions().items()]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(outdir / "schakel2010_strict_flux_coefficients.csv", index=False)

    summary_rows = []
    for scan, part in df.groupby("scan"):
        for metric in ["RE_EE", "minus_TE_TM_TM", "L_abs"]:
            idx = part[metric].idxmax()
            min_val = part[metric].min()
            max_val = part[metric].max()
            summary_rows.append({
                "scan": scan,
                "metric": metric,
                "min": min_val,
                "max": max_val,
                "max_over_min": max_val / min_val if min_val > 0 else np.nan,
                "scan_value_at_max": df.loc[idx, "scan_value"],
            })
    pd.DataFrame(summary_rows).to_csv(outdir / "schakel2010_strict_flux_summary.csv", index=False)

    table_rows = []
    for theta in [30.0, 45.0]:
        result = reference_energy_flux_coefficients(omega=1.0e6, theta_deg=theta)
        table_rows.append({
            "omega_rad_s": 1.0e6,
            "theta_deg": theta,
            "RE_EE": result["RE_EE"],
            "TE_TM_TM": result["TE_TM_TM"],
            "minus_TE_TM_TM": result["minus_TE_TM_TM"],
            "RE_Pr_Pr": result["RE_Pr_Pr"],
        })
    pd.DataFrame(table_rows).to_csv(outdir / "schakel2010_tableIII_validation.csv", index=False)

    xlabels = {
        "electrolyte_concentration_molL": "C (mol/L)",
        "viscosity_Pa_s": r"$\eta$ (Pa s)",
        "permeability_m2": r"$k_0$ (m$^2$)",
        "tortuosity": r"$\alpha_\infty$",
        "pH": "pH",
        "pore_fluid_bulk_modulus_Pa": r"$K_f$ (Pa)",
    }
    xscale = {
        "electrolyte_concentration_molL": "log",
        "viscosity_Pa_s": "log",
        "permeability_m2": "log",
        "pore_fluid_bulk_modulus_Pa": "log",
    }
    plot_two_by_two(
        df,
        ["electrolyte_concentration_molL", "viscosity_Pa_s"],
        outdir / "fig4_strict_concentration_viscosity.png",
        "Strict Schakel Fig. 4 reproduction",
        xlabels,
        xscale,
        yscale="linear",
    )
    plot_two_by_two(
        df,
        ["permeability_m2", "tortuosity"],
        outdir / "fig5_strict_permeability_tortuosity.png",
        "Strict Schakel Fig. 5 reproduction",
        xlabels,
        xscale,
        yscale="linear",
    )
    plot_two_by_two(
        df,
        ["pH", "pore_fluid_bulk_modulus_Pa"],
        outdir / "fig6_strict_ph_bulk_modulus.png",
        "Strict Schakel Fig. 6 reproduction",
        xlabels,
        xscale,
        yscale="linear",
    )
    plot_porosity(df, outdir / "fig7_strict_porosity.png")
    return df


def write_readme(outdir: Path) -> None:
    text = """# Strict Schakel & Smeulders (2010) Sensitivity Reproduction

This directory contains a strict reproduction of the two energy-flux
coefficients plotted in Schakel & Smeulders Fig. 4-7:

- `RE_EE`: reflected EM orthodox vertical flux coefficient.
- `TE_TM_TM`: transmitted TM orthodox vertical flux coefficient.
- `minus_TE_TM_TM`: the plotted positive value `-TE_TM_TM`.

The validation file compares the implementation with Schakel Table III at
`omega = 1e6 rad/s` and `theta = 30, 45 deg`.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="sensitivity_results/schakel2010_strict_reproduction")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    run_strict_sensitivity(outdir)
    write_readme(outdir)


if __name__ == "__main__":
    main()
