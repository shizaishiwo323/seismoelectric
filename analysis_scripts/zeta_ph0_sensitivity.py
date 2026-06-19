#!/usr/bin/env python3
"""Sensitivity test for the empirical zeta-potential zero-crossing pH.

This analysis is intentionally non-invasive: it reads existing Schakel 2011
Sommerfeld RTM result tables, recomputes only the dynamic electrokinetic
properties under alternative pH0 values in the semi-empirical zeta relation,
and writes all outputs to a separate results directory.
"""

from __future__ import annotations

import argparse
import math
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import seismoelectric_offset_liu2018_spectral as se


RUNS = {
    "Pe0.1": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe0p1"),
    "Pe1": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe1"),
    "Pe10": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe10"),
}

DEFAULT_PH0 = [1.5, 2.0, 2.5, 3.0]


def _electrochemistry_with_ph0(
    cH: float,
    cfg: se.SEConfig,
    *,
    ph0: float,
    C_override_molL: float | None = None,
) -> dict[str, float]:
    cH_molL = se.h_conc_to_molL(cH, cfg)
    pH = -math.log10(cH_molL)
    if C_override_molL is not None and np.isfinite(C_override_molL) and C_override_molL > 0:
        C = float(C_override_molL)
    else:
        C = cfg.C_background_molL + cH_molL
    zeta = (0.010 + 0.025 * math.log10(C)) * (pH - ph0) / 5.0
    return {"cH_molL": cH_molL, "pH": pH, "C_molL": C, "zeta": zeta}


@contextmanager
def _patched_zeta_zero(ph0: float) -> Iterator[None]:
    original = se.electrochemistry_from_h

    def replacement(
        cH: float,
        cfg: se.SEConfig,
        C_override_molL: float | None = None,
    ) -> dict[str, float]:
        return _electrochemistry_with_ph0(cH, cfg, ph0=ph0, C_override_molL=C_override_molL)

    se.electrochemistry_from_h = replacement
    try:
        yield
    finally:
        se.electrochemistry_from_h = original


def _read_run(label: str, path: Path) -> pd.DataFrame:
    table_path = path / "seismoelectric_timeseries_results.csv"
    df = pd.read_csv(table_path)
    df.insert(0, "Pe_label", label)
    df.insert(1, "source_table", str(table_path))
    return df


def _optional_float(row: pd.Series, name: str) -> float | None:
    if name not in row:
        return None
    value = row[name]
    if pd.isna(value):
        return None
    return float(value)


def _recompute_for_ph0(df: pd.DataFrame, ph0: float, cfg: se.SEConfig) -> pd.DataFrame:
    records = []
    with _patched_zeta_zero(ph0):
        for _, row in df.iterrows():
            coeff = se.dynamic_coefficients(
                float(row["Porosity_used"]),
                float(row["k0_m2"]),
                float(row["Tortuosity"]),
                float(row["OutletHConc_raw"]),
                float(row["omega0_rad_s"]),
                cfg,
                C_override_molL=_optional_float(row, "ElectrolyteConcentration_input_molL"),
                sigma_f_override=_optional_float(row, "FluidConductivity_input_S_m"),
            )
            records.append(
                {
                    "Pe_label": row["Pe_label"],
                    "pH0": ph0,
                    "Time_s": float(row["Time_s"]),
                    "Time_min": float(row["Time_min"]),
                    "Porosity_used": float(row["Porosity_used"]),
                    "Permeability_mD": float(row["Permeability_mD"]),
                    "k0_m2": float(row["k0_m2"]),
                    "Tortuosity": float(row["Tortuosity"]),
                    "OutletHConc_raw": float(row["OutletHConc_raw"]),
                    "cH_molL": float(coeff["cH_molL"]),
                    "pH": float(coeff["pH"]),
                    "C_molL": float(coeff["C_molL"]),
                    "zeta": float(coeff["zeta"]),
                    "omega_t": float(coeff["omega_t"]),
                    "Lambda_m": float(coeff["Lambda"]),
                    "debye_d": float(coeff["debye_d"]),
                    "L_abs": abs(coeff["L"]),
                    "sigma_abs": abs(coeff["sigma"]),
                    "L_over_sigma": abs(coeff["L"]) / abs(coeff["sigma"]) if abs(coeff["sigma"]) > 0 else np.nan,
                    "Amax_waveform_schakel2011_RE": float(row["Amax_waveform_schakel2011_RE"]),
                    "Amax_waveform_schakel2011_TE": float(row["Amax_waveform_schakel2011_TE"]),
                }
            )
    out = pd.DataFrame(records)
    first = out.sort_values("Time_s").groupby("Pe_label", sort=False).first()
    out["L_abs_norm"] = out.apply(lambda r: r["L_abs"] / first.loc[r["Pe_label"], "L_abs"], axis=1)
    out["sigma_abs_norm"] = out.apply(
        lambda r: r["sigma_abs"] / first.loc[r["Pe_label"], "sigma_abs"], axis=1
    )
    out["L_over_sigma_norm"] = out.apply(
        lambda r: r["L_over_sigma"] / first.loc[r["Pe_label"], "L_over_sigma"], axis=1
    )
    return out


def _crossing_time(group: pd.DataFrame, ph0: float) -> float:
    ordered = group.sort_values("Time_s")
    x = ordered["pH"].to_numpy(dtype=float) - ph0
    t = ordered["Time_min"].to_numpy(dtype=float)
    for i in range(len(x) - 1):
        if x[i] == 0:
            return float(t[i])
        if x[i] * x[i + 1] < 0:
            return float(t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]))
    return float("nan")


def _build_summary(sensitivity: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (pe_label, ph0), group in sensitivity.groupby(["Pe_label", "pH0"], sort=False):
        ordered = group.sort_values("Time_s")
        min_row = ordered.loc[ordered["L_abs"].idxmin()]
        nearest_ph0_row = ordered.loc[(ordered["pH"] - ph0).abs().idxmin()]
        base = ordered.iloc[0]
        rows.append(
            {
                "Pe_label": pe_label,
                "pH0": ph0,
                "pH_crossing_time_min": _crossing_time(ordered, ph0),
                "nearest_sample_time_min_to_pH0": float(nearest_ph0_row["Time_min"]),
                "nearest_sample_abs_pH_error_to_pH0": abs(float(nearest_ph0_row["pH"]) - ph0),
                "L_min_time_min": float(min_row["Time_min"]),
                "L_min_pH": float(min_row["pH"]),
                "L_min_abs_pH_error_to_pH0": abs(float(min_row["pH"]) - ph0),
                "L_min_zeta": float(min_row["zeta"]),
                "L_min_abs": float(min_row["L_abs"]),
                "L_min_norm": float(min_row["L_abs"] / base["L_abs"]),
                "sigma_at_L_min_abs": float(min_row["sigma_abs"]),
                "sigma_at_L_min_norm": float(min_row["sigma_abs"] / base["sigma_abs"]),
                "Porosity_at_L_min": float(min_row["Porosity_used"]),
                "Permeability_factor_at_L_min": float(min_row["Permeability_mD"] / base["Permeability_mD"]),
            }
        )
    return pd.DataFrame(rows)


def _plot_pe_timeseries(sensitivity: pd.DataFrame, outdir: Path) -> None:
    for pe_label, group in sensitivity.groupby("Pe_label", sort=False):
        safe_label = pe_label.replace(".", "p")
        fig, axes = plt.subplots(3, 1, figsize=(8.0, 8.8), sharex=True)
        for ph0, ph_group in group.groupby("pH0", sort=True):
            ordered = ph_group.sort_values("Time_min")
            axes[0].plot(ordered["Time_min"], ordered["zeta"], marker="o", ms=2.6, lw=1.2, label=f"pH0={ph0:g}")
            axes[1].plot(
                ordered["Time_min"],
                ordered["L_abs_norm"].where(ordered["L_abs_norm"] > 0),
                marker="o",
                ms=2.6,
                lw=1.2,
            )
            axes[2].plot(
                ordered["Time_min"],
                ordered["L_over_sigma_norm"].where(ordered["L_over_sigma_norm"] > 0),
                marker="o",
                ms=2.6,
                lw=1.2,
            )
        axes[0].axhline(0, color="0.25", lw=0.9, ls="--")
        axes[0].set_ylabel("zeta (V)")
        axes[1].set_ylabel("normalized |L(omega)|")
        axes[2].set_ylabel("normalized |L|/|sigma|")
        axes[2].set_xlabel("dissolution time (min)")
        axes[1].set_yscale("log")
        axes[2].set_yscale("log")
        axes[0].legend(frameon=False, ncol=4)
        for ax in axes:
            ax.grid(True, which="major", alpha=0.25)
            ax.grid(True, which="minor", alpha=0.12)
        fig.suptitle(f"{pe_label}: zeta zero-crossing sensitivity", fontsize=13, fontweight="bold")
        fig.tight_layout()
        fig.savefig(outdir / f"{safe_label}_zeta_ph0_sensitivity_timeseries.png", dpi=260)
        plt.close(fig)


def _plot_summary(summary: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0))
    for pe_label, group in summary.groupby("Pe_label", sort=False):
        ordered = group.sort_values("pH0")
        axes[0].plot(ordered["pH0"], ordered["pH_crossing_time_min"], marker="o", lw=1.5, label=pe_label)
        axes[1].plot(ordered["pH0"], ordered["L_min_norm"], marker="o", lw=1.5, label=pe_label)
    axes[0].set_xlabel("assumed zeta zero-crossing pH0")
    axes[0].set_ylabel("pH crossing time (min)")
    axes[1].set_xlabel("assumed zeta zero-crossing pH0")
    axes[1].set_ylabel("minimum normalized |L(omega)|")
    axes[1].set_yscale("log")
    for ax in axes:
        ax.grid(True, which="major", alpha=0.25)
        ax.grid(True, which="minor", alpha=0.12)
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / "zeta_ph0_sensitivity_summary.png", dpi=260)
    plt.close(fig)


def _plot_ph0_waveform_style(sensitivity: pd.DataFrame, outdir: Path) -> None:
    """Make one log-time, two-panel figure per pH0 with pH on twin y-axes."""
    colors = {"Pe0.1": "tab:blue", "Pe1": "tab:orange", "Pe10": "tab:green"}
    for ph0, group in sensitivity.groupby("pH0", sort=True):
        fig, axes = plt.subplots(2, 1, figsize=(11.0, 7.0), sharex=True)
        twin_axes = [ax.twinx() for ax in axes]

        for pe_label, pe_group in group.groupby("Pe_label", sort=False):
            ordered = pe_group.sort_values("Time_s")
            color = colors.get(pe_label, None)
            t = ordered["Time_s"].astype(float)
            valid_t = t > 0

            axes[0].plot(
                t[valid_t],
                ordered.loc[valid_t, "L_abs_norm"].where(ordered.loc[valid_t, "L_abs_norm"] > 0),
                marker="o",
                ms=3.0,
                lw=1.6,
                color=color,
                label=pe_label,
            )
            axes[1].plot(
                t[valid_t],
                ordered.loc[valid_t, "L_over_sigma_norm"].where(
                    ordered.loc[valid_t, "L_over_sigma_norm"] > 0
                ),
                marker="o",
                ms=3.0,
                lw=1.6,
                color=color,
                label=pe_label,
            )

            for twin in twin_axes:
                twin.plot(
                    t[valid_t],
                    ordered.loc[valid_t, "pH"],
                    ls=":",
                    lw=1.2,
                    color=color,
                    alpha=0.45,
                )

        for ax in axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="major", alpha=0.25)
            ax.grid(True, which="minor", alpha=0.10)
            ax.legend(frameon=False, ncol=3, loc="lower left")

        for twin in twin_axes:
            twin.axhline(ph0, color="0.25", lw=1.1, ls="--", alpha=0.75)
            twin.set_ylabel("pH", color="0.25")
            twin.tick_params(axis="y", colors="0.25")
            twin.grid(False)

        axes[0].set_title(f"Dynamic coupling with pH trajectory, zeta zero pH0 = {ph0:g}")
        axes[0].set_ylabel("normalized |L(omega)|")
        axes[1].set_title("Electrokinetic efficiency index")
        axes[1].set_ylabel("normalized |L|/|sigma|")
        axes[1].set_xlabel("dissolution time (s)")
        fig.tight_layout()

        ph0_label = f"{ph0:g}".replace(".", "p")
        fig.savefig(outdir / f"ph0_{ph0_label}_waveform_style_dynamic_with_pH_logx.png", dpi=260)
        plt.close(fig)


def _plot_ph0_waveform_peak_with_ph(sensitivity: pd.DataFrame, outdir: Path) -> None:
    """Match the existing waveform-peak log-time style, adding pH as twin axes."""
    colors = {"Pe0.1": "tab:blue", "Pe1": "tab:orange", "Pe10": "tab:green"}
    for ph0, group in sensitivity.groupby("pH0", sort=True):
        fig, axes = plt.subplots(2, 1, figsize=(11.0, 7.0), sharex=True)
        twin_axes = [ax.twinx() for ax in axes]

        for pe_label, pe_group in group.groupby("Pe_label", sort=False):
            ordered = pe_group.sort_values("Time_s")
            color = colors.get(pe_label, None)
            t = ordered["Time_s"].astype(float)
            valid_t = t > 0
            re_mv = ordered.loc[valid_t, "Amax_waveform_schakel2011_RE"].abs() * 1.0e3
            te_mv = ordered.loc[valid_t, "Amax_waveform_schakel2011_TE"].abs() * 1.0e3

            axes[0].plot(t[valid_t], re_mv.where(re_mv > 0), marker="o", ms=3.0, lw=1.6, color=color, label=pe_label)
            axes[1].plot(t[valid_t], te_mv.where(te_mv > 0), marker="o", ms=3.0, lw=1.6, color=color, label=pe_label)

            for twin in twin_axes:
                twin.plot(
                    t[valid_t],
                    ordered.loc[valid_t, "pH"],
                    ls=":",
                    lw=1.2,
                    color=color,
                    alpha=0.45,
                )

        for ax in axes:
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.grid(True, which="major", alpha=0.25)
            ax.grid(True, which="minor", alpha=0.10)
            ax.legend(frameon=False, ncol=3, loc="lower left")

        for twin in twin_axes:
            twin.axhline(ph0, color="0.25", lw=1.1, ls="--", alpha=0.75)
            twin.set_ylabel("pH", color="0.25")
            twin.tick_params(axis="y", colors="0.25")
            twin.grid(False)

        axes[0].set_title(f"Fluid-side R_E peak with pH trajectory, zeta zero pH0 = {ph0:g}")
        axes[0].set_ylabel("electric potential (mV)")
        axes[1].set_title("Porous-side T_E peak")
        axes[1].set_ylabel("electric potential (mV)")
        axes[1].set_xlabel("dissolution time (s)")
        fig.tight_layout()

        ph0_label = f"{ph0:g}".replace(".", "p")
        fig.savefig(outdir / f"ph0_{ph0_label}_waveform_peak_mV_with_pH_logx.png", dpi=260)
        plt.close(fig)


def run_analysis(run_root: Path, outdir: Path, ph0_values: list[float]) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = se.SEConfig()
    frames = []
    for label, rel_path in RUNS.items():
        frames.append(_read_run(label, run_root / rel_path))
    source = pd.concat(frames, ignore_index=True)

    sensitivity = pd.concat(
        [_recompute_for_ph0(source, ph0, cfg) for ph0 in ph0_values],
        ignore_index=True,
    )
    summary = _build_summary(sensitivity)

    outdir.mkdir(parents=True, exist_ok=True)
    sensitivity.to_csv(outdir / "zeta_ph0_sensitivity_timeseries.csv", index=False)
    summary.to_csv(outdir / "zeta_ph0_sensitivity_summary.csv", index=False)
    _plot_pe_timeseries(sensitivity, outdir)
    _plot_summary(summary, outdir)
    _plot_ph0_waveform_style(sensitivity, outdir)
    _plot_ph0_waveform_peak_with_ph(sensitivity, outdir)
    return sensitivity, summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, default=ROOT, help="Repository root containing the results directory.")
    parser.add_argument(
        "--outdir",
        type=Path,
        default=ROOT / "results" / "zeta_ph0_sensitivity_20260617",
        help="Separate output directory for sensitivity-test products.",
    )
    parser.add_argument(
        "--ph0",
        type=float,
        nargs="+",
        default=DEFAULT_PH0,
        help="Alternative zero-crossing pH values in zeta = factor * (pH - pH0) / 5.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _, summary = run_analysis(args.run_root, args.outdir, args.ph0)
    safe_outdir = str(args.outdir).encode("ascii", "backslashreplace").decode("ascii")
    print(f"Wrote sensitivity outputs to: {safe_outdir}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
