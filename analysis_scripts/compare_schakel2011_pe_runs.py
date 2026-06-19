#!/usr/bin/env python3
"""Compare Schakel 2011 Sommerfeld RT-SE outputs across Pe runs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RUNS = {
    "Pe0.1": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe0p1"),
    "Pe1": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe1"),
    "Pe10": Path("results/schakel2011_sommerfeld_rtm_20260617/Pe10"),
}


def _read_run(label: str, run_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(run_dir / "seismoelectric_timeseries_results.csv")
    df.insert(0, "Pe_label", label)
    df.insert(1, "run_dir", str(run_dir))
    return df


def _load_all() -> pd.DataFrame:
    frames = [_read_run(label, run_dir) for label, run_dir in RUNS.items()]
    return pd.concat(frames, ignore_index=True)


def _nearest_summary(df: pd.DataFrame, phi_targets: list[float]) -> pd.DataFrame:
    rows = []
    for label, group in df.groupby("Pe_label", sort=False):
        ordered = group.sort_values("Time_s")
        for target in phi_targets:
            idx = (ordered["Porosity_used"] - target).abs().idxmin()
            row = ordered.loc[idx].copy()
            row["target_Porosity"] = target
            row["abs_phi_error"] = abs(float(row["Porosity_used"]) - target)
            rows.append(row)
    return pd.DataFrame(rows)


def _endpoint_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, group in df.groupby("Pe_label", sort=False):
        ordered = group.sort_values("Time_s")
        start = ordered.iloc[0]
        end = ordered.iloc[-1]
        peak_re = ordered.loc[ordered["Amax_waveform_schakel2011_RE"].abs().idxmax()]
        peak_te = ordered.loc[ordered["Amax_waveform_schakel2011_TE"].abs().idxmax()]
        rows.append(
            {
                "Pe_label": label,
                "n_steps": len(ordered),
                "Time_start_s": float(start["Time_s"]),
                "Time_end_s": float(end["Time_s"]),
                "Porosity_start": float(start["Porosity_used"]),
                "Porosity_end": float(end["Porosity_used"]),
                "Permeability_start_mD": float(start["Permeability_mD"]),
                "Permeability_end_mD": float(end["Permeability_mD"]),
                "Tortuosity_start": float(start["Tortuosity"]),
                "Tortuosity_end": float(end["Tortuosity"]),
                "OutletHConc_start": float(start["OutletHConc_raw"]),
                "OutletHConc_end": float(end["OutletHConc_raw"]),
                "RE_peak_abs_max_mV": abs(float(peak_re["Amax_waveform_schakel2011_RE"])) * 1.0e3,
                "RE_peak_abs_max_Time_s": float(peak_re["Time_s"]),
                "RE_peak_abs_max_Porosity": float(peak_re["Porosity_used"]),
                "TE_peak_abs_max_mV": abs(float(peak_te["Amax_waveform_schakel2011_TE"])) * 1.0e3,
                "TE_peak_abs_max_Time_s": float(peak_te["Time_s"]),
                "TE_peak_abs_max_Porosity": float(peak_te["Porosity_used"]),
            }
        )
    return pd.DataFrame(rows)


def _plot_lines(
    df: pd.DataFrame,
    outdir: Path,
    x_col: str,
    y_cols: list[tuple[str, str, float]],
    filename: str,
    xlabel: str,
    ylabel: str,
    logy: bool = False,
    logx: bool = False,
) -> None:
    fig, axes = plt.subplots(len(y_cols), 1, figsize=(7.2, 2.7 * len(y_cols)), sharex=True)
    if len(y_cols) == 1:
        axes = [axes]
    for ax, (col, title, scale) in zip(axes, y_cols):
        for label, group in df.groupby("Pe_label", sort=False):
            ordered = group.sort_values(x_col)
            y = ordered[col].astype(float).to_numpy() * scale
            if logy:
                y = np.abs(y)
                y[y <= 0] = np.nan
            ax.plot(ordered[x_col], y, marker="o", ms=2.5, lw=1.2, label=label)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
        if logx:
            ax.set_xscale("log")
        if logy:
            ax.set_yscale("log")
        ax.legend(frameon=False, ncol=3)
    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(outdir / filename, dpi=220)
    plt.close(fig)


def main() -> None:
    outdir = Path("results/schakel2011_sommerfeld_rtm_20260617/comparison")
    outdir.mkdir(parents=True, exist_ok=True)

    df = _load_all()
    df.to_csv(outdir / "combined_timeseries_results.csv", index=False)
    _endpoint_summary(df).to_csv(outdir / "endpoint_and_peak_summary.csv", index=False)
    _nearest_summary(df, [0.60, 0.70, 0.80, 0.90, 0.99]).to_csv(
        outdir / "nearest_porosity_milestone_summary.csv", index=False
    )

    _plot_lines(
        df,
        outdir,
        "Time_s",
        [
            ("Porosity_used", "Porosity", 1.0),
            ("Permeability_mD", "Permeability", 1.0),
            ("Tortuosity", "Tortuosity", 1.0),
            ("OutletHConc_raw", "Outlet H+ concentration", 1.0),
        ],
        "rt_parameters_vs_time.png",
        "dissolution time (s)",
        "value",
        logy=False,
    )
    _plot_lines(
        df,
        outdir,
        "Porosity_used",
        [
            ("sigma_abs", "Dynamic conductivity |sigma|", 1.0),
            ("L_abs", "Dynamic coupling |L|", 1.0),
            ("zeta", "Zeta potential", 1.0),
        ],
        "dynamic_coefficients_vs_porosity.png",
        "porosity",
        "value",
        logy=True,
    )
    _plot_lines(
        df,
        outdir,
        "Porosity_used",
        [
            ("RE_abs", "Reflected interface coefficient |R_E|", 1.0),
            ("TTM_abs", "Transmitted TM coefficient |T_TM|", 1.0),
        ],
        "interface_coefficients_vs_porosity.png",
        "porosity",
        "coefficient",
        logy=True,
    )
    _plot_lines(
        df,
        outdir,
        "Porosity_used",
        [
            ("Amax_waveform_schakel2011_RE", "Fluid-side R_E peak", 1.0e3),
            ("Amax_waveform_schakel2011_TE", "Porous-side T_E peak", 1.0e3),
        ],
        "waveform_peak_mV_vs_porosity.png",
        "porosity",
        "electric potential (mV)",
        logy=True,
    )
    _plot_lines(
        df,
        outdir,
        "Time_s",
        [
            ("Amax_waveform_schakel2011_RE", "Fluid-side R_E peak", 1.0e3),
            ("Amax_waveform_schakel2011_TE", "Porous-side T_E peak", 1.0e3),
        ],
        "waveform_peak_mV_vs_time.png",
        "dissolution time (s)",
        "electric potential (mV)",
        logy=True,
    )
    _plot_lines(
        df,
        outdir,
        "Time_s",
        [
            ("Amax_waveform_schakel2011_RE", "Fluid-side R_E peak", 1.0e3),
            ("Amax_waveform_schakel2011_TE", "Porous-side T_E peak", 1.0e3),
        ],
        "waveform_peak_mV_vs_time_logx.png",
        "dissolution time (s)",
        "electric potential (mV)",
        logy=True,
        logx=True,
    )


if __name__ == "__main__":
    main()
