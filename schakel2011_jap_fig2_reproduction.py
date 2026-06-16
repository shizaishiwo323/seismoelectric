#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce Schakel et al. (2011, JAP) Figure 2.

The implementation follows the JAP theory-paper geometry and Eq. (5)/(8)
Sommerfeld integrals. Fluid-side positions 1-5 use the reflected EM field of
Eq. (5). Position 6 uses the porous-medium expression of Eq. (8), including
the first transmitted TM field, the Pf coseismic field, and the Eq. (7)
back/front interface-response terms.
"""

from __future__ import annotations

import argparse
import cmath
import math
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

import schakel2010_strict_sensitivity as strict
from schakel2011_geophysics_reproduction import _j1, _source_fft_table

try:
    from scipy.integrate import quad as scipy_quad
    from scipy.special import j0 as scipy_j0
except Exception:  # pragma: no cover
    scipy_quad = None
    scipy_j0 = None


@dataclass(frozen=True)
class JAPParameters:
    K_s: float = 49.9e9
    K_f: float = 2.2e9
    eta: float = 1.0e-3
    rho_f: float = 998.0
    eps_f: float = 80.1
    eps_s: float = 4.0
    mu0: float = 4.0 * math.pi * 1.0e-7
    Lambda: float = 9.43e-6
    K_b: float = 0.93e9
    G: float = 0.88e9
    phi: float = 0.52
    rho_s: float = 2570.0
    k0_m2: float = 3.4e-12
    alpha_inf: float = 1.7
    sigma_f_s_m: float = 4.8e-2
    zeta_v: float = -4.0e-2
    sample_width_m: float = 3.21e-2
    source_z_m: float = -0.15
    transducer_radius_m: float = 1.125 * 0.0254 / 2.0
    temperature_k: float = 293.15
    bandpass_low_hz: float = 144_000.0
    bandpass_high_hz: float = 896_000.0

    @property
    def cP_fluid_m_s(self) -> float:
        return math.sqrt(self.K_f / self.rho_f)

    @property
    def source_z_to_interface_m(self) -> float:
        return abs(self.source_z_m)


JAP_POSITIONS = pd.DataFrame(
    [
        {"position": 1, "r_m": 0.0, "z_m": -0.023, "medium": "fluid"},
        {"position": 2, "r_m": 0.0, "z_m": -0.018, "medium": "fluid"},
        {"position": 3, "r_m": 0.0, "z_m": -0.013, "medium": "fluid"},
        {"position": 4, "r_m": 0.0, "z_m": -0.008, "medium": "fluid"},
        {"position": 5, "r_m": 0.0, "z_m": -0.003, "medium": "fluid"},
        {"position": 6, "r_m": 0.0, "z_m": 0.010, "medium": "porous"},
    ]
)


FIG2_PAGE3_AXES = {
    "predicted": {"x0": 1489, "x1": 1829, "ylim_mv": 0.5},
    "observed": {"x0": 1936, "x1": 2275, "ylim_mv": 0.2},
}

FIG2_PAGE3_ROWS = [
    {"position": 1, "y_top": 1772, "y_zero": 1820, "y_bottom": 1868},
    {"position": 2, "y_top": 1990, "y_zero": 2038, "y_bottom": 2086},
    {"position": 3, "y_top": 2208, "y_zero": 2256, "y_bottom": 2305},
    {"position": 4, "y_top": 2430, "y_zero": 2478, "y_bottom": 2527},
    {"position": 5, "y_top": 2648, "y_zero": 2697, "y_bottom": 2745},
    {"position": 6, "y_top": 2867, "y_zero": 2915, "y_bottom": 2963},
]

FIG2_TIME_MIN_MS = 0.055
FIG2_TIME_MAX_MS = 0.145


def _project_root() -> Path:
    return Path(__file__).resolve().parent


def _jap_pdf_path() -> Path:
    return (
        _project_root()
        / "uploaded_files"
        / "01_classical_theory_foundational_papers"
        / "02_Schakel_2011_JAP_laboratory_theory.pdf"
    )


def _render_jap_page3(outdir: Path | None = None) -> Path:
    """Render the JAP paper page containing Fig. 2 at 300 dpi."""
    if outdir is None:
        outdir = Path(tempfile.mkdtemp(prefix="schakel2011_jap_fig2_"))
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / "jap_page"
    target = outdir / "jap_page-3.png"
    if target.exists():
        return target
    pdf_relative = Path("uploaded_files") / "01_classical_theory_foundational_papers" / _jap_pdf_path().name
    completed = subprocess.run(
        [
            "pdftoppm",
            "-f",
            "3",
            "-l",
            "3",
            "-r",
            "300",
            "-png",
            str(pdf_relative),
            str(prefix),
        ],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=str(_project_root()),
    )
    if not target.exists() and completed.returncode != 0:
        message = completed.stderr.decode(errors="replace")
        raise RuntimeError(f"pdftoppm did not render JAP page 3: {message}")
    return target


def _digitize_axis_trace(
    page_gray: np.ndarray,
    x0: int,
    x1: int,
    y_top: int,
    y_zero: int,
    y_bottom: int,
    ylim_mv: float,
    n_time: int,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Digitize one compact Fig. 2 trace from the 300 dpi page raster.

    Frame ticks are removed before sampling pixels that deviate from the zero
    baseline. Missing columns are filled by the zero baseline rather than by a
    fitted waveform, preserving the fact that weak position-1/2 pulses are
    barely resolved in the printed raster.
    """
    crop = page_gray[y_top + 2 : y_bottom - 2, x0 + 2 : x1 - 2]
    dark = crop < 150
    margin = 12
    dark[:margin, :] = False
    dark[-margin:, :] = False
    col_count = dark.sum(axis=0)
    remove_cols = col_count > crop.shape[0] * 0.45
    clean = dark & (~remove_cols[None, :])
    local_y_zero = y_zero - (y_top + 2)

    xs = []
    ys = []
    for local_x in range(clean.shape[1]):
        local_y = np.where(clean[:, local_x])[0]
        if len(local_y):
            deviating = np.abs(local_y - local_y_zero) > 1.0
            if np.any(deviating):
                xs.append(x0 + 2 + local_x)
                ys.append(y_top + 2 + float(np.median(local_y[deviating])))

    time_ms = np.linspace(FIG2_TIME_MIN_MS, FIG2_TIME_MAX_MS, int(n_time))
    if len(xs) < 2:
        return time_ms, np.zeros_like(time_ms), len(xs)

    xs_arr = np.asarray(xs, dtype=float)
    ys_arr = np.asarray(ys, dtype=float)
    source_time = FIG2_TIME_MIN_MS + (xs_arr - x0) / (x1 - x0) * (FIG2_TIME_MAX_MS - FIG2_TIME_MIN_MS)
    pixel_half_height = (y_bottom - y_top) / 2.0
    source_mv = (y_zero - ys_arr) / pixel_half_height * float(ylim_mv)
    digitized_mv = np.interp(time_ms, source_time, source_mv, left=0.0, right=0.0)
    return time_ms, digitized_mv, len(xs)


def digitized_published_fig2_traces(
    n_time: int = 361,
    page_png: str | Path | None = None,
) -> pd.DataFrame:
    """Return page-image digitizations of the published JAP Fig. 2 traces."""
    page = Path(page_png) if page_png is not None else _render_jap_page3()
    gray = np.asarray(Image.open(page).convert("L"))
    rows = []
    for column, axis in FIG2_PAGE3_AXES.items():
        for row in FIG2_PAGE3_ROWS:
            time_ms, digitized_mv, n_pixels = _digitize_axis_trace(
                gray,
                axis["x0"],
                axis["x1"],
                row["y_top"],
                row["y_zero"],
                row["y_bottom"],
                axis["ylim_mv"],
                n_time,
            )
            for t, v in zip(time_ms, digitized_mv):
                rows.append(
                    {
                        "column": column,
                        "position": int(row["position"]),
                        "time_ms": float(t),
                        "digitized_mv": float(v),
                        "source": "Schakel_2011_JAP_Fig2_300dpi_page_raster",
                        "axis_ylim_mv": float(axis["ylim_mv"]),
                        "dark_curve_pixels_used": int(n_pixels),
                    }
                )
    return pd.DataFrame(rows)


def jap_parameters() -> JAPParameters:
    return JAPParameters()


def _j0(x):
    if scipy_j0 is not None:
        return scipy_j0(x)
    x = np.asarray(x, dtype=complex)
    t = np.linspace(0.0, 2.0 * math.pi, 512, endpoint=False)
    return np.mean(np.exp(1j * x[..., None] * np.cos(t)), axis=-1)


def _jap_config(params: JAPParameters) -> strict.SchakelConfig:
    cfg = strict.SchakelConfig()
    cfg.K_b = params.K_b
    cfg.G = params.G
    cfg.K_s = params.K_s
    cfg.K_f = params.K_f
    cfg.K_fl = params.K_f
    cfg.eta = params.eta
    cfg.rho_f = params.rho_f
    cfg.rho_fl = params.rho_f
    cfg.rho_s = params.rho_s
    cfg.alpha_inf = params.alpha_inf
    cfg.phi = params.phi
    cfg.k0_m2 = params.k0_m2
    cfg.temperature = params.temperature_k
    cfg.eps_f = params.eps_f
    cfg.eps_s = params.eps_s
    cfg.eps_fl = params.eps_f
    cfg.sigma_fl = params.sigma_f_s_m
    cfg.mu0 = params.mu0
    # Table II gives Lambda explicitly. In the Schakel dynamic-permeability
    # notation Lambda=sqrt(8 alpha_inf k0/(phi M)), so choose M accordingly.
    cfg.M_similarity = 8.0 * params.alpha_inf * params.k0_m2 / (params.phi * params.Lambda**2)
    # Table II says pH=6 and zeta=-0.04 V. The Schakel 2010 zeta relation is
    # retained by choosing the equivalent NaCl concentration for that zeta.
    cfg.pH = 6.0
    prefactor = params.zeta_v / ((cfg.pH - 2.0) / 5.0)
    cfg.C_molL = 10.0 ** ((prefactor - 0.010) / 0.025)
    return cfg


def _cH(cfg: strict.SchakelConfig) -> float:
    return strict.h_concentration_for_ph(cfg.pH)


@lru_cache(maxsize=65536)
def _fluid_incident_coefficients_cached(frequency_hz: float, kx_key: float | None) -> dict[str, complex]:
    params = jap_parameters()
    cfg = _jap_config(params)
    omega = 2.0 * math.pi * frequency_hz
    kx = None if kx_key is None else float(kx_key)
    theta_deg = 0.0 if kx is None else 90.0
    coeff = strict.se_coefficients(
        params.phi,
        params.k0_m2,
        params.alpha_inf,
        _cH(cfg),
        omega,
        theta_deg,
        cfg,
        kx_override=kx,
        C_override_molL=cfg.C_molL,
    )
    return coeff


def _fluid_incident_coefficients(frequency_hz: float, kx: float | None) -> dict[str, complex]:
    key = None if kx is None else round(float(kx), 6)
    return _fluid_incident_coefficients_cached(round(float(frequency_hz), 3), key)


def _source_A_spectrum(frequency_hz: np.ndarray, params: JAPParameters) -> np.ndarray:
    """Eq. (1)-consistent source amplitude from the recorded on-axis pressure."""
    fft_freqs, fft_amp = _source_fft_table()
    freq = np.asarray(frequency_hz, dtype=float)
    real = np.interp(freq, fft_freqs, np.real(fft_amp), left=0.0, right=0.0)
    imag = np.interp(freq, fft_freqs, np.imag(fft_amp), left=0.0, right=0.0)
    # Eq. (1): p_hat = A(omega) / R * exp(-i k R) * D(theta).
    # The recorded source waveform is at (r,z)=(0,0), on the source axis:
    # R=|z_s| and D(0)=lim J1(x)/x=0.5. The local record is therefore
    # converted to the source amplitude A(omega) by multiplying by R/D(0).
    return (real + 1j * imag) * (params.source_z_to_interface_m / 0.5) * _bandpass_taper(freq, params)


def _bandpass_taper(frequency_hz: np.ndarray, params: JAPParameters, taper_hz: float = 40_000.0) -> np.ndarray:
    """Smooth implementation of the paper's 144-896 kHz numerical band-pass."""
    f = np.asarray(frequency_hz, dtype=float)
    w = np.ones_like(f)
    low = params.bandpass_low_hz
    high = params.bandpass_high_hz
    w[(f <= low) | (f >= high)] = 0.0
    lo = (f > low) & (f < low + taper_hz)
    hi = (f < high) & (f > high - taper_hz)
    w[lo] = 0.5 - 0.5 * np.cos(math.pi * (f[lo] - low) / taper_hz)
    w[hi] = 0.5 - 0.5 * np.cos(math.pi * (high - f[hi]) / taper_hz)
    return w


def _em_kz_fluid(omega: float, kx: complex, params: JAPParameters) -> complex:
    s2_e = params.mu0 * strict.SchakelConfig.eps0 * params.eps_f - 1j * params.mu0 * params.sigma_f_s_m / omega
    return complex(strict.complex_sqrt_branch(omega**2 * s2_e - kx**2))


def _boundary_matrix_for_fluid_incidence(
    state: dict[str, complex],
    omega: float,
    k1: float,
    cfg: strict.SchakelConfig,
) -> np.ndarray:
    """Rebuild the Schakel 2010 Appendix-B matrix used by strict.se_coefficients."""
    phi = cfg.phi
    c_fl = math.sqrt(cfg.K_fl / cfg.rho_fl)
    k3_fl = strict.complex_sqrt_branch((omega / c_fl) ** 2 - k1**2)
    s2_E = cfg.mu0 * cfg.eps0 * cfg.eps_fl - 1j * cfg.mu0 * cfg.sigma_fl / omega
    k3_E = strict.complex_sqrt_branch(omega**2 * s2_E - k1**2)
    k3_Pf = strict.complex_sqrt_branch(omega**2 * state["s2_Pf"] - k1**2)
    k3_Ps = strict.complex_sqrt_branch(omega**2 * state["s2_Ps"] - k1**2)
    k3_TM = strict.complex_sqrt_branch(omega**2 * state["s2_TM"] - k1**2)
    k3_SV = strict.complex_sqrt_branch(omega**2 * state["s2_SV"] - k1**2)

    Q, R, P = state["Q"], state["R"], state["P"]
    beta_Pf, beta_Ps = state["beta_Pf"], state["beta_Ps"]
    beta_TM, beta_SV = state["beta_TM"], state["beta_SV"]
    alpha_Pf, alpha_Ps = state["alpha_Pf"], state["alpha_Ps"]
    alpha_TM, alpha_SV = state["alpha_TM"], state["alpha_SV"]
    s2_Pf, s2_Ps, s2_TM, s2_SV = state["s2_Pf"], state["s2_Ps"], state["s2_TM"], state["s2_SV"]

    N1 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Pf
    N2 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Ps

    A = np.zeros((6, 6), dtype=complex)
    A[0, 1] = k3_fl
    A[0, 2] = k3_Pf * (1.0 - phi + phi * beta_Pf)
    A[0, 3] = k3_Ps * (1.0 - phi + phi * beta_Ps)
    A[0, 4] = k1 * (1.0 - phi + phi * beta_TM)
    A[0, 5] = k1 * (1.0 - phi + phi * beta_SV)

    A[1, 1] = -phi * cfg.rho_fl
    A[1, 2] = (Q + R * beta_Pf) * s2_Pf
    A[1, 3] = (Q + R * beta_Ps) * s2_Ps

    A[2, 2] = k1 * k3_Pf
    A[2, 3] = k1 * k3_Ps
    A[2, 4] = k1**2 - 0.5 * omega**2 * s2_TM
    A[2, 5] = k1**2 - 0.5 * omega**2 * s2_SV

    A[3, 2] = k1**2 - omega**2 * s2_Pf * N1 / (2.0 * cfg.G)
    A[3, 3] = k1**2 - omega**2 * s2_Ps * N2 / (2.0 * cfg.G)
    A[3, 4] = -k1 * k3_TM
    A[3, 5] = -k1 * k3_SV

    A[4, 0] = -s2_E / cfg.mu0
    A[4, 4] = alpha_TM * s2_TM / cfg.mu0
    A[4, 5] = alpha_SV * s2_SV / cfg.mu0

    A[5, 0] = -k3_E
    A[5, 2] = k1 * alpha_Pf
    A[5, 3] = k1 * alpha_Ps
    A[5, 4] = -k3_TM * alpha_TM
    A[5, 5] = -k3_SV * alpha_SV
    return A


def _pf_incident_column(state: dict[str, complex], omega: float, k1: float, cfg: strict.SchakelConfig, upward: bool) -> np.ndarray:
    """Boundary-equation contribution for a unit Pf wave incident from porous medium."""
    phi = cfg.phi
    k3 = strict.complex_sqrt_branch(omega**2 * state["s2_Pf"] - k1**2)
    sign = -1.0 if upward else 1.0
    beta = state["beta_Pf"]
    alpha = state["alpha_Pf"]
    s2 = state["s2_Pf"]
    Q, R, P = state["Q"], state["R"], state["P"]
    N = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta
    col = np.zeros(6, dtype=complex)
    col[0] = sign * k3 * (1.0 - phi + phi * beta)
    col[1] = (Q + R * beta) * s2
    col[2] = sign * k1 * k3
    col[3] = k1**2 - omega**2 * s2 * N / (2.0 * cfg.G)
    col[5] = k1 * alpha
    return col


@lru_cache(maxsize=65536)
def _pf_interface_coefficients_cached(frequency_hz: float, kx_key: float) -> dict[str, complex]:
    """Pf incidence on a porous/fluid boundary using the same Appendix-B rows.

    The returned reflected coefficients are for waves propagating back into the
    porous medium. `R_TM_potential` includes the alpha_TM electric-potential
    ratio, while `R_Pf` is the acoustic Pf reflection coefficient.
    """
    params = jap_parameters()
    cfg = _jap_config(params)
    omega = 2.0 * math.pi * frequency_hz
    k1 = float(kx_key)
    state = strict.wave_slownesses(
        params.phi,
        params.k0_m2,
        params.alpha_inf,
        _cH(cfg),
        omega,
        cfg,
        C_override_molL=cfg.C_molL,
    )
    A = _boundary_matrix_for_fluid_incidence(state, omega, k1, cfg)
    incident = _pf_incident_column(state, omega, k1, cfg, upward=True)
    try:
        x = np.linalg.solve(A, -incident)
    except np.linalg.LinAlgError:
        x = np.full(6, np.nan + 1j * np.nan)
    return {
        **state,
        "R_Pf": x[2],
        # JAP Eq. (8) is written for electric scalar potential. The reused
        # Schakel 2010 solver stores alpha with the field-potential sign used
        # in its Appendix-B boundary rows; converting it to the JAP scalar
        # potential gives the polarity reversal across the interface noted
        # below Fig. 2.
        "R_TM_potential": -state["alpha_TM"] * x[4],
        "R_TM_raw": x[4],
    }


def _pf_interface_coefficients(frequency_hz: float, kx: float) -> dict[str, complex]:
    return _pf_interface_coefficients_cached(round(float(frequency_hz), 3), round(float(kx), 6))


def _pressure_normalized_re(frequency_hz: float, kx: float) -> complex:
    params = jap_parameters()
    omega = 2.0 * math.pi * frequency_hz
    coeff = _fluid_incident_coefficients(frequency_hz, kx)
    return coeff["R_E"] / (params.rho_f * omega**2)


def _front_terms(coeff: dict[str, complex], omega: float, params: JAPParameters) -> tuple[complex, complex, complex]:
    """JAP Table I pressure-normalized front-interface terms."""
    scale = params.rho_f * omega**2
    tf_tm = -coeff["alpha_TM"] * coeff["T_TM"] / scale
    tf_pf = coeff["T_Pf"] / scale
    alpha_pf = -coeff["alpha_Pf"]
    return tf_tm, tf_pf, alpha_pf


def _porous_eq8_terms(frequency_hz: float, kx: float, z_receiver: float) -> dict[str, complex]:
    params = jap_parameters()
    omega = 2.0 * math.pi * frequency_hz
    ws = params.sample_width_m
    coeff = _fluid_incident_coefficients(frequency_hz, kx)
    pf_coeff = _pf_interface_coefficients(frequency_hz, kx)
    tf_tm, tf_pf, alpha_pf = _front_terms(coeff, omega, params)
    k3_tm = coeff["k3_TM"]
    k3_pf = coeff["k3_Pf"]
    rb_tm = pf_coeff["R_TM_potential"]
    rb_pf = pf_coeff["R_Pf"]
    rf_tm = pf_coeff["R_TM_potential"]

    se_back_tm = rb_tm * cmath.exp(-1j * k3_tm * (ws - z_receiver))
    se_back_pf_coseismic = rb_pf * alpha_pf * cmath.exp(-1j * k3_pf * (ws - z_receiver))
    se_front_return_tm = (
        rb_pf
        * cmath.exp(-1j * k3_pf * ws)
        * rf_tm
        * cmath.exp(-1j * k3_tm * z_receiver)
    )
    se_total = se_back_tm + se_back_pf_coseismic + se_front_return_tm
    front_tm = tf_tm * cmath.exp(-1j * k3_tm * z_receiver)
    front_pf_coseismic = tf_pf * alpha_pf * cmath.exp(-1j * k3_pf * z_receiver)
    back_interface_sequence = tf_pf * cmath.exp(-1j * k3_pf * ws) * se_total
    return {
        "front_tm": front_tm,
        "front_pf_coseismic": front_pf_coseismic,
        "se_back_tm": se_back_tm,
        "se_back_pf_coseismic": se_back_pf_coseismic,
        "se_front_return_tm": se_front_return_tm,
        "se_total": se_total,
        "back_interface_sequence": back_interface_sequence,
        "eq8_total": front_tm + front_pf_coseismic + back_interface_sequence,
    }


def _integrate_complex(func, lo: float, hi: float) -> complex:
    if scipy_quad is None:
        grid = np.linspace(lo, hi, 80)
        return complex(np.trapezoid(np.asarray([func(float(x)) for x in grid]), grid))
    return complex(
        scipy_quad(lambda x: func(x).real, lo, hi, epsrel=1.0e-6, limit=80)[0],
        scipy_quad(lambda x: func(x).imag, lo, hi, epsrel=1.0e-6, limit=80)[0],
    )


def _fluid_response(frequency_hz: float, r_m: float, z_m: float, n_theta: int, method: str) -> complex:
    params = jap_parameters()
    omega = 2.0 * math.pi * frequency_hz
    k = omega / params.cP_fluid_m_s
    a = params.transducer_radius_m

    def real_integrand(theta: float) -> complex:
        sin_h = math.sin(theta)
        kx = k * sin_h
        re = _pressure_normalized_re(frequency_hz, kx)
        kz_e = _em_kz_fluid(omega, kx, params)
        term = _j0(k * r_m * sin_h) * _j1(k * a * sin_h)
        term *= cmath.exp(1j * k * params.source_z_m * math.cos(theta))
        term *= cmath.exp(1j * kz_e * z_m) * re
        return complex(term)

    def gamma_integrand(gamma: float) -> complex:
        root = math.sqrt(gamma**2 + 1.0)
        kx = k * root
        re = _pressure_normalized_re(frequency_hz, kx)
        kz_e = _em_kz_fluid(omega, kx, params)
        term = _j0(k * r_m * root) * (_j1(k * a * root) / root)
        term *= math.exp(k * params.source_z_m * gamma)
        term *= cmath.exp(1j * kz_e * z_m) * re
        return complex(term)

    if method == "adaptive":
        first = _integrate_complex(real_integrand, 0.0, 0.5 * math.pi)
        second = _integrate_complex(gamma_integrand, 0.0, 8.0)
    else:
        h = np.linspace(0.0, 0.5 * math.pi, n_theta)
        g = np.linspace(0.0, 8.0, max(12, n_theta // 2))
        first = complex(np.trapezoid(np.asarray([real_integrand(float(x)) for x in h]), h))
        second = complex(np.trapezoid(np.asarray([gamma_integrand(float(x)) for x in g]), g))
    A = _source_A_spectrum(np.array([frequency_hz]), params)[0]
    return -(1j * A / a) * first + (A / a) * second


def _porous_response(frequency_hz: float, r_m: float, z_m: float, n_theta: int, method: str) -> complex:
    params = jap_parameters()
    omega = 2.0 * math.pi * frequency_hz
    k = omega / params.cP_fluid_m_s
    a = params.transducer_radius_m

    def bracket(kx: float, z_receiver: float) -> complex:
        return _porous_eq8_terms(frequency_hz, kx, z_receiver)["eq8_total"]

    def real_integrand(theta: float) -> complex:
        sin_h = math.sin(theta)
        kx = k * sin_h
        term = _j0(k * r_m * sin_h) * _j1(k * a * sin_h)
        term *= cmath.exp(1j * k * params.source_z_m * math.cos(theta))
        return complex(term * bracket(kx, z_m))

    def gamma_integrand(gamma: float) -> complex:
        root = math.sqrt(gamma**2 + 1.0)
        kx = k * root
        term = _j0(k * r_m * root) * (_j1(k * a * root) / root)
        term *= math.exp(k * params.source_z_m * gamma)
        return complex(term * bracket(kx, z_m))

    if method == "adaptive":
        first = _integrate_complex(real_integrand, 0.0, 0.5 * math.pi)
        second = _integrate_complex(gamma_integrand, 0.0, 8.0)
    else:
        h = np.linspace(0.0, 0.5 * math.pi, n_theta)
        g = np.linspace(0.0, 8.0, max(12, n_theta // 2))
        first = complex(np.trapezoid(np.asarray([real_integrand(float(x)) for x in h]), h))
        second = complex(np.trapezoid(np.asarray([gamma_integrand(float(x)) for x in g]), g))
    A = _source_A_spectrum(np.array([frequency_hz]), params)[0]
    return -(1j * A / a) * first + (A / a) * second


def _frequency_grid(n_frequencies: int) -> np.ndarray:
    params = jap_parameters()
    return np.linspace(params.bandpass_low_hz, params.bandpass_high_hz, int(n_frequencies))


def _fig2_wave_packet_trace(time_ms: np.ndarray, raw_trace_mv: np.ndarray) -> np.ndarray:
    """Keep resolved arrival packets and suppress out-of-arrival side lobes.

    The continuous band-limited spectrum produces low-amplitude side lobes over
    the whole plotted time range. JAP Fig. 2 displays resolved arrival packets,
    not those numerical side lobes. This deterministic mask keeps samples near
    wave packets whose envelope exceeds the printed-trace visibility threshold.
    """
    time_ms = np.asarray(time_ms, dtype=float)
    raw = np.asarray(raw_trace_mv, dtype=float)
    peak = float(np.max(np.abs(raw))) if len(raw) else 0.0
    if peak <= 0.0:
        return raw.copy()
    visible = np.abs(raw) >= 0.30 * peak
    half_width_ms = 0.0018
    keep = np.zeros_like(visible, dtype=bool)
    for idx in np.where(visible)[0]:
        keep |= np.abs(time_ms - time_ms[idx]) <= half_width_ms
    return np.where(keep, raw, 0.0)


def synthesize_jap_fig2(
    n_frequencies: int = 81,
    n_theta: int = 64,
    integration_method: str = "adaptive",
) -> pd.DataFrame:
    frequencies = _frequency_grid(n_frequencies)
    time_s = np.linspace(0.055e-3, 0.145e-3, 720)
    rows = []
    for _, pos in JAP_POSITIONS.iterrows():
        response_fn = _fluid_response if pos["medium"] == "fluid" else _porous_response
        spec = np.array(
            [
                response_fn(
                    float(f),
                    float(pos["r_m"]),
                    float(pos["z_m"]),
                    n_theta,
                    integration_method,
                )
                for f in frequencies
            ],
            dtype=complex,
        )
        phase = np.exp(1j * 2.0 * math.pi * frequencies[:, None] * time_s[None, :])
        raw_trace_mv = 2.0 * np.real(np.trapezoid(spec[:, None] * phase, frequencies, axis=0)) * 1.0e3
        time_ms = time_s * 1.0e3
        trace_mv = _fig2_wave_packet_trace(time_ms, raw_trace_mv)
        for t, v, raw_v in zip(time_ms, trace_mv, raw_trace_mv):
            rows.append(
                {
                    "position": int(pos["position"]),
                    "medium": pos["medium"],
                    "z_m": float(pos["z_m"]),
                    "r_m": float(pos["r_m"]),
                    "time_ms": float(t),
                    "model_mv": float(v),
                    "raw_model_mv": float(raw_v),
                }
            )
    return pd.DataFrame(rows)


def _arrival_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pos, part in df.groupby("position"):
        t = part["time_ms"].to_numpy()
        y = part["model_mv"].to_numpy()
        peak_idx = int(np.argmax(np.abs(y)))
        early = part[(part["time_ms"] >= 0.096) & (part["time_ms"] <= 0.102)]
        second = part[(part["time_ms"] >= 0.102) & (part["time_ms"] <= 0.111)]
        early_time = np.nan
        early_peak = np.nan
        if len(early):
            early_idx = int(np.argmax(np.abs(early["model_mv"].to_numpy(dtype=float))))
            early_time = float(early["time_ms"].iloc[early_idx])
            early_peak = float(early["model_mv"].iloc[early_idx])
        second_time = np.nan
        second_peak = np.nan
        if len(second):
            second_idx = int(np.argmax(np.abs(second["model_mv"].to_numpy(dtype=float))))
            second_time = float(second["time_ms"].iloc[second_idx])
            second_peak = float(second["model_mv"].iloc[second_idx])
        rows.append(
            {
                "position": int(pos),
                "z_cm": float(part["z_m"].iloc[0] * 100.0),
                "medium": str(part["medium"].iloc[0]),
                "max_abs_mv": float(np.max(np.abs(y))),
                "global_peak_time_ms": float(t[peak_idx]),
                "early_0p10ms_peak_time_ms": early_time,
                "early_0p10ms_peak_mv": early_peak,
                "second_0p105ms_peak_time_ms": second_time,
                "second_0p105ms_peak_mv": second_peak,
                "second_minus_early_us": float((second_time - early_time) * 1000.0)
                if np.isfinite(second_time) and np.isfinite(early_time)
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _plot_jap_fig2(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(6, 2, figsize=(7.8, 8.8), sharex=True)
    for idx, pos in enumerate(range(1, 7)):
        part = df[df["position"] == pos]
        for col in range(2):
            ax = axes[idx, col]
            ax.axhline(0.0, color="0.3", linewidth=0.7)
            ax.set_xlim(0.055, 0.145)
            ax.set_ylim(-0.5 if col == 0 else -0.2, 0.5 if col == 0 else 0.2)
            ax.set_yticks([-0.5, 0.0, 0.5] if col == 0 else [-0.2, 0.0, 0.2])
            ax.text(0.50, 0.86, f"z\n{pos}", transform=ax.transAxes, ha="center", va="top", fontsize=9)
            if col == 0:
                ax.plot(part["time_ms"], part["model_mv"], color="black", linewidth=0.8)
            else:
                # The paper reports observed panels scaled down by a factor 2.5.
                # This panel is a diagnostic overlay using only that published
                # global scale, not an independent fit to the trace shape.
                ax.plot(part["time_ms"], part["model_mv"] / 2.5, color="black", linewidth=0.8)
                ax.text(0.02, 0.78, "model/2.5", transform=ax.transAxes, fontsize=6, color="0.35")
            if idx < 5:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Time (ms)")
    axes[0, 0].set_title("Predicted")
    axes[0, 1].set_title("Observed scale")
    axes[2, 0].set_ylabel("Electric potential (mV)")
    fig.tight_layout(h_pad=0.75, w_pad=0.8)
    fig.savefig(outdir / "jap_fig2_reproduction.png", dpi=300)
    plt.close(fig)


def _plot_model_vs_digitized(model_df: pd.DataFrame, digitized_df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(6, 2, figsize=(7.8, 8.8), sharex=True)
    for idx, pos in enumerate(range(1, 7)):
        model_part = model_df[model_df["position"] == pos]
        for col_idx, column in enumerate(["predicted", "observed"]):
            ax = axes[idx, col_idx]
            published = digitized_df[
                (digitized_df["position"] == pos) & (digitized_df["column"] == column)
            ]
            ax.axhline(0.0, color="0.35", linewidth=0.7)
            ax.plot(
                published["time_ms"],
                published["digitized_mv"],
                color="0.65",
                linewidth=1.0,
                label="published raster",
            )
            if column == "predicted":
                ax.plot(model_part["time_ms"], model_part["model_mv"], color="black", linewidth=0.75, label="model")
                ax.set_ylim(-0.5, 0.5)
                ax.set_yticks([-0.5, 0.0, 0.5])
            else:
                ax.plot(
                    model_part["time_ms"],
                    model_part["model_mv"] / 2.5,
                    color="black",
                    linewidth=0.75,
                    label="model/2.5",
                )
                ax.set_ylim(-0.2, 0.2)
                ax.set_yticks([-0.2, 0.0, 0.2])
            ax.set_xlim(FIG2_TIME_MIN_MS, FIG2_TIME_MAX_MS)
            ax.text(0.50, 0.86, f"z\n{pos}", transform=ax.transAxes, ha="center", va="top", fontsize=9)
            if idx < 5:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Time (ms)")
    axes[0, 0].set_title("Predicted")
    axes[0, 1].set_title("Observed")
    axes[2, 0].set_ylabel("Electric potential (mV)")
    axes[0, 0].legend(frameon=False, fontsize=6, loc="upper left")
    axes[0, 1].legend(frameon=False, fontsize=6, loc="upper left")
    fig.tight_layout(h_pad=0.75, w_pad=0.8)
    fig.savefig(outdir / "jap_fig2_model_vs_digitized.png", dpi=300)
    plt.close(fig)


def _write_fig2_digitization_artifacts(outdir: Path) -> pd.DataFrame:
    source_dir = outdir / "source_pages"
    page_png = _render_jap_page3(source_dir)
    page = Image.open(page_png)
    crop_box = (1370, 1580, 2310, 3095)
    page.crop(crop_box).save(outdir / "jap_fig2_source_crop.png")

    digitized = digitized_published_fig2_traces(page_png=page_png)
    digitized.to_csv(outdir / "digitized_published_fig2_waveforms.csv", index=False)
    pd.DataFrame(
        [
            {
                "figure": "Schakel_2011_JAP_Fig2",
                "source_pdf": str(_jap_pdf_path()),
                "rendered_page": "source_pages/jap_page-3.png",
                "source_crop": "jap_fig2_source_crop.png",
                "rendering": "pdftoppm page 3 at 300 dpi",
                "time_axis_ms": f"{FIG2_TIME_MIN_MS}..{FIG2_TIME_MAX_MS}",
                "predicted_axis_ylim_mv": FIG2_PAGE3_AXES["predicted"]["ylim_mv"],
                "observed_axis_ylim_mv": FIG2_PAGE3_AXES["observed"]["ylim_mv"],
                "method": (
                    "Raster trace extraction after removing long axis/tick lines; weak "
                    "unresolved columns are set to the zero baseline. The digitized data "
                    "are used only for visual correspondence checks, not for model tuning."
                ),
            }
        ]
    ).to_csv(outdir / "digitization_metadata_fig2.csv", index=False)
    return digitized


def _write_parameters(outdir: Path) -> None:
    params = jap_parameters()
    rows = [{"parameter": k, "value": v} for k, v in asdict(params).items()]
    cfg = _jap_config(params)
    rows.extend(
        [
            {"parameter": "M_similarity_computed_from_TableII_Lambda", "value": cfg.M_similarity},
            {"parameter": "C_molL_equivalent_for_pH6_zeta_minus_0p04", "value": cfg.C_molL},
        ]
    )
    pd.DataFrame(rows).to_csv(outdir / "parameters_used.csv", index=False)
    JAP_POSITIONS.to_csv(outdir / "jap_fig2_positions.csv", index=False)


def _write_audit(outdir: Path) -> None:
    params = jap_parameters()
    cfg = _jap_config(params)
    text = f"""# JAP Figure 2 Formula Audit

Implemented target:

- JAP Eq. (5) for fluid positions 1-5, using the same Sommerfeld path as the paper: real angles `0..pi/2` plus the `gamma` substitution `theta = pi/2 + i ln(sqrt(gamma^2+1)+gamma)`.
- JAP Eq. (8) for position 6 inside the porous sample. The bracket contains the front-interface TM term, the Pf coseismic term, and the Eq. (7) back/front interface-response terms.
- JAP Table II parameters are used. The Schakel dynamic-permeability helper computes `Lambda=sqrt(8 alpha_inf k0/(phi M))`; therefore `M={cfg.M_similarity:.8g}` is chosen to reproduce the Table II value `Lambda={params.Lambda:.6g} m`.
- The source is the digitized pressure trace already used for the Geophysics reproduction. Because Eq. (1) defines `A(omega)` through `p_hat=A/R exp(-ikR) D(theta)`, the on-axis pressure record at `(r,z)=(0,0)` is converted to `A(omega)` by multiplying by `R/D(0)=|z_s|/0.5`. No extra amplitude fit is applied.
- The 144-896 kHz numerical band-pass stated below Eq. (8) is used. Because
  the paper does not specify the digital filter shape, a 40 kHz raised-cosine
  edge is used to avoid nonphysical ringing from a rectangular spectral cut.

Coefficient normalization:

- `R^E`, `T_f^TM`, and `alpha T_f^Pf` use JAP Table I pressure normalization by `rho_f omega^2`.
- The porous-side electric scalar potential uses the opposite sign of the
  Schakel Appendix-B field-potential `alpha` storage convention. This sign
  conversion is applied only to the JAP Eq. (8) porous-potential terms and
  gives the position-6 polarity reversal explicitly described below Fig. 2.
- The fluid-incident coefficients come directly from `schakel2010_strict_sensitivity.se_coefficients`, the Schakel and Smeulders (2010) Appendix-B 6x6 system.
- The Pf-incident coefficients `R_b^TM`, `R_b^Pf`, and `R_f^TM` are obtained by reusing the same Appendix-B boundary rows and replacing the right-hand side by the boundary contribution of a unit upward Pf wave in the porous medium. `R_TM_potential` includes the `alpha_TM` electric-potential ratio. Because the JAP paper states only that these Pf-incident coefficients are derived by a procedure similar to Ref. 25 and does not print the matrix, this part must be treated as a transparent reconstruction of an omitted derivation, not as a verbatim printed-equation implementation.
- `R_b` and `R_f` use the same local porous/fluid-interface coefficient because the experiment has the same fluid and porous material on both sides of the sample; the code differs only through the propagation phases from Eq. (7).

Published Fig. 2 digitization:

- `digitized_published_fig2_waveforms.csv` is extracted from a 300 dpi render of the published Fig. 2. The predicted column is digitized on the `+/-0.5 mV` axes and the observed column is digitized on the `+/-0.2 mV` axes.
- These digitized traces are used only in `jap_fig2_model_vs_digitized.png` to compare timing, polarity, and approximate printed-trace amplitude. They are not fed back into the theoretical model and are not used as a fitting target.
- `jap_fig2_model_waveforms.csv` keeps both `raw_model_mv` and `model_mv`. `raw_model_mv` is the direct finite-band spectral integral and contains low-amplitude side lobes outside physical arrivals. `model_mv` is the Fig. 2 display trace after deterministic wave-packet visibility masking, matching the paper's presentation of resolved arrival packets rather than continuous finite-band side lobes.

Important limitations:

- The JAP paper states that the Pf-incident coefficients are derived "in a procedure similar" to Ref. 25 but does not print the corresponding matrix. The implementation records the matrix extension explicitly in code; it is not an empirical fit, but it should still be treated as a reproducible reconstruction rather than a verbatim printed formula.
- The original measured pressure waveform is not machine-readable. The right column in `jap_fig2_reproduction.png` shows the model divided by the paper's published global observed/predicted scale factor 2.5; the separate `jap_fig2_model_vs_digitized.png` overlays the page-digitized observed traces.
- Later multipulse timing in position 6 follows Eq. (7)-(8). No pulse-specific delay or amplitude tuning is introduced.
- The generated observed-scale column is `model/2.5`, using the global scale factor stated in the paper. It is not a digitized observed waveform and should not be used for measured-trace residuals.
- The wave-packet visibility mask is a plotting/reproduction step, not a replacement for the frequency-domain solution. Use `raw_model_mv` for diagnosing finite-band spectral leakage.
"""
    (outdir / "formula_audit_jap_fig2.md").write_text(text, encoding="utf-8")


def run_reproduction(
    outdir: str | Path = "results/Schakel2011_JAP_Fig2",
    n_frequencies: int = 81,
    n_theta: int = 64,
    integration_method: str = "adaptive",
) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    _write_parameters(out)
    df = synthesize_jap_fig2(n_frequencies=n_frequencies, n_theta=n_theta, integration_method=integration_method)
    df.to_csv(out / "jap_fig2_model_waveforms.csv", index=False)
    _arrival_summary(df).to_csv(out / "jap_fig2_arrival_summary.csv", index=False)
    digitized = _write_fig2_digitization_artifacts(out)
    _plot_jap_fig2(df, out)
    _plot_model_vs_digitized(df, digitized, out)
    _write_audit(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/Schakel2011_JAP_Fig2")
    parser.add_argument("--n-frequencies", type=int, default=81)
    parser.add_argument("--n-theta", type=int, default=64)
    parser.add_argument("--integration-method", choices=["adaptive", "fixed"], default="adaptive")
    args = parser.parse_args()
    run_reproduction(args.outdir, args.n_frequencies, args.n_theta, args.integration_method)


if __name__ == "__main__":
    main()
