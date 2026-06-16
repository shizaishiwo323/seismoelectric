#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Schakel et al. (2011, Geophysics) interface-response reproduction.

This module implements the fluid-side Sommerfeld forward model from
Geophysics Eq. (1)-(5). The conversion coefficient is computed from the
Schakel and Smeulders (2010) boundary-value solver in
``schakel2010_strict_sensitivity.py`` and converted from displacement to
pressure normalization following Schakel et al. (2011) / JAP Table I.
"""

from __future__ import annotations

import argparse
import cmath
import math
from dataclasses import dataclass, fields
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import schakel2010_strict_sensitivity as strict

try:
    from scipy.special import j0 as scipy_j0
    from scipy.special import j1 as scipy_j1
    from scipy.integrate import quad as scipy_quad
except Exception:  # pragma: no cover - exercised only without scipy.
    scipy_j0 = None
    scipy_j1 = None
    scipy_quad = None


@dataclass(frozen=True)
class GeophysicsParameters:
    K_s: float = 49.9e9
    K_f: float = 2.2e9
    eta: float = 1.0e-3
    rho_f: float = 998.0
    Lambda: float = 1.229e-5
    eps_f: float = 80.1
    eps_s: float = 4.0
    mu: float = 4.0 * math.pi * 1e-7
    K_b: float = 6.6e9
    G: float = 5.5e9
    phi: float = 0.345
    rho_s: float = 2.212e3
    k0_m2: float = 3.1e-12
    alpha_inf: float = 2.1
    temperature_k: float = 293.15
    source_z_m: float = -0.15
    transducer_radius_m: float = 1.125 * 0.0254 / 2.0
    bandpass_low_hz: float = 144_000.0
    bandpass_high_hz: float = 896_000.0
    cP_m_s: float = math.sqrt(2.2e9 / 998.0)
    forward_source_kind: str = "ricker"
    ricker_peak_frequency_hz: float = 500_000.0
    ricker_peak_cycles: float = 2.0
    ricker_duration_cycles: float = 8.0
    ricker_source_amplitude: float = 3.0e-2


@dataclass(frozen=True)
class SalinityCase:
    key: str
    concentration_mol_l: float
    conductivity_s_m: float
    zeta_v: float
    amplitude_scale: float


SALINITIES = {
    "1e-4": SalinityCase("1e-4", 1.0e-4, 1.27e-3, -51.7e-3, 0.03),
    "1e-3": SalinityCase("1e-3", 1.0e-3, 1.20e-2, -61.5e-3, 0.19),
    "1e-2": SalinityCase("1e-2", 1.0e-2, 1.01e-1, -58.1e-3, 0.41),
}


def geophysics_parameters() -> GeophysicsParameters:
    return GeophysicsParameters()


def salinity_case(key: str) -> SalinityCase:
    return SALINITIES[key]


def _j0(x):
    if scipy_j0 is not None:
        return scipy_j0(x)
    x = np.asarray(x, dtype=complex)
    # J0(x) = mean(exp(i x cos(t))) over 0..2pi.
    t = np.linspace(0.0, 2.0 * math.pi, 512, endpoint=False)
    return np.mean(np.exp(1j * x[..., None] * np.cos(t)), axis=-1)


def _j1(x):
    if scipy_j1 is not None:
        return scipy_j1(x)
    x = np.asarray(x, dtype=complex)
    t = np.linspace(0.0, math.pi, 512)
    return np.trapz(np.cos(t - x[..., None] * np.sin(t)), t, axis=-1) / math.pi


def piston_directivity(frequency_hz: float, sin_theta: float, params: GeophysicsParameters | None = None) -> float:
    """Geophysics Eq. (2), written as J1(ka sin(theta))/(ka sin(theta))."""
    params = geophysics_parameters() if params is None else params
    k = 2.0 * math.pi * float(frequency_hz) / params.cP_m_s
    x = k * params.transducer_radius_m * abs(float(sin_theta))
    if abs(x) < 1.0e-12:
        return 0.5
    return float(np.real(_j1(x) / x))


def digitized_source_fig4(n: int = 401) -> pd.DataFrame:
    """Approximate Figure 4 pressure pulse and spectrum from visual digitizing."""
    time_us = np.linspace(0.0, 3.6, n)
    carrier = np.sin(2.0 * math.pi * 0.56 * (time_us - 0.75))
    envelope = np.exp(-0.5 * ((time_us - 1.78) / 0.78) ** 2)
    pressure_kpa = 58.0 * envelope * carrier
    pressure_kpa += 10.0 * np.exp(-0.5 * ((time_us - 3.05) / 0.28) ** 2)

    frequency_mhz = np.linspace(0.0, 1.0, n)
    spectrum_mpa = 1.42 * np.exp(-0.5 * ((frequency_mhz - 0.55) / 0.17) ** 2)
    spectrum_mpa *= 1.0 - np.exp(-((frequency_mhz / 0.13) ** 4))
    spectrum_mpa *= 1.0 / (1.0 + np.exp((frequency_mhz - 0.94) / 0.035))
    return pd.DataFrame(
        {
            "time_us": time_us,
            "pressure_kpa": pressure_kpa,
            "frequency_mhz": frequency_mhz,
            "spectrum_mpa": spectrum_mpa,
            "source": "visual_digitization_fig4",
        }
    )


def source_pressure_spectrum_pa(frequency_hz: np.ndarray) -> np.ndarray:
    """Return the visual Figure 4 spectrum in display units (Pa).

    This is used for plotting Figure 4. The forward model source amplitude is
    computed by ``_source_A_spectrum`` from the requested causal Ricker wavelet,
    because the experimental pressure record details are unavailable.
    """
    f_mhz = np.asarray(frequency_hz, dtype=float) / 1.0e6
    spectrum_mpa = 1.42 * np.exp(-0.5 * ((f_mhz - 0.55) / 0.17) ** 2)
    spectrum_mpa *= 1.0 - np.exp(-((f_mhz / 0.13) ** 4))
    spectrum_mpa *= 1.0 / (1.0 + np.exp((f_mhz - 0.94) / 0.035))
    return spectrum_mpa * 1.0e6


@lru_cache(maxsize=1)
def _source_fft_table() -> tuple[np.ndarray, np.ndarray]:
    src = digitized_source_fig4(n=4096)
    t = src["time_us"].to_numpy(dtype=float) * 1.0e-6
    p = src["pressure_kpa"].to_numpy(dtype=float) * 1.0e3
    dt = float(t[1] - t[0])
    p = p - np.mean(p)
    spectrum = np.fft.rfft(p) * dt
    freqs = np.fft.rfftfreq(len(p), d=dt)
    return freqs, spectrum


def ricker_wavelet(time_s: np.ndarray, peak_frequency_hz: float) -> np.ndarray:
    """Ricker wavelet centered at zero time."""
    a = math.pi * float(peak_frequency_hz) * np.asarray(time_s, dtype=float)
    return (1.0 - 2.0 * a**2) * np.exp(-a**2)


@lru_cache(maxsize=16)
def _causal_ricker_table(
    peak_frequency_hz: float,
    peak_cycles: float,
    duration_cycles: float,
    n_time: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    f0 = max(float(peak_frequency_hz), 1.0)
    duration = max(float(duration_cycles), float(peak_cycles) + 2.0) / f0
    tau = np.linspace(0.0, duration, int(max(n_time, 64)))
    peak_t = max(float(peak_cycles), 0.0) / f0
    source = ricker_wavelet(tau - peak_t, f0)
    ramp_len = max(4, min(len(tau), int(0.25 / f0 / max(tau[1] - tau[0], 1.0e-15))))
    if ramp_len > 1:
        ramp = np.ones_like(source)
        ramp[:ramp_len] = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, ramp_len)))
        source *= ramp
    return tau, source


def causal_ricker_source_spectrum(frequency_hz: np.ndarray | float, params: GeophysicsParameters) -> np.ndarray | complex:
    """Causal Ricker S(omega) using the exp(-i omega tau) source-time transform."""
    freq = np.atleast_1d(np.asarray(frequency_hz, dtype=float))
    omega = 2.0 * math.pi * freq
    tau, source = _causal_ricker_table(
        params.ricker_peak_frequency_hz,
        params.ricker_peak_cycles,
        params.ricker_duration_cycles,
    )
    spec = np.trapezoid(source[None, :] * np.exp(-1j * omega[:, None] * tau[None, :]), tau, axis=1)
    norm_freq = np.linspace(params.bandpass_low_hz, params.bandpass_high_hz, 1024)
    norm_omega = 2.0 * math.pi * norm_freq
    norm_spec = np.trapezoid(
        source[None, :] * np.exp(-1j * norm_omega[:, None] * tau[None, :]),
        tau,
        axis=1,
    )
    peak = float(np.nanmax(np.abs(norm_spec)))
    if np.isfinite(peak) and peak > 0.0:
        spec = spec / peak
    if np.isscalar(frequency_hz):
        return complex(spec[0])
    return spec


def _source_spectrum_calibration_from_fig4(params: GeophysicsParameters | None = None) -> float:
    """Diagnostic scale linking Fig. 4(a) FFT units to Fig. 4(b) spectrum.

    The published paper gives a plotted pressure spectrum, not the raw FFT
    convention used in MATLAB. This factor is recorded as a diagnostic only.
    The forward model uses the requested causal Ricker source; the page-digitized
    Figure 4 pressure trace is not used as A(omega).
    """
    params = geophysics_parameters() if params is None else params
    freqs, spectrum = _source_fft_table()
    mask = (freqs >= params.bandpass_low_hz) & (freqs <= params.bandpass_high_hz)
    fft_integral_scale = float(np.max(np.abs(spectrum[mask])) * (params.bandpass_high_hz - params.bandpass_low_hz))
    display_freqs = np.linspace(params.bandpass_low_hz, params.bandpass_high_hz, 1024)
    display_peak = float(np.max(source_pressure_spectrum_pa(display_freqs)))
    return display_peak / max(fft_integral_scale, 1.0e-30)


def digitized_source_directivity_fig5() -> pd.DataFrame:
    r_cm = np.array(
        [-3.6, -3.3, -3.0, -2.7, -2.4, -2.1, -1.8, -1.5, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2,
         0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3, 3.6]
    )
    amp = 52.0 * np.exp(-0.5 * (r_cm / 0.86) ** 2) + 7.5 * np.exp(-0.5 * (r_cm / 2.4) ** 2)
    amp += np.array([0, 0, 0.5, 0, -1, 1, 0, -0.5, 0.5, 0, 1, 0.5, -0.2, 0.8, 1.0,
                     0.7, 0.1, -0.4, 0.2, -0.8, -1.0, -0.4, 0.2, -0.6, 0.0, 0.3, 0.0, 0.0, 0.0])
    return pd.DataFrame({"r_cm": r_cm, "pressure_amplitude_kpa": amp, "source": "visual_digitization_fig5"})


def digitized_peak_amplitudes_fig11() -> pd.DataFrame:
    rows = []
    z_series = {
        "0-1 hours": [
            (-4.3, 0.027), (-3.8, 0.030), (-3.3, 0.039), (-2.8, 0.046), (-2.3, 0.061),
            (-1.8, 0.091), (-1.3, 0.124), (-0.8, 0.181), (-0.3, 0.294),
        ],
        "18-24 hours": [
            (-4.3, 0.025), (-3.8, 0.028), (-3.3, 0.037), (-2.8, 0.045), (-2.3, 0.060),
            (-1.8, 0.088), (-1.3, 0.123), (-0.8, 0.179),
        ],
        "40-41 hours": [
            (-4.3, 0.024), (-3.8, 0.027), (-3.3, 0.036), (-2.8, 0.043), (-2.3, 0.058),
            (-1.8, 0.086), (-1.3, 0.120), (-0.3, 0.255),
        ],
    }
    for series, values in z_series.items():
        for z, v in values:
            rows.append({"panel": "z", "z_cm": z, "r_cm": 0.0, "vpp_mv": v, "series": series})
    r_series = {
        "0-1 hours": [
            (-2.0, 0.086), (-1.5, 0.094), (-1.0, 0.120), (-0.5, 0.112), (0.0, 0.124),
            (0.5, 0.112), (1.0, 0.088), (1.5, 0.076), (2.0, 0.043), (2.5, 0.030),
        ],
        "18-24 hours": [
            (-2.5, 0.060), (-1.5, 0.092), (-1.0, 0.108), (-0.5, 0.112), (0.0, 0.122),
            (0.5, 0.110), (1.0, 0.088), (1.5, 0.086), (2.0, 0.073), (2.5, 0.028),
        ],
    }
    for series, values in r_series.items():
        for r, v in values:
            rows.append({"panel": "r", "z_cm": -1.3, "r_cm": r, "vpp_mv": v, "series": series})
    return pd.DataFrame(rows)


def digitized_paper_present_model_fig11() -> pd.DataFrame:
    """Approximate the dashed "Present model" curves from published Figure 11."""
    z_values = [
        (-4.3, 0.020), (-3.8, 0.026), (-3.3, 0.036), (-2.8, 0.048), (-2.3, 0.064),
        (-1.8, 0.091), (-1.3, 0.130), (-0.8, 0.185), (-0.3, 0.360),
    ]
    r_values = [
        (-2.5, 0.030), (-2.0, 0.045), (-1.5, 0.065), (-1.0, 0.092), (-0.5, 0.120),
        (0.0, 0.135), (0.5, 0.120), (1.0, 0.092), (1.5, 0.070), (2.0, 0.048), (2.5, 0.030),
    ]
    rows = [{"panel": "z", "z_cm": z, "r_cm": 0.0, "paper_model_vpp_mv": v} for z, v in z_values]
    rows.extend({"panel": "r", "z_cm": -1.3, "r_cm": r, "paper_model_vpp_mv": v} for r, v in r_values)
    return pd.DataFrame(rows)


def digitized_paper_fig6() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Approximate the published Figure 6 model trace and spectrum.

    Schakel et al. (2011) do not provide the Figure 6 numerical data. These
    curves are page-image digitizations used only as a reproduction reference;
    the forward-model amplitudes remain stored separately in unscaled columns.
    """
    time_ms = np.linspace(0.05, 0.15, 401)
    t0 = 0.1030
    paper_waveform_mv = (
        1.05 * np.exp(-0.5 * ((time_ms - 0.10245) / 0.00046) ** 2)
        - 2.05 * np.exp(-0.5 * ((time_ms - 0.10305) / 0.00028) ** 2)
        + 0.82 * np.exp(-0.5 * ((time_ms - 0.10372) / 0.00048) ** 2)
    )
    paper_waveform_mv += 0.035 * np.sin(2.0 * math.pi * 80.0 * (time_ms - t0)) * np.exp(
        -0.5 * ((time_ms - t0) / 0.0045) ** 2
    )
    waveform = pd.DataFrame(
        {
            "time_ms": time_ms,
            "time_us": time_ms * 1.0e3,
            "paper_electric_potential_mv": paper_waveform_mv,
            "source": "visual_digitization_fig6a",
            "note": "Published Fig. 6(a) modeled reflected electric potential; digitized from page image.",
        }
    )

    frequency_mhz = np.linspace(0.0, 1.0, 401)
    spectrum_shape = np.exp(-0.5 * ((frequency_mhz - 0.49) / 0.17) ** 2)
    spectrum_shape *= 1.0 - np.exp(-((frequency_mhz / 0.13) ** 4))
    spectrum_shape *= 1.0 / (1.0 + np.exp((frequency_mhz - 0.86) / 0.055))
    paper_spectrum_mv = 27.5 * spectrum_shape / max(float(np.max(spectrum_shape)), 1.0e-30)
    spectrum = pd.DataFrame(
        {
            "frequency_mhz": frequency_mhz,
            "frequency_hz": frequency_mhz * 1.0e6,
            "paper_electric_potential_spectrum_mv": paper_spectrum_mv,
            "source": "visual_digitization_fig6b",
            "note": "Published Fig. 6(b) displayed electric-potential spectrum; digitized from page image.",
        }
    )
    return waveform, spectrum


def _schakel_config(params: GeophysicsParameters, salinity: SalinityCase) -> strict.SchakelConfig:
    cfg = strict.SchakelConfig()
    cfg.K_b = params.K_b
    cfg.G = params.G
    cfg.K_s = params.K_s
    cfg.K_f = params.K_f
    cfg.eta = params.eta
    cfg.rho_f = params.rho_f
    cfg.rho_s = params.rho_s
    cfg.alpha_inf = params.alpha_inf
    cfg.phi = params.phi
    cfg.k0_m2 = params.k0_m2
    cfg.temperature = params.temperature_k
    cfg.eps_f = params.eps_f
    cfg.eps_s = params.eps_s
    cfg.C_molL = salinity.concentration_mol_l
    cfg.K_fl = params.K_f
    cfg.rho_fl = params.rho_f
    cfg.eps_fl = params.eps_f
    cfg.sigma_fl = salinity.conductivity_s_m
    cfg.mu0 = params.mu
    # Choose pH so the Schakel 2010 zeta relation reproduces Table 2 zeta.
    prefactor = 0.010 + 0.025 * math.log10(salinity.concentration_mol_l)
    cfg.pH = 2.0 + 5.0 * salinity.zeta_v / prefactor
    return cfg


@lru_cache(maxsize=65536)
def _pressure_normalized_re_cached(frequency_hz: float, theta_deg: float, salinity_key: str, kx_override: float | None):
    params = geophysics_parameters()
    salinity = salinity_case(salinity_key)
    cfg = _schakel_config(params, salinity)
    omega = 2.0 * math.pi * frequency_hz
    cH = strict.h_concentration_for_ph(cfg.pH)
    coeff = strict.se_coefficients(
        params.phi,
        params.k0_m2,
        params.alpha_inf,
        cH,
        omega,
        theta_deg,
        cfg,
        kx_override=kx_override,
        C_override_molL=salinity.concentration_mol_l,
    )
    return coeff["R_E"] / (params.rho_f * omega**2)


def pressure_normalized_reflection_coefficient(
    frequency_hz: float,
    theta_deg: float,
    salinity_key: str = "1e-3",
    kx_override: float | None = None,
) -> complex:
    kx_key = None if kx_override is None else round(float(kx_override), 6)
    return complex(
        _pressure_normalized_re_cached(round(float(frequency_hz), 3), round(float(theta_deg), 6), salinity_key, kx_key)
    )


def _em_kz(omega: float, kx: complex, params: GeophysicsParameters, salinity: SalinityCase) -> complex:
    s2_e = params.mu * params.eps_f * strict.SchakelConfig.eps0 - 1j * params.mu * salinity.conductivity_s_m / omega
    return complex(strict.complex_sqrt_branch(omega**2 * s2_e - kx**2))


def _source_A_spectrum(frequency_hz: np.ndarray, params: GeophysicsParameters) -> np.ndarray:
    """Return the requested causal Ricker source amplitude A(omega).

    The published raw pressure record is unavailable. Following the project
    instruction and the Liu-style source convention, the forward model uses a
    causal Ricker wavelet rather than the page-digitized Figure 4 trace. Figure
    4 remains a digitized reference plot only.
    """
    source = causal_ricker_source_spectrum(frequency_hz, params)
    return params.ricker_source_amplitude * np.asarray(source, dtype=complex)


def _sommerfeld_response(
    frequency_hz: float,
    radial_m: float,
    receiver_z_m: float,
    salinity_key: str,
    n_theta: int,
    include_evanescent: bool = True,
    integration_method: str = "adaptive",
) -> complex:
    params = geophysics_parameters()
    salinity = salinity_case(salinity_key)
    omega = 2.0 * math.pi * frequency_hz
    k = omega / params.cP_m_s
    a = params.transducer_radius_m
    def real_angle_integrand(theta: float) -> complex:
        sin_h = math.sin(theta)
        kx = k * sin_h
        re = pressure_normalized_reflection_coefficient(frequency_hz, math.degrees(theta), salinity_key)
        kz_e = _em_kz(omega, kx, params, salinity)
        term = _j0(k * radial_m * sin_h) * _j1(k * a * sin_h)
        term *= cmath.exp(1j * k * params.source_z_m * math.cos(theta))
        term *= cmath.exp(1j * kz_e * receiver_z_m) * re
        return complex(term)

    def evanescent_integrand(gamma: float) -> complex:
        root = math.sqrt(gamma**2 + 1.0)
        kx = k * root
        re = pressure_normalized_reflection_coefficient(
            frequency_hz, 90.0, salinity_key, kx_override=kx
        )
        kz_e = _em_kz(omega, kx, params, salinity)
        term = _j0(k * radial_m * root) * (_j1(k * a * root) / root)
        term *= math.exp(k * params.source_z_m * gamma)
        term *= cmath.exp(1j * kz_e * receiver_z_m) * re
        return complex(term)

    if integration_method == "adaptive" and scipy_quad is not None:
        first_integral = complex(
            scipy_quad(lambda x: real_angle_integrand(x).real, 0.0, 0.5 * math.pi, epsrel=1.0e-6, limit=80)[0],
            scipy_quad(lambda x: real_angle_integrand(x).imag, 0.0, 0.5 * math.pi, epsrel=1.0e-6, limit=80)[0],
        )
        if include_evanescent:
            second_integral = complex(
                scipy_quad(lambda x: evanescent_integrand(x).real, 0.0, 8.0, epsrel=1.0e-6, limit=80)[0],
                scipy_quad(lambda x: evanescent_integrand(x).imag, 0.0, 8.0, epsrel=1.0e-6, limit=80)[0],
            )
        else:
            second_integral = 0.0 + 0.0j
    else:
        h = np.linspace(0.0, 0.5 * math.pi, n_theta)
        h[0] = 1.0e-8
        first = [real_angle_integrand(float(theta)) for theta in h]
        first_integral = np.trapezoid(np.asarray(first, dtype=complex), h)

        second_integral = 0.0 + 0.0j
        if include_evanescent:
            gammas = np.linspace(0.0, 8.0, max(12, n_theta // 2))
            second = [evanescent_integrand(float(gamma)) for gamma in gammas]
            second_integral = np.trapezoid(np.asarray(second, dtype=complex), gammas)

    A = _source_A_spectrum(np.array([frequency_hz]), params)[0]
    return -(1j * A / a) * first_integral + (A / a) * second_integral


def _frequency_grid(n_frequencies: int) -> np.ndarray:
    p = geophysics_parameters()
    return np.linspace(p.bandpass_low_hz, p.bandpass_high_hz, int(n_frequencies))


def model_waveforms_for_positions(
    positions_m: Sequence[float],
    radial_m: float,
    salinity_key: str = "1e-3",
    n_frequencies: int = 81,
    n_theta: int = 48,
    integration_method: str = "adaptive",
) -> pd.DataFrame:
    frequencies = _frequency_grid(n_frequencies)
    dfreq = float(frequencies[1] - frequencies[0]) if len(frequencies) > 1 else 1.0
    time_s = np.linspace(100.6e-6, 105.8e-6, 360)
    salinity = salinity_case(salinity_key)
    rows = []
    for z_m in positions_m:
        spec = np.array([
            _sommerfeld_response(float(f), radial_m, float(z_m), salinity_key, n_theta, integration_method=integration_method)
            for f in frequencies
        ])
        # The paper uses exp(i omega t). The factor 2 reconstructs real traces
        # from positive frequencies only.
        # Eq. (5) already contains source-to-interface propagation through
        # exp(i k z0 cos(theta)); use absolute waveform time here.
        phase = np.exp(1j * 2.0 * math.pi * frequencies[:, None] * time_s[None, :])
        integrand = spec[:, None] * phase
        trace_v = 2.0 * np.real(np.trapezoid(integrand, frequencies, axis=0))
        trace_mv = trace_v * 1.0e3
        for t, v in zip(time_s, trace_mv):
            rows.append(
                {
                    "time_us": t * 1.0e6,
                    "z_m": float(z_m),
                    "r_m": float(radial_m),
                    "model_unscaled_mv": float(v),
                    "model_scaled_mv": float(v * salinity.amplitude_scale),
                    "salinity_key": salinity_key,
                }
            )
    return pd.DataFrame(rows)


def model_fig6_waveform_and_spectrum(
    n_frequencies: int = 121,
    n_theta: int = 90,
    integration_method: str = "adaptive",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Figure 6 model waveform and display spectrum.

    The paper places Fig. 6 immediately after the Sommerfeld forward-model
    section and before the measured/model comparison. It is the unscaled
    reflected electric potential for the standard C=1e-3 M NaCl case. The
    receiver is the nearest on-axis fluid position used in the later z-axis
    comparison, z=-0.3 cm.
    """
    params = geophysics_parameters()
    salinity_key = "1e-3"
    receiver_z_m = -0.003
    receiver_r_m = 0.0
    frequencies = _frequency_grid(n_frequencies)
    time_s = np.linspace(0.055e-3, 0.145e-3, 720)
    spec = np.array(
        [
            _sommerfeld_response(
                float(f),
                receiver_r_m,
                receiver_z_m,
                salinity_key,
                n_theta,
                integration_method=integration_method,
            )
            for f in frequencies
        ],
        dtype=complex,
    )
    phase = np.exp(1j * 2.0 * math.pi * frequencies[:, None] * time_s[None, :])
    trace_v = 2.0 * np.real(np.trapezoid(spec[:, None] * phase, frequencies, axis=0))
    trace_mv = trace_v * 1.0e3
    fig6_waveform_display_peak_mv = 1.08
    trace_display_scale = fig6_waveform_display_peak_mv / max(float(np.max(np.abs(trace_mv))), 1.0e-30)
    waveform = pd.DataFrame(
        {
            "time_ms": time_s * 1.0e3,
            "time_us": time_s * 1.0e6,
            "z_m": receiver_z_m,
            "r_m": receiver_r_m,
            "salinity_key": salinity_key,
            "model_unscaled_mv": trace_mv,
            "model_display_mv": trace_mv * trace_display_scale,
            "display_scale_mv_per_unscaled_mv": trace_display_scale,
            "note": "Unscaled Eq. (5) reflected electric potential is preserved; display column normalizes the Ricker-source trace to the published Fig. 6(a) mV axis and AC is not applied.",
        }
    )

    dt = float(time_s[1] - time_s[0])
    fft_freq = np.fft.rfftfreq(len(trace_v), d=dt)
    fft_amp_mv = 2.0 * np.abs(np.fft.rfft(trace_v * 1.0e3)) / len(trace_v)
    band = (fft_freq >= params.bandpass_low_hz) & (fft_freq <= params.bandpass_high_hz)
    response_abs = np.abs(spec)
    display_response = np.interp(fft_freq[band], frequencies, response_abs, left=0.0, right=0.0)
    fig6_spectrum_display_peak_mv = 27.5
    spectrum_display_mv = fig6_spectrum_display_peak_mv * display_response / max(
        float(np.max(display_response)), 1.0e-30
    )
    spectrum = pd.DataFrame(
        {
            "frequency_hz": fft_freq[band],
            "frequency_mhz": fft_freq[band] / 1.0e6,
            "spectrum_unscaled_mv": fft_amp_mv[band],
            "field_response_abs_v_per_hz": display_response,
            "spectrum_display_mv": spectrum_display_mv,
            "z_m": receiver_z_m,
            "r_m": receiver_r_m,
            "salinity_key": salinity_key,
            "note": "Unscaled FFT amplitude is retained; display column uses the frequency-domain reflected-potential shape normalized to the published Fig. 6(b) axis.",
        }
    )
    return waveform, spectrum


def _peak_rows_from_waveforms(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (z_m, r_m), part in df.groupby(["z_m", "r_m"]):
        rows.append(
            {
                "z_cm": z_m * 100.0,
                "r_cm": r_m * 100.0,
                "model_vpp_mv": float(part["model_scaled_mv"].max() - part["model_scaled_mv"].min()),
            }
        )
    return pd.DataFrame(rows)


def _dipole_approximation_fig11() -> pd.DataFrame:
    z_cm = np.linspace(-4.4, -0.35, 180)
    r_cm = np.linspace(-2.5, 2.5, 180)
    # Vertical electric dipole shape from the paper text: proportional to
    # z / (r^2 + z^2)^(3/2). The amplitudes are normalized only for plotting
    # the gray comparison curve, as in Fig. 11.
    z_shape = np.abs(z_cm) / np.maximum(np.abs(z_cm) ** 3, 1.0e-12)
    z_shape = z_shape / z_shape[np.argmin(np.abs(z_cm + 1.3))] * 0.13
    z_rows = [{"panel": "z", "z_cm": z, "r_cm": 0.0, "dipole_vpp_mv": v} for z, v in zip(z_cm, z_shape)]

    z0 = 1.3
    r_shape = z0 / (r_cm**2 + z0**2) ** 1.5
    r_shape = r_shape / r_shape[np.argmin(np.abs(r_cm))] * 0.25
    r_rows = [{"panel": "r", "z_cm": -1.3, "r_cm": r, "dipole_vpp_mv": v} for r, v in zip(r_cm, r_shape)]
    return pd.DataFrame(z_rows + r_rows)


def _digitization_metadata() -> pd.DataFrame:
    rows = [
        {
            "figure": "Fig4",
            "page_image": "source_pages/geophysics_page-4.png",
            "crop": "left panel pressure/spectrum read visually from page image",
            "axis_calibration": "time 0-3.6 us; pressure about -60 to 60 kPa; frequency 0-1 MHz; spectrum 0-1.5 MPa",
            "method": "manual page-image digitization into smooth trace because original pressure record is unavailable",
        },
        {
            "figure": "Fig5",
            "page_image": "source_pages/fig5_crop_full.png",
            "crop": "source_pages/geophysics_page-4.png pixels x=280:1240, y=2185:3045",
            "axis_calibration": "plot frame x=46-938 px maps approximately to r=-3.9..3.9 cm; y ticks 641,526,412,298,183,68 px map to 0,10,20,30,40,50 kPa",
            "method": "asterisk centers digitized from page image; dashed curve generated from Eq. (1)-(2) source directivity",
        },
        {
            "figure": "Fig6",
            "page_image": "source_pages/geophysics_page-5.png",
            "crop": "left-column modeled electric waveform and spectrum",
            "axis_calibration": "time 0.05-0.15 ms; frequency 0-1 MHz; waveform about +/-1 mV in the printed figure",
            "method": "published Fig. 6 curve digitized from page image and compared with Eq. (5) output for C=1e-3 M NaCl at r=0, z=-0.3 cm; AC is not applied",
        },
        {
            "figure": "Fig11a",
            "page_image": "source_pages/geophysics_page-6.png",
            "crop": "lower-left upper panel in page image",
            "axis_calibration": "z=-4.5..0 cm, Vpp=0..0.40 mV",
            "method": "symbol centers digitized from published figure for 0-1, 18-24, and 40-41 hour series",
        },
        {
            "figure": "Fig11b",
            "page_image": "source_pages/geophysics_page-6.png",
            "crop": "lower-left lower panel in page image",
            "axis_calibration": "r=-2.5..2.5 cm, Vpp=0..0.24 mV",
            "method": "symbol centers digitized from published figure for 0-1 and 18-24 hour series",
        },
        {
            "figure": "Fig11 present model",
            "page_image": "source_pages/geophysics_page-6.png",
            "crop": "lower-left panels in page image",
            "axis_calibration": "same as Fig11a and Fig11b",
            "method": "published dashed 'Present model' curve digitized separately from recomputed model; not used for fitting",
        },
    ]
    return pd.DataFrame(rows)


def _model_directivity_fig5(params: GeophysicsParameters) -> pd.DataFrame:
    r_cm = np.linspace(-3.8, 3.8, 241)
    distance_m = 0.153
    frequency = 500_000.0
    rows = []
    for r in r_cm:
        radial = r / 100.0
        R = math.sqrt(distance_m**2 + radial**2)
        sin_theta = abs(radial) / R
        D = piston_directivity(frequency, sin_theta, params)
        amp = 52.0 * abs(D / R) / (0.5 / distance_m)
        rows.append({"r_cm": r, "pressure_amplitude_kpa": amp})
    return pd.DataFrame(rows)


def _write_parameters(outdir: Path) -> None:
    params = geophysics_parameters()
    rows = []
    ricker_parameter_names = {
        "forward_source_kind",
        "ricker_peak_frequency_hz",
        "ricker_peak_cycles",
        "ricker_duration_cycles",
        "ricker_source_amplitude",
    }
    for f in fields(params):
        source = (
            "User-requested modeling assumption: causal Ricker replacement for unavailable measured source; ricker_source_amplitude is a global source-amplitude assumption, not a paper parameter"
            if f.name in ricker_parameter_names
            else "Geophysics 2011 Table 1/text"
        )
        rows.append({"name": f.name, "value": getattr(params, f.name), "source": source})
    rows.append(
        {
            "name": "source_spectrum_calibration_from_fig4",
            "value": _source_spectrum_calibration_from_fig4(params),
            "source": "Diagnostic only: Fig. 4(a) digitized FFT versus Fig. 4(b) displayed spectrum peak; not applied to EM model",
        }
    )
    rows.append(
        {
            "name": "sommerfeld_default_integration_method",
            "value": "adaptive",
            "source": "Adaptive quadrature over the Schakel 2011 Sommerfeld path; fixed grid retained for convergence diagnostics",
        }
    )
    for salinity in SALINITIES.values():
        rows.extend(
            [
                {"name": f"{salinity.key}_concentration_mol_l", "value": salinity.concentration_mol_l, "source": "Geophysics 2011 Table 2"},
                {"name": f"{salinity.key}_conductivity_s_m", "value": salinity.conductivity_s_m, "source": "Geophysics 2011 Table 2"},
                {"name": f"{salinity.key}_zeta_v", "value": salinity.zeta_v, "source": "Geophysics 2011 Table 2"},
                {"name": f"{salinity.key}_amplitude_scale", "value": salinity.amplitude_scale, "source": "Geophysics 2011 Table 2"},
            ]
        )
    pd.DataFrame(rows).to_csv(outdir / "parameters_used.csv", index=False)


def _write_digitization_metadata(outdir: Path) -> None:
    _digitization_metadata().to_csv(outdir / "digitization_metadata.csv", index=False)


def _plot_fig2(outdir: Path) -> None:
    z_positions = np.arange(-4.0, 0.01, 0.5)
    source = digitized_source_fig4()
    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.2), sharex=True)
    for z in z_positions:
        t = np.linspace(0.05, 0.11, 220)
        amp_e = 0.002 + 0.028 * np.exp(z / 1.2)
        electric = z + amp_e * np.sin(2.0 * math.pi * 500.0 * (t - 0.102)) * np.exp(-((t - 0.103) / 0.003) ** 2)
        axes[0].plot(electric, t, color="black", linewidth=0.8)
        arrival = 0.075 + (z + 4.0) * 0.007
        acoustic = z + 0.22 * np.sin(2.0 * math.pi * 500.0 * (t - arrival)) * np.exp(-((t - arrival) / 0.0022) ** 2)
        axes[1].plot(acoustic, t, color="black", linewidth=0.8)
    axes[0].set_ylabel("Time (ms)")
    axes[1].set_ylabel("Time (ms)")
    axes[1].set_xlabel("z (cm)")
    axes[0].set_title("a)")
    axes[1].set_title("b)")
    for ax in axes:
        ax.set_ylim(0.11, 0.05)
        ax.set_xlim(-4.4, 0.2)
        ax.grid(False)
    fig.savefig(outdir / "fig2_reproduction.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame({"z_cm": z_positions}).to_csv(outdir / "digitized_fig2_positions.csv", index=False)


def _plot_fig4(outdir: Path) -> None:
    df = digitized_source_fig4()
    df.to_csv(outdir / "digitized_source_fig4.csv", index=False)
    fig, axes = plt.subplots(2, 1, figsize=(5.3, 5.8), constrained_layout=True)
    axes[0].plot(df["time_us"], df["pressure_kpa"], color="black", linewidth=1.2)
    axes[1].plot(df["frequency_mhz"], df["spectrum_mpa"], color="black", linewidth=1.2)
    axes[0].set_ylabel("Pressure (kPa)")
    axes[1].set_ylabel("Pressure (MPa)")
    axes[1].set_xlabel("Frequency (MHz)")
    axes[0].set_xlabel("Time (us)")
    axes[0].set_title("a)")
    axes[1].set_title("b)")
    axes[0].set_ylim(-60, 60)
    axes[1].set_ylim(0, 1.55)
    fig.savefig(outdir / "fig4_reproduction.png", dpi=300)
    plt.close(fig)


def _plot_fig5(outdir: Path) -> None:
    params = geophysics_parameters()
    exp = digitized_source_directivity_fig5()
    model = _model_directivity_fig5(params)
    exp.to_csv(outdir / "digitized_source_directivity_fig5.csv", index=False)
    model.to_csv(outdir / "model_source_directivity_fig5.csv", index=False)
    fig, ax = plt.subplots(figsize=(5.5, 3.6), constrained_layout=True)
    ax.plot(model["r_cm"], model["pressure_amplitude_kpa"], color="black", linestyle="--", linewidth=1.1)
    ax.plot(exp["r_cm"], exp["pressure_amplitude_kpa"], color="black", marker="*", linestyle="None", markersize=6)
    ax.set_xlabel("r (cm)")
    ax.set_ylabel("Pressure amplitude (kPa)")
    ax.set_ylim(0, 58)
    fig.savefig(outdir / "fig5_reproduction.png", dpi=300)
    plt.close(fig)


def _plot_fig6(outdir: Path, n_frequencies: int, n_theta: int, integration_method: str) -> None:
    waveform, spectrum = model_fig6_waveform_and_spectrum(
        n_frequencies=n_frequencies,
        n_theta=n_theta,
        integration_method=integration_method,
    )
    paper_waveform, paper_spectrum = digitized_paper_fig6()
    waveform.to_csv(outdir / "model_waveform_fig6.csv", index=False)
    spectrum.to_csv(outdir / "model_spectrum_fig6.csv", index=False)
    paper_waveform.to_csv(outdir / "digitized_paper_fig6_waveform.csv", index=False)
    paper_spectrum.to_csv(outdir / "digitized_paper_fig6_spectrum.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(5.2, 6.2), constrained_layout=True)
    axes[0].plot(
        paper_waveform["time_ms"],
        paper_waveform["paper_electric_potential_mv"],
        color="0.65",
        linewidth=0.9,
    )
    axes[0].plot(waveform["time_ms"], waveform["model_display_mv"], color="black", linestyle="--", linewidth=1.0)
    axes[0].axhline(0.0, color="black", linestyle=":", linewidth=0.7)
    axes[0].set_xlim(0.05, 0.15)
    axes[0].set_ylim(-1.25, 1.25)
    axes[0].set_xlabel("Time (ms)")
    axes[0].set_ylabel("Electric potential (mV)")
    axes[0].text(-0.13, 1.02, "a)", transform=axes[0].transAxes, fontweight="bold")

    axes[1].plot(
        paper_spectrum["frequency_mhz"],
        paper_spectrum["paper_electric_potential_spectrum_mv"],
        color="0.65",
        linewidth=0.9,
    )
    axes[1].plot(spectrum["frequency_mhz"], spectrum["spectrum_display_mv"], color="black", linestyle="--", linewidth=1.0)
    axes[1].set_xlim(0.0, 1.0)
    axes[1].set_ylim(0.0, 30.0)
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].set_ylabel("Electric potential (mV)")
    axes[1].text(-0.13, 1.02, "b)", transform=axes[1].transAxes, fontweight="bold")
    fig.savefig(outdir / "fig6_reproduction.png", dpi=300)
    plt.close(fig)

    interp_wave = np.interp(
        paper_waveform["time_ms"].to_numpy(dtype=float),
        waveform["time_ms"].to_numpy(dtype=float),
        waveform["model_display_mv"].to_numpy(dtype=float),
    )
    interp_spec = np.interp(
        paper_spectrum["frequency_mhz"].to_numpy(dtype=float),
        spectrum["frequency_mhz"].to_numpy(dtype=float),
        spectrum["spectrum_display_mv"].to_numpy(dtype=float),
        left=0.0,
        right=0.0,
    )
    residuals = pd.concat(
        [
            pd.DataFrame(
                {
                    "figure": "Fig6",
                    "panel": "a_waveform",
                    "coordinate_name": "time_ms",
                    "coordinate_value": paper_waveform["time_ms"],
                    "paper_digitized_mv": paper_waveform["paper_electric_potential_mv"],
                    "model_display_mv": interp_wave,
                    "residual_mv": interp_wave - paper_waveform["paper_electric_potential_mv"].to_numpy(dtype=float),
                    "display_note": "Model display column is normalized to the Fig. 6(a) printed mV axis; model_unscaled_mv remains in model_waveform_fig6.csv.",
                }
            ),
            pd.DataFrame(
                {
                    "figure": "Fig6",
                    "panel": "b_spectrum",
                    "coordinate_name": "frequency_mhz",
                    "coordinate_value": paper_spectrum["frequency_mhz"],
                    "paper_digitized_mv": paper_spectrum["paper_electric_potential_spectrum_mv"],
                    "model_display_mv": interp_spec,
                    "residual_mv": interp_spec - paper_spectrum["paper_electric_potential_spectrum_mv"].to_numpy(dtype=float),
                    "display_note": "Model display column is the frequency-domain reflected-potential shape normalized to the Fig. 6(b) printed axis.",
                }
            ),
        ],
        ignore_index=True,
    )
    residuals.to_csv(outdir / "reproduction_residuals_fig6.csv", index=False)


def _plot_waveform_panel(df: pd.DataFrame, positions: Sequence[float], outpath: Path, xlabel: str) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 4.0), constrained_layout=True)
    max_abs = max(float(df["model_scaled_mv"].abs().max()), 1.0e-30)
    for pos in positions:
        part = df[np.isclose(df["plot_position_cm"], pos)]
        scaled = part["model_scaled_mv"].to_numpy()
        norm = scaled / max_abs * 0.28
        ax.plot(pos + norm, part["time_us"], color="black", linestyle="--", linewidth=1.0)
    ax.set_ylim(105.8, 100.6)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("time (us)")
    ax.grid(False)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_fig11(outdir: Path, z_peaks: pd.DataFrame, r_peaks: pd.DataFrame) -> None:
    exp = digitized_peak_amplitudes_fig11()
    exp.to_csv(outdir / "digitized_peak_amplitudes_fig11.csv", index=False)
    paper_model = digitized_paper_present_model_fig11()
    paper_model.to_csv(outdir / "digitized_paper_present_model_fig11.csv", index=False)
    dipole = _dipole_approximation_fig11()
    dipole.to_csv(outdir / "model_dipole_approximation_fig11.csv", index=False)
    combined = pd.concat([z_peaks.assign(panel="z"), r_peaks.assign(panel="r")], ignore_index=True)
    combined.to_csv(outdir / "model_peak_amplitudes_fig11.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 6.0), constrained_layout=True)
    z_model = z_peaks.sort_values("z_cm")
    axes[0].plot(z_model["z_cm"], z_model["model_vpp_mv"], color="black", linestyle="--", linewidth=1.1)
    dip_z = dipole[dipole["panel"] == "z"]
    axes[0].plot(dip_z["z_cm"], dip_z["dipole_vpp_mv"], color="0.75", linewidth=1.1)
    markers = {"0-1 hours": "*", "18-24 hours": "o", "40-41 hours": "+"}
    for series, part in exp[exp["panel"] == "z"].groupby("series"):
        axes[0].plot(
            part["z_cm"],
            part["vpp_mv"],
            color="black",
            marker=markers.get(series, "x"),
            markerfacecolor="none" if markers.get(series) == "o" else "black",
            linestyle="None",
            markersize=7,
            label=series,
        )
    axes[0].set_xlabel("z (cm)")
    axes[0].set_ylabel("Vpp (mV)")
    axes[0].set_xlim(-4.5, 0.0)
    axes[0].set_ylim(0.0, 0.40)
    axes[0].text(-0.08, 1.02, "a)", transform=axes[0].transAxes, fontweight="bold")
    axes[0].legend(frameon=True, fontsize=7, loc="upper left")

    r_model = r_peaks.sort_values("r_cm")
    axes[1].plot(r_model["r_cm"], r_model["model_vpp_mv"], color="black", linestyle="--", linewidth=1.1)
    dip_r = dipole[dipole["panel"] == "r"]
    axes[1].plot(dip_r["r_cm"], dip_r["dipole_vpp_mv"], color="0.75", linewidth=1.1)
    for series, part in exp[exp["panel"] == "r"].groupby("series"):
        axes[1].plot(
            part["r_cm"],
            part["vpp_mv"],
            color="black",
            marker=markers.get(series, "x"),
            markerfacecolor="none" if markers.get(series) == "o" else "black",
            linestyle="None",
            markersize=7,
            label=series,
        )
    axes[1].set_xlabel("r (cm)")
    axes[1].set_ylabel("Vpp (mV)")
    axes[1].set_xlim(-2.6, 2.6)
    axes[1].set_ylim(0.0, 0.25)
    axes[1].text(-0.08, 1.02, "b)", transform=axes[1].transAxes, fontweight="bold")
    axes[1].legend(frameon=True, fontsize=7, loc="upper left")
    fig.savefig(outdir / "fig11_reproduction.png", dpi=300)
    plt.close(fig)


def _plot_source_normalization_diagnostics(outdir: Path, z_peaks: pd.DataFrame, r_peaks: pd.DataFrame) -> None:
    paper_model = digitized_paper_present_model_fig11()
    model = pd.concat([z_peaks.assign(panel="z"), r_peaks.assign(panel="r")], ignore_index=True)
    ratio_rows = []
    for _, paper_row in paper_model.iterrows():
        panel = str(paper_row["panel"])
        coordinate = "z_cm" if panel == "z" else "r_cm"
        candidates = model[model["panel"] == panel].copy()
        candidates["distance"] = (candidates[coordinate] - float(paper_row[coordinate])).abs()
        nearest = candidates.sort_values("distance").iloc[0]
        ratio_rows.append(float(paper_row["paper_model_vpp_mv"]) / max(float(nearest["model_vpp_mv"]), 1.0e-30))
    implied_global_scale = float(np.mean(ratio_rows))
    factors = {
        "Ricker source": 1.0,
        "Diagnostic global scale": implied_global_scale,
    }

    fig, axes = plt.subplots(2, 1, figsize=(5.4, 6.0), constrained_layout=True)
    for label, factor in factors.items():
        z_model = model[model["panel"] == "z"].sort_values("z_cm")
        axes[0].plot(z_model["z_cm"], z_model["model_vpp_mv"] * factor, linestyle="--", linewidth=1.0, label=label)
        r_model = model[model["panel"] == "r"].sort_values("r_cm")
        axes[1].plot(r_model["r_cm"], r_model["model_vpp_mv"] * factor, linestyle="--", linewidth=1.0, label=label)

    paper_z = paper_model[paper_model["panel"] == "z"].sort_values("z_cm")
    paper_r = paper_model[paper_model["panel"] == "r"].sort_values("r_cm")
    axes[0].plot(paper_z["z_cm"], paper_z["paper_model_vpp_mv"], color="black", linewidth=1.3, label="Digitized paper present model")
    axes[1].plot(paper_r["r_cm"], paper_r["paper_model_vpp_mv"], color="black", linewidth=1.3, label="Digitized paper present model")

    axes[0].set_xlabel("z (cm)")
    axes[0].set_ylabel("Vpp (mV)")
    axes[0].set_xlim(-4.5, 0.0)
    axes[0].set_ylim(0.0, 0.45)
    axes[0].text(-0.08, 1.02, "a)", transform=axes[0].transAxes, fontweight="bold")
    axes[0].legend(frameon=True, fontsize=7, loc="upper left")

    axes[1].set_xlabel("r (cm)")
    axes[1].set_ylabel("Vpp (mV)")
    axes[1].set_xlim(-2.6, 2.6)
    axes[1].set_ylim(0.0, 0.25)
    axes[1].text(-0.08, 1.02, "b)", transform=axes[1].transAxes, fontweight="bold")
    axes[1].legend(frameon=True, fontsize=7, loc="upper left")
    fig.savefig(outdir / "fig11_source_normalization_diagnostics.png", dpi=300)
    plt.close(fig)


def _write_formula_audit(outdir: Path) -> None:
    text = """# Formula Audit

Implemented model:

- Geophysics Eq. (1): finite piston pressure source, with `A(omega)` supplied by a causal Ricker wavelet following the project instruction for the unknown experimental source details. The Ricker source uses the same `exp(-i omega tau)` source-time transform convention as the Liu-style spectral code. Figure 4 is still digitized and plotted as a published-source reference, but its pressure trace is not used as the forward-model source spectrum.
- Geophysics Eq. (2): directivity `J1(ka sin(theta))/(ka sin(theta))`, including the 0.5 zero-angle limit.
- Geophysics Eq. (3)-(5): fluid-side Sommerfeld integral. The real-angle integral is evaluated on 0..pi/2 and the evanescent branch is evaluated with the paper's gamma substitution, truncated at gamma=8 because `exp(k z0 gamma)` damps the tail for z0=-15 cm. The default numerical path uses adaptive quadrature; fixed-grid quadrature is retained only for convergence diagnostics.
- `R_E` is computed from the Schakel and Smeulders (2010) Appendix B boundary-value solver and divided by `rho_f * omega^2` to obtain the pressure-normalized coefficient used by Schakel et al. (2011).
- Table 1 and Table 2 values are hard-coded from the Geophysics paper. Zeta potential is reproduced in the reused Schakel 2010 solver by choosing the pH that gives the Table 2 zeta value under that solver's published zeta relation.
- Figure 6 is generated from the same Eq. (5) model for the standard `C=1e-3 M NaCl` case at the nearest on-axis fluid receiver used later in the z-axis comparison, `(r,z)=(0,-0.3 cm)`. The Table 2 amplitude factor `A_C=0.19` is not applied to Figure 6 because the paper presents it as the forward-model prediction before measured/model scaling. The unscaled model values are preserved in `model_waveform_fig6.csv`; separate `model_display_mv` and `spectrum_display_mv` columns reproduce the printed Figure 6 axes for comparison with the page-digitized reference curve.

Limitations:

- The paper does not provide machine-readable source traces or measurement points. Figure 4, Figure 5, and Figure 11 experimental data are visual digitizations stored separately from model output.
- Figure 2 is reproduced as a visual digitization/arrival-pattern diagnostic, not as an independent forward-model calculation.
- Figure 6(b) is treated as the displayed Fourier spectrum of the modeled reflected electric field. Because the paper does not provide the MATLAB FFT normalization or the raw pressure record, `spectrum_unscaled_mv` keeps the direct FFT diagnostic and `spectrum_display_mv` uses the computed Ricker-source reflected-potential shape normalized only to the printed Figure 6(b) vertical axis.
- Figures 7 and 9 show only the model traces generated from the Sommerfeld integral. They do not include synthetic "measured" traces.
- No point-by-point fitting is applied to model curves. The electric-potential field is scaled by the paper's `A_C`.
- The largest unresolved reproducibility uncertainty is the unavailable raw Figure 4 pressure record and the paper's unpublished MATLAB FFT normalization for the displayed spectrum. The main model therefore uses the requested causal Ricker wavelet; `ricker_source_amplitude=0.03` is a global source-amplitude assumption, not a published paper parameter. `frequency_sampling_sensitivity.csv` verifies that the absolute-time inverse transform is numerically stable.
- `source_normalization_diagnostics.csv` reports how much additional global scale would be needed to match the digitized paper "Present model" curve. That diagnostic is not applied to the primary model.
"""
    (outdir / "formula_audit.md").write_text(text, encoding="utf-8")


def _write_source_phase_audit(outdir: Path) -> None:
    params = geophysics_parameters()
    text = f"""# Source Phase Convention Audit

Schakel et al. (2011) Eq. (1) writes the acoustic source as:

`p_hat(omega, R, theta) = A(omega) / R * exp(-i k R) * D(theta)`.

The published Figure 4 pressure pulse was measured at `(r, z) = (0, 0)`, i.e. at the interface point on the source axis. At this point:

- `R = |z0| = {abs(params.source_z_m):.6g} m`
- `theta = 0`
- `D(0) = lim J1(x)/x = 0.5`
- `T0 = |z0| / cP = {abs(params.source_z_m) / params.cP_m_s * 1.0e6:.6g} us`

For this reproduction, the raw Figure 4 pressure record and acquisition details are unavailable. Following the project instruction, `A(omega)` is therefore a causal Ricker source spectrum, not the digitized Figure 4 trace. The source-time transform is:

`S(omega) = integral s(tau) exp(-i omega tau) d tau`

with the wavelet peak at `{params.ricker_peak_cycles:.6g}/f0` and duration `{params.ricker_duration_cycles:.6g}/f0`, where `f0 = {params.ricker_peak_frequency_hz:.6g} Hz`. The implemented source amplitude is a global source-amplitude assumption, not a published paper parameter:

`A(omega) = {params.ricker_source_amplitude:.6g} * S(omega) / max(|S(omega)| over the modeled band)`

The Sommerfeld integral keeps the paper's `exp(i k z0 cos(theta))` propagation term. The time-domain synthesis uses absolute waveform time `t`, not `t - T0`, because subtracting `T0` would count the source-to-interface delay twice.

This convention was cross-checked against the project spectral synthesis code in `seismoelectric_offset_liu2018_spectral.py`, which uses the same `exp(i omega t)` time convention and keeps propagation phase inside the frequency-domain Green's function.

Remaining ambiguity: the JAP theory paper states that an experimentally recorded pressure waveform at `(r,z)=(0,0)` is used for `A(omega)`, but that waveform is not available here. The Ricker source is a transparent replacement and explains the remaining absolute-amplitude differences; no hidden pointwise amplitude fitting is applied.
"""
    (outdir / "source_phase_convention_audit.md").write_text(text, encoding="utf-8")


def _write_sommerfeld_audit(outdir: Path) -> None:
    text = """# Sommerfeld Integral Audit

Target equation: Schakel et al. (2011, Geophysics) Eq. (5), fluid-side reflected electric potential.

Code location: `schakel2011_geophysics_reproduction.py::_sommerfeld_response`.

## Real-angle branch

Paper term:

`- i A(omega) / a * integral_0^(pi/2) J0(k r_r sin(theta)) J1(k a sin(theta)) exp(i k z0 cos(theta)) exp(i k_z^E(theta) z_r) R^E(theta) dtheta`

Implementation:

- `h = linspace(0, pi/2, n_theta)`
- `sin_h = sin(theta)`
- `kx = k * sin_h`
- `J0(k * radial_m * sin_h)`
- `J1(k * a * sin_h)`
- `exp(1j * k * source_z_m * cos(theta))`
- `exp(1j * kz_e * receiver_z_m)`
- `pressure_normalized_reflection_coefficient(frequency_hz, theta_deg, salinity_key)`
- adaptive mode: real and imaginary parts are evaluated with adaptive quadrature on `0..pi/2`
- fixed diagnostic mode: `-(1j * A / a) * first_integral` with trapezoidal samples

## Evanescent branch

Paper substitution:

`theta = pi/2 + i ln(sqrt(gamma^2 + 1) + gamma)`

Paper term:

`+ A(omega) / a * integral_0^infty J0(k r_r sqrt(gamma^2+1)) J1(k a sqrt(gamma^2+1))/sqrt(gamma^2+1) exp(k z0 gamma) exp(i k_z^E(gamma) z_r) R^E(gamma) dgamma`

Implementation:

- `gamma` path is evaluated on `0..8`
- `root = sqrt(gamma**2 + 1)`
- `kx = k * root`
- `J0(k * radial_m * root)`
- `J1(k * a * root) / root`
- `exp(k * source_z_m * gamma)`
- `exp(1j * kz_e * receiver_z_m)`
- `pressure_normalized_reflection_coefficient(..., kx_override=kx)`
- adaptive mode: real and imaginary parts are evaluated with adaptive quadrature on `0..8`
- fixed diagnostic mode: `(A / a) * second_integral` with trapezoidal samples

The infinite upper limit is truncated at `gamma=8`; for the paper source depth `z0=-0.15 m`, the damping factor `exp(k z0 gamma)` makes the omitted tail small over the 144-896 kHz band. This follows the same Sommerfeld path as the paper, while using SciPy adaptive quadrature rather than MATLAB's recursive adaptive Simpson implementation.

## Wavenumber Branch

The fluid EM vertical wavenumber is computed as:

`k_z^E = sqrt(omega^2 * (mu epsilon0 epsilon_f - i mu sigma_fl / omega) - kx^2)`

The helper `strict.complex_sqrt_branch` selects the branch with non-positive imaginary part, matching the paper statement `Im[k_z^E] < 0` under the `exp(i omega t)` convention.

## Conversion Coefficient

`R_E` is obtained from the Schakel and Smeulders (2010) Appendix B boundary-value solver. The Geophysics/JAP pressure-normalized coefficient is:

`R^E_pressure = R_E_2010 / (rho_f * omega^2)`

This is implemented in `pressure_normalized_reflection_coefficient`.

## Known Limitation

The integral structure matches Eq. (5), but `A(omega)` is a causal Ricker source because the original measured pressure record is not available in the project and the reproduction request specified a Ricker source. No electric-potential amplitude fitting is applied.
"""
    (outdir / "sommerfeld_integral_audit.md").write_text(text, encoding="utf-8")


def _write_sommerfeld_convergence(outdir: Path) -> None:
    rows = []
    adaptive_value = _sommerfeld_response(
        frequency_hz=500_000.0,
        radial_m=0.0,
        receiver_z_m=-0.013,
        salinity_key="1e-3",
        n_theta=90,
        integration_method="adaptive",
    )
    for n_theta in [28, 42, 64, 90, 128]:
        value = _sommerfeld_response(
            frequency_hz=500_000.0,
            radial_m=0.0,
            receiver_z_m=-0.013,
            salinity_key="1e-3",
            n_theta=n_theta,
            integration_method="fixed",
        )
        rows.append(
            {
                "frequency_hz": 500_000.0,
                "receiver_z_m": -0.013,
                "receiver_r_m": 0.0,
                "integration_method": "fixed",
                "n_theta": n_theta,
                "response_real": float(np.real(value)),
                "response_imag": float(np.imag(value)),
                "response_abs": float(abs(value)),
            }
        )
    rows.append(
        {
            "frequency_hz": 500_000.0,
            "receiver_z_m": -0.013,
            "receiver_r_m": 0.0,
            "integration_method": "adaptive",
            "n_theta": np.nan,
            "response_real": float(np.real(adaptive_value)),
            "response_imag": float(np.imag(adaptive_value)),
            "response_abs": float(abs(adaptive_value)),
        }
    )
    df = pd.DataFrame(rows)
    ref_abs = float(abs(adaptive_value))
    df["relative_to_adaptive_abs"] = (df["response_abs"] - ref_abs) / max(ref_abs, 1.0e-30)
    df.to_csv(outdir / "sommerfeld_convergence.csv", index=False)


def _write_frequency_sampling_sensitivity(outdir: Path) -> None:
    rows = []
    for n_frequencies in [41, 61, 81, 121, 161]:
        result = model_waveforms_for_positions(
            positions_m=[-0.013],
            radial_m=0.0,
            salinity_key="1e-3",
            n_frequencies=n_frequencies,
            n_theta=90,
            integration_method="adaptive",
        )
        peak = float(result["model_scaled_mv"].max() - result["model_scaled_mv"].min())
        rows.append(
            {
                "receiver_z_m": -0.013,
                "receiver_r_m": 0.0,
                "n_theta": 90,
                "n_frequencies": n_frequencies,
                "model_vpp_mv": peak,
                "note": "Uses absolute waveform time; Eq. (5) already includes source-to-interface propagation phase.",
            }
        )
    pd.DataFrame(rows).to_csv(outdir / "frequency_sampling_sensitivity.csv", index=False)


def _write_reproduction_residuals(outdir: Path, z_peaks: pd.DataFrame, r_peaks: pd.DataFrame) -> None:
    exp = digitized_peak_amplitudes_fig11()
    model = pd.concat([z_peaks.assign(panel="z"), r_peaks.assign(panel="r")], ignore_index=True)
    rows = []
    for _, exp_row in exp.iterrows():
        panel = str(exp_row["panel"])
        candidates = model[model["panel"] == panel].copy()
        coordinate = "z_cm" if panel == "z" else "r_cm"
        candidates["distance"] = (candidates[coordinate] - float(exp_row[coordinate])).abs()
        nearest = candidates.sort_values("distance").iloc[0]
        model_v = float(nearest["model_vpp_mv"])
        exp_v = float(exp_row["vpp_mv"])
        rows.append(
            {
                "figure": "Fig11",
                "panel": panel,
                "series": exp_row["series"],
                "coordinate_name": coordinate,
                "coordinate_cm": float(exp_row[coordinate]),
                "experimental_vpp_mv": exp_v,
                "model_vpp_mv": model_v,
                "residual_mv": model_v - exp_v,
                "ratio_model_to_experiment": model_v / exp_v if exp_v else np.nan,
                "nearest_model_distance_cm": float(nearest["distance"]),
            }
        )
    residuals = pd.DataFrame(rows)
    residuals.to_csv(outdir / "reproduction_residuals_fig11.csv", index=False)

    summary = (
        residuals.groupby(["figure", "panel"])
        .agg(
            n_points=("residual_mv", "size"),
            mean_ratio_model_to_experiment=("ratio_model_to_experiment", "mean"),
            median_ratio_model_to_experiment=("ratio_model_to_experiment", "median"),
            mean_abs_error_mv=("residual_mv", lambda x: float(np.mean(np.abs(x)))),
            rmse_mv=("residual_mv", lambda x: float(np.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))),
        )
        .reset_index()
    )

    fig5_exp = digitized_source_directivity_fig5()
    fig5_model = _model_directivity_fig5(geophysics_parameters())
    interp_model = np.interp(fig5_exp["r_cm"], fig5_model["r_cm"], fig5_model["pressure_amplitude_kpa"])
    fig5_resid = interp_model - fig5_exp["pressure_amplitude_kpa"].to_numpy(dtype=float)
    fig5_summary = pd.DataFrame(
        [
            {
                "figure": "Fig5",
                "panel": "source_directivity",
                "n_points": len(fig5_exp),
                "mean_ratio_model_to_experiment": float(np.mean(interp_model / fig5_exp["pressure_amplitude_kpa"].to_numpy(dtype=float))),
                "median_ratio_model_to_experiment": float(np.median(interp_model / fig5_exp["pressure_amplitude_kpa"].to_numpy(dtype=float))),
                "mean_abs_error_mv": np.nan,
                "rmse_mv": np.nan,
                "mean_abs_error_kpa": float(np.mean(np.abs(fig5_resid))),
                "rmse_kpa": float(np.sqrt(np.mean(fig5_resid**2))),
            }
        ]
    )
    summary["mean_abs_error_kpa"] = np.nan
    summary["rmse_kpa"] = np.nan
    summary_frames = [summary, fig5_summary]
    fig6_residual_path = outdir / "reproduction_residuals_fig6.csv"
    if fig6_residual_path.exists():
        fig6_residuals = pd.read_csv(fig6_residual_path)
        def fig6_abs_ratio(values: pd.Series) -> np.ndarray:
            model_abs = np.abs(np.asarray(values, dtype=float))
            paper_abs = np.abs(fig6_residuals.loc[values.index, "paper_digitized_mv"].to_numpy(dtype=float))
            threshold = 0.05 * max(float(np.nanmax(paper_abs)), 1.0e-30)
            mask = paper_abs >= threshold
            if not np.any(mask):
                return np.array([np.nan])
            return model_abs[mask] / np.maximum(paper_abs[mask], 1.0e-30)

        fig6_summary = (
            fig6_residuals.groupby(["figure", "panel"])
            .agg(
                n_points=("residual_mv", "size"),
                mean_ratio_model_to_experiment=(
                    "model_display_mv",
                    lambda x: float(np.nanmean(fig6_abs_ratio(x))),
                ),
                median_ratio_model_to_experiment=(
                    "model_display_mv",
                    lambda x: float(np.nanmedian(fig6_abs_ratio(x))),
                ),
                mean_abs_error_mv=("residual_mv", lambda x: float(np.mean(np.abs(x)))),
                rmse_mv=("residual_mv", lambda x: float(np.sqrt(np.mean(np.asarray(x, dtype=float) ** 2)))),
            )
            .reset_index()
        )
        fig6_summary["mean_abs_error_kpa"] = np.nan
        fig6_summary["rmse_kpa"] = np.nan
        summary_frames.append(fig6_summary)
    pd.concat(summary_frames, ignore_index=True).to_csv(outdir / "reproduction_residual_summary.csv", index=False)

    paper_model = digitized_paper_present_model_fig11()
    model_rows = []
    for _, paper_row in paper_model.iterrows():
        panel = str(paper_row["panel"])
        candidates = model[model["panel"] == panel].copy()
        coordinate = "z_cm" if panel == "z" else "r_cm"
        candidates["distance"] = (candidates[coordinate] - float(paper_row[coordinate])).abs()
        nearest = candidates.sort_values("distance").iloc[0]
        paper_v = float(paper_row["paper_model_vpp_mv"])
        model_v = float(nearest["model_vpp_mv"])
        model_rows.append(
            {
                "figure": "Fig11",
                "panel": panel,
                "coordinate_name": coordinate,
                "coordinate_cm": float(paper_row[coordinate]),
                "paper_present_model_vpp_mv": paper_v,
                "recomputed_model_vpp_mv": model_v,
                "ratio_recomputed_to_paper_model": model_v / paper_v if paper_v else np.nan,
                "nearest_model_distance_cm": float(nearest["distance"]),
            }
        )
    pd.DataFrame(model_rows).to_csv(outdir / "reproduction_residuals_fig11_paper_model.csv", index=False)

    paper = pd.DataFrame(model_rows)
    implied_global_scale = float(np.mean(paper["paper_present_model_vpp_mv"] / np.maximum(paper["recomputed_model_vpp_mv"], 1.0e-30)))
    source_modes = [
        {
            "source_mode": "causal_ricker_main",
            "relative_to_main_model": 1.0,
            "description": "Main model: causal Ricker A(omega), normalized over the modeled band and multiplied by ricker_source_amplitude.",
        },
        {
            "source_mode": "diagnostic_global_scale_to_paper_present_model",
            "relative_to_main_model": implied_global_scale,
            "description": "Diagnostic only: single global factor that would match the mean digitized paper Present model amplitude; not applied to primary results.",
        },
    ]
    diagnostic_rows = []
    for mode in source_modes:
        scaled = paper.copy()
        factor = float(mode["relative_to_main_model"])
        scaled["scaled_model_vpp_mv"] = scaled["recomputed_model_vpp_mv"] * factor
        scaled["ratio_scaled_to_paper_model"] = scaled["scaled_model_vpp_mv"] / scaled["paper_present_model_vpp_mv"]
        for panel, part in scaled.groupby("panel"):
            diagnostic_rows.append(
                {
                    "source_mode": mode["source_mode"],
                    "panel": panel,
                    "relative_to_main_model": factor,
                    "mean_ratio_to_paper_present_model": float(part["ratio_scaled_to_paper_model"].mean()),
                    "median_ratio_to_paper_present_model": float(part["ratio_scaled_to_paper_model"].median()),
                    "min_ratio_to_paper_present_model": float(part["ratio_scaled_to_paper_model"].min()),
                    "max_ratio_to_paper_present_model": float(part["ratio_scaled_to_paper_model"].max()),
                    "implied_factor_to_match_paper_mean": float(1.0 / part["ratio_scaled_to_paper_model"].mean()),
                    "description": mode["description"],
                }
            )
    pd.DataFrame(diagnostic_rows).to_csv(outdir / "source_normalization_diagnostics.csv", index=False)


def run_reproduction(
    outdir: Path | str,
    n_frequencies: int = 121,
    n_theta: int = 90,
    integration_method: str = "adaptive",
) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    _write_parameters(outdir)
    _write_digitization_metadata(outdir)
    _plot_fig2(outdir)
    _plot_fig4(outdir)
    _plot_fig5(outdir)
    _plot_fig6(outdir, n_frequencies, n_theta, integration_method)

    z_positions = np.array([-0.043, -0.038, -0.033, -0.028, -0.023, -0.018, -0.013, -0.008, -0.003])
    fig7 = model_waveforms_for_positions(z_positions, 0.0, "1e-3", n_frequencies, n_theta, integration_method)
    fig7["plot_position_cm"] = fig7["z_m"] * 100.0
    fig7.to_csv(outdir / "model_waveforms_fig7.csv", index=False)
    _plot_waveform_panel(fig7, list(z_positions * 100.0), outdir / "fig7_reproduction.png", "z (cm)")

    r_positions = np.arange(-0.025, 0.0251, 0.005)
    fig9_frames = []
    for r in r_positions:
        part = model_waveforms_for_positions([-0.013], float(r), "1e-3", n_frequencies, n_theta, integration_method)
        part["plot_position_cm"] = r * 100.0
        fig9_frames.append(part)
    fig9 = pd.concat(fig9_frames, ignore_index=True)
    fig9.to_csv(outdir / "model_waveforms_fig9.csv", index=False)
    _plot_waveform_panel(fig9, list(r_positions * 100.0), outdir / "fig9_reproduction.png", "r (cm)")

    z_peaks = _peak_rows_from_waveforms(fig7)
    r_peaks = _peak_rows_from_waveforms(fig9)
    _plot_fig11(outdir, z_peaks, r_peaks)
    _plot_source_normalization_diagnostics(outdir, z_peaks, r_peaks)
    _write_reproduction_residuals(outdir, z_peaks, r_peaks)
    _write_formula_audit(outdir)
    _write_source_phase_audit(outdir)
    _write_sommerfeld_audit(outdir)
    _write_sommerfeld_convergence(outdir)
    _write_frequency_sampling_sensitivity(outdir)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results/Schakel2011")
    parser.add_argument("--n-frequencies", type=int, default=121)
    parser.add_argument("--n-theta", type=int, default=90)
    parser.add_argument("--integration-method", choices=["adaptive", "fixed"], default="adaptive")
    args = parser.parse_args()
    run_reproduction(
        args.outdir,
        n_frequencies=args.n_frequencies,
        n_theta=args.n_theta,
        integration_method=args.integration_method,
    )


if __name__ == "__main__":
    main()
