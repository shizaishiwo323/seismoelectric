#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reactive-transport seismoelectric model with a Schakel 2011 waveform core.

This script intentionally creates a new workflow rather than modifying
``seismoelectric_offset_liu2018_spectral.py``.  It keeps that script's
reactive-transport parameter mapping, receiver geometry, plotting style, and
result filenames, but forces the finite offset to zero and replaces the Liu
frequency-wavenumber waveform synthesis with the Schakel et al. (2011)
Sommerfeld integral used in ``schakel2011_geophysics_reproduction.py``.

Scope note
----------
Fluid-side receivers use the Schakel et al. (2011) Eq. (5) reflected
electric-potential response.  Porous-side receivers use the front-interface
transmitted terms from the JAP Eq. (8) structure as implemented in
``schakel2011_jap_fig2_reproduction.py``: the transmitted TM electric-potential
term and the Pf coseismic electric-potential term.  Back-interface/sample-width
reflections from the laboratory slab are not included because the present RT-SE
geometry is a single interface.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seismoelectric_offset_liu2018_spectral as base

try:
    from scipy.special import j0 as scipy_j0
    from scipy.special import j1 as scipy_j1
except Exception:  # pragma: no cover - only used if scipy is unavailable.
    scipy_j0 = None
    scipy_j1 = None


@dataclass
class ZeroOffsetSchakelConfig:
    """Configuration for the zero-offset Schakel 2011 RT-SE simulator.

    All default parameters are declared in this file.  We still reuse numerical
    helper functions from ``seismoelectric_offset_liu2018_spectral.py``, but the
    source/receiver geometry, source wavelet, material properties, and numerical
    controls below are not inherited from that module's ``SEConfig``.

    - ``z_s``: positive source-to-interface distance in the upper fluid.  The
      Sommerfeld integral uses the signed source coordinate ``source_z_m=-z_s``,
      so the acoustic transmitter is above the interface on the fluid side.
    - ``receiver_z_min`` and ``receiver_z_max``: lower/upper receiver-electrode
      positions relative to the interface.  Negative values are fluid-side
      electrodes; positive values are porous-side electrodes; ``z=0`` is the
      interface line.
    - ``receiver_spacing``: receiver-electrode interval.  The number of
      receivers is not typed manually; it is computed from
      ``receiver_z_min:receiver_spacing:receiver_z_max`` in
      ``synthesize_waveforms_schakel2011``.
    - ``offset_D``: horizontal offset between the source axis and receiver line.
      This zero-offset Schakel workflow forces it to ``0.0``.
    """

    # Pore-scale RTM geometry retained for parameter bookkeeping and plotting.
    Lx: float = 600e-6
    Ly: float = 400e-6
    grain_radius: float = 30e-6
    pore_throat_width: float = 10e-6

    # Source and receiver geometry used by the waveform gather.
    # z_s is the positive source-interface distance in the upper fluid; the
    # Sommerfeld integral converts it to the signed source coordinate -z_s.
    z_s: float = 15.0e-3

    # Receiver/electrode line positions relative to the interface.  Negative
    # values are fluid-side receivers, positive values are porous-side receivers.
    # Receiver count is computed later from min/max/spacing, not set separately.
    receiver_z_min: float = -20.0e-3
    receiver_z_max: float = 20.0e-3
    receiver_spacing: float = 1.0e-3

    # Horizontal source-receiver offset.  This script is the zero-offset
    # Schakel version, so even if CLI accepts an offset for compatibility,
    # the value is reset to 0 m before waveform synthesis.
    offset_D: float = 0.0

    # Acoustic source and waveform sampling.
    f0: float = 500.0e3
    source_pressure_amp: float = 1.0
    source_kb_m_inv: float = 80.0
    source_beam_theta_deg: float = 0.0
    source_peak_cycles: float = 2.0
    source_duration_cycles: float = 8.0
    spectral_f_min_factor: float = 0.25
    spectral_f_max_factor: float = 2.5
    spectral_n_omega: int = 48
    spectral_n_k: int = 401
    spectral_k_limit_factor: float = 0.98
    # Retained for parameter-file compatibility; the Schakel gather now starts
    # at 0 s by construction and marks T0 explicitly in the plot.
    waveform_t_before: float = 0.0
    waveform_t_after: float = 18.0e-6
    waveform_nt: int = 1200

    # Plane-wave coefficient diagnostic angle.
    coeff_theta_deg: float = 45.0

    # Baseline porous frame / mineral / fluid parameters.
    K_b: float = 5.8e9
    G: float = 3.4e9
    K_s: float = 76.8e9
    K_f: float = 2.22e9
    K_fl: float = 2.22e9
    eta: float = 1.0e-3
    rho_f: float = 1000.0
    rho_fl: float = 1000.0
    rho_s: float = 2710.0
    eps_f: float = 80.0
    eps_s: float = 8.0
    eps_fl: float = 80.0
    sigma_fl_default: float = 5.0e-3

    # Electrochemistry and reactive-transport mapping controls.
    temperature: float = 295.0
    C_background_molL: float = 1.0e-3
    H_min_molL: float = 1.0e-7
    outlet_h_unit: str = "mol_cm3"
    upper_fluid_conductivity_mode: str = "constant"
    M_similarity: float = 1.0

    # Numerical stability / model-validity bounds.
    phi_min: float = 0.05
    phi_max_valid: float = 0.95
    k0_min: float = 1.0e-18
    eps_complex: float = 1.0e-30

    # Physical constants.
    eps0: float = 8.854187817e-12
    mu0: float = 4.0 * math.pi * 1e-7
    e_charge: float = 1.602176634e-19
    k_B: float = 1.380649e-23
    N_A: float = 6.02214076e23
    z1: int = 1
    z2: int = -1
    b1: float = 3.246e11
    b2: float = 4.931e11

    # Schakel et al. (2011) pressure-source/Sommerfeld settings.
    # Physical radius of the circular piston transducer.  This affects source
    # directivity through J1(k a sin(theta)) in the Sommerfeld integral, but it
    # does not change the source center location; that location is controlled by
    # z_s above.
    transducer_radius_m: float = 1.125 * 0.0254 / 2.0
    schakel_bandpass_low_hz: float = 1.0
    schakel_bandpass_high_hz: float = 1.25e6
    schakel_bandpass_taper_hz: float = 40_000.0
    schakel_source_taper_us: float = 0.35
    schakel_source_mode: str = "causal_ricker"  # "fig4_digitized" is a legacy key for the Fig. 4 visual source approximation.
    include_porous_pf_coseismic: bool = False
    schakel_gamma_max: float = 8.0
    schakel_include_evanescent: bool = True


def _j0(x):
    if scipy_j0 is not None:
        return scipy_j0(x)
    x = np.asarray(x, dtype=complex)
    theta = np.linspace(0.0, 2.0 * math.pi, 512, endpoint=False)
    return np.mean(np.exp(1j * x[..., None] * np.cos(theta)), axis=-1)


def _j1(x):
    if scipy_j1 is not None:
        return scipy_j1(x)
    x = np.asarray(x, dtype=complex)
    theta = np.linspace(0.0, math.pi, 512)
    return np.trapz(np.cos(theta - x[..., None] * np.sin(theta)), theta, axis=-1) / math.pi


def _trapz(values: np.ndarray, x: np.ndarray, axis: int = -1) -> np.ndarray:
    integrate = getattr(np, "trapezoid", np.trapz)
    return integrate(values, x, axis=axis)


def _frequency_grid(cfg: ZeroOffsetSchakelConfig, n_frequencies: int) -> np.ndarray:
    if int(n_frequencies) < 2:
        raise ValueError("n_frequencies must be >= 2")
    return np.linspace(cfg.schakel_bandpass_low_hz, cfg.schakel_bandpass_high_hz, int(n_frequencies))


@lru_cache(maxsize=1)
def _source_time_table() -> tuple[np.ndarray, np.ndarray]:
    """Approximate the causal Schakel 2011 Fig. 4 pressure pulse.

    The original pressure record is not distributed with the paper.  This uses
    the same visual Fig. 4 approximation as ``schakel2011_geophysics_reproduction``.
    """
    time_us = np.linspace(0.0, 3.6, 4096)
    carrier = np.sin(2.0 * math.pi * 0.56 * (time_us - 0.75))
    envelope = np.exp(-0.5 * ((time_us - 1.78) / 0.78) ** 2)
    pressure_kpa = 58.0 * envelope * carrier
    pressure_kpa += 10.0 * np.exp(-0.5 * ((time_us - 3.05) / 0.28) ** 2)

    t = time_us * 1.0e-6
    p = pressure_kpa * 1.0e3
    return t, p - np.mean(p)


def schakel_source_A_spectrum(frequency_hz: np.ndarray, cfg: ZeroOffsetSchakelConfig) -> np.ndarray:
    """Schakel 2011 Eq. (1) source amplitude using the causal Ricker source.

    The default RT-SE workflow inserts the causal Ricker spectrum from
    ``seismoelectric_offset_liu2018_spectral.py`` directly as ``A(omega)``.
    The optional ``fig4_digitized`` branch is retained only as an explicit
    legacy diagnostic for the visual Schakel Fig. 4 pressure approximation.
    """
    freq = np.asarray(frequency_hz, dtype=float)
    if cfg.schakel_source_mode == "causal_ricker":
        omega = 2.0 * math.pi * freq
        return cfg.source_pressure_amp * base.causal_ricker_source_spectrum(omega, cfg)
    if cfg.schakel_source_mode != "fig4_digitized":
        raise ValueError("schakel_source_mode must be 'causal_ricker' or 'fig4_digitized'.")
    tau, pressure = _source_time_table()
    pressure = pressure * _source_time_taper(tau, cfg)
    omega = 2.0 * math.pi * freq
    spectrum = _trapz(pressure[None, :] * np.exp(-1j * omega[:, None] * tau[None, :]), tau, axis=1)
    return spectrum * _bandpass_taper(freq, cfg) * abs(float(cfg.z_s)) / 0.5


def _source_time_taper(tau: np.ndarray, cfg: ZeroOffsetSchakelConfig) -> np.ndarray:
    """Smooth the finite visual Fig. 4 pressure-record approximation endpoints without receiver gating."""
    tau = np.asarray(tau, dtype=float)
    taper_duration = max(0.0, float(cfg.schakel_source_taper_us)) * 1.0e-6
    if taper_duration <= 0.0 or len(tau) < 3:
        return np.ones_like(tau)
    dt = float(tau[1] - tau[0])
    n_taper = max(1, min(len(tau) // 2, int(round(taper_duration / max(dt, 1.0e-15)))))
    taper = np.ones_like(tau)
    edge = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, n_taper)))
    taper[:n_taper] = edge
    taper[-n_taper:] = edge[::-1]
    return taper


def _bandpass_taper(frequency_hz: np.ndarray, cfg: ZeroOffsetSchakelConfig) -> np.ndarray:
    """Smooth Schakel 2011 band-pass edges to avoid rectangular-window ringing."""
    f = np.asarray(frequency_hz, dtype=float)
    w = np.ones_like(f)
    low = float(cfg.schakel_bandpass_low_hz)
    high = float(cfg.schakel_bandpass_high_hz)
    taper = max(0.0, float(cfg.schakel_bandpass_taper_hz))
    w[(f <= low) | (f >= high)] = 0.0
    if taper <= 0.0:
        return w
    taper = min(taper, max((high - low) / 2.0, 1.0))
    lo = (f > low) & (f < low + taper)
    hi = (f < high) & (f > high - taper)
    w[lo] = 0.5 - 0.5 * np.cos(math.pi * (f[lo] - low) / taper)
    w[hi] = 0.5 - 0.5 * np.cos(math.pi * (high - f[hi]) / taper)
    return w


def pressure_normalized_re_from_coeff(coeff: Dict[str, complex], rho_f: float, omega: float) -> complex:
    """Convert Schakel 2010 displacement-normalized ``R_E`` to Schakel 2011 pressure normalization."""
    return complex(coeff["R_E"]) / (float(rho_f) * float(omega) ** 2)


def pressure_normalized_porous_terms_from_coeff(
    coeff: Dict[str, complex], rho_f: float, omega: float
) -> tuple[complex, complex]:
    """Return Schakel/JAP front-interface porous electric-potential terms.

    The signs follow the scalar-potential conversion used in
    ``schakel2011_jap_fig2_reproduction.py::_front_terms``.
    """
    scale = float(rho_f) * float(omega) ** 2
    tm_term = -complex(coeff["alpha_TM"]) * complex(coeff["T_TM"]) / scale
    pf_term = -complex(coeff["alpha_Pf"]) * complex(coeff["T_Pf"]) / scale
    return tm_term, pf_term


def _row_material_values(row: pd.Series, cfg: ZeroOffsetSchakelConfig) -> Dict[str, float | None]:
    phi = float(np.clip(float(row["Porosity"]), cfg.phi_min, cfg.phi_max_valid))
    return {
        "phi": phi,
        "k0_m2": max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min),
        "tau": max(float(row["Tortuosity"]), 1.0 + 1.0e-6),
        "cH": float(row["OutletHConc"]),
        "C_override": base.optional_float(row, "ElectrolyteConcentration_molL"),
        "sigma_f_override": base.optional_float(row, "FluidConductivity_S_m"),
    }


def _se_coeff_for_row(
    row_values: Dict[str, float | None],
    omega: float,
    cfg: ZeroOffsetSchakelConfig,
    theta_deg: float | None = None,
    kx_override: float | None = None,
) -> Dict[str, complex]:
    return base.se_coefficients(
        float(row_values["phi"]),
        float(row_values["k0_m2"]),
        float(row_values["tau"]),
        float(row_values["cH"]),
        omega,
        theta_deg,
        cfg,
        kx_override=kx_override,
        C_override_molL=row_values["C_override"],
        sigma_f_override=row_values["sigma_f_override"],
    )


def _safe_interface_terms(
    row_values: Dict[str, float | None],
    omega: float,
    cfg: ZeroOffsetSchakelConfig,
    theta_deg: float | None = None,
    kx_override: float | None = None,
) -> tuple[complex, complex, complex, complex, complex, complex]:
    try:
        coeff = _se_coeff_for_row(row_values, omega, cfg, theta_deg=theta_deg, kx_override=kx_override)
        re = pressure_normalized_re_from_coeff(coeff, cfg.rho_fl, omega)
        tm_term, pf_term = pressure_normalized_porous_terms_from_coeff(coeff, cfg.rho_fl, omega)
        kz_e = complex(coeff["k3_E"])
        k3_tm = complex(coeff["k3_TM"])
        k3_pf = complex(coeff["k3_Pf"])
        vals = [re, kz_e, tm_term, pf_term, k3_tm, k3_pf]
        if not all(np.isfinite(v.real) and np.isfinite(v.imag) for v in vals):
            return (0.0 + 0.0j,) * 6
        return re, kz_e, tm_term, pf_term, k3_tm, k3_pf
    except Exception:
        return (0.0 + 0.0j,) * 6


def _frequency_response_for_receivers(
    row: pd.Series,
    cfg: ZeroOffsetSchakelConfig,
    frequency_hz: float,
    z_receivers: np.ndarray,
    n_theta: int,
    integration_method: str = "fixed",
) -> np.ndarray:
    """Evaluate Schakel 2011 Eq. (5) for all receiver depths at one frequency."""
    if integration_method != "fixed":
        # The row-dependent coefficient is expensive and receiver-vectorized fixed
        # quadrature is more transparent here.  The original reproduction script
        # keeps adaptive quadrature for paper figures.
        raise ValueError("This RT-SE workflow currently supports integration_method='fixed'.")
    if int(n_theta) < 3:
        raise ValueError("n_theta must be >= 3")

    omega = 2.0 * math.pi * float(frequency_hz)
    k = omega / math.sqrt(cfg.K_fl / cfg.rho_fl)
    a = float(cfg.transducer_radius_m)

    # Receiver-line horizontal offset from the source axis.  In this file it is
    # forced to zero, so all receiver electrodes are vertically below/above the
    # source axis in a VSP-style line.
    radial_m = abs(float(cfg.offset_D))

    # Signed transmitter position used in Schakel Eq. (5).  The user-facing
    # config value ``z_s`` is a positive source-interface distance; the physical
    # source coordinate is negative because the upper fluid side is z < 0.
    source_z_m = -abs(float(cfg.z_s))
    row_values = _row_material_values(row, cfg)
    response = np.zeros(len(z_receivers), dtype=complex)

    # Receiver electrodes are split by side of the interface.  The y-axis in the
    # plotted gather is exactly these z_receivers values converted to millimetres:
    # z < 0 fluid-side/reflected EM, z > 0 porous-side/transmitted EM.
    interface_mask = np.isclose(z_receivers, 0.0, atol=1.0e-12, rtol=0.0)
    fluid_mask = z_receivers < 0.0
    porous_mask = z_receivers > 0.0
    z_fluid = z_receivers[fluid_mask]
    z_porous = z_receivers[porous_mask]

    theta = np.linspace(0.0, 0.5 * math.pi, int(n_theta))
    theta[0] = 1.0e-8
    real_fluid_terms = []
    real_fluid_kz = []
    real_porous_tm_terms = []
    real_porous_pf_terms = []
    real_porous_tm_kz = []
    real_porous_pf_kz = []
    for h in theta:
        sin_h = math.sin(float(h))
        coeff_re, kz_e, tm_term, pf_term, k3_tm, k3_pf = _safe_interface_terms(
            row_values, omega, cfg, theta_deg=math.degrees(float(h))
        )
        common = _j0(k * radial_m * sin_h) * _j1(k * a * sin_h)
        common *= np.exp(1j * k * source_z_m * math.cos(float(h)))
        real_fluid_terms.append(complex(common * coeff_re))
        real_fluid_kz.append(kz_e)
        real_porous_tm_terms.append(complex(common * tm_term))
        real_porous_pf_terms.append(complex(common * pf_term if cfg.include_porous_pf_coseismic else 0.0j))
        real_porous_tm_kz.append(k3_tm)
        real_porous_pf_kz.append(k3_pf)

    first_fluid = np.zeros(len(z_fluid), dtype=complex)
    if np.any(fluid_mask):
        terms = np.asarray(real_fluid_terms, dtype=complex)
        kz = np.asarray(real_fluid_kz, dtype=complex)
        integrand = terms[None, :] * np.exp(1j * z_fluid[:, None] * kz[None, :])
        first_fluid = _trapz(integrand, theta, axis=1)

    first_porous = np.zeros(len(z_porous), dtype=complex)
    if np.any(porous_mask):
        tm_terms = np.asarray(real_porous_tm_terms, dtype=complex)
        pf_terms = np.asarray(real_porous_pf_terms, dtype=complex)
        k3_tm_arr = np.asarray(real_porous_tm_kz, dtype=complex)
        k3_pf_arr = np.asarray(real_porous_pf_kz, dtype=complex)
        integrand = (
            tm_terms[None, :] * np.exp(-1j * z_porous[:, None] * k3_tm_arr[None, :])
            + pf_terms[None, :] * np.exp(-1j * z_porous[:, None] * k3_pf_arr[None, :])
        )
        first_porous = _trapz(integrand, theta, axis=1)

    second_fluid = np.zeros_like(first_fluid)
    second_porous = np.zeros_like(first_porous)
    if cfg.schakel_include_evanescent:
        gammas = np.linspace(0.0, float(cfg.schakel_gamma_max), max(4, int(n_theta) // 2))
        ev_fluid_terms = []
        ev_fluid_kz = []
        ev_porous_tm_terms = []
        ev_porous_pf_terms = []
        ev_porous_tm_kz = []
        ev_porous_pf_kz = []
        for gamma in gammas:
            root = math.sqrt(float(gamma) ** 2 + 1.0)
            kx = k * root
            coeff_re, kz_e, tm_term, pf_term, k3_tm, k3_pf = _safe_interface_terms(
                row_values, omega, cfg, theta_deg=90.0, kx_override=kx
            )
            common = _j0(k * radial_m * root) * (_j1(k * a * root) / root)
            common *= math.exp(k * source_z_m * float(gamma))
            ev_fluid_terms.append(complex(common * coeff_re))
            ev_fluid_kz.append(kz_e)
            ev_porous_tm_terms.append(complex(common * tm_term))
            ev_porous_pf_terms.append(complex(common * pf_term if cfg.include_porous_pf_coseismic else 0.0j))
            ev_porous_tm_kz.append(k3_tm)
            ev_porous_pf_kz.append(k3_pf)

        if np.any(fluid_mask):
            terms = np.asarray(ev_fluid_terms, dtype=complex)
            kz = np.asarray(ev_fluid_kz, dtype=complex)
            integrand = terms[None, :] * np.exp(1j * z_fluid[:, None] * kz[None, :])
            second_fluid = _trapz(integrand, gammas, axis=1)

        if np.any(porous_mask):
            tm_terms = np.asarray(ev_porous_tm_terms, dtype=complex)
            pf_terms = np.asarray(ev_porous_pf_terms, dtype=complex)
            k3_tm_arr = np.asarray(ev_porous_tm_kz, dtype=complex)
            k3_pf_arr = np.asarray(ev_porous_pf_kz, dtype=complex)
            integrand = (
                tm_terms[None, :] * np.exp(-1j * z_porous[:, None] * k3_tm_arr[None, :])
                + pf_terms[None, :] * np.exp(-1j * z_porous[:, None] * k3_pf_arr[None, :])
            )
            second_porous = _trapz(integrand, gammas, axis=1)

    A_omega = schakel_source_A_spectrum(np.array([frequency_hz]), cfg)[0]
    if np.any(fluid_mask):
        response[fluid_mask] = -(1j * A_omega / a) * first_fluid + (A_omega / a) * second_fluid
    if np.any(porous_mask):
        response[porous_mask] = -(1j * A_omega / a) * first_porous + (A_omega / a) * second_porous
    if np.any(interface_mask) and np.any(fluid_mask) and np.any(porous_mask):
        i_fluid = np.where(fluid_mask)[0][np.argmin(np.abs(z_receivers[fluid_mask]))]
        i_porous = np.where(porous_mask)[0][np.argmin(np.abs(z_receivers[porous_mask]))]
        response[interface_mask] = 0.5 * (response[i_fluid] + response[i_porous])
    return response


def synthesize_waveforms_schakel2011(
    row: pd.Series,
    cfg: ZeroOffsetSchakelConfig,
    n_frequencies: int | None = None,
    n_theta: int = 48,
    integration_method: str = "fixed",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Synthesize a zero-offset Schakel 2011 interface EM waveform gather."""
    cfg.offset_D = 0.0
    n_frequencies = int(n_frequencies or cfg.spectral_n_omega)

    # Receiver/electrode array definition.  These three configuration values
    # control the vertical positions and the receiver count in the figure:
    #
    #   receiver_z_min   first electrode position relative to interface (m)
    #   receiver_z_max   last electrode position relative to interface (m)
    #   receiver_spacing electrode spacing (m)
    #
    # Example default inherited from base.SEConfig:
    #   -0.100 m to +0.100 m every 0.001 m -> about 201 receiver traces.
    # The plotted y-axis multiplies these values by 1000, giving -100..100 mm.
    z_receivers = np.arange(
        cfg.receiver_z_min,
        cfg.receiver_z_max + 0.5 * cfg.receiver_spacing,
        cfg.receiver_spacing,
    )
    z_receivers[np.isclose(z_receivers, 0.0, atol=1.0e-12, rtol=0.0)] = 0.0
    vf = math.sqrt(cfg.K_fl / cfg.rho_fl)

    # Interface-arrival time from source to interface.  The source is located
    # one ``z_s`` above the interface, so T0 = z_s / upper-fluid acoustic speed.
    # Keep the output/gather time axis starting at 0 s for visual context; the
    # T0 arrival is still marked by plot_waveform_gather().
    T0 = cfg.z_s / vf
    t = np.linspace(0.0, T0 + cfg.waveform_t_after, cfg.waveform_nt)
    frequencies = _frequency_grid(cfg, n_frequencies)

    spectrum = np.zeros((len(z_receivers), len(frequencies)), dtype=complex)
    for i, freq in enumerate(frequencies):
        spectrum[:, i] = _frequency_response_for_receivers(
            row,
            cfg,
            float(freq),
            z_receivers,
            n_theta=int(n_theta),
            integration_method=integration_method,
        )

    phase = np.exp(1j * 2.0 * math.pi * frequencies[None, :] * t[:, None])
    traces = 2.0 * np.real(_trapz(spectrum[:, None, :] * phase[None, :, :], frequencies, axis=2))
    return z_receivers, t, traces


def compute_peak_amplitude_schakel2011(
    ts: pd.DataFrame,
    df_raw: pd.DataFrame,
    cfg: ZeroOffsetSchakelConfig,
    n_frequencies: int | None = None,
    n_theta: int = 48,
    integration_method: str = "fixed",
) -> pd.DataFrame:
    """Compute waveform peak columns while preserving the original output schema."""
    out = ts.copy()
    vals_all = []
    vals_re = []
    vals_te = []
    for idx, row in df_raw.reset_index(drop=True).iterrows():
        if idx >= len(out) or not bool(out.loc[idx, "valid_poroelastic"]):
            vals_all.append(np.nan)
            vals_re.append(np.nan)
            vals_te.append(np.nan)
            continue
        try:
            z, _, u = synthesize_waveforms_schakel2011(
                row,
                cfg,
                n_frequencies=n_frequencies,
                n_theta=n_theta,
                integration_method=integration_method,
            )
            vals_all.append(float(np.nanmax(np.abs(u))))
            re_mask = z <= 0.0
            te_mask = z > 0.0
            vals_re.append(float(np.nanmax(np.abs(u[re_mask, :]))) if np.any(re_mask) else np.nan)
            vals_te.append(float(np.nanmax(np.abs(u[te_mask, :]))) if np.any(te_mask) else np.nan)
        except Exception:
            vals_all.append(np.nan)
            vals_re.append(np.nan)
            vals_te.append(np.nan)

    out["Amax_waveform_schakel2011"] = vals_all
    out["Amax_waveform_schakel2011_RE"] = vals_re
    out["Amax_waveform_schakel2011_TE"] = vals_te

    # Compatibility with the plotting functions from the Liu-style script.
    out["Amax_waveform_spectral"] = out["Amax_waveform_schakel2011"]
    out["Amax_waveform_spectral_RE"] = out["Amax_waveform_schakel2011_RE"]
    out["Amax_waveform_spectral_TE"] = out["Amax_waveform_schakel2011_TE"]

    for col in [
        "Amax_waveform_schakel2011",
        "Amax_waveform_schakel2011_RE",
        "Amax_waveform_schakel2011_TE",
        "Amax_waveform_spectral",
        "Amax_waveform_spectral_RE",
        "Amax_waveform_spectral_TE",
    ]:
        valid_vals = out.loc[out["valid_poroelastic"] & np.isfinite(out[col]), col]
        ref = valid_vals.iloc[0] if len(valid_vals) else np.nan
        out[f"{col}_norm"] = out[col] / ref if ref and np.isfinite(ref) and ref != 0 else np.nan
    return out


def save_parameter_table(cfg: ZeroOffsetSchakelConfig, outdir: Path) -> None:
    base.save_parameter_table(cfg, outdir)
    path = outdir / "parameters_used.csv"
    df = pd.read_csv(path)
    extra = pd.DataFrame(
        [
            {
                "parameter": "transducer_radius_m",
                "meaning": "Schakel 2011 finite piston transducer radius",
                "value": cfg.transducer_radius_m,
                "unit": "m",
            },
            {
                "parameter": "schakel_bandpass_low_hz",
                "meaning": "Schakel 2011 lower frequency bound",
                "value": cfg.schakel_bandpass_low_hz,
                "unit": "Hz",
            },
            {
                "parameter": "schakel_bandpass_high_hz",
                "meaning": "Schakel 2011 upper frequency bound",
                "value": cfg.schakel_bandpass_high_hz,
                "unit": "Hz",
            },
            {
                "parameter": "schakel_bandpass_taper_hz",
                "meaning": "Smooth taper width applied at Schakel source-spectrum band edges",
                "value": cfg.schakel_bandpass_taper_hz,
                "unit": "Hz",
            },
            {
                "parameter": "schakel_source_mode",
                "meaning": "causal Ricker source spectrum used inside the Schakel 2011 Sommerfeld integral by default",
                "value": cfg.schakel_source_mode,
                "unit": "-",
            },
            {
                "parameter": "include_porous_pf_coseismic",
                "meaning": "Include JAP Eq. (8) Pf coseismic term in addition to interface TM term",
                "value": cfg.include_porous_pf_coseismic,
                "unit": "-",
            },
            {
                "parameter": "schakel_gamma_max",
                "meaning": "Upper gamma limit for Schakel 2011 evanescent Sommerfeld branch",
                "value": cfg.schakel_gamma_max,
                "unit": "-",
            },
            {
                "parameter": "schakel_include_evanescent",
                "meaning": "Whether Eq. (5) evanescent branch is included",
                "value": cfg.schakel_include_evanescent,
                "unit": "-",
            },
        ]
    )
    pd.concat([df, extra], ignore_index=True).to_csv(path, index=False)


def write_formula_audit(outdir: Path) -> None:
    text = """# Formula Audit

Implemented model:

- The reactive-transport mapping from porosity, permeability, tortuosity, and H+ concentration to dynamic permeability, electrokinetic coupling, dynamic conductivity, and Schakel interface coefficients is reused from `seismoelectric_offset_liu2018_spectral.py`.
- Receiver geometry is reused from `seismoelectric_offset_liu2018_spectral.py`, but `offset_D` is forced to `0.0 m`.
- Waveform synthesis replaces the Liu 2018 frequency-wavenumber integral with a Schakel et al. (2011) Sommerfeld integral. Fluid-side receivers use the reflected electric-potential structure of Schakel/JAP Eq. (5). Porous-side receivers use the front-interface transmitted TM electric-potential term from the JAP Eq. (8) structure, `-alpha_TM*T_TM/(rho_fl*omega**2)`, propagated with its Schakel vertical wavenumber. The Pf coseismic term `-alpha_Pf*T_Pf/(rho_fl*omega**2)` is available as an explicit diagnostic option but is off by default so the waveform emphasizes the interface EM response rather than a porous acoustic-coseismic arrival.
- The default RT-SE source mode is the existing causal Ricker source spectrum from `seismoelectric_offset_liu2018_spectral.py`, inserted as `A(omega)` in the Schakel 2011 pressure-source integral. Its default frequency band starts near zero and extends to 2.5 times the 500 kHz source frequency to reduce noncausal-looking band-pass side lobes. The optional legacy `fig4_digitized` mode uses a visual approximation to the Schakel 2011 laboratory Fig. 4 source and evaluates it by direct causal Fourier integration with `exp(-i omega tau)`, not by periodic FFT interpolation.
- Schakel 2011 Sommerfeld integration is evaluated with a finite piston term `J1(k a sin(theta))`, a real-angle integral over `0..pi/2`, and the evanescent gamma branch.
- The interface conversion coefficient is computed from the Schakel and Smeulders (2010) Appendix B boundary-value solver and converted to the Schakel 2011 pressure-normalized coefficient as `R_E / (rho_fl * omega**2)`.
- Time synthesis uses positive frequencies and `exp(i omega t)`, consistent with the Schakel convention.
- The default waveform output window starts at 0 s and marks T0, the acoustic arrival time at the interface, so the saved interface-EM gather displays the full pre-interface-arrival interval without receiver-side trace gating.

Important limitation:

- This is a single-interface RT-SE forward model, not the finite-thickness Schakel 2011 laboratory slab. The porous side includes the front-interface Eq. (8) TM term but omits back-interface/sample-width multiple reflections. If `include_porous_pf_coseismic=True`, the Pf coseismic term is added for diagnostics and should not be interpreted as pure interface EM.
- Finite frequency bands can leave T0-before side lobes because the default RT-SE window now starts at 0 s. The model does not gate receiver traces; these pre-T0 samples are retained for explicit side-lobe audits. The default near-zero low-frequency limit is chosen for the RT-SE causal-source simulation, while Schakel laboratory band-pass settings should be treated as a separate reproduction/diagnostic mode.
"""
    (outdir / "formula_audit.md").write_text(text, encoding="utf-8")


def interface_em_polarity_diagnostics(
    z: np.ndarray,
    t: np.ndarray,
    u: np.ndarray,
    cfg: ZeroOffsetSchakelConfig,
    distances_m: Iterable[float] = (1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3),
) -> pd.DataFrame:
    """Diagnose polarity at a common interface-EM arrival time for symmetric receivers."""
    T0 = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)
    rows = []
    for distance in distances_m:
        i_fluid = int(np.argmin(np.abs(z + distance)))
        i_porous = int(np.argmin(np.abs(z - distance)))
        window = (t >= T0) & (t <= T0 + 8.0e-6)
        if not np.any(window):
            continue
        combined = np.abs(u[i_fluid, window]) + np.abs(u[i_porous, window])
        local_idx = int(np.nanargmax(combined))
        global_idx = np.where(window)[0][local_idx]
        fluid_val = float(u[i_fluid, global_idx])
        porous_val = float(u[i_porous, global_idx])
        rows.append(
            {
                "distance_from_interface_m": float(distance),
                "distance_from_interface_mm": float(distance * 1.0e3),
                "common_arrival_time_s": float(t[global_idx]),
                "common_arrival_time_us": float(t[global_idx] * 1.0e6),
                "time_after_T0_us": float((t[global_idx] - T0) * 1.0e6),
                "fluid_z_m": float(z[i_fluid]),
                "porous_z_m": float(z[i_porous]),
                "fluid_signed": fluid_val,
                "porous_signed": porous_val,
                "signed_product": fluid_val * porous_val,
                "polarity_reversed": bool(fluid_val * porous_val < 0.0),
                "fluid_abs": abs(fluid_val),
                "porous_abs": abs(porous_val),
            }
        )
    return pd.DataFrame(rows)


def t0_causality_diagnostics(z: np.ndarray, t: np.ndarray, u: np.ndarray,
                             cfg: ZeroOffsetSchakelConfig) -> pd.DataFrame:
    """Summarize the saved T0 window and T0-after main response without trace gating."""
    T0 = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)
    receiver_mask = ~np.isclose(z, 0.0, atol=1.0e-12, rtol=0.0)
    z_eval = z[receiver_mask]
    u_eval = u[receiver_mask, :]
    pre_mask = t < T0
    post_mask = t >= T0
    rows = []
    if np.any(pre_mask):
        pre_abs = np.abs(u_eval[:, pre_mask])
        pre_idx = np.unravel_index(int(np.nanargmax(pre_abs)), pre_abs.shape)
        pre_global_time = np.where(pre_mask)[0][pre_idx[1]]
        pre_max = float(pre_abs[pre_idx])
        pre_z = float(z_eval[pre_idx[0]])
        pre_time = float(t[pre_global_time])
    else:
        pre_max = 0.0
        pre_z = np.nan
        pre_time = np.nan
    if np.any(post_mask):
        post_abs = np.abs(u_eval[:, post_mask])
        post_idx = np.unravel_index(int(np.nanargmax(post_abs)), post_abs.shape)
        post_global_time = np.where(post_mask)[0][post_idx[1]]
        post_max = float(post_abs[post_idx])
        post_z = float(z_eval[post_idx[0]])
        post_time = float(t[post_global_time])
    else:
        post_max = np.nan
        post_z = np.nan
        post_time = np.nan
    rows.append(
        {
            "T0_s": float(T0),
            "T0_us": float(T0 * 1.0e6),
            "pre_T0_max_abs": pre_max,
            "pre_T0_max_time_s": pre_time,
            "pre_T0_max_time_us": pre_time * 1.0e6 if np.isfinite(pre_time) else np.nan,
            "pre_T0_max_time_before_T0_us": (T0 - pre_time) * 1.0e6 if np.isfinite(pre_time) else np.nan,
            "pre_T0_max_z_m": pre_z,
            "pre_T0_max_z_mm": pre_z * 1.0e3 if np.isfinite(pre_z) else np.nan,
            "post_T0_max_abs": post_max,
            "post_T0_max_time_s": post_time,
            "post_T0_max_time_us": post_time * 1.0e6 if np.isfinite(post_time) else np.nan,
            "post_T0_max_time_after_T0_us": (post_time - T0) * 1.0e6 if np.isfinite(post_time) else np.nan,
            "post_T0_max_z_m": post_z,
            "post_T0_max_z_mm": post_z * 1.0e3 if np.isfinite(post_z) else np.nan,
            "pre_to_post_abs_ratio": pre_max / post_max if post_max and np.isfinite(post_max) else np.nan,
            "main_peak_after_T0": bool(np.isfinite(post_time) and post_max >= pre_max),
            "receiver_traces_gated": False,
            "pre_T0_samples_present": bool(np.any(pre_mask)),
            "interpretation": (
                "output window starts at or after T0; pre_T0_max_abs=0 is a window statement, not trace gating"
                if not np.any(pre_mask)
                else "pre-T0 samples are present; any pre-T0 energy is finite-band inverse-transform side lobe, not receiver-side gating"
            ),
        }
    )
    return pd.DataFrame(rows)


def run_simulation(
    input_path: str | Path,
    outdir: str | Path,
    snapshot_target_phi: float = 0.75,
    n_frequencies: int | None = None,
    n_theta: int = 48,
    peak_n_frequencies: int | None = None,
    peak_n_theta: int = 32,
    integration_method: str = "fixed",
    cfg: ZeroOffsetSchakelConfig | None = None,
) -> Path:
    cfg = ZeroOffsetSchakelConfig() if cfg is None else cfg
    cfg.offset_D = 0.0
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = base.load_reactive_transport_table(input_path)
    ts = base.compute_time_series(df, cfg)
    ts = compute_peak_amplitude_schakel2011(
        ts,
        df,
        cfg,
        n_frequencies=peak_n_frequencies or max(3, cfg.spectral_n_omega // 2),
        n_theta=peak_n_theta,
        integration_method=integration_method,
    )
    ts.to_csv(outdir / "seismoelectric_timeseries_results.csv", index=False)
    save_parameter_table(cfg, outdir)
    base.plot_coefficients(ts, outdir)
    base.plot_peak_amplitude(ts, outdir)

    idx = base.choose_snapshot(df, target_phi=snapshot_target_phi)
    row = df.iloc[idx]
    z, t, u = synthesize_waveforms_schakel2011(
        row,
        cfg,
        n_frequencies=n_frequencies or cfg.spectral_n_omega,
        n_theta=n_theta,
        integration_method=integration_method,
    )
    plot_name = "waveform_snapshot_schakel2011"
    base.save_waveform_arrays(z, t, u, outdir, plot_name)
    diag = base.save_waveform_spatial_peak_diagnostics(z, t, u, cfg, outdir)
    polarity_diag = interface_em_polarity_diagnostics(z, t, u, cfg)
    polarity_diag.to_csv(outdir / "waveform_interface_em_polarity_diagnostics.csv", index=False)
    causality_diag = t0_causality_diagnostics(z, t, u, cfg)
    causality_diag.to_csv(outdir / "waveform_t0_causality_diagnostics.csv", index=False)
    amax_snapshot = base.plot_waveform_gather(z, t, u, row, cfg, outdir, plot_name)

    re_diag = diag[diag["side"].isin(["R_E", "interface"])]
    te_diag = diag[diag["side"] == "T_E"]
    re_peak_distance_mm = (
        float(re_diag.loc[re_diag["peak_abs"].idxmax(), "distance_from_interface_mm"]) if not re_diag.empty else np.nan
    )
    te_peak_distance_mm = (
        float(te_diag.loc[te_diag["peak_abs"].idxmax(), "distance_from_interface_mm"]) if not te_diag.empty else np.nan
    )
    T0 = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)
    summary = {
        "input": str(input_path),
        "outdir": str(outdir),
        "snapshot_index": int(idx),
        "snapshot_Time_s": float(row["Time_s"]),
        "snapshot_Porosity": float(row["Porosity"]),
        "snapshot_Amax": amax_snapshot,
        "T0_us": T0 * 1.0e6,
        "waveform_mode": "schakel2011_sommerfeld_zerooffset",
        "offset_D_m": float(cfg.offset_D),
        "pre_T0_max_abs": float(np.nanmax(np.abs(u[:, t < T0]))) if np.any(t < T0) else 0.0,
        "post_T0_max_abs": float(np.nanmax(np.abs(u[:, t >= T0]))) if np.any(t >= T0) else np.nan,
        "noninterface_pre_to_post_abs_ratio": float(causality_diag["pre_to_post_abs_ratio"].iloc[0]) if not causality_diag.empty else np.nan,
        "noninterface_pre_T0_max_time_before_T0_us": float(causality_diag["pre_T0_max_time_before_T0_us"].iloc[0]) if not causality_diag.empty else np.nan,
        "noninterface_main_peak_after_T0": bool(causality_diag["main_peak_after_T0"].iloc[0]) if not causality_diag.empty else False,
        "RE_peak_distance_from_interface_mm": re_peak_distance_mm,
        "TE_peak_distance_from_interface_mm": te_peak_distance_mm,
        "schakel_n_frequencies": int(n_frequencies or cfg.spectral_n_omega),
        "schakel_n_theta": int(n_theta),
        "peak_schakel_n_frequencies": int(peak_n_frequencies or max(3, cfg.spectral_n_omega // 2)),
        "peak_schakel_n_theta": int(peak_n_theta),
        "schakel_bandpass_low_hz": float(cfg.schakel_bandpass_low_hz),
        "schakel_bandpass_high_hz": float(cfg.schakel_bandpass_high_hz),
        "schakel_gamma_max": float(cfg.schakel_gamma_max),
        "source_mode": cfg.schakel_source_mode,
        "include_porous_pf_coseismic": bool(cfg.include_porous_pf_coseismic),
        "interface_em_polarity_pairs_checked": int(len(polarity_diag)),
        "interface_em_polarity_pairs_reversed": int(polarity_diag["polarity_reversed"].sum()) if not polarity_diag.empty else 0,
        "interface_em_mean_arrival_after_T0_us": float(polarity_diag["time_after_T0_us"].mean()) if not polarity_diag.empty else np.nan,
        "porous_side_waveform_policy": "front_interface_TM_term_from_JAP_Eq8; Pf coseismic optional; no back-interface slab multiples",
    }
    pd.Series(summary).to_csv(outdir / "run_summary.csv")
    write_formula_audit(outdir)
    return outdir


def _apply_common_overrides(args: argparse.Namespace, cfg: ZeroOffsetSchakelConfig) -> None:
    if args.z_s is not None:
        cfg.z_s = args.z_s
    if args.z_s_mm is not None:
        cfg.z_s = args.z_s_mm * 1.0e-3
    if args.receiver_z_min_mm is not None:
        cfg.receiver_z_min = args.receiver_z_min_mm * 1.0e-3
    if args.receiver_z_max_mm is not None:
        cfg.receiver_z_max = args.receiver_z_max_mm * 1.0e-3
    if args.receiver_spacing_mm is not None:
        cfg.receiver_spacing = args.receiver_spacing_mm * 1.0e-3
    if args.f0 is not None:
        cfg.f0 = args.f0
    if args.upper_fluid_conductivity_mode is not None:
        cfg.upper_fluid_conductivity_mode = args.upper_fluid_conductivity_mode
    if args.transducer_radius_mm is not None:
        cfg.transducer_radius_m = args.transducer_radius_mm * 1.0e-3
    if args.schakel_bandpass_low_hz is not None:
        cfg.schakel_bandpass_low_hz = args.schakel_bandpass_low_hz
    if args.schakel_bandpass_high_hz is not None:
        cfg.schakel_bandpass_high_hz = args.schakel_bandpass_high_hz
    if args.schakel_gamma_max is not None:
        cfg.schakel_gamma_max = args.schakel_gamma_max
    if args.source_mode is not None:
        cfg.schakel_source_mode = args.source_mode


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="global_evolution.xlsx")
    parser.add_argument("--outdir", type=str, default="se_results_zerooffset_schakel2011")
    parser.add_argument("--snapshot-target-phi", type=float, default=0.75)
    parser.add_argument("--n-frequencies", type=int, default=None)
    parser.add_argument("--n-theta", type=int, default=48)
    parser.add_argument("--peak-n-frequencies", type=int, default=None)
    parser.add_argument("--peak-n-theta", type=int, default=32)
    parser.add_argument("--integration-method", choices=["fixed"], default="fixed")
    parser.add_argument("--z_s", type=float, default=None)
    parser.add_argument("--z_s_mm", type=float, default=None)
    parser.add_argument("--receiver-z-min-mm", type=float, default=None)
    parser.add_argument("--receiver-z-max-mm", type=float, default=None)
    parser.add_argument("--receiver-spacing-mm", type=float, default=None)
    parser.add_argument("--offset-D-mm", type=float, default=None, help="Accepted for compatibility; always forced to zero.")
    parser.add_argument("--f0", type=float, default=None)
    parser.add_argument("--upper-fluid-conductivity-mode", choices=["constant", "dynamic_pore_fluid"], default=None)
    parser.add_argument("--transducer-radius-mm", type=float, default=None)
    parser.add_argument("--schakel-bandpass-low-hz", type=float, default=None)
    parser.add_argument("--schakel-bandpass-high-hz", type=float, default=None)
    parser.add_argument("--schakel-gamma-max", type=float, default=None)
    parser.add_argument("--source-mode", choices=["causal_ricker", "fig4_digitized"], default=None)
    args = parser.parse_args()

    cfg = ZeroOffsetSchakelConfig()
    _apply_common_overrides(args, cfg)
    if args.offset_D_mm not in (None, 0, 0.0):
        print("Warning: --offset-D-mm is ignored; this Schakel 2011 workflow forces offset_D=0.")
    cfg.offset_D = 0.0
    outdir = run_simulation(
        args.input,
        args.outdir,
        snapshot_target_phi=args.snapshot_target_phi,
        n_frequencies=args.n_frequencies,
        n_theta=args.n_theta,
        peak_n_frequencies=args.peak_n_frequencies,
        peak_n_theta=args.peak_n_theta,
        integration_method=args.integration_method,
        cfg=cfg,
    )
    print("Done. Outputs written to:", outdir)
    for path in sorted(Path(outdir).glob("*")):
        print(" -", path.name)


if __name__ == "__main__":
    main()
