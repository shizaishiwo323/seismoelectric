#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reactive-transport-driven seismoelectric prototype.

Purpose
-------
Read pore-scale reactive transport outputs, map porosity/permeability/tortuosity/H+ concentration
into the Schakel & Smeulders (2010, JASA) fluid--porous-medium interface model, compute:
  1) seismoelectric reflection coefficient R_E versus dissolution time,
  2) TM-mode seismoelectric transmission coefficient T_TM versus dissolution time,
  3) VSEP-style interface EM waveform gather for one selected dissolution snapshot,
  4) maximum waveform peak amplitude versus dissolution time.

Important scope
---------------
- R_E and T_TM are computed from the Schakel & Smeulders 6x6 boundary-value system.
- The default waveform gather is a Liu et al. (2018)-style frequency--wavenumber spectral
  synthesis. In that mode, the receiver waveforms are
  synthesized by summing positive-frequency and horizontal-wavenumber components, each
  weighted by the acoustic source spectrum A(k,omega), the Schakel conversion coefficient
  converted to Liu's electrical-potential coefficient, and the corresponding propagation phase.
- Liu Fig. 2(a) is the frequency-wavenumber modeling result. Interpretation
  models used for separate amplitude-directivity checks are intentionally kept
  outside this main spectral forward script.

Usage
-----
python seismoelectric_reactive_transport.py \
    --input b1ce6e99-a989-4127-86f5-92aaef7851f7.xlsx \
    --outdir se_results

The input file should contain columns similar to:
Time_s, Porosity, Permeability_mD, Tortuosity, OutletHConc.
"""

from __future__ import annotations

import argparse
import math
import os
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Tuple, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# 1. Configuration parameters
# -----------------------------

@dataclass
class SEConfig:
    # Geometry: user's pore-scale RTM domain
    Lx: float = 600e-6                 # m, domain width in x direction
    Ly: float = 400e-6                 # m, porous-domain height in z direction
    grain_radius: float = 30e-6        # m
    pore_throat_width: float = 10e-6   # m

    # VSEP geometry
    z_s: float = 80.0e-3               # m, source-to-interface distance in fluid, Liu-style default
    receiver_z_min: float = -100.0e-3  # m, fluid side negative; Liu Fig. 2 uses -100..100 mm
    receiver_z_max: float = 100.0e-3   # m, porous side positive
    receiver_spacing: float = 1.0e-3   # m, Liu Fig. 2 receiver trace interval
    offset_D: float = 45e-3           # m, receiver-line horizontal offset

    # Acoustic source
    f0: float = 500.0e3                # Hz, Liu-style ultrasonic source central frequency
    source_pressure_amp: float = 1.0   # arbitrary unit for waveform visualization
    source_kb_m_inv: float = 80.0      # 1/m, Liu finite-width source k_b for A(k,omega)
    source_beam_theta_deg: float = 0.0 # deg, central acoustic beam angle in Liu A(k,omega)
    source_peak_cycles: float = 2.0    # Ricker peak time after source onset, in f0 cycles
    source_duration_cycles: float = 8.0 # finite causal source duration used to compute S(omega)
    spectral_f_min_factor: float = 0.25 # f_min = factor * f0 for positive-frequency integration
    spectral_f_max_factor: float = 2.5  # f_max = factor * f0 for positive-frequency integration
    spectral_n_omega: int = 48          # default number of frequency samples for spectral mode
    spectral_n_k: int = 401             # k integration needs fine sampling for finite-offset cancellation
    spectral_k_limit_factor: float = 0.98 # incident acoustic branch limit |k|<=factor*omega/Vf
    waveform_t_before: float = 15.0e-6 # s, time window before T0
    waveform_t_after: float = 18.0e-6  # s, time window after T0
    waveform_nt: int = 1200            # number of time samples

    # Coefficient evaluation setting
    coeff_theta_deg: float = 45.0      # deg, because plane-wave R_E is zero at theta=0

    # Baseline porous frame / mineral / fluid parameters
    # K_s and rho_s are set for calcite; K_b and G are framework-scale effective values and should
    # be calibrated if independent mechanical data become available.
    K_b: float = 5.8e9                 # Pa, drained framework bulk modulus
    G: float = 3.4e9                   # Pa, framework shear modulus
    K_s: float = 76.8e9                # Pa, calcite grain bulk modulus
    K_f: float = 2.22e9                # Pa, pore fluid bulk modulus, water-like acid solution
    K_fl: float = 2.22e9               # Pa, upper fluid bulk modulus
    eta: float = 1.0e-3                # Pa s, dynamic viscosity
    rho_f: float = 1000.0              # kg/m^3, pore fluid density
    rho_fl: float = 1000.0             # kg/m^3, upper fluid density
    rho_s: float = 2710.0              # kg/m^3, calcite density
    eps_f: float = 80.0                # relative permittivity of pore fluid
    eps_s: float = 8.0                 # relative permittivity of calcite/solid
    eps_fl: float = 80.0               # relative permittivity of upper fluid
    sigma_fl_default: float = 5.0e-3   # S/m, used if upper_fluid_conductivity_mode='constant'

    # Electrochemistry
    temperature: float = 295.0         # K
    C_background_molL: float = 1.0e-3  # mol/L, background electrolyte concentration
    H_min_molL: float = 1.0e-7         # mol/L, neutral-water floor for pH calculation
    outlet_h_unit: str = "mol_cm3"     # "mol_cm3" or "mol_L"
    upper_fluid_conductivity_mode: str = "constant"  # "constant" strictly follows Schakel Table I; "dynamic_pore_fluid" is optional
    M_similarity: float = 1.0          # Schakel Table I uses M=1
    # Numerical stability / validity
    phi_min: float = 0.05
    phi_max_valid: float = 0.95        # beyond this, poroelastic framework is nearly gone
    k0_min: float = 1e-18              # m^2
    eps_complex: float = 1e-30

    # Physical constants
    eps0: float = 8.854187817e-12      # F/m
    mu0: float = 4.0 * math.pi * 1e-7  # H/m
    e_charge: float = 1.602176634e-19  # C
    k_B: float = 1.380649e-23          # J/K
    N_A: float = 6.02214076e23         # 1/mol
    z1: int = 1
    z2: int = -1
    b1: float = 3.246e11               # m/(N s), Schakel Table I
    b2: float = 4.931e11               # m/(N s), Schakel Table I


# -----------------------------
# 2. Utilities
# -----------------------------

def complex_sqrt_branch(x: complex | np.ndarray) -> complex | np.ndarray:
    """Square-root branch with non-negative imaginary part; if imag≈0, positive real part."""
    y = np.sqrt(x + 0j)
    if np.isscalar(y):
        if np.imag(y) < 0 or (abs(np.imag(y)) < 1e-18 and np.real(y) < 0):
            y = -y
        return y
    mask = (np.imag(y) < 0) | ((np.abs(np.imag(y)) < 1e-18) & (np.real(y) < 0))
    y[mask] = -y[mask]
    return y


def ricker(t: np.ndarray, f0: float) -> np.ndarray:
    """Ricker wavelet centered at t=0."""
    a = math.pi * f0 * t
    return (1.0 - 2.0 * a**2) * np.exp(-a**2)


def load_reactive_transport_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(path)
    elif path.suffix.lower() in [".csv", ".txt"]:
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported input format: {path.suffix}. Use xlsx/csv/txt.")

    # Standardize a few likely column names.
    rename = {}
    for c in df.columns:
        key = c.strip().lower()
        if key in ["time_s", "time", "times", "t_s"]:
            rename[c] = "Time_s"
        elif key in ["porosity", "phi", "孔隙率"]:
            rename[c] = "Porosity"
        elif key in ["permeability_md", "k_md", "permeability"]:
            rename[c] = "Permeability_mD"
        elif key in ["tortuosity", "tau", "迂曲度", "曲折度"]:
            rename[c] = "Tortuosity"
        elif key in ["outlethconc", "outlet_h_conc", "h_conc", "hplus", "h+"]:
            rename[c] = "OutletHConc"
        elif key in ["surfacearea_cm2", "surface_area_cm2"]:
            rename[c] = "SurfaceArea_cm2"
        elif key in ["fluidconductivity_s_m", "fluid_conductivity_s_m", "sigma_f", "sigma_fluid"]:
            rename[c] = "FluidConductivity_S_m"
        elif key in ["electrolyteconcentration_moll", "electrolyte_concentration_moll", "c_moll", "c_mol_l"]:
            rename[c] = "ElectrolyteConcentration_molL"
    df = df.rename(columns=rename)

    required = ["Time_s", "Porosity", "Permeability_mD", "Tortuosity", "OutletHConc"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Existing columns: {list(df.columns)}")

    df = df.copy()
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=required).reset_index(drop=True)
    return df


# -----------------------------
# 3. Schakel & Smeulders dynamic coefficients
# -----------------------------

def h_conc_to_molL(cH: float, cfg: SEConfig) -> float:
    """Convert OutletHConc to mol/L, with a lower floor for pH calculation."""
    if not np.isfinite(cH):
        cH = 0.0
    cH = max(float(cH), 0.0)
    if cfg.outlet_h_unit.lower() in ["mol_cm3", "mol/cm3", "mol cm-3", "molcm3"]:
        cH_molL = cH * 1000.0  # 1 mol/cm^3 = 1000 mol/L
    elif cfg.outlet_h_unit.lower() in ["mol_l", "mol/l", "moll", "mol_liter"]:
        cH_molL = cH
    else:
        raise ValueError("outlet_h_unit must be 'mol_cm3' or 'mol_L'.")
    return max(cH_molL, cfg.H_min_molL)


def electrochemistry_from_h(cH: float, cfg: SEConfig, C_override_molL: float | None = None) -> Dict[str, float]:
    cH_molL = h_conc_to_molL(cH, cfg)
    pH = -math.log10(cH_molL)
    # If measured electrolyte concentration is available, use it. Otherwise, if only outlet H+
    # is available, approximate the electrolyte as background salt plus added HCl.
    # The latter is a model assumption, not an extra formula from Schakel & Smeulders.
    if C_override_molL is not None and np.isfinite(C_override_molL) and C_override_molL > 0:
        C = float(C_override_molL)
    else:
        C = cfg.C_background_molL + cH_molL
    zeta = (0.010 + 0.025 * math.log10(C)) * (pH - 2.0) / 5.0
    return {"cH_molL": cH_molL, "pH": pH, "C_molL": C, "zeta": zeta}


def dynamic_coefficients(phi: float, k0: float, alpha_inf: float, cH: float,
                         omega: float, cfg: SEConfig,
                         C_override_molL: float | None = None,
                         sigma_f_override: float | None = None) -> Dict[str, complex]:
    """Compute k(omega), L(omega), sigma(omega), epsbar following Appendix A."""
    phi = float(np.clip(phi, cfg.phi_min, cfg.phi_max_valid))
    k0 = max(float(k0), cfg.k0_min)
    alpha_inf = max(float(alpha_inf), 1.0 + 1e-6)
    M = cfg.M_similarity
    eta = cfg.eta
    rho_f = cfg.rho_f

    ec = electrochemistry_from_h(cH, cfg, C_override_molL=C_override_molL)
    C = ec["C_molL"]
    zeta = ec["zeta"]
    T = cfg.temperature
    eps0 = cfg.eps0
    eps_f = cfg.eps_f
    e = cfg.e_charge
    kB = cfg.k_B
    NA = cfg.N_A
    z_vals = np.array([cfg.z1, cfg.z2], dtype=float)
    b_vals = np.array([cfg.b1, cfg.b2], dtype=float)

    # Convert mol/L to ion number density, m^-3. Binary symmetric electrolyte: N1=N2.
    N_each = C * 1000.0 * NA
    N_vals = np.array([N_each, N_each], dtype=float)

    # Appendix A, Eq. A2-A3.
    omega_t = phi * eta / (alpha_inf * k0 * rho_f)
    Lambda = math.sqrt(max(8.0 * alpha_inf * k0 / (phi * M), cfg.eps_complex))

    # Appendix A, Eq. A1.
    k_dyn = k0 / (np.sqrt(1.0 + 1j * (omega / omega_t) * M / 2.0) + 1j * (omega / omega_t))

    # Appendix A, Eq. A11.
    denom_debye = np.sum((e * z_vals)**2 * N_vals) / (eps0 * eps_f * kB * T)
    d = math.sqrt(1.0 / max(float(np.real(denom_debye)), cfg.eps_complex))

    # Appendix A, Eq. A4.
    one_minus = 1.0 - 2.0 * d / Lambda
    L = - (phi / alpha_inf) * (eps0 * eps_f * zeta / eta) * one_minus
    L *= (1.0 + 2j * omega / (M * omega_t) * one_minus**2 *
          (1.0 + d * np.sqrt(1j * omega * rho_f / eta))**2) ** (-0.5)

    # Appendix A, Eq. A7-A10.
    sigma_f = float(np.sum((e * z_vals)**2 * b_vals * N_vals))
    if sigma_f_override is not None and np.isfinite(sigma_f_override) and sigma_f_override > 0:
        # Prefer measured/calibrated pore-fluid conductivity when available.
        sigma_f = float(sigma_f_override)

    if abs(zeta) < 1e-12 or sigma_f <= 0:
        C_em = 0.0
        C_os = 0.0 + 0j
        P_os = 0.0
    else:
        exp_terms = np.exp(-(e * z_vals * zeta) / (2.0 * kB * T)) - 1.0
        C_em = float(2.0 * d * np.sum((e * z_vals)**2 * b_vals * N_vals * exp_terms))
        P_os = float((8.0 * kB * T * d**2) / (eps0 * eps_f * zeta**2) * np.sum(N_vals * exp_terms))
        if abs(P_os) < 1e-30:
            C_os = 0.0 + 0j
        else:
            C_os = ((eps0 * eps_f)**2 * zeta**2 / (2.0 * d * eta)) * P_os * \
                   (1.0 + (2.0 / P_os) * d * np.sqrt(1j * omega * rho_f / eta)) ** (-1)

    sigma = (phi * sigma_f / alpha_inf) * (1.0 + 2.0 * (C_em + C_os) / max(sigma_f * Lambda, cfg.eps_complex))

    # Bulk electrical permittivity in text below Eq. (7): eps = eps0[phi(eps_f-eps_s)/alpha_inf + eps_s]
    eps_bulk = eps0 * (phi * (eps_f - cfg.eps_s) / alpha_inf + cfg.eps_s)

    # Eq. (20): effective electrical permittivity.
    eps_bar = eps_bulk - 1j * sigma / omega + 1j * eta * L**2 / (omega * k_dyn)

    return {
        **ec,
        "omega_t": omega_t,
        "Lambda": Lambda,
        "debye_d": d,
        "k_dyn": k_dyn,
        "L": L,
        "sigma_f": sigma_f,
        "C_em": C_em,
        "C_os": C_os,
        "sigma": sigma,
        "eps_bulk": eps_bulk,
        "eps_bar": eps_bar,
    }


# -----------------------------
# 4. Schakel boundary-value solver
# -----------------------------

def biot_elastic_coefficients(phi: float, cfg: SEConfig) -> Tuple[float, float, float, float]:
    """Schakel Eq. (8)-(10): A, Q, R, P=A+2G."""
    phi = float(np.clip(phi, cfg.phi_min, cfg.phi_max_valid))
    Ks, Kf, Kb, G = cfg.K_s, cfg.K_f, cfg.K_b, cfg.G
    denom = Kf * (1.0 - phi - Kb / Ks) + phi * Ks
    A = (((1.0 - phi)**2 * Ks * Kf - (1.0 - phi) * Kb * Kf + phi * Ks * Kb) / denom
         - 2.0 * G / 3.0)
    Q = phi * (Ks * (1.0 - phi) - Kb) * Kf / denom
    R = phi**2 * Ks * Kf / denom
    P = A + 2.0 * G
    return A, Q, R, P


def wave_slownesses(phi: float, k0: float, alpha_inf: float, cH: float,
                    omega: float, cfg: SEConfig,
                    C_override_molL: float | None = None,
                    sigma_f_override: float | None = None) -> Dict[str, complex]:
    """Compute complex slownesses and amplitude ratios Eq. (24)-(39)."""
    phi = float(np.clip(phi, cfg.phi_min, cfg.phi_max_valid))
    dyn = dynamic_coefficients(phi, k0, alpha_inf, cH, omega, cfg,
                               C_override_molL=C_override_molL,
                               sigma_f_override=sigma_f_override)
    k_dyn = dyn["k_dyn"]
    L = dyn["L"]
    eps_bar = dyn["eps_bar"]
    eta = cfg.eta

    A, Q, R, P = biot_elastic_coefficients(phi, cfg)
    G = cfg.G
    rho_f, rho_s = cfg.rho_f, cfg.rho_s

    # Eq. (11)-(13).
    # Schakel Eq. (12): rho12 = phi*rho_f [1 + i phi*eta/(omega*rho_f*k(omega))].
    rho12 = phi * rho_f * (1.0 + 1j * phi * eta / (omega * rho_f * k_dyn))
    rho11 = (1.0 - phi) * rho_s - rho12
    rho22 = phi * rho_f - rho12

    # Eq. (26)-(29).
    # Schakel Eq. (29): E_K includes phi^2.
    E_K = eta**2 * phi**2 * L**2 / (k_dyn**2 * eps_bar * omega**2)
    rb11 = rho11 - E_K
    rb12 = rho12 + E_K
    rb22 = rho22 - E_K

    # Longitudinal roots Eq. (24)-(25).
    d0 = rb11 * rb22 - rb12**2
    d1 = -(P * rb22 + R * rb11 - 2.0 * Q * rb12)
    d2 = P * R - Q**2
    disc = (d1 / d2)**2 - 4.0 * d0 / d2
    roots_long = [(-d1 / d2 + np.sqrt(disc)) / 2.0, (-d1 / d2 - np.sqrt(disc)) / 2.0]
    roots_long = sorted(roots_long, key=lambda x: abs(x))
    s2_Pf, s2_Ps = roots_long[0], roots_long[1]

    # Transversal roots Eq. (24), with Eq. (30).
    mu = cfg.mu0
    d0_t = mu * eps_bar * (rb11 * rb22 - rb12**2) / G
    d1_t = -mu * eps_bar * rb22 - (rho11 * rho22 - rho12**2) / G
    d2_t = rho22
    disc_t = (d1_t / d2_t)**2 - 4.0 * d0_t / d2_t
    roots_trans = [(-d1_t / d2_t + np.sqrt(disc_t)) / 2.0,
                   (-d1_t / d2_t - np.sqrt(disc_t)) / 2.0]
    roots_trans = sorted(roots_trans, key=lambda x: abs(x))
    s2_TM, s2_SV = roots_trans[0], roots_trans[1]

    # Eq. (36)-(39).
    beta_Pf = (rb11 - P * s2_Pf) / (Q * s2_Pf - rb12)
    beta_Ps = (rb11 - P * s2_Ps) / (Q * s2_Ps - rb12)
    beta_TM = (G * s2_TM - (1.0 - phi) * rho_s) / (phi * rho_f)
    beta_SV = (G * s2_SV - (1.0 - phi) * rho_s) / (phi * rho_f)

    # Schakel Eqs. (38)-(39): alpha includes phi in eta*phi*L.
    alpha_Pf = eta * phi * L / (k_dyn * eps_bar) * (1.0 - beta_Pf)
    alpha_Ps = eta * phi * L / (k_dyn * eps_bar) * (1.0 - beta_Ps)
    alpha_TM = (mu * eta * phi * L) / (k_dyn * (mu * eps_bar - s2_TM)) * (1.0 - beta_TM)
    alpha_SV = (mu * eta * phi * L) / (k_dyn * (mu * eps_bar - s2_SV)) * (1.0 - beta_SV)

    return {
        **dyn,
        "A": A, "Q": Q, "R": R, "P": P,
        "rho11": rho11, "rho12": rho12, "rho22": rho22,
        "rb11": rb11, "rb12": rb12, "rb22": rb22,
        "s2_Pf": s2_Pf, "s2_Ps": s2_Ps, "s2_TM": s2_TM, "s2_SV": s2_SV,
        "beta_Pf": beta_Pf, "beta_Ps": beta_Ps, "beta_TM": beta_TM, "beta_SV": beta_SV,
        "alpha_Pf": alpha_Pf, "alpha_Ps": alpha_Ps, "alpha_TM": alpha_TM, "alpha_SV": alpha_SV,
    }


def se_coefficients(phi: float, k0_m2: float, alpha_inf: float, cH: float,
                    omega: float, theta_deg: float | None,
                    cfg: SEConfig, kx_override: float | None = None,
                    C_override_molL: float | None = None,
                    sigma_f_override: float | None = None) -> Dict[str, complex]:
    """Solve Schakel Appendix B 6x6 system and return R_E, T_TM, etc."""
    phi = float(np.clip(phi, cfg.phi_min, cfg.phi_max_valid))
    state = wave_slownesses(phi, k0_m2, alpha_inf, cH, omega, cfg,
                            C_override_molL=C_override_molL,
                            sigma_f_override=sigma_f_override)

    c_fl = math.sqrt(cfg.K_fl / cfg.rho_fl)
    if kx_override is None:
        theta = math.radians(theta_deg if theta_deg is not None else cfg.coeff_theta_deg)
        k1 = omega / c_fl * math.sin(theta)
    else:
        k1 = float(kx_override)

    k3_fl = complex_sqrt_branch((omega / c_fl)**2 - k1**2)

    s2_E = cfg.mu0 * cfg.eps0 * cfg.eps_fl - 1j * cfg.mu0 * _upper_fluid_sigma(state, cfg) / omega
    k3_E = complex_sqrt_branch(omega**2 * s2_E - k1**2)
    k3_Pf = complex_sqrt_branch(omega**2 * state["s2_Pf"] - k1**2)
    k3_Ps = complex_sqrt_branch(omega**2 * state["s2_Ps"] - k1**2)
    k3_TM = complex_sqrt_branch(omega**2 * state["s2_TM"] - k1**2)
    k3_SV = complex_sqrt_branch(omega**2 * state["s2_SV"] - k1**2)

    Q, R, P, G = state["Q"], state["R"], state["P"], cfg.G
    beta_Pf, beta_Ps = state["beta_Pf"], state["beta_Ps"]
    beta_TM, beta_SV = state["beta_TM"], state["beta_SV"]
    alpha_Pf, alpha_Ps = state["alpha_Pf"], state["alpha_Ps"]
    alpha_TM, alpha_SV = state["alpha_TM"], state["alpha_SV"]
    s2_Pf, s2_Ps, s2_TM, s2_SV = state["s2_Pf"], state["s2_Ps"], state["s2_TM"], state["s2_SV"]

    N1 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Pf
    N2 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Ps

    A_mat = np.zeros((6, 6), dtype=complex)

    # Appendix B, Eq. B2.
    A_mat[0, 0] = 0.0
    A_mat[0, 1] = k3_fl
    A_mat[0, 2] = k3_Pf * (1.0 - phi + phi * beta_Pf)
    A_mat[0, 3] = k3_Ps * (1.0 - phi + phi * beta_Ps)
    A_mat[0, 4] = k1 * (1.0 - phi + phi * beta_TM)
    A_mat[0, 5] = k1 * (1.0 - phi + phi * beta_SV)

    # Eq. B3.
    A_mat[1, 0] = 0.0
    A_mat[1, 1] = -phi * cfg.rho_fl
    A_mat[1, 2] = (Q + R * beta_Pf) * s2_Pf
    A_mat[1, 3] = (Q + R * beta_Ps) * s2_Ps
    A_mat[1, 4] = 0.0
    A_mat[1, 5] = 0.0

    # Eq. B4.
    A_mat[2, 0] = 0.0
    A_mat[2, 1] = 0.0
    A_mat[2, 2] = k1 * k3_Pf
    A_mat[2, 3] = k1 * k3_Ps
    A_mat[2, 4] = k1**2 - 0.5 * omega**2 * s2_TM
    A_mat[2, 5] = k1**2 - 0.5 * omega**2 * s2_SV

    # Eq. B5.
    A_mat[3, 0] = 0.0
    A_mat[3, 1] = 0.0
    A_mat[3, 2] = k1**2 - omega**2 * s2_Pf * N1 / (2.0 * G)
    A_mat[3, 3] = k1**2 - omega**2 * s2_Ps * N2 / (2.0 * G)
    A_mat[3, 4] = -k1 * k3_TM
    A_mat[3, 5] = -k1 * k3_SV

    # Eq. B6.
    A_mat[4, 0] = -s2_E / cfg.mu0
    A_mat[4, 1] = 0.0
    A_mat[4, 2] = 0.0
    A_mat[4, 3] = 0.0
    A_mat[4, 4] = alpha_TM * s2_TM / cfg.mu0
    A_mat[4, 5] = alpha_SV * s2_SV / cfg.mu0

    # Eq. B7.
    A_mat[5, 0] = -k3_E
    A_mat[5, 1] = 0.0
    A_mat[5, 2] = k1 * alpha_Pf
    A_mat[5, 3] = k1 * alpha_Ps
    A_mat[5, 4] = -k3_TM * alpha_TM
    A_mat[5, 5] = -k3_SV * alpha_SV

    b_vec = np.array([k3_fl, phi * cfg.rho_fl, 0.0, 0.0, 0.0, 0.0], dtype=complex)

    try:
        x = np.linalg.solve(A_mat, b_vec)
    except np.linalg.LinAlgError:
        x = np.full(6, np.nan + 1j * np.nan, dtype=complex)

    return {
        **state,
        "k1": k1,
        "k3_fl": k3_fl,
        "k3_E": k3_E,
        "k3_Pf": k3_Pf,
        "k3_Ps": k3_Ps,
        "k3_TM": k3_TM,
        "k3_SV": k3_SV,
        "s2_E": s2_E,
        "R_E": x[0],
        "R_M": x[1],
        "T_Pf": x[2],
        "T_Ps": x[3],
        "T_TM": x[4],
        "T_SV": x[5],
        "matrix_cond": np.linalg.cond(A_mat) if np.all(np.isfinite(A_mat)) else np.nan,
    }


def _upper_fluid_sigma(state: Dict[str, complex], cfg: SEConfig) -> float:
    if cfg.upper_fluid_conductivity_mode == "dynamic_pore_fluid":
        # Use real part of pore-fluid conductivity estimate as first-order upper fluid conductivity.
        return max(float(np.real(state.get("sigma_f", cfg.sigma_fl_default))), cfg.sigma_fl_default)
    return cfg.sigma_fl_default


# -----------------------------
# 5. Time-series processing
# -----------------------------

def optional_float(row: pd.Series, name: str) -> float | None:
    if name in row.index:
        try:
            val = float(row[name])
            return val if np.isfinite(val) else None
        except Exception:
            return None
    return None


def compute_time_series(df: pd.DataFrame, cfg: SEConfig) -> pd.DataFrame:
    omega0 = 2.0 * math.pi * cfg.f0
    rows = []
    for _, r in df.iterrows():
        phi_raw = float(r["Porosity"])
        valid = (cfg.phi_min < phi_raw < cfg.phi_max_valid)
        phi = float(np.clip(phi_raw, cfg.phi_min, cfg.phi_max_valid))
        k0_m2 = max(float(r["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
        tau = max(float(r["Tortuosity"]), 1.0 + 1e-6)
        cH = float(r["OutletHConc"])
        C_override = optional_float(r, "ElectrolyteConcentration_molL")
        sigma_f_override = optional_float(r, "FluidConductivity_S_m")
        ec = electrochemistry_from_h(cH, cfg, C_override_molL=C_override)
        try:
            coeff = se_coefficients(phi, k0_m2, tau, cH, omega0, cfg.coeff_theta_deg, cfg,
                                    C_override_molL=C_override,
                                    sigma_f_override=sigma_f_override)
            RE = coeff["R_E"]
            TTM = coeff["T_TM"]
            Lc = coeff["L"]
            sig = coeff["sigma"]
            wt = coeff["omega_t"]
            lam = coeff["Lambda"]
            cond = coeff["matrix_cond"]
        except Exception as e:
            RE = TTM = Lc = sig = np.nan + 1j * np.nan
            wt = lam = cond = np.nan
        rows.append({
            "Time_s": r["Time_s"],
            "Time_min": r["Time_s"] / 60.0,
            "Porosity_raw": phi_raw,
            "Porosity_used": phi,
            "valid_poroelastic": valid,
            "Permeability_mD": r["Permeability_mD"],
            "k0_m2": k0_m2,
            "Tortuosity": tau,
            "OutletHConc_raw": cH,
            "ElectrolyteConcentration_input_molL": C_override,
            "FluidConductivity_input_S_m": sigma_f_override,
            **ec,
            "omega0_rad_s": omega0,
            "theta_deg": cfg.coeff_theta_deg,
            "omega_t": wt,
            "Lambda_m": lam,
            "L_abs": abs(Lc),
            "sigma_abs": abs(sig),
            "RE_real": np.real(RE),
            "RE_imag": np.imag(RE),
            "RE_abs": abs(RE),
            "TTM_real": np.real(TTM),
            "TTM_imag": np.imag(TTM),
            "TTM_abs": abs(TTM),
            "matrix_cond": cond,
        })
    out = pd.DataFrame(rows)

    # Normalize coefficients by the first valid value for clearer monitoring curves.
    for col in ["RE_abs", "TTM_abs", "L_abs", "sigma_abs"]:
        valid_vals = out.loc[out["valid_poroelastic"] & np.isfinite(out[col]), col]
        ref = valid_vals.iloc[0] if len(valid_vals) else np.nan
        out[col + "_norm"] = out[col] / ref if ref and np.isfinite(ref) and ref != 0 else np.nan
    return out


# -----------------------------
# 6. Waveform synthesis
# -----------------------------

def choose_snapshot(df: pd.DataFrame, target_phi: float = 0.75) -> int:
    vals = np.abs(df["Porosity"].to_numpy(float) - target_phi)
    return int(np.nanargmin(vals))


def causal_ricker_source_spectrum(omega: np.ndarray | float, cfg: SEConfig,
                                  n_time: int = 2048) -> np.ndarray | complex:
    """Complex S(omega) for Liu Eq. (1)-(2) from a causal source wavelet.

    Liu's A(k,omega) contains the source spectrum S(omega).  Using only the
    real amplitude spectrum creates a zero-phase pulse centered around T0.  Here
    S(omega) is computed from a source-time wavelet that starts at tau=0, so the
    propagation phase exp(i k_z z_s) naturally moves the response to T0 without
    output-side clipping.
    """
    omega_arr = np.atleast_1d(np.asarray(omega, dtype=float))
    f0 = max(float(cfg.f0), 1.0)
    duration = max(float(cfg.source_duration_cycles), float(cfg.source_peak_cycles) + 2.0) / f0
    tau = np.linspace(0.0, duration, int(max(n_time, 64)))
    peak_t = max(0.0, float(cfg.source_peak_cycles)) / f0
    src = ricker(tau - peak_t, f0)
    # Smoothly start the emitted pulse at zero; this is source definition, not
    # receiver-side gating.  It avoids injecting a finite value at tau=0.
    ramp_len = max(4, min(len(tau), int(0.25 / f0 / max(tau[1] - tau[0], 1e-15))))
    ramp = np.ones_like(src)
    if ramp_len > 1:
        ramp[:ramp_len] = 0.5 * (1.0 - np.cos(np.linspace(0.0, math.pi, ramp_len)))
    src *= ramp
    integrate = getattr(np, "trapezoid", np.trapz)
    spec = integrate(src[None, :] * np.exp(1j * omega_arr[:, None] * tau[None, :]), tau, axis=1)
    peak = np.nanmax(np.abs(spec))
    if np.isfinite(peak) and peak > 0:
        spec = spec / peak
    if np.isscalar(omega):
        return complex(spec[0])
    return spec


def _trapz_weights(n: int, dx: float) -> np.ndarray:
    """Trapezoidal-rule weights for an equally spaced grid."""
    w = np.full(n, dx, dtype=float)
    if n > 1:
        w[0] *= 0.5
        w[-1] *= 0.5
    return w


def liu_interface_coefficient_kx(kx: float) -> float:
    """Horizontal wavenumber used to evaluate Liu's interface coefficients.

    Liu Eq. (1)-(2) use signed horizontal wavenumber in the spatial phase kx,
    but the isotropic plane-interface conversion coefficients are functions of
    incidence-angle magnitude. Keeping the sign inside Schakel's vector-potential
    convention makes R_E and T_TM odd in k and incorrectly turns the finite-offset
    scalar-potential response into a sin(kD)-type kernel.
    """
    return abs(float(kx))


def liu_electrical_potential_coefficients(coeff: Dict[str, complex],
                                          kx: float) -> Tuple[complex, complex]:
    """Convert Schakel potentials to Liu's reflected/transmitted electrical potentials.

    Schakel Eq. (45) defines R_E as the reflected EM vector-potential coefficient
    and T_TM as the transmitted solid TM vector-potential coefficient. Schakel
    Eq. (39) gives the transmitted electric vector potential as alpha_TM*T_TM.
    For a TM vector potential Psi_E=(0, psi_E, 0), Eq. (23) gives
    E=curl(Psi_E). Relating its horizontal component to an electrical potential
    through E_x=-du/dx produces opposite propagation-direction signs:

        R_u = -(k3_E/|kx|) R_E
        T_u = +(k3_TM/|kx|) alpha_TM T_TM

    At exact normal incidence Schakel's raw coefficients vanish while these
    ratios have finite limits. The caller must evaluate that limit at a small
    positive wavenumber rather than pass kx=0.
    """
    k_abs = abs(float(kx))
    if k_abs <= 1e-14:
        raise ValueError("kx must be nonzero; evaluate the normal-incidence potential limit first")
    re_potential = -(coeff["k3_E"] / k_abs) * coeff["R_E"]
    te_potential = (coeff["k3_TM"] / k_abs) * coeff["alpha_TM"] * coeff["T_TM"]
    return re_potential, te_potential


def synthesize_waveforms_spectral(row: pd.Series, cfg: SEConfig,
                                  n_omega: int | None = None,
                                  n_k: int | None = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Liu et al. (2018)-style frequency--wavenumber synthesis of interface EM waveforms.

    The implementation follows the structure of Liu et al. Eq. (1)--(2):

        u_r(x,z,t) = ∬ A(k,ω) R_u(k,ω)
                       exp[i(k_z^i z_s - k_z^{Er} z + k x - ωt)] dk dω

        u_t(x,z,t) = ∬ A(k,ω) T_u(k,ω)
                       exp[i(k_z^i z_s + k_z^{Et} z + k x - ωt)] dk dω

    where R_u and T_u are Liu electrical-potential coefficients obtained from the
    Schakel & Smeulders 6x6 solution for each (k,ω). It preserves frequency
    dependence, horizontal-wavenumber/incident-angle effects, and propagation
    phase.

    Notes
    -----
    - Only positive frequencies are explicitly integrated; the real waveform is
      obtained by adding the complex-conjugate contribution via a factor of 2.
    - The acoustic k integration is restricted to propagating incident-fluid
      components |k| <= spectral_k_limit_factor * ω/V_f. This avoids making the
      source itself strongly evanescent; the EM fields can still be diffusive or
      evanescent through their complex vertical wavenumbers.
    - The returned amplitudes are in arbitrary source-scaled units. The plotting
      routine normalizes traces for display but keeps the actual maximum amplitude
      as Amax.
    """
    n_omega = int(n_omega or cfg.spectral_n_omega)
    n_k = int(n_k or cfg.spectral_n_k)
    if n_omega < 2:
        raise ValueError("n_omega must be >= 2 for spectral synthesis")
    if n_k < 3:
        raise ValueError("n_k must be >= 3 for spectral synthesis")

    z_receivers = np.arange(cfg.receiver_z_min, cfg.receiver_z_max + 0.5 * cfg.receiver_spacing, cfg.receiver_spacing)
    Vf = math.sqrt(cfg.K_fl / cfg.rho_fl)
    T0 = cfg.z_s / Vf
    t = np.linspace(max(0.0, T0 - cfg.waveform_t_before), T0 + cfg.waveform_t_after, cfg.waveform_nt)
    x = float(cfg.offset_D)

    phi = float(np.clip(row["Porosity"], cfg.phi_min, cfg.phi_max_valid))
    k0_m2 = max(float(row["Permeability_mD"]) * 9.869233e-16, cfg.k0_min)
    tau = max(float(row["Tortuosity"]), 1.0 + 1e-6)
    cH = float(row["OutletHConc"])
    C_override = optional_float(row, "ElectrolyteConcentration_molL")
    sigma_f_override = optional_float(row, "FluidConductivity_S_m")

    f_min = max(1.0, cfg.spectral_f_min_factor * cfg.f0)
    f_max = max(f_min * 1.01, cfg.spectral_f_max_factor * cfg.f0)
    omega_grid = 2.0 * math.pi * np.linspace(f_min, f_max, n_omega)
    domega = float(omega_grid[1] - omega_grid[0])
    omega_weights = _trapz_weights(n_omega, domega)

    theta0 = math.radians(cfg.source_beam_theta_deg)
    kb = max(float(cfg.source_kb_m_inv), 1e-12)
    U = np.zeros((len(z_receivers), len(t)), dtype=complex)

    # Integral normalization. Liu's published formula omits explicit 2π convention details;
    # this factor gives a standard inverse-transform scaling but the source amplitude remains arbitrary.
    inv_2pi2 = 1.0 / (2.0 * math.pi)**2

    source_spectrum = causal_ricker_source_spectrum(omega_grid, cfg)
    fluid_mask = z_receivers < 0
    porous_mask = ~fluid_mask

    for iw, omega in enumerate(omega_grid):
        k_ac = omega / Vf
        k_center = k_ac * math.sin(theta0)
        k_lim = max(1e-12, cfg.spectral_k_limit_factor * k_ac)
        k_grid = np.linspace(-k_lim, k_lim, n_k)
        dk = float(k_grid[1] - k_grid[0])
        k_weights = _trapz_weights(n_k, dk)
        S_omega = source_spectrum[iw]
        phase_t = np.exp(-1j * omega * t)
        coeff_cache: Dict[float, Dict[str, complex]] = {}

        for ik, kx in enumerate(k_grid):
            # Liu Eq. (1)-(2) finite-width source spectrum:
            # A(k,omega)=exp[-(k-k0)^2/k_b^2] S(omega).
            Akw = cfg.source_pressure_amp * S_omega * np.exp(-((kx - k_center) / kb)**2)
            if not np.isfinite(Akw) or abs(Akw) < 1e-14:
                continue
            coeff_kx = liu_interface_coefficient_kx(float(kx))
            # R_E and T_TM vanish at exact normal incidence, while the Liu
            # electrical-potential ratios R_E/k and T_TM/k have finite limits.
            # Evaluate that single quadrature node just off normal incidence;
            # the source spectrum and Fourier phase still use the true kx=0.
            potential_kx = coeff_kx if coeff_kx > 0.0 else k_lim * 1e-8
            cache_key = round(potential_kx, 12)
            try:
                if cache_key not in coeff_cache:
                    coeff_cache[cache_key] = se_coefficients(
                        phi, k0_m2, tau, cH, omega, None, cfg,
                        kx_override=potential_kx,
                        C_override_molL=C_override,
                        sigma_f_override=sigma_f_override,
                    )
                coeff = coeff_cache[cache_key]
            except Exception:
                continue
            RE_potential, TE_potential = liu_electrical_potential_coefficients(coeff, potential_kx)
            if not (
                np.isfinite(RE_potential.real)
                and np.isfinite(RE_potential.imag)
                and np.isfinite(TE_potential.real)
                and np.isfinite(TE_potential.imag)
            ):
                continue

            kzi = coeff["k3_fl"]
            kEr = coeff["k3_E"]
            kEt = coeff["k3_TM"]
            weight = 2.0 * Akw * omega_weights[iw] * k_weights[ik] * inv_2pi2
            common_x = np.exp(1j * kx * x)

            if np.any(fluid_mask):
                z_fluid = z_receivers[fluid_mask]
                amp_fluid = RE_potential * np.exp(1j * (kzi * cfg.z_s - kEr * z_fluid)) * common_x
                U[fluid_mask, :] += weight * amp_fluid[:, None] * phase_t[None, :]
            if np.any(porous_mask):
                z_porous = z_receivers[porous_mask]
                amp_porous = TE_potential * np.exp(1j * (kzi * cfg.z_s + kEt * z_porous)) * common_x
                U[porous_mask, :] += weight * amp_porous[:, None] * phase_t[None, :]

    return z_receivers, t, np.real(U)


# -----------------------------
# 7. Plotting
# -----------------------------

def plot_coefficients(ts: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    valid = ts["valid_poroelastic"].to_numpy(bool)
    x = ts["Time_min"].to_numpy(float)
    ax.plot(x[valid], ts.loc[valid, "RE_abs_norm"], marker="o", linewidth=1.6, label=r"$|R_E|/|R_E(t_0)|$")
    ax.plot(x[valid], ts.loc[valid, "TTM_abs_norm"], marker="s", linewidth=1.6, label=r"$|T_{TM}|/|T_{TM}(t_0)|$")
    ax.set_xlabel("Dissolution time (min)")
    ax.set_ylabel("Normalized coefficient magnitude")
    ax.set_title("Seismoelectric coefficients during calcite dissolution")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "coefficients_vs_dissolution_time.png", dpi=300)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(7.6, 4.8))
    ax1.plot(x[valid], ts.loc[valid, "L_abs_norm"], marker="o", linewidth=1.6, label=r"$|L(\omega)|$ norm")
    ax1.plot(x[valid], ts.loc[valid, "sigma_abs_norm"], marker="s", linewidth=1.6, label=r"$|\sigma(\omega)|$ norm")
    ax1.set_xlabel("Dissolution time (min)")
    ax1.set_ylabel("Normalized dynamic coefficient")
    ax1.set_title("Dynamic electrokinetic coefficients")
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    fig.tight_layout()
    fig.savefig(outdir / "dynamic_coefficients_vs_dissolution_time.png", dpi=300)
    plt.close(fig)


def save_waveform_arrays(z: np.ndarray, t: np.ndarray, U: np.ndarray, outdir: Path, name: str) -> None:
    """Save waveform gather both as compressed NumPy arrays and as a wide CSV table."""
    np.savez_compressed(outdir / f"{name}.npz", z_m=z, t_s=t, U=U)
    data = {"z_m": z}
    # A wide CSV is convenient for quick inspection in Excel; keep time labels compact.
    for j, tj in enumerate(t):
        data[f"t_{tj*1e6:.4f}_us"] = U[:, j]
    pd.DataFrame(data).to_csv(outdir / f"{name}.csv", index=False)


def waveform_spatial_peak_diagnostics(z: np.ndarray, t: np.ndarray, U: np.ndarray,
                                      offset_D: float) -> pd.DataFrame:
    """Return per-receiver peak diagnostics for auditing Liu Fig. 2-style gathers."""
    z = np.asarray(z, dtype=float)
    t = np.asarray(t, dtype=float)
    U = np.asarray(U, dtype=float)
    if U.shape != (len(z), len(t)):
        raise ValueError("U must have shape (len(z), len(t)).")

    peak_idx = np.nanargmax(np.abs(U), axis=1)
    peak_signed = U[np.arange(len(z)), peak_idx]
    peak_abs = np.abs(peak_signed)
    interface_mask = np.isclose(z, 0.0, atol=1e-12, rtol=0.0)
    side = np.where(interface_mask, "interface", np.where(z < 0.0, "R_E", "T_E"))
    ref_distance = abs(float(offset_D)) / math.sqrt(2.0) if np.isfinite(offset_D) else np.nan
    return pd.DataFrame({
        "z_m": z,
        "z_mm": z * 1e3,
        "side": side,
        "distance_from_interface_m": np.abs(z),
        "distance_from_interface_mm": np.abs(z) * 1e3,
        "peak_abs": peak_abs,
        "peak_signed": peak_signed,
        "peak_time_s": t[peak_idx],
        "peak_time_us": t[peak_idx] * 1e6,
        "liu_dipole_peak_distance_m": ref_distance,
        "liu_dipole_peak_distance_mm": ref_distance * 1e3 if np.isfinite(ref_distance) else np.nan,
    })


def save_waveform_spatial_peak_diagnostics(z: np.ndarray, t: np.ndarray, U: np.ndarray,
                                           cfg: SEConfig, outdir: Path,
                                           name: str = "waveform_spatial_peak_diagnostics") -> pd.DataFrame:
    diag = waveform_spatial_peak_diagnostics(z, t, U, cfg.offset_D)
    diag.to_csv(outdir / f"{name}.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for side, color, label in [("R_E", "tab:blue", r"$R_E$ side"), ("T_E", "tab:red", r"$T_E$ side")]:
        part = diag[diag["side"] == side]
        if part.empty:
            continue
        vals = part["peak_abs"].to_numpy(float)
        vmax = np.nanmax(vals)
        vals_norm = vals / vmax if np.isfinite(vmax) and vmax > 0 else vals
        ax.plot(part["distance_from_interface_mm"], vals_norm, marker=".", linestyle="-",
                linewidth=1.0, markersize=3.0, color=color, label=label)
    ref = abs(cfg.offset_D) / math.sqrt(2.0) * 1e3
    ax.axvline(ref, color="0.35", linestyle=":", linewidth=1.2, label=r"$D/\sqrt{2}$")
    ax.set_xlabel("Distance from interface (mm)")
    ax.set_ylabel("Side-normalized peak amplitude")
    ax.set_title("Spatial peak audit of spectral waveform gather")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300)
    plt.close(fig)
    return diag


def waveform_display_indices(z: np.ndarray, cfg: SEConfig) -> np.ndarray:
    """Return receiver indices drawn in a waveform gather."""
    z = np.asarray(z, dtype=float)
    display_stride = max(1, int(round(5.0e-3 / max(float(cfg.receiver_spacing), 1e-12))))
    display_idx = np.arange(0, len(z), display_stride)
    if len(z) - 1 not in display_idx:
        display_idx = np.append(display_idx, len(z) - 1)
    if abs(float(cfg.offset_D)) > 1e-12:
        interface_mask = np.isclose(z[display_idx], 0.0, atol=1e-12, rtol=0.0)
        display_idx = display_idx[~interface_mask]
    return display_idx


def plot_waveform_gather(z: np.ndarray, t: np.ndarray, U: np.ndarray,
                         row: pd.Series, cfg: SEConfig, outdir: Path,
                         name: str = "waveform_snapshot") -> float:
    # Normalize traces for display; keep actual max separately.
    # The reflected fluid-side R_E and transmitted porous-side T_E electrical
    # potentials can have different numerical scales, so
    # normalize each side independently for a readable Liu-style gather.
    Amax = float(np.nanmax(np.abs(U)))
    U_plot = np.zeros_like(U, dtype=float)
    for mask in (z < 0, z > 0):
        if not np.any(mask):
            continue
        side_max = float(np.nanmax(np.abs(U[mask, :])))
        if side_max > 0 and np.isfinite(side_max):
            U_plot[mask, :] = U[mask, :] / side_max

    # Use mm for Liu-style laboratory-scale receiver positions; use µm for pore-scale views.
    if np.nanmax(np.abs(z)) >= 5e-3:
        z_plot = z * 1e3
        z_unit = "mm"
    else:
        z_plot = z * 1e6
        z_unit = "µm"
    t_us = t * 1e6
    display_idx = waveform_display_indices(z, cfg)
    z_display = z_plot[display_idx]
    U_display = U_plot[display_idx, :]
    trace_spacing = abs(np.median(np.diff(z_display))) if len(z_display) > 1 else 1.0
    scale = 0.42 * trace_spacing

    fig, ax = plt.subplots(figsize=(7.8, 6.0))
    for zi, tr in zip(z_display, U_display):
        color = "tab:red" if zi > 0 else ("tab:blue" if abs(zi) < 1e-12 else "0.25")
        ax.plot(t_us, zi + scale * tr, color=color, linewidth=0.9)
    T0_us = cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl) * 1e6
    ax.axvline(T0_us, color="tab:blue", linestyle=":", linewidth=1.4)
    ax.axhline(0.0, color="tab:blue", linewidth=1.0)
    ax.text(T0_us, np.nanmax(z_plot), r" $T_0$", color="tab:blue", va="top")
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel(f"Electrode position relative to interface ({z_unit})")
    ax.set_title(f"Interface EM waveforms at t_d={row['Time_s']:.1f} s, phi={row['Porosity']:.3f} [R_E/T_E norm]")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / f"{name}.png", dpi=300)
    plt.close(fig)
    return Amax


def compute_peak_amplitude_spectral(ts: pd.DataFrame, df_raw: pd.DataFrame, cfg: SEConfig,
                                    n_omega: int | None = None,
                                    n_k: int | None = None) -> pd.DataFrame:
    """Compute Amax from the same Liu-style spectral waveform model used for plots."""
    out = ts.copy()
    vals = []
    for idx, r in df_raw.reset_index(drop=True).iterrows():
        if idx >= len(out) or not bool(out.loc[idx, "valid_poroelastic"]):
            vals.append(np.nan)
            continue
        try:
            _, _, U = synthesize_waveforms_spectral(r, cfg, n_omega=n_omega, n_k=n_k)
            vals.append(float(np.nanmax(np.abs(U))))
        except Exception:
            vals.append(np.nan)
    out["Amax_waveform_spectral"] = vals
    valid_vals = out.loc[out["valid_poroelastic"] & np.isfinite(out["Amax_waveform_spectral"]), "Amax_waveform_spectral"]
    ref = valid_vals.iloc[0] if len(valid_vals) else np.nan
    out["Amax_waveform_spectral_norm"] = out["Amax_waveform_spectral"] / ref if ref and np.isfinite(ref) and ref != 0 else np.nan
    return out


def plot_peak_amplitude(ts: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    valid = ts["valid_poroelastic"].to_numpy(bool)
    x = ts["Time_min"].to_numpy(float)
    ax.plot(x[valid], ts.loc[valid, "Amax_waveform_spectral_norm"], marker="o", linewidth=1.6)
    ax.set_xlabel("Dissolution time (min)")
    ax.set_ylabel("Normalized maximum waveform peak")
    ax.set_title("Maximum interface EM waveform peak during dissolution")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "peak_amplitude_vs_dissolution_time.png", dpi=300)
    plt.close(fig)


def save_parameter_table(cfg: SEConfig, outdir: Path) -> None:
    rows = []
    labels = {
        "Lx": "RTM domain length in x", "Ly": "RTM domain height in z",
        "grain_radius": "CaCO3 grain radius", "pore_throat_width": "designed pore-throat width",
        "z_s": "acoustic source distance to interface", "receiver_z_min": "minimum receiver/electrode position relative to interface",
        "receiver_z_max": "maximum receiver/electrode position relative to interface", "receiver_spacing": "receiver/electrode trace interval",
        "offset_D": "receiver-line offset", "f0": "source central frequency",
        "source_kb_m_inv": "Liu frequency-wavenumber beam-width k_b",
        "source_peak_cycles": "Ricker peak time after source onset",
        "source_duration_cycles": "finite causal source duration",
        "K_s": "calcite grain bulk modulus",
        "K_b": "drained framework bulk modulus", "G": "framework shear modulus",
        "K_f": "pore fluid bulk modulus", "K_fl": "upper fluid bulk modulus",
        "eta": "pore fluid viscosity", "rho_f": "pore fluid density", "rho_fl": "upper fluid density",
        "rho_s": "calcite density", "eps_f": "pore fluid relative permittivity",
        "eps_s": "solid relative permittivity", "eps_fl": "upper fluid relative permittivity",
        "C_background_molL": "background electrolyte concentration", "H_min_molL": "minimum H+ concentration for pH",
        "M_similarity": "similarity parameter M", "phi_max_valid": "maximum porosity for valid poroelastic model",
    }
    units = {
        "Lx": "m", "Ly": "m", "grain_radius": "m", "pore_throat_width": "m", "z_s": "m",
        "receiver_z_min": "m", "receiver_z_max": "m", "receiver_spacing": "m",
        "offset_D": "m", "f0": "Hz", "source_beam_theta_deg": "deg",
        "source_kb_m_inv": "1/m", "source_peak_cycles": "cycles",
        "source_duration_cycles": "cycles",
        "spectral_f_min_factor": "-", "spectral_f_max_factor": "-", "spectral_n_omega": "-",
        "spectral_n_k": "-", "spectral_k_limit_factor": "-", "K_s": "Pa", "K_b": "Pa", "G": "Pa", "K_f": "Pa", "K_fl": "Pa",
        "eta": "Pa s", "rho_f": "kg/m3", "rho_fl": "kg/m3", "rho_s": "kg/m3", "eps_f": "-",
        "eps_s": "-", "eps_fl": "-", "C_background_molL": "mol/L", "H_min_molL": "mol/L",
        "M_similarity": "-", "phi_max_valid": "-",
    }
    for k, v in asdict(cfg).items():
        if k in labels:
            rows.append({"parameter": k, "meaning": labels[k], "value": v, "unit": units.get(k, "")})
    pd.DataFrame(rows).to_csv(outdir / "parameters_used.csv", index=False)


# -----------------------------
# 8. Main
# -----------------------------

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=str,
        default=str(Path(__file__).with_name("global_evolution.xlsx")),
        help="Input reactive transport table (xlsx/csv/txt). Defaults to global_evolution.xlsx next to this script.",
    )
    parser.add_argument("--outdir", type=str, default="se_results_offset", help="Output directory for results and plots")
    parser.add_argument("--snapshot-target-phi", type=float, default=0.75)
    parser.add_argument("--spectral-n-omega", type=int, default=None, help="Number of positive-frequency samples in spectral mode")
    parser.add_argument("--spectral-n-k", type=int, default=None, help="Number of horizontal-wavenumber samples per frequency in spectral mode")
    parser.add_argument("--peak-spectral-n-omega", type=int, default=24, help="Coarser spectral frequency samples for Amax-vs-dissolution curve")
    parser.add_argument("--peak-spectral-n-k", type=int, default=101, help="Coarser spectral wavenumber samples for Amax-vs-dissolution curve")
    parser.add_argument("--spectral-f-min-factor", type=float, default=None, help="Minimum spectral frequency as a factor of f0")
    parser.add_argument("--spectral-f-max-factor", type=float, default=None, help="Maximum spectral frequency as a factor of f0")
    parser.add_argument("--source-beam-theta-deg", type=float, default=None, help="Central acoustic beam angle used in A(k,omega), degrees")
    parser.add_argument("--source-kb-m-inv", type=float, default=None, help="Liu finite-width acoustic source k_b, in 1/m")
    parser.add_argument("--source-peak-cycles", type=float, default=None, help="Ricker source peak time after source onset, in f0 cycles")
    parser.add_argument("--source-duration-cycles", type=float, default=None, help="Finite causal source duration used to compute S(omega), in f0 cycles")
    parser.add_argument("--z_s", type=float, default=None, help="Override source-interface distance in meters")
    parser.add_argument("--z_s_mm", type=float, default=None, help="Override source-interface distance in millimeters")
    parser.add_argument("--receiver-z-min-mm", type=float, default=None, help="Minimum receiver/electrode position relative to interface, in mm")
    parser.add_argument("--receiver-z-max-mm", type=float, default=None, help="Maximum receiver/electrode position relative to interface, in mm")
    parser.add_argument("--receiver-spacing-mm", type=float, default=None, help="Receiver/electrode trace interval, in mm")
    parser.add_argument("--offset-D-mm", type=float, default=None, help="Horizontal offset between source line and receiver line, in mm")
    parser.add_argument("--f0", type=float, default=None, help="Override central frequency in Hz")
    parser.add_argument("--upper-fluid-conductivity-mode", type=str, default=None, choices=["constant", "dynamic_pore_fluid"],
                        help="constant follows Schakel Table I sigma_fl; dynamic_pore_fluid is a model assumption")
    args = parser.parse_args()

    cfg = SEConfig()
    if args.z_s is not None:
        cfg.z_s = args.z_s
    if args.z_s_mm is not None:
        cfg.z_s = args.z_s_mm * 1e-3
    if args.receiver_z_min_mm is not None:
        cfg.receiver_z_min = args.receiver_z_min_mm * 1e-3
    if args.receiver_z_max_mm is not None:
        cfg.receiver_z_max = args.receiver_z_max_mm * 1e-3
    if args.receiver_spacing_mm is not None:
        cfg.receiver_spacing = args.receiver_spacing_mm * 1e-3
    if args.offset_D_mm is not None:
        cfg.offset_D = args.offset_D_mm * 1e-3
    if args.f0 is not None:
        cfg.f0 = args.f0
    if args.spectral_n_omega is not None:
        cfg.spectral_n_omega = args.spectral_n_omega
    if args.spectral_n_k is not None:
        cfg.spectral_n_k = args.spectral_n_k
    if args.spectral_f_min_factor is not None:
        cfg.spectral_f_min_factor = args.spectral_f_min_factor
    if args.spectral_f_max_factor is not None:
        cfg.spectral_f_max_factor = args.spectral_f_max_factor
    if args.source_beam_theta_deg is not None:
        cfg.source_beam_theta_deg = args.source_beam_theta_deg
    if args.source_kb_m_inv is not None:
        cfg.source_kb_m_inv = args.source_kb_m_inv
    if args.source_peak_cycles is not None:
        cfg.source_peak_cycles = args.source_peak_cycles
    if args.source_duration_cycles is not None:
        cfg.source_duration_cycles = args.source_duration_cycles
    if args.upper_fluid_conductivity_mode is not None:
        cfg.upper_fluid_conductivity_mode = args.upper_fluid_conductivity_mode

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_reactive_transport_table(args.input)
    ts = compute_time_series(df, cfg)
    ts = compute_peak_amplitude_spectral(ts, df, cfg,
                                         n_omega=args.peak_spectral_n_omega,
                                         n_k=args.peak_spectral_n_k)
    ts.to_csv(outdir / "seismoelectric_timeseries_results.csv", index=False)
    save_parameter_table(cfg, outdir)

    plot_coefficients(ts, outdir)
    plot_peak_amplitude(ts, outdir)

    idx = choose_snapshot(df, target_phi=args.snapshot_target_phi)
    row = df.iloc[idx]
    z, t, U = synthesize_waveforms_spectral(row, cfg, n_omega=cfg.spectral_n_omega, n_k=cfg.spectral_n_k)
    plot_name = "waveform_snapshot_spectral"
    save_waveform_arrays(z, t, U, outdir, plot_name)
    diag = save_waveform_spatial_peak_diagnostics(z, t, U, cfg, outdir)
    Amax_snapshot = plot_waveform_gather(z, t, U, row, cfg, outdir, plot_name)
    re_diag = diag[diag["side"] == "R_E"]
    te_diag = diag[diag["side"] == "T_E"]
    re_peak_distance_mm = float(re_diag.loc[re_diag["peak_abs"].idxmax(), "distance_from_interface_mm"]) if not re_diag.empty else np.nan
    te_peak_distance_mm = float(te_diag.loc[te_diag["peak_abs"].idxmax(), "distance_from_interface_mm"]) if not te_diag.empty else np.nan

    summary = {
        "input": str(args.input),
        "outdir": str(outdir),
        "snapshot_index": int(idx),
        "snapshot_Time_s": float(row["Time_s"]),
        "snapshot_Porosity": float(row["Porosity"]),
        "snapshot_Amax": Amax_snapshot,
        "T0_us": cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl) * 1e6,
        "waveform_mode": "spectral",
        "pre_T0_max_abs": float(np.nanmax(np.abs(U[:, t < cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)]))) if np.any(t < cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)) else 0.0,
        "post_T0_max_abs": float(np.nanmax(np.abs(U[:, t >= cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)]))) if np.any(t >= cfg.z_s / math.sqrt(cfg.K_fl / cfg.rho_fl)) else np.nan,
        "RE_peak_distance_from_interface_mm": re_peak_distance_mm,
        "TE_peak_distance_from_interface_mm": te_peak_distance_mm,
        "liu_dipole_reference_peak_distance_mm": abs(cfg.offset_D) / math.sqrt(2.0) * 1e3,
        "spectral_n_omega": cfg.spectral_n_omega,
        "spectral_n_k": cfg.spectral_n_k,
        "peak_spectral_n_omega": args.peak_spectral_n_omega,
        "peak_spectral_n_k": args.peak_spectral_n_k,
        "spectral_f_min_factor": cfg.spectral_f_min_factor,
        "spectral_f_max_factor": cfg.spectral_f_max_factor,
        "source_beam_theta_deg": cfg.source_beam_theta_deg,
        "source_kb_m_inv": cfg.source_kb_m_inv,
        "source_peak_cycles": cfg.source_peak_cycles,
        "source_duration_cycles": cfg.source_duration_cycles,
    }
    pd.Series(summary).to_csv(outdir / "run_summary.csv")

    print("Done. Outputs written to:", outdir)
    for f in sorted(outdir.glob("*")):
        print(" -", f.name)


if __name__ == "__main__":
    main()
