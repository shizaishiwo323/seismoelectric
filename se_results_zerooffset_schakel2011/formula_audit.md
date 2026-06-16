# Formula Audit

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
