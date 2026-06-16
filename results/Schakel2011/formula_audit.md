# Formula Audit

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
