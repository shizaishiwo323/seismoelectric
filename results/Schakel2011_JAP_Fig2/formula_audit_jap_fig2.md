# JAP Figure 2 Formula Audit

Implemented target:

- JAP Eq. (5) for fluid positions 1-5, using the same Sommerfeld path as the paper: real angles `0..pi/2` plus the `gamma` substitution `theta = pi/2 + i ln(sqrt(gamma^2+1)+gamma)`.
- JAP Eq. (8) for position 6 inside the porous sample. The bracket contains the front-interface TM term, the Pf coseismic term, and the Eq. (7) back/front interface-response terms.
- JAP Table II parameters are used. The Schakel dynamic-permeability helper computes `Lambda=sqrt(8 alpha_inf k0/(phi M))`; therefore `M=0.9999795` is chosen to reproduce the Table II value `Lambda=9.43e-06 m`.
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
