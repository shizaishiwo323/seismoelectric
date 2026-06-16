# Subagent Review Summary

## Post-Review Correction

After the first review, receiver position 6 was found to be visibly wrong: the synthetic trace clipped the paper-style axes and showed unrealistically large ringing. A systematic decomposition showed that the mismatch was not caused only by the back-interface Pf-incidence reconstruction; the front TM and front Pf coseismic terms were already too large, indicating a source-amplitude normalization error.

The implemented correction is Eq. (1)-based, not a trace fit:

- Eq. (1) defines the acoustic pressure as `p_hat=A(omega)/R exp(-ikR) D(theta)`.
- The recorded pressure waveform used as source input is on axis at `(r,z)=(0,0)`, so `R=|z_s|` and `D(0)=lim J1(x)/x=0.5`.
- Therefore the pressure record is converted to the source amplitude with `A(omega)=P_record(omega) |z_s|/0.5`.

After this correction, the predicted/model peak-to-peak ratios relative to the page-digitized published predicted panels are approximately:

- z1: 0.84
- z2: 0.93
- z3: 0.91
- z4: 0.91
- z5: 0.96
- z6: 1.30

The z6 maximum absolute amplitude is now about `0.315 mV`, compared with about `1.05 mV` before the correction. The z6 timing and polarity checks remain consistent with the paper.

## Formula Fidelity Audit

Verdict: `PASS_WITH_RESERVATIONS`.

The implementation does not show evidence of per-trace hard fitting. The fluid-side response follows JAP Eq. (5), and the porous-side response follows JAP Eq. (7)-(8) as a deterministic forward model. The published Fig. 2 digitization is written only as diagnostic data and is not fed back into the model.

Evidence:

- Fluid positions 1-5 are synthesized by `schakel2011_jap_fig2_reproduction.py::_fluid_response`, using the real-angle Sommerfeld integral plus the gamma-path evanescent integral, `J0`, `J1`, source-depth phase, pressure-normalized `R^E`, and `exp(i k_z^E z_r)`.
- Position 6 is synthesized by `_porous_response`, including the front-interface TM wave, transmitted Pf coseismic field, and the three Eq. (7) interface-response terms.
- Table II parameters are encoded in `JAPParameters`; receiver geometry is encoded in `JAP_POSITIONS`.
- Table I pressure normalization is applied through `_pressure_normalized_re` and `_front_terms`.
- Fluid-incident conversion coefficients come from `schakel2010_strict_sensitivity.se_coefficients`, the Schakel and Smeulders (2010) Appendix-B boundary matrix implementation.

Remaining reservations:

- Pf-incident coefficients are reconstructed by reusing the Schakel 2010 boundary rows with a hand-built unit Pf incident column. This is transparent and theory-based, but it is not a verbatim printed formula from the JAP paper.
- `R_b` and `R_f` are treated as the same local coefficient because the front and back media are the same fluid/porous pair; the reversed interface orientation is not separately derived.
- The sign bridge from Schakel 2010 field-potential notation to JAP scalar electric potential is plausible but not explicitly printed in JAP Table I.
- The pressure source waveform is a visual/synthetic approximation inherited from the existing Schakel 2011 Geophysics reproduction, not the unavailable raw JAP pressure record. It is now converted to `A(omega)` using the Eq. (1) source geometry.
- The paper states a 144-896 kHz band-pass but not the filter shape; the code uses a 40 kHz raised-cosine taper for numerical smoothness.
- The gamma integral is truncated at `gamma=8`; the path form matches the paper, but this output does not yet include a tail-convergence table.

## Figure 2 Result Audit

Verdict: `PASS_WITH_RESERVATIONS`.

The generated results correspond to the main Fig. 2 kinematic and polarity features, but they should not be described as an exact absolute-amplitude match to the printed/digitized traces.

Evidence:

- First pulse timing is `0.09956-0.09969 ms` for positions 1-6, matching the paper statement that the pulse appears near `0.10 ms`.
- Fluid positions 1-5 show increasing amplitude toward the interface: `0.1156, 0.1592, 0.2265, 0.3462, 0.5554 mV`.
- Position 6 has reversed first-pulse polarity relative to positions 1-5.
- Position 6 second-pulse delay is `5.633 us`, close to the paper value of `5.7 us`.
- The observed-scale plot uses the paper's published global scale factor, `model / 2.5`.

Remaining reservations:

- Absolute amplitudes are now close for positions 1-5 and still somewhat high for position 6. The z6 predicted peak-to-peak amplitude is about `1.30x` the page-digitized published predicted trace.
- Position 6 still contains more late oscillatory energy than the page-digitized trace. This is most likely tied to the reconstructed Pf-incident/back-interface coefficients and the unavailable raw source waveform, not to a per-trace fit.
- The raster digitization is useful for timing/polarity/shape comparison only; weak printed traces and baseline pixels limit pointwise residual accuracy.

Recommended interpretation:

Describe the result as a transparent, theory-based qualitative/kinematic reproduction of JAP Fig. 2: timing, spatial decay trend, polarity reversal, and the 5.7 us position-6 delay are reproduced. Do not claim calibrated absolute-amplitude agreement with the printed measured traces.
