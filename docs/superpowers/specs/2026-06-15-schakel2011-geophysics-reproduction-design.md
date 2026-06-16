# Schakel 2011 Geophysics Reproduction Design

## Goal

Reproduce the Schakel et al. (2011, Geophysics) fluid/porous-medium interface forward-model figures 2, 4, 5, 6, 7, 9, and 11 under `results/Schakel2011`.

## Scope

The reproduction implements the paper's fluid-side reflected electric-potential model:

- Source pressure: Geophysics Eq. (1), using a page-digitized approximation of Figure 4 because the raw pressure record is not provided.
- Piston directivity: Geophysics Eq. (2).
- Sommerfeld integral: Geophysics Eq. (3)-(5).
- Pressure-normalized seismoelectric conversion coefficient `R_E`, computed from Schakel and Smeulders (2010) Appendix B boundary-value system using Pride dynamic coefficients.
- Experimental parameters from Geophysics Table 1 and Table 2.
- Source location `(r0, z0) = (0, -15 cm)`.
- Band-pass filter 144-896 kHz.
- Figure 6 target: modeled reflected electric potential and displayed frequency spectrum for `1e-3 M NaCl`; `A_C` is not applied.
- Target NaCl case for figures 7, 9, and 11: `1e-3 M`, with `A_C = 0.19`.

The implementation is intentionally separate from the Liu 2018 finite-offset workflow. It may reuse the already validated Schakel 2010 coefficient solver, but not the Liu spectral waveform model.

## Digitized Experimental Data

The paper does not provide machine-readable data. The reproduction stores digitized or visually approximated data in CSV files:

- `digitized_source_fig4.csv`: source pressure pulse and spectrum approximated from Figure 4.
- `digitized_source_directivity_fig5.csv`: asterisk points approximated from Figure 5.
- `digitized_paper_fig6_waveform.csv` and `digitized_paper_fig6_spectrum.csv`: published Figure 6 model curve digitizations used only to check display reproduction.
- `digitized_peak_amplitudes_fig11.csv`: symbol points approximated from Figure 11.

These data are used for overlay and error summaries only. Model curves are generated independently from the equations and are not fitted point-by-point to those digitized data. The only paper-prescribed amplitude scaling applied to electric-potential model outputs is `A_C`.

The waveform panels do not synthesize measured traces. Where the paper's full measured traces are not digitized, the reproduction shows model traces only and documents that limitation.

## Output Contract

The pipeline writes:

- Reproduction figures: `fig2_reproduction.png`, `fig4_reproduction.png`, `fig5_reproduction.png`, `fig6_reproduction.png`, `fig7_reproduction.png`, `fig9_reproduction.png`, `fig11_reproduction.png`.
- Model and digitized tables for each figure.
- `parameters_used.csv` with provenance and values from Geophysics Table 1 and Table 2.
- `formula_audit.md` explaining the formula mapping and limitations.

## Tests

Tests verify:

- Table 1 and Table 2 constants match the paper.
- The source directivity obeys the finite piston limit at zero angle.
- The source pulse spectrum peaks near 0.5 MHz and is band-limited by 144-896 kHz.
- The Schakel 2011 config produces finite pressure-normalized `R_E` values through the Schakel 2010 boundary solver.
- Figure 6 writes unscaled model values, paper-style display columns, and digitized published reference curves.
- The model's peak amplitudes decay with increasing distance from the interface along the negative z-axis.
- The pipeline writes the expected figures and CSV files.
