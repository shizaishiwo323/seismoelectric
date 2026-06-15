# Recommended Figure Sequence Analysis Manifest

This directory contains the current paper-facing analysis sequence derived from `doc/文献调研_反应输运耦合震电.md`.

## Inputs

- `Pe=0.1`: `C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI\dissolution_results-Da_0.0369_Pe_0.1000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`
- `Pe=1`: `C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI\dissolution_results-Da_0.0369_Pe_1.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`
- `Pe=10`: `C:\Users\imgw\Documents\MATLAB\RTSPHEM-main\T2single-RI\dissolution_results-Da_0.0369_Pe_10.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`

## Core Scripts

- Figure 1: `analysis_scripts/figure1_conceptual_workflow.py`
- Figure 2: `analysis_scripts/figure2_rt_evolution.py`
- Figure 3: `analysis_scripts/figure3_dynamic_bridge.py`
- Figure 4: `analysis_scripts/figure4_interface_coefficients.py`
- Figure 5: `analysis_scripts/figure5_waveform_panels.py`
- Figure 6: `analysis_scripts/figure6_spatial_peak_dipole.py`
- Figure 7: `analysis_scripts/figure7_mechanism_decomposition.py`
- Figure 8: `analysis_scripts/figure8_monitoring_metrics.py`
- Full run: `analysis_scripts/run_recommended_figure_sequence.py`
- T0 leakage diagnostic: `analysis_scripts/diagnose_waveform_t0_leakage.py`
- Frequency convergence diagnostic: `analysis_scripts/check_waveform_frequency_convergence.py`

## Figure Outputs

- `figures/figure1_imagegen_concept.png`: imagegen-generated conceptual asset.
- `figures/figure1_conceptual_workflow_timescale_separation.png`: deterministic manuscript workflow diagram.
- `figures/figure2_reactive_transport_evolution.png`
- `figures/figure3_dynamic_electrokinetic_bridge.png`
- `figures/figure4a_conversion_coefficients_time.png`
- `figures/figure4b_conversion_frequency_angle_heatmaps.png`
- `figures/figure5_finite_offset_waveform_panels.png`
- `figures/figure6_spatial_peak_dipole_interpretation.png`
- `figures/figure7_mechanism_decomposition.png`
- `figures/figure8_normalized_monitoring_metrics.png`
- `figures/waveform_t0_leakage_diagnostics.png`
- `figures/waveform_frequency_convergence.png`

## Key Numerical Finding

The cross-case summary is saved in `tables/cross_case_key_findings.csv`. Within the valid poroelastic window, permeability increases by about 20 times for Pe = 0.1 and by more than 200 times for Pe = 1 and Pe = 10. In contrast, the higher Pe cases show strong conductivity growth: $|\sigma|$ increases by about 114-158 times, while $|L|$ decreases to about 1% of its initial value and $R_E$ decreases to about $10^{-4}$ of its initial value.

## Audit Response

An independent subagent audit found that the first analysis pass was a real mechanism-chain analysis, not a superficial figure list, but flagged a serious Figure 5 issue: coarse frequency sampling produced strong pre-$T_0$ aliases. The waveform scripts were revised to use $n_\omega=128$ for Figure 5. The new `waveform_t0_leakage_diagnostics.csv` shows `pre_to_post_T0_ratio` of about 0.011-0.012 across all nine waveform panels.

## Remaining Scientific Caution

The current Figure 7 chemistry pathway groups `OutletHConc`, electrolyte concentration, and fluid conductivity together. This is appropriate for the current RT table because independent measured $\sigma_f$ is not available, but a future mechanism-decomposition version should split $\mathrm{H}^+/\zeta$ and $\sigma_f$ if independent conductivity data or a better electrochemical model is added.
