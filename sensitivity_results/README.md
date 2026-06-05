# Parameter Sensitivity Outputs

This directory separates two tasks.

- `schakel_reproduction/`: Schakel-style parameter sweeps using Table I reference settings. Figures use the current boundary-value solver's `R_E` and `T_TM` potentials plus squared-magnitude energy proxies, not the full orthodox/interference vertical-flux decomposition.
- `research_data/`: Sensitivity analysis for `global_evolution.xlsx`. One-at-a-time curves isolate model-input groups along the dissolution path. The waveform contribution bar uses signed changes in `log10(Amax)` between the first valid and last valid poroelastic snapshots; the residual is the nonlinear interaction not assigned to one parameter.
