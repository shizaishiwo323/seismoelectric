# Strict Schakel & Smeulders (2010) Sensitivity Reproduction

This directory contains a strict reproduction of the two energy-flux
coefficients plotted in Schakel & Smeulders Fig. 4-7:

- `RE_EE`: reflected EM orthodox vertical flux coefficient.
- `TE_TM_TM`: transmitted TM orthodox vertical flux coefficient.
- `minus_TE_TM_TM`: the plotted positive value `-TE_TM_TM`.

The validation file compares the implementation with Schakel Table III at
`omega = 1e6 rad/s` and `theta = 30, 45 deg`.
