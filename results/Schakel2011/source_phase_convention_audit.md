# Source Phase Convention Audit

Schakel et al. (2011) Eq. (1) writes the acoustic source as:

`p_hat(omega, R, theta) = A(omega) / R * exp(-i k R) * D(theta)`.

The published Figure 4 pressure pulse was measured at `(r, z) = (0, 0)`, i.e. at the interface point on the source axis. At this point:

- `R = |z0| = 0.15 m`
- `theta = 0`
- `D(0) = lim J1(x)/x = 0.5`
- `T0 = |z0| / cP = 101.029 us`

For this reproduction, the raw Figure 4 pressure record and acquisition details are unavailable. Following the project instruction, `A(omega)` is therefore a causal Ricker source spectrum, not the digitized Figure 4 trace. The source-time transform is:

`S(omega) = integral s(tau) exp(-i omega tau) d tau`

with the wavelet peak at `2/f0` and duration `8/f0`, where `f0 = 500000 Hz`. The implemented source amplitude is a global source-amplitude assumption, not a published paper parameter:

`A(omega) = 0.03 * S(omega) / max(|S(omega)| over the modeled band)`

The Sommerfeld integral keeps the paper's `exp(i k z0 cos(theta))` propagation term. The time-domain synthesis uses absolute waveform time `t`, not `t - T0`, because subtracting `T0` would count the source-to-interface delay twice.

This convention was cross-checked against the project spectral synthesis code in `seismoelectric_offset_liu2018_spectral.py`, which uses the same `exp(i omega t)` time convention and keeps propagation phase inside the frequency-domain Green's function.

Remaining ambiguity: the JAP theory paper states that an experimentally recorded pressure waveform at `(r,z)=(0,0)` is used for `A(omega)`, but that waveform is not available here. The Ricker source is a transparent replacement and explains the remaining absolute-amplitude differences; no hidden pointwise amplitude fitting is applied.
