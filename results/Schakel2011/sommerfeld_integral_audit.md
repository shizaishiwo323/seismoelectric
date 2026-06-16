# Sommerfeld Integral Audit

Target equation: Schakel et al. (2011, Geophysics) Eq. (5), fluid-side reflected electric potential.

Code location: `schakel2011_geophysics_reproduction.py::_sommerfeld_response`.

## Real-angle branch

Paper term:

`- i A(omega) / a * integral_0^(pi/2) J0(k r_r sin(theta)) J1(k a sin(theta)) exp(i k z0 cos(theta)) exp(i k_z^E(theta) z_r) R^E(theta) dtheta`

Implementation:

- `h = linspace(0, pi/2, n_theta)`
- `sin_h = sin(theta)`
- `kx = k * sin_h`
- `J0(k * radial_m * sin_h)`
- `J1(k * a * sin_h)`
- `exp(1j * k * source_z_m * cos(theta))`
- `exp(1j * kz_e * receiver_z_m)`
- `pressure_normalized_reflection_coefficient(frequency_hz, theta_deg, salinity_key)`
- adaptive mode: real and imaginary parts are evaluated with adaptive quadrature on `0..pi/2`
- fixed diagnostic mode: `-(1j * A / a) * first_integral` with trapezoidal samples

## Evanescent branch

Paper substitution:

`theta = pi/2 + i ln(sqrt(gamma^2 + 1) + gamma)`

Paper term:

`+ A(omega) / a * integral_0^infty J0(k r_r sqrt(gamma^2+1)) J1(k a sqrt(gamma^2+1))/sqrt(gamma^2+1) exp(k z0 gamma) exp(i k_z^E(gamma) z_r) R^E(gamma) dgamma`

Implementation:

- `gamma` path is evaluated on `0..8`
- `root = sqrt(gamma**2 + 1)`
- `kx = k * root`
- `J0(k * radial_m * root)`
- `J1(k * a * root) / root`
- `exp(k * source_z_m * gamma)`
- `exp(1j * kz_e * receiver_z_m)`
- `pressure_normalized_reflection_coefficient(..., kx_override=kx)`
- adaptive mode: real and imaginary parts are evaluated with adaptive quadrature on `0..8`
- fixed diagnostic mode: `(A / a) * second_integral` with trapezoidal samples

The infinite upper limit is truncated at `gamma=8`; for the paper source depth `z0=-0.15 m`, the damping factor `exp(k z0 gamma)` makes the omitted tail small over the 144-896 kHz band. This follows the same Sommerfeld path as the paper, while using SciPy adaptive quadrature rather than MATLAB's recursive adaptive Simpson implementation.

## Wavenumber Branch

The fluid EM vertical wavenumber is computed as:

`k_z^E = sqrt(omega^2 * (mu epsilon0 epsilon_f - i mu sigma_fl / omega) - kx^2)`

The helper `strict.complex_sqrt_branch` selects the branch with non-positive imaginary part, matching the paper statement `Im[k_z^E] < 0` under the `exp(i omega t)` convention.

## Conversion Coefficient

`R_E` is obtained from the Schakel and Smeulders (2010) Appendix B boundary-value solver. The Geophysics/JAP pressure-normalized coefficient is:

`R^E_pressure = R_E_2010 / (rho_f * omega^2)`

This is implemented in `pressure_normalized_reflection_coefficient`.

## Known Limitation

The integral structure matches Eq. (5), but `A(omega)` is a causal Ricker source because the original measured pressure record is not available in the project and the reproduction request specified a Ricker source. No electric-potential amplitude fitting is applied.
