# R1 action-native gradient / dispersion audit v1

**Status:** `PASS_NUMERICAL_PRINCIPAL_GRADIENT_DISPERSION_AUDIT`

This audit uses the existing action-native ADM + Schutz--Sorkin quadratic system. It performs no observational fit, no parameter retuning, and no Gate-7 operation.

## Owner-reproduction gate

- Worst n=2 action-block difference: `0.000000e+00`
- Worst n=2 reduced-block difference: `0.000000e+00`

## Kinetic gate

- Minimum reduced kinetic eigenvalue over the scan: `1.940993223685e-07`
- All scanned kinetic matrices positive: `True`

## Full fluid-coupled principal dispersion

- High-n sequence: `[80, 160, 320, 640, 1280]`
- Maximum extrapolated Re(s)/sqrt(k^2) intercept: `0.00038453528688022127`
- Minimum extrapolated non-radiation propagating c^2: `0.9103733185808105`
- Maximum radiation-branch |c^2-1/3|: `0.0003349432769985672`
- All high-k growth tails nonincreasing: `True`

A pressureless-dust growing mode is not classified as a gradient instability unless its growth scales as sqrt(k^2) at high k.

## Independent pure gravity/scalar owner diagnostic

- min G_S,2: `1.471470144828e-03`
- min F_S,2 (gravity/scalar only): `5.037267318952e-01`
- min F_S,2/G_S,2: `2.590795631561e+00`

The F_S/G_S result is reported as an independent owner diagnostic; the manuscript-level stability classification comes from the full matter+radiation dispersion pencil.

## Claim boundary

- This is a linear/quadratic stability audit on frozen R1.
- It is not nonlinear fluid closure.
- The n>=3 continuation is used to expose the principal spatial symbol; it is not a new high-n observational fit.
- A numerical PASS supports the paper-level statement that no ghost or high-k gradient instability was detected over the audited redshift/harmonic domain. It is not a microscopic BHSM theorem.
