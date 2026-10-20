# Finite conventional refit of the prepared R1 encounter

This is an inverse-model application to the **existing** R1 three-image
injection. Its purpose is to determine what its light observables identify
after conventional lens freedoms are admitted. The frozen state, physical
incident marker, source, redshifts, clock and distance maps that generated
the injection are unchanged. Fitted source coordinates are inverse-model
nuisances; they do not replace the source used in the physical comparison.

The original leading offsets, 25.959258149566722 s and
8.649019375383205 s, compare the R1 response to its prescribed fixed-body
reference. They are not residuals after a conventional model has been fit.

## The original applications fitted

Use the same finite Tolman VII conventional profile, now with mass M,
radius R, angular center theta_c, source angle beta and one physical sheet
as inverse-model coordinates. For each observed image theta_i, evaluate

\[
 r_i=\theta_i-\beta-q\hat\alpha_{M,R}
       [D_{ol}(\theta_i-\theta_c)]
       -\kappa(\theta_i-\theta_c),\qquad q=D_{ls}/D_{os}.
\]

The physical encounter potential is

\[
 \hat\psi(x)=\hat\psi_{M,R}(x-D_{ol}\theta_c)
  +\frac{\kappa}{2D_{ol}q}|x-D_{ol}\theta_c|^2.
\]

The corresponding original Fermat row is

\[
 ct_i=\frac{1+z_l}{2}\frac{D_{ol}}q
                (\theta_i-\beta)^2
       -(1+z_l)\hat\psi(D_{ol}\theta_i).
\]

Fit the three r_i=0 rows and the two relative arrival-time rows. The five
coordinates are (delta ln M, delta ln R, delta beta in arcsec, theta_c in
arcsec, kappa). The sheet is in the encounter once and is not also placed
in the propagation maps. No free H0 or arbitrary independent delay scale
is introduced. Fixing the center leaves four coordinates for five rows.

The numerical scales 0.001 arcsec and 0.1 day condition the least-squares
problem. They do not define an observing likelihood, uncertainty forecast,
state population or quantum contraction. All other axial coordinates
remain zero. The fits are local, bounded and noncritical, and roots are
independently re-solved after fitting; observed image positions are not
simply copied into the final prediction.

## Actual result

| Coordinate | Fixed center | Free center |
|---|---:|---:|
| delta ln M | 3.9959695571e-6 | 4.0205625603e-6 |
| delta ln R | 3.9827138748e-7 | 4.1234000691e-7 |
| delta beta (arcsec) | 4.3495339029e-9 | 6.8147778325e-9 |
| theta_c (arcsec) | 0 | 1.7743086205e-9 |
| kappa | 9.9468749779e-6 | 9.9228734337e-6 |
| Maximum re-solved image mismatch (microarcsec) | 1.0143e-4 | 8.74e-11 |
| First delay mismatch (s) | +0.00330640 | +9.21e-10 |
| Second delay mismatch (s) | -0.00503119 | -1.53e-10 |

The free-center final lens-equation residual is 4.24e-22 radians and the
stationarity residual is 1.30e-18 Mpc. Its scaled five-by-five Jacobian
has singular values approximately (3641.6,1875.9,1564.7,415.4,19.37).
Those numbers concern the declared numerical coordinate scales.

This finite nonlinear application confirms the earlier first-order
absorption result for this particular three-image preparation. Five
adjustable coordinates fitting five scalar observables is intentionally
a saturated local example, not a universal theorem about arbitrary lenses
or a successful observational inference. The constructive positive-halo
result is the stronger retained-field equivalence for all incident rays.
Additional source systems constrain a shared restricted profile, while
a sufficiently flexible positive effective mass still reproduces the
retained R1 metric. No distinctive GR departure is required or claimed.

## Reproduction and checks

```text
python scripts/run_strong_lensing_refit.py
python -m pytest -q tests/test_strong_lensing_refit.py
```

Four focused tests pass: the exact finite mass-sheet image/time symmetry;
recovery of an independently generated finite conventional injection;
the original-row/independent-root R1 fit; and the remaining mismatch when
the center is fixed. This does not rerun the original nine-command
calculation or alter its artifacts. `nonlinear_refit.json` records the
source hashes and original residual applications.
