# Constructed local encounter and manuscript corrections, 20 October 2026

This continues draft PR #23 from reviewed commit
`97bf5571335cf60e15e2e486e1197ca226198ed2` in the existing worktree and branch.
That commit is a provenance anchor. No reset, broad cosmology rerun, source
import, CI repair, muon preparation or backup was performed. The original
SIS and smooth-core numerical inputs/arrays are unchanged.

## Results evaluated

The corrected endpoint derivation distinguishes fixed-affine and
fixed-radial frequency variations, includes the ray-coordinate shift in
the actual radial event shift, and uses total frequency variation at a
moving vertex. The numerical event-root implementation was preserved.

The scalar-map SIS/smooth-core response factorizes analytically as
`B = N_min R`. Application to the saved arrays gives relative whitened
discrepancies `2.452944603865752e-12` and `6.377724303772925e-12`.
This explains the covariance-independent nuisance null in these families.
It is an observable degeneracy between different source/galaxy preparations.
The fixed-physical-source SIS delay and its qualified H0 cancellation remain.
The old smooth-core central image has magnification ratio `2.79565e-5`;
equal image precision remains an idealized control.

The new local reference is the physical Path-B S6 meridian
`eta_q=(cos q, sin q e1/sqrt(2))`, with flat color connection and regular,
source-free target boundary data. Its radial equation is

```
partial_r[r² exp(nu) sqrt(f_r) w(kappa_eta+X³) qprime] = 0,
X = f_r qprime².
```

Regularity forces the conserved flux to zero; positive `w,kappa_eta` then
forces `qprime=0`. Actual sections at `q=0,0.2,0.6` are evaluated. Target,
unit, color Gauss/Euler, momentum and stress rows vanish exactly. This is
a continuous conditional state family, rather than normalization of an
arbitrary texture. Nonconstant phase fluctuations still have positive
Dirichlet stiffness and a positive original Lorentzian kinetic matrix.

Ordinary finite matter supplies the nonzero gravitational field:
`M=1e11 Msun`, areal radius `R=3 kpc`, density
`15M/(8piR³)(1-r²/R²)` inside and zero outside. Pressure and lapse are
solved from Einstein/TOV, with `p(R)=0` and Schwarzschild lapse at R.
The resulting barotropic `p(rho)` is positive and has maximum
`c_s²/c²=1.9939296285e-6`. Compactness is `1.5951386139e-6`.
Induced metric/extrinsic curvature and opposite-side target flux match
without an Israel surface layer. Pressure/lapse/lapse-derivative surface
residuals are zero at recorded precision. Independent TOV and lapse
derivative residuals are `4.27e-13` and `2.38e-12`; Hamiltonian/radial
Einstein residuals scaled by compactness are below `4.3e-15`.
The angular residual scaled by compactness squared is `2.23e-9`.

An exact physical log-lapse/spatial-conformal metric Hessian is evaluated
on the vacuum isotropic annulus outside this same body, with homogeneous
Dirichlet variations and the same normalized radial profile in both metric
columns and the phase column. In static energy convention, defining
`K=c^4/(8piG)` and `k_eta=w*kappa_eta`,

```
(a,b,c)/(K R) = (-1.27e-21, 19.739240289, -19.739240289),
u=v=0,
c_phi/(k_eta R) = Schur/(k_eta R) = 9.86962014448,
det(A_g)/(K R)^2 = -389.637607185,
condition(A_g)=2.618033989.
```

The exact simplified block is `2K I_E [[0,1],[1,-1]]`; the original
unsimplified second variations independently reproduce it. The profile
norm is `1+2.2e-16`; the original vacuum metric Euler residual is
`3.27e-21` in R units. The projected metric inverse is justified only
on this two-mode Dirichlet domain; the static scalar energy signature is
not used as a propagating ghost test. The normalized target-phase source
operator is applied with homogeneous Dirichlet data and zero derived
phase forcing/flux. Positivity proves its unique zero solution.

For the same finite matter, source and propagation maps, the weak-field
line-of-sight integral gives projected mass fraction
`F2(y)=1-(1-y²)^(5/2)` for `y<1` and 1 outside. Its potential,
deflection, Hessian and fixed-impact log-mass/log-radius derivatives are
computed. This encounter replaces the SIS in its own control; it is not
added to the SIS or to the smooth propagation potential.

| Application | Evaluated value |
|---|---|
| Image x positions, arcsec | 0.83547603394, -0.08216193978, -0.65749552508 |
| Relative delays from positive image, days | 21.7160765661, 14.2714523869 |
| d delays / d ln M, days | 31.5546082822, 7.10166443362 |
| d delays / d ln R, days | -28.4907205308, 0 |
| d images/dq and d delays/dq | exactly zero |
| Original image residual, radians | 8.47e-22 |
| Mass/radius finite-change whitened errors | 1.08e-10 / 1.03e-9 |
| New finite-profile nuisance factorization error | 1.3152e-16 |
| Six cosmic state columns, fixed/free nuisance rank | 3 / 0 |

The constant-state BHSM target sector gives the conventional-equivalent
encounter for the same matter. This is a derived zero additional response
with nonzero gravity. The mass/radius columns measure the shared body;
they are not an extra BHSM effect. The finite-profile nuisance application
was recomputed rather than borrowed from the SIS theorem. No independent
kinematic likelihood or degeneracy-breaking information was invented.

## Exact scope and remaining dependency

The conditional effective gravity coefficient is
`C_partial=c^4/(16piG)`, with no additional local cap traction or vacuum
energy. Target fields add no eta-dependent sigma/cap source, but the
independent local sigma/cap equations are outside this truncation.
Changing the local metric could affect those sectors. Full stratified
admissibility and matching to the rolling closed-cosmology background are
not established. Pure background-curvature corrections scale as
`(H R/c)^2`; this does not bound a separate rolling-scalar/cap mixed response.
Optical values use leading weak-field/thin-plane accuracy, with relative
omitted scales compactness and R/D. Exact TOV does not make this an
all-orders optical calculation.

The recovered target action defines the canonical target-phase operator,
source and tangent field normalization. It does **not** identify the
separately exported topographic scalar psi with that phase. Its defining
consumed relation remains a common-action map from exported forcing to
physical phase/metric, including the gravitational source in
`L_psi psi=J_grav+E_Sigma[Q_Sigma]+xi` if that branch is prescribed.
Zero target-phase forcing does not set this separate J_grav to zero.
The solved normalized target-phase problem and the nonzero conventional
gravitational encounter remain valid within the declared sector without
claiming that missing constitutive map. A nonzero seam source would
require specified charge/traction data and a new coupled solution.
No moving-seam continuation is assumed.

The light remains the common probe. No statistical state kernel is
inferred from a retarded array or native heat functional. Exclusive
chronological two-point, inclusive contour four-point and complete matched
occupation calculations remain distinct downstream tasks. No Wick
replacement, calibration change, accounting duplication or change to the
inherited `2e-10` physical guard is made; this conditional local calculation
is not physical admission of an unfinished preparation.

## Reproduce and review

From the worktree root:

```powershell
python code/strong_lensing_degeneracy.py
python scripts/run_strong_lensing_local_family.py
python -m pytest tests/test_strong_lensing_comparison.py tests/test_strong_lensing_transport.py tests/test_strong_lensing_envelope.py tests/test_bhsm_optical_transfer.py tests/test_bhsm_optical_constitutive.py tests/test_bhsm_optical_hessian_projection.py tests/test_bhsm_optical_current_gate.py tests/test_bhsm_optical_low_rank_covariance.py tests/test_strong_lensing_degeneracy.py tests/test_strong_lensing_intrinsic_encounter.py tests/test_strong_lensing_local_family.py tests/test_strong_lensing_local_hessian.py -q
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error -file-line-error -outdir=build/strong-lensing strong_lensing_main.tex
```

Seventy-five focused tests pass. Independent scientific review records
24 tests plus source/action, signs, constraints, domain, source ownership,
matching, once-only accounting and nuisance audits in
`artifacts/strong_lensing_local_family_20261020/review.json`.
Original BHSM source commit `524ed90689bd5923c249bba2e699abf627e703cd`
and eight exact source hashes are pinned in the source receipt. The
retained baseline files are hashed and checked against the reviewed
checkpoint. Full transport is reused, not rerun as a dependency.

The standalone manuscript is 17 pages. Its final LaTeX log has no
overfull/underfull boxes or undefined references. All pages were rendered
and visually inspected as a contact sheet, with full-page inspection of
the endpoint/proposition/action/Hessian/source/result pages. Publication
verification hashes the PDF and sources and retains the contact sheet.
The broader Universe manuscript and other worktrees remain unchanged.
