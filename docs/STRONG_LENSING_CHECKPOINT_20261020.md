# First matched strong-lensing checkpoint - 20 October 2026

> The [acceptance checkpoint](STRONG_LENSING_ACCEPTANCE_20261020.md) is the
> current reader entry point. It follows the retained source-to-light transfer
> into finite inverse modeling and controlled profile/accuracy applications.

> Historical checkpoint. Its calculations and validation belong to the first
> comparison recorded below. The [current scope](STRONG_LENSING_SCOPE_20261020.md)
> and [R1 research review](strong_lensing_r1_screened_20261020/GALAXIES_RESEARCHER_REVIEW.md)
> supersede its project framing and earlier submission schedule. The active
> BHSM paper is for Galaxies; it treats leading conventional effective-mass
> equivalence as a constructive representation of the retained envelope response.
> The old local-source inventory is followed by the evaluated constant-target
> family and the separate R1 sourced branch. No historical absence statement
> overrides those later calculations.

> The [current integration checkpoint](STRONG_LENSING_BHSM_INTEGRATION_20261020.md)
> additionally records the evaluated coupled R1 source/flux application,
> positive-lift nonuniqueness proof and exact field-continuation derivation.
> Those later equations and source residuals supersede this receipt's old
> unworked-source inventory. They do not claim an evolved higher-order body or
> global matching solution; their consumed material/environment/outer-state
> prescriptions are explicit in the continuation record. Author declarations
> and contact have since been confirmed in the active lensing materials.

Historical baseline receipt. The subsequent constructed local encounter,
corrected endpoint derivation and analytic nuisance result are recorded in
[the local checkpoint](STRONG_LENSING_LOCAL_CHECKPOINT_20261020.md).
The original numerical baseline below remains unchanged.

The existing checkout was inspected at revision
`3427c5297fac71faa501685bd38e5d420619c761`. Work is isolated on
`codex/strong-lensing-first-comparison-20261020`; the original checkout's
modified manuscript PDF/scalar bridge and unrelated local research are
preserved. The recovery revision `5fd65ea18a6055dbc6028839c2918415e3c00f71`
is provenance, not a reset target. The action-native files used here agree
with that later recovery's relevant source versions; the actual base above
is recorded in the input manifest.

The [scope record](STRONG_LENSING_SCOPE_20261020.md) preserves the accepted
invitation and 25 February 2027 submission target. No journal submission,
correspondence or fee commitment is part of this checkpoint. The local
handoff text was recovered; the named ZIP was not available. Its essential
current-owner definitions were recovered from the checkout and verified
directly rather than treated as missing theoretical equations.

## Executed experiment

Two realizations use the same transient narrowband light preparation,
observer proper clock, measured redshifts `z_l=0.3,z_s=1.0`, physical galaxy
and transverse physical source displacement:

1. The conventional SIS galaxy on the frozen R1 background, with velocity
   dispersion220km/s. The source is away from caustics at reference
   `beta=0.2 theta_E` and displacement1.457663367kpc.
2. The same encounter with the smooth optical response of a selected
   normalized cosmological n=2 profile. Six physical anchor amplitudes are
   propagated from the existing coupled action. The additional established
   local static spatial phase channel has zero mixed response under its
   stated current gate.

The second realization is a conditional metric correction experiment.
A microscopic envelope replacing the conventional galaxy is **not**
established. No conventional potential, smooth metric response or residual
line-of-sight uncertainty is counted twice. The incoming light, radiation
fluid perturbations and muon/reference covariance are distinct; no muon
preparation or inclusive four-point calculation was restarted.

The selected spatial profile obeys
`Q=sqrt(3)/(pi R^(3/2)) sin(2 chi/R) mu`,
`integral Q^2 dV_gamma=1`, with conditional axial ray`mu=1`.
Its normalization is0.02709350135 in H0=1 units. The nine-component n=2
spatial eigenspace is not claimed to have only six coefficients: six remain
after selecting this one profile. The reference-normalized Pi2 is mapped to
the full canonical state; it is not silently identified with p_zeta.

## Numerical results

| Quantity | Executed result |
| --- | ---: |
| D_ol [Mpc] | 945.16283937 |
| D_os [Mpc] | 1689.31223969 |
| D_ls [Mpc] | 1076.99434806 |
| SIS image angles [arcsec] | +1.0678830532, -0.7119220355 |
| Relative arrival delay [days] | 17.0945320516 |
| State rank with nuisance parameters fixed | 3 |
| State rank with declared nuisances unrestricted | 0 |
| State rank with synthetic independent nuisance priors | 3 |

The six fractional delay derivatives, in anchor order
`(zeta,Pi2,delta_m_ADM,delta_r_ADM,v_m_ADM,v_r_ADM)`, are approximately
`(-1.16307e-5,2.83592e-3,2.79299e-4,7.83676e-8,-2.69536e-3,5.31244e-7)`.
They are unit sensitivities, not derived physical amplitudes.

A specified `zeta=1e-3` numerical probe changes the delay by
`-1.98831e-7 days = -0.0171790 seconds`. Image locations are independently
re-solved. The primary matched experiment fixes transverse physical source
position; inferred angular source position changes. Fixed observed redshift
sets the radial endpoint. A fully fixed three-dimensional source worldline
is not additionally imposed.

The0.05kpc smooth-core control has three noncaustic images and two relative
delays30.259543 and17.093193days. Its delay ratio0.56488603 removes common
mass-sheet scaling, but source/galaxy freedoms still leave state rank0.
More independent directions alone cannot remove one free mass-sheet scale
per lens. Synthetic prior-assisted information is identified as extra
assumed information, not inferred from these illustrative observables.

## Equations and validation

`code/strong_lensing_transport.py` integrates the original metric null and
angular geodesic equations, matter tetrads, redshift roots and lens-vertex
Jacobi reset. It excludes the thin galaxy. All three maps share their
screen and affine/frequency conventions; D_ls is not a difference of
observer distances. Reverse propagation checks reciprocity without
symmetrizing maps. The inherited seven-state optical equation supplies
an independent conformal Jacobi control, with measured-redshift shifts
restored explicitly. Its unused frequency slot is not a claim of zero
physical frequency response.

`code/strong_lensing_comparison.py` then evaluates

```text
beta = theta - Dos^-1 Dls alpha_hat(Dol theta),
T = (1+zl) Dol^T Dls^-1 Dos,
ct = (theta-beta)^T T (theta-beta)/2 -(1+zl)psi_hat(Dol theta)+f(beta).
```

The derivative includes dT, physical-source remapping, physical impact and
any intrinsic encounter term. Image stationarity removes dtheta only from
the first-order arrival variation. Relative delays cancel f(beta); the
stored stationary arrival values are not absolute cosmological travel times.

The available native constrained ADM projection is recorded in
[the envelope derivation](STRONG_LENSING_ENVELOPE_DERIVATION_20261020.md),
with its original equations and units. At z0.3 its lapse/shift block has
`(a,b,c)=(-7.96725028,0.56954480,0.48267978)`, determinant-4.17001193,
original constraint residual5.27e-16 and Schur agreement4.15e-16.
These are not relabelled as local Newtonian metric/phase coefficients.
The current-null/constant-texture projection gives u=v=0 in its specified
static spatial phase channel. The full phase self-curvature is not
assigned zero from this fact.

Targeted verification covers curved-FLRW recovery, all-six-coordinate
constraint reconstruction, complete noncommuting matrix responses,
independently re-solved finite images, original source/impact terms,
MST/H0 convention controls, reciprocity, half-step convergence and reuse
of the existing optical modules. The new derivative check reaches
1.45e-7 maximum covariance-whitened relative difference; the independent
seven-state map control agrees to9.72e-11 in c/H0 per anchor unit.
These discretization checks are not physical uncertainty bars.

## Reproduce and inspect

From the worktree root:

```powershell
python scripts/run_strong_lensing_comparison.py
python -m pytest tests/test_strong_lensing_comparison.py tests/test_strong_lensing_transport.py tests/test_strong_lensing_envelope.py tests/test_bhsm_optical_transfer.py tests/test_bhsm_optical_constitutive.py tests/test_bhsm_optical_hessian_projection.py tests/test_bhsm_optical_current_gate.py tests/test_bhsm_optical_low_rank_covariance.py -q
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=build/strong-lensing strong_lensing_main.tex
```

Results and input hashes are in
`artifacts/strong_lensing_first_comparison_20261020/`. The machine-readable
results retain every map, endpoint shift, response, covariance, nuisance
column, singular value and state null direction. The posterior of a null
direction is not represented by a zero-uncertainty pseudoinverse estimate.
The focused paper is `manuscript/strong_lensing_main.tex`; it preserves the
existing Universe manuscript and bibliography. The numerical figure and
table are generated by the same executable as the results.

## Remaining dependency and next calculation

The local geometric-envelope amplitude needs a normalized physical
reference texture/profile and the common-domain constrained metric/phase
Hessian plus its source, normal and seam matching. The original
`L_psi psi = J_grav+E_Sigma[Q_Sigma]+xi` and
`h_topo=-A_g^-1 B_gpsi psi` identify where these inputs are consumed.
The static seam export remains distinct from a proposed moving-seam
kinematic continuation. Missing normalizations are not fitted to lens
data or inferred from unavailable arrays.

The result supports a separately specified envelope state and independent
lens/environment information, or a general ray calculation for the same
declared profile. It does not support adding a new angular channel to
manufacture a detection. A real observed-lens likelihood can follow once
its source/lens/environment contracts and independent degeneracy-breaking
information are declared.
