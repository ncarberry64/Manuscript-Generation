# BHSM strong-lensing integration checkpoint

> Historical integration checkpoint. Its science pin, 18-page PDF audit,
> numerical replay and remote outcomes belong to the completed integration
> snapshot below. The [acceptance checkpoint](STRONG_LENSING_ACCEPTANCE_20261020.md)
> supersedes its active manuscript status with the subsequent coherent
> transfer/inference revision, finite nonlinear re-fit and profile/accuracy
> controls. The older receipts are preserved and do not certify a later PDF.

This continues draft PR #23 on
`codex/strong-lensing-first-comparison-20261020`, from public provenance
anchor `49a1dabbe148f8ba9da4ac079a8d93fe97c5de76`. The existing worktree
was clean at that revision. The separate original checkout's cosmology
and other research edits were inspected and preserved. No reset was used.

The integrated calculations, inputs, outputs, source/continuation derivations
and checks are pinned at science commit
`0a329fa13f6c90bb2e3b5999b0155aa1e81fe8bd` (68 new files). The focused
paper's bibliography and data-availability statement cite that concrete
calculation revision. The manuscript/publication commit follows it on the
same draft-PR branch; the older `49a1dab` remains its provenance anchor.

The active paper is **Strong-Lensing Time Delays in the Berger–Hopf
Framework: An R1 Effective Calculation**, for **Galaxies**, special issue
**Strong Gravitational Lens Modeling**. No submission target date is set.
No PR merge, journal submission or cover-letter transmission was performed.

## Governing question and scientific order

**How do BHSM surface contributions determine gravitational lensing and
time delays for light traversing an envelope event?**

BHSM determines the question and physical interpretation. The central
object is the prepared incoming-to-outgoing light transfer through the
envelope. The manuscript follows envelope state and environmental
preparation; action, normalization, source and boundary conditions;
metric/optical response; image directions and observer arrival-time
differences; then conventional gravitational representation and checks.

Conventional matter and scalar descriptions are developed as
representations of the same surface contributions. That is the governing
BHSM premise. The actual calculation derives a conditional response from
the independently normalized R1 effective action, and proves that a
positive effective-mass profile reproduces the retained leading metric
and complete optical field. That constructive equivalence does not
complete the microscopic common-source reduction. A distinguishable
departure from conventional gravity is not the paper's success criterion.

The source-free constant-target result remains valid in its declared
boundary problem. The independent rolling scalar T, target phase and
metric conformal coordinate remain distinct. A scalar contribution and
its equivalent effective-mass contribution are alternative descriptions,
counted once. The local encounter and smooth propagation are also each
counted once.

## Supplied research package and integration

The supplied package was recovered at
`C:/Users/carbe/Downloads/Galaxies_Lensing_Research_Package.zip`.
Its SHA256 is
`5a66dae72f678baf7d5ef5a03919a8670e79edd91db03092b71ce80881458256`.
All 57 manifest entries and all 12 pinned dependencies were checked.
The package's 52-file overlay was reconciled with the actual worktree.
The four existing destinations were inspected rather than overwritten
as a group: README, local-family derivation, focused TeX main and PDF.
The active README and documents were merged; the source adopts the real
MDPI main and the PDF is rebuilt from the revised sources.

The new R1 screening, finite-body GR2, independent projection, coefficient,
composition and nuisance calculations, technical derivations, supplied
outputs, official MDPI dependencies and focused build workflow are now
integrated. Earlier package publication/validation claims remain
historical records under `artifacts/strong_lensing_bhsm_integration_20261020/`.
Active research review, scope, README, cover letter and handoff are
superseded consistently by the user's BHSM scope. Historical checkpoints
and review documents have explicit superseding notes; their numerical
evidence is retained. The separate Universe manuscript is unchanged.

## Preserved calculations and meaning

| Quantity | Supplied value retained in the paper |
|---|---:|
| Screening radius | 212.02922497712947 kpc |
| Leading fixed-body-reference delay offsets | 25.959258149566722 s; 8.649019375383205 s |
| Common finite-body GR2 reference-relative offsets | 25.959327139730703 s; 8.649009238058625 s |
| GR2 changes of these offsets | +68.99 microseconds; -10.14 microseconds |
| Effective asymptotic extra mass | 1.50788624e8 solar masses |

These seconds compare the sourced response with the specified fixed-body
reference. They do not represent an unavoidable disagreement between
fully matched scalar and effective-mass representations, which agree
where the leading equivalence is demonstrated. No parameter was retuned
and no new amplitude was introduced. The supplied scientific equations and
numerical results are preserved; explicit source applications and continuation
equations have now been added as described below. Output routing and UTF-8/LF serialization
were added to reproduction scripts so replays preserve historical evidence.
The local-family text correction limits zero export to the specified
material-surface seam contribution; it does not set every gravitational
source to zero.

The nuisance ranks retain their declared scope: single source with the
centre fixed or free, multiple sources with the finite radial profile,
and separate smooth-map columns. Weighting covariances are numerical
choices, not demonstrated observing precision. Mass-sheet, source,
centering, profile and line-of-sight freedoms remain explicit.

## Common-source application and the remaining physical map

The [common-source derivation](strong_lensing_r1_screened_20261020/BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
now writes the consumed source after eliminating the same metric constraint:

```text
2 Mpl² Delta Phi = rho + rho_Sigma - b_T Delta pi_T,
J_eff = J_T - b_T(rho + rho_Sigma)/(2 Mpl²),
F_R1 = C grad pi_T
       - (2/lambda)[(Delta pi_T) grad pi_T - Hess(pi_T) grad pi_T],
r² F_R1.n = C r² pi_T' - (4/lambda) r pi_T'^2,
[r² F_R1.n] = R²[J_T,Sigma - b_T sigma_Sigma/(2 Mpl²)].
```

The original coupled metric/scalar application is evaluated on the existing
finite-body section, with no new fitted scalar charge. At the material radius,
`pi_T' = -2.119825778845613e-4 Mpl H0` and
`integral_0^R r² J_eff dr = -1.216417408185084e-13` in the declared natural
units. Across the evaluated radii, the scalar-flux relative residual is at
most `4.15e-16`, the original metric Gauss residual is zero at recorded
precision, and the material flux-jump relative residual is `6.23e-16`.
The exact-GR reference material junction and the separately evaluated leading
R1 material junction each give `Delta K_uu = 0`; their extrinsic-curvature
values are separately labeled. Thus the specified material seam charge is
zero while the R1 volume matter-curvature source is nonzero. This is not an
exact modified-gravity TOV matching claim.
This ordinary material surface is distinct from a cosmological cap equator.
The exact applications, normalization and residuals are retained in
`artifacts/strong_lensing_bhsm_integration_20261020/source_application.json`.

The common-source work also proves a remaining nonuniqueness rather than
inferring a missing relation from an unavailable array. Two explicit positive
quadratic lifts on the same normalized Dirichlet domain have the same frozen
leading `C,D` but different normalized fourth-order coefficient
`B/ell0² = 0.2359125661` or `0.05897814152`, and expansion scales
`p ell0² = 1` or `4`. These are mathematical witnesses in the inherited Schur
family, not selected galaxy states or new physical parameters. The physical
embedding, complement spectrum and length normalization are not fixed by the
retained leading R1 coefficients. A direct constant rescaling of the inherited
nonzero `Theta(r-r_star)/r²` source is also excluded by its `1/r` flux tail,
whereas the isolated R1 canonical flux has a `1/r²` tail. This excludes that
direct identification; an action-owned finite-support export or nonlocal map
remains possible. Seven focused checks passed for this new source application
and its mathematical counterexample.

## Derived field continuation and its physical data

The [continuation derivation](strong_lensing_r1_screened_20261020/R1_REMAINING_FIELD_DERIVATION_20261020.md)
supplies the exact scalar current, cubic stress and scalar/Einstein component
applications from the same retained action. Its spherical unknowns are mass
length `m(t,r)`, log lapse `nu(t,r)`, rolling scalar `T(t,r)`, and material
density/velocity when the body is evolved. It includes the material continuity
and Euler rows, nonlinear environmental cross-Hessian, regular-branch overlap
coefficients, matching conditions and the subsequent null-ray Hamiltonian.
Twelve independent symbolic contraction/Einstein/fluid identities reduce to
zero in `remaining_field_algebra_check.json`.

This is substantive equation and domain work, not an evolved higher-order
physical solution. The exact current balance shows why a static metric with
`T = v t + pi(r)` and discarded background time-current cannot be the isolated
continuation of the rolling R1 state. The evolved material rows consume an EOS
or declared support relation; the prescribed common density alone does not
supply it. The environmental row consumes the three-dimensional scalar Hessian,
which projected convergence/shear does not determine. The overlap must join an
outer constrained Cauchy state in the same closed geometry. Extending only the
flat isolated Poisson equation to a compact slice fails its integral
compatibility condition for a positive uncompensated excess; this is a test
of that naive extension, not a prohibition of local positive mass in the
full closed Einstein problem.

These equations have been derived and algebraically checked. A physical
continuation still requires a declared material EOS/support, environmental
density and Cauchy data, and outer homogeneous/compensation and matching
prescription before coupled evolution and finite-distance ray propagation.
No new vacuum, stochastic factorization, moving-seam law or lensing amplitude
has been selected to fill those prescriptions.

## Local validation and reproducibility

Run from the repository root:

```powershell
python scripts/validate_strong_lensing_bhsm_integration.py
python scripts/validate_strong_lensing_bhsm_integration.py --assess-only
python scripts/build_strong_lensing_manuscript.py
```

The first command executes the nine supplied numerical/review/plot
commands once, routing results into the integration replay directory.
The second reuses those receipts without repeating the calculations.
The actual local environment is Python 3.14.3, NumPy 2.4.6, SciPy 1.17.1,
and Matplotlib 3.11.0. All nine commands succeeded; 1,161 numerical values
were compared and all 16 focused equation/rank checks passed. All 12
historical numerical and figure files remain byte-identical. The two
scientific numerical modules remain byte-identical to the supplied package.
An AST audit distinguishes output-routing/serialization changes from the
unchanged scientific applications.

The replay screening radius is 212.0292249771777 kpc; the largest
headline delay difference from the supplied results is 3.07e-10 seconds.
This is platform roundoff, not an adjustment of the scientific values.
The integrated scalar flux residual is 8.30e-16; the independent
projection discrepancy is 4.41e-13; the lens-equation residual is
8.47e-22 radians. The finite-body GR2 optical remainder is 1.71e-16
radians and the mixed stationary-action drift discrepancy is 6.15e-10
seconds. Near-zero diagnostics are checked using explicit absolute
tolerances; their large relative changes do not imply changed physics.

Numerical checks concern retained equations. They do not bound omitted
scalar higher orders, modified-fluid corrections, environmental mixed
terms, finite-distance cosmological effects or global matching. The
common GR2 result tests one specified ordinary-gravity correction; it
is not an all-orders scalar or environmental error budget. Historical
75-test evidence is reused for unchanged baseline implementations.
The seven new source checks and twelve symbolic continuation identities are
separate from that reused evidence and from the nine-command package replay;
they do not certify an evolved coupled continuation.

## Focused PDF and CI

The focused build targets `manuscript/strong_lensing_main.tex`, using
the actual `Definitions/mdpi.cls` with its `galaxies` option and real
dependencies. The verified class, original logos and supplied converted
logos are unchanged. When Ghostscript is available the build converts
the original assets; otherwise it uses the supplied converted official
assets. It never substitutes a class stub. The focused GitHub workflow
builds this target explicitly, independently of the other manuscript.

The final page count, visual audit, source hashes and build diagnostics
are recorded in the regenerated integration publication manifest. The author
confirmed reuse of his prior declarations: Independent Researcher,
Oconomowoc, Wisconsin, USA; ORCID `0009-0000-6650-3485`; no external funding,
no conflicts and the existing sole-author CRediT roles. He supplied
`carberry8878130@hotmail.com` and confirmed that the manuscript is neither
published nor under consideration elsewhere. These declarations are now
written in the manuscript and unsent cover letter. Historical AI model labels
are attributed to their records; substantive ChatGPT/Codex assistance is
disclosed without an unsupported current-model assertion. The calculation is
identified by science commit `0a329fa13f6c90bb2e3b5999b0155aa1e81fe8bd`;
no release DOI or unmade journal submission is invented.

The final focused PDF has 18 pages, including the source application and
continuation appendices. It uses the real Galaxies class and dependencies;
the final log has no overfull/underfull boxes,
undefined citations/references or duplicate equation destinations.
MiKTeX's package installation was enabled for this build command only,
without changing global settings. All 18 pages were rendered and
inspected as a contact sheet; pages 1, 3, 5, 6, 7, 8, 10, 14, 15, 16, 17
and 18 were also inspected at full-page size. The source table, full action/source
equations, image/delay and GR2 tables, figure, equivalence proof,
appendices, declarations, journal metadata and references were also
inspected. The visible correspondence email and science-commit citation were
also verified from the final PDF text. All 32 original displayed equation blocks
and every original numerical literal remain present. The added inherited
BHSM source definition is explicitly distinguished from a completed
constitutive reduction. Official MDPI asset bytes are preserved by a
narrow Git attribute so their published hashes remain valid on checkout.

The previously recorded broad numerical CI `background` import failure
is separate from this integration. Its old failing job was not used as
evidence about the focused manuscript or the new scoped calculations,
and no unrelated CI repair was undertaken. Actual remote outcomes are
reported in the PR and final delivery rather than inferred from another
target's build.

## Next specific physical application

The consumed R1 source, nonlinear flux and shell relation have now been
derived and applied. The independent remaining BHSM relation is an
action-owned physical embedding and source pullback from the cap/surface
sector into canonical `T`, including its cubic vertex, complement spectrum,
length normalization and common boundary domain. The explicit positive-lift
witnesses show which data the leading effective coefficients do not select.
The target phase, rolling scalar and metric coordinate are not identified by
name. For a field continuation of the retained effective action, prescribe
the EOS/support and environmental/outer Cauchy data identified above, solve
the coupled component equations, and then propagate the prepared light through
that matched metric. The effective-mass equivalence represents the computed
leading response; neither the derived continuation equations nor the present
source checks claim that evolved physical solution has already been obtained.
