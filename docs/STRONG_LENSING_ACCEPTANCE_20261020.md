# BHSM strong-lensing transfer and identifiability checkpoint

> **Historical 19-page synthetic/control checkpoint.** The observed HE 0435-1223 application and new manuscript title now supersede this document as the active checkpoint. See the [observed study](STRONG_LENSING_OBSERVED_20261020.md) and [current release artifacts](../artifacts/strong_lensing_observed_20261020/README.md), science `14ad777`. All values, source pin and PDF evidence below retain their original historical scope; they are not evidence for the later 24-page observed-light revision.

This checkpoint continues draft PR
[#23](https://github.com/ncarberry64/Manuscript-Generation/pull/23) in the
existing `codex/strong-lensing-first-comparison-20261020` worktree. The active
paper is **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1
Effective Calculation**, for **Galaxies**, special issue **Strong
Gravitational Lens Modeling**. No submission target date is set. The cover
letter remains unsent; no journal submission or PR merge is authorized here.

## Contribution and manuscript order

The paper develops one conditional forward and inverse argument:
**specified envelope state and source -> normalized metric/light transfer ->
image and delay observables -> identifiability under ordinary lens freedoms**.
BHSM supplies the organizing physical interpretation. The retained cubic R1
action supplies the effective source, normalization, regular screened branch
and evaluated optical response. Conventional gravity supplies both the
declared benchmark and a positive effective-mass representation of the same
retained response.

The contribution is the controlled application and its explicit inference
limit. Cubic screening, Fermat optics and mass-sheet freedom are established
methods; no universal novelty claim is made. The
[literature and call assessment](strong_lensing_r1_screened_20261020/LITERATURE_AND_CALL_POSITIONING_20261020.md)
verifies primary references, the issue identity and an actual restricted
forward-model methods article in that issue. The full invitation body was
not retrievable, so no reconstructed call wording or current submission-window
status is asserted.

The focused main now follows source, regular branch, metric/optics,
state-coefficient partials, images/delays, positive effective mass, finite
inverse modeling, multi-source scope, and profile/accuracy controls. Long
cap/complement/source and exact-continuation material is retained in the
[technical supplement source](../manuscript/sections/strong_lensing_technical_supplement_20261020.tex)
and readable supplementary records. All original 32 package equation blocks
remain in the main text and appendices A–C. This reorganization does not
discard those equations or promote supplementary derivations into new optical
predictions.

## Executed applications

The original preparation is retained: the injected pulse, finite body,
physical source coordinate, measured redshifts, observer clock and compatible
distance maps. The compact reference has `M=1e11` solar masses and `R=3 kpc`.
The screening radius is `212.02922497712947 kpc`. Its leading sourced-minus-
fixed-body delay offsets remain `25.959258149566722 s` and
`8.649019375383205 s`. Those are offsets from the specified body alone, not
unavoidable differences between fully matched conventional and scalar
representations.

| Application | Actual result | Meaning and domain |
|---|---|---|
| Finite conventional inverse model | With mass, scale, source position, center and one physical sheet free, the original-row norm is `1.38947e-13`; independently re-solved branch delays differ by at most `9.21e-10 s`. | Five coordinates fit five axial observables: a saturated local inverse example, using the same clock and maps. Not a likelihood, global root census or population result. |
| Fixed-center inverse control | Delay remainders are `+0.00330640 s` and `-0.00503119 s`. | A four-coordinate restricted-profile fit under declared numerical row scales; not an irreducible signal or detection forecast. |
| Smooth Plummer control | Same mass and three-dimensional half-mass radius `1.929416346490062 kpc`; scale `a=1.478745083254196 kpc`. Reference-relative offsets are `22.450967893860252 s` and `9.620034449068271 s`. | Smooth infinite-support finite-mass control, with its own tail integration. No artificial edge at the compact body's radius. |
| Ordinary profile sensitivity | Baseline conventional delays change by `-167124.924460958 s` and `-124500.39008295165 s`; the scalar offsets change by `-3.50829025570647 s` and `+0.9710150736850665 s`. | Controlled sensitivity at fixed mass and half-mass radius, not observational error bars, a prior or a posterior. |
| Eight matched profile/input cases | Both profiles, plus one-at-a-time `+20%` mass, half-mass-radius and source-displacement controls, retain three images. | Each scalar/body comparison shares its own preparation. No parameter is fitted to strengthen an offset. |
| Direct effective-mass projection | Independent additional-deflection agreement within `5.27e-13` relative and matched delays within `6.14e-10 s`. | Positive effective-mass equivalence survives both profiles at the retained leading order, and describes the same field once. |

The [finite-refit derivation](strong_lensing_r1_screened_20261020/NONLINEAR_LENS_REFIT_20261020.md)
gives the original lens and delay rows, bounds and conditioning scales. The
physical sheet potential is `s |x-D_ol theta_c|²/2`, with
`s=kappa/(D_ol q)` and `q=D_ls/D_os`. Its gradient supplies the one
`kappa(theta-theta_c)` angular lens term, and its potential supplies the
same observer-clock delay. No independent time scale is fitted.

The [profile and coefficient record](strong_lensing_r1_screened_20261020/LENS_MODEL_ROBUSTNESS_20261020.md)
gives the profile formulas, pressure-support scope and projected derivatives.
For `S=sqrt(1+r_V³ f_M/r³)`, at fixed body, impact and `lambda`, the scalar
force partials are `d ln(delta Phi')/d ln b_T=3/2+1/(2S)` and
`d ln(delta Phi')/d ln C=-1/S`. Projection is a positive line-of-sight weighted
mean, verified by finite coefficient changes. These are coefficient
sensitivities, not independently selected physical envelope states. A state
derivative still composes them with the physical state-to-coefficient
Jacobian.

The all-ray effective-mass equivalence follows from the field construction;
it is stronger than what the saturated five-row fit alone proves. The
multi-source finite-profile remainder remains restricted to its declared
nuisance family. None of these inverse applications establishes unique
microscopic state recovery.

## Practical endpoint: a predictive optical map for each envelope

Existing lens observations can provide an empirical baseline for the same
forward/inverse problem. A prospective application would calibrate the
conventional benchmark and conditional BHSM response jointly, stating the
admitted mass-profile, source, sheet, environment and cosmological freedoms
and which constraints come from independent data or priors. Freeze that
calibration, its uncertainty model and its prediction protocol before using
held-out observations or new monitoring to test image, delay or other
actually modeled light-transfer predictions. New observations can then
refine the envelope's optical response map under those stated assumptions.

This is a proposed empirical workflow, not an application performed in the
present synthetic study. A separately inferred optical map for each envelope
does not imply a unique microscopic envelope state: observationally
equivalent mass, source, geometry or state combinations must remain visible.
The positive effective-mass construction and finite conventional re-fit show
why that distinction is necessary. State uniqueness would require additional
independent observables and an action-owned state-to-response map, followed
by an identifiability demonstration in the admitted physical domain.

## Accuracy, source and physical boundaries

The common finite-body ordinary-GR second-order control is preserved,
including the central ray. It changes the two fixed-reference scalar offsets
by approximately `+68.99` and `-10.14` microseconds. This tests a common
conventional term; it is not a bound on every omitted scalar, finite-distance,
material or environmental contribution.

The [distributed-encounter application](strong_lensing_r1_screened_20261020/ENCOUNTER_ACCURACY_CONTROL_20261020.md)
solves `x''+kappa x=-2 Phi'(sqrt(x²+z²)) x/sqrt(x²+z²)` between separately
derived closed-FLRW exterior wavefronts. It fixes the same physical source
and observer clock, and compares a thin and distributed encounter in the
same field window. At a `20 Mpc` half-window, the preferred stationary
optical-homotopy application gives delay-offset changes
`-0.000619043479865411 s` and `+0.000001006405318193515 s`. Its integration
coordinate interpolates two fixed optical fields; it is not a fitted
amplitude or a family of selected physical states. Direct finite-arrival
subtraction agrees within `7 ns`. The isolated additional-field,
fixed-impact Born tail beyond the window has delay bounds of approximately
`35.04` and `13.47` microseconds. This narrow bound does not cover arbitrary
environmental fields, a nonlinear re-solved tail or global matching.

The distributed calculation uses the retained leading weak-field paraxial
action. It does not include GR2, scalar-dependent higher-order metric terms,
time evolution or a completed finite-distance null solution. Its first
delay-pair width correction exceeds the separate common-GR2 drift, making
the distinction between numerical digits and physical accuracy concrete.
The two controls are not added as a complete higher-order prediction.

The supplementary [source application](strong_lensing_r1_screened_20261020/BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
already evaluates `J_eff=J_T-b_T(rho+rho_Sigma)/(2 Mpl²)` and its nonlinear
flux on the computed finite-body section. Its source integral is nonzero;
zero smooth material seam charge does not remove the volume source. Two
positive same-domain lifts explicitly demonstrate what the leading
coefficients leave free. The [exact continuation](strong_lensing_r1_screened_20261020/R1_REMAINING_FIELD_DERIVATION_20261020.md)
derives scalar/metric/fluid, environmental cross-Hessian and matching
relations, with twelve checked identities. Neither record is an evolved
higher-order physical solution.

The microscopic surface embedding/source pullback, complement domain and
physical state-to-coefficient map remain unselected. A physical continuation
also consumes material EOS/support, environmental density/Cauchy state and
outer homogeneous/compensation/matching data. The illustrative profiles and
frozen effective coefficients do not secretly supply those relations. The
constant-target phase null remains a separate variable/boundary control.
The scalar response and its effective mass, and the local encounter and
reference propagation, are each counted once. No statistical state kernel,
Wick closure or later inclusive pair-production object is inferred from this
deterministic light calculation.

## Verification and provenance

Four finite-refit tests passed, covering the finite mass-sheet identity,
conventional injection recovery, original-row/independent-branch application
and fixed-center remainder. Seven profile/coefficient tests passed, covering
the profiles and projections, matched scale, preserved baseline, effective
mass, pressure/rank scope and finite derivative checks.
Seven distributed-optics tests also passed: individual physical wavefronts,
zero-field propagation, independent potential/deflection integrals, the
isolated tail bound, computed-section force/boundary balance and frozen
force reproduction. The maximum integrated momentum residual is below
`3.5e-16 rad`. Recorded collocation and radial-quadrature refinement changes
the direct stationary width response by less than `0.65 ns`; this is a
numerical convergence diagnostic for the retained optical action.

Independent AI-assisted code/equation and researcher reviews are saved in the
[acceptance artifact index](../artifacts/strong_lensing_acceptance_20261020/README.md).
They are not external human peer review. Original numerical replays and
earlier focused test evidence are reused; they were not repeated for these
document changes.

None of these provisional effective applications is physical admission. The
existing `2e-10` physical guard is unchanged, as are inherited information,
calibration and exact-statistics requirements outside the deterministic
calculation actually executed here.

The [previous integration checkpoint](STRONG_LENSING_BHSM_INTEGRATION_20261020.md)
pins science commit `0a329fa13f6c90bb2e3b5999b0155aa1e81fe8bd` and its
18-page PDF audit. Those records are historical to that completed snapshot,
and do not include the later acceptance code or certify the rewritten PDF.
The new applications retain source hashes and saved outputs at their actual
paths. The new immutable science checkpoint is
[`87503ae9d7672405e07b573fd140c784b8ce6cbd`](https://github.com/ncarberry64/Manuscript-Generation/tree/87503ae9d7672405e07b573fd140c784b8ce6cbd).
It contains the three new applications, eighteen new focused tests, four
output/review records, three method notes and the technical supplement.
The manuscript's bibliography and data-availability statement pin this
calculation revision. Subsequent manuscript and publication-record commits
are distinct from that immutable science checkpoint.

The current focused PDF has **19 pages** and **361639 bytes**. It was built
with the real Galaxies/MDPI class and dependencies. The build log has no
overfull/underfull boxes, undefined references/citations or duplicate
destinations. The release's PDF hash and rendered visual review belong to
the [publication manifest](../artifacts/strong_lensing_acceptance_20261020/publication_manifest.json)
and [contact sheet](../artifacts/strong_lensing_acceptance_20261020/pdf_contact_sheet.png).
The PDF SHA256 is
`d0ad684149c931e515731324d8771e0ff77b844fa7e8dfb8bec1dbcfc4a99051`.
All nineteen rendered pages were inspected in the contact sheet; full pages
1, 5, 6, 11, 12, 14, 16, 17 and 19 were checked without layout defects.
The author email, ORCID and current science hyperlinks were verified, and
the PDF metadata contains no replacement characters. The preservation audit
finds all 32 original equation/align blocks unchanged in the focused main
and no missing original numerical tokens. The manifest retains the twelve
official journal assets and twelve original numerical/figure-file checks.
The earlier 18-page audit is not reused as evidence for this PDF.

Author identity, correspondence, no-external-funding, no-conflicts,
sole-author roles and no-concurrent-publication facts are confirmed and
written in the manuscript and [unsent cover letter](strong_lensing_r1_screened_20261020/GALAXIES_COVER_LETTER_DRAFT.md).
No factual author placeholder or invented DOI remains. The separate Universe
manuscript and unrelated work are preserved.
