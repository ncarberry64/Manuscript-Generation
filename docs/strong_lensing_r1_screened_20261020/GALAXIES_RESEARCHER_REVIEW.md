# Researcher Review and Corrections for Galaxies

> **Later observed-data revision:** the [HE 0435-1223 study](../STRONG_LENSING_OBSERVED_20261020.md) and [observed release records](../../artifacts/strong_lensing_observed_20261020/README.md), science `14ad777`, now provide the active empirical checkpoint. This review and its subsequent 19-page control links remain historical. They do not establish the later observed-test result or review of the 24-page manuscript.

> This report records the earlier package/integration review. Its concrete
> science pin and 18-page PDF inspection below apply to that snapshot. The
> [acceptance checkpoint](../STRONG_LENSING_ACCEPTANCE_20261020.md) and
> [current researcher receipt](../../artifacts/strong_lensing_acceptance_20261020/researcher_review.json)
> review the subsequent transfer/inference revision, finite nonlinear re-fit,
> matched smooth-profile and distributed-encounter controls. Earlier
> validation is reused within its scope; it is not reattributed to a later PDF.

## Review decision

The manuscript is a **BHSM strong-lensing paper** for **Galaxies**, associated with **“Strong Gravitational Lens Modeling.”** Its working title is **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Effective Calculation**. It asks: **How do BHSM surface contributions determine gravitational lensing and time delays for light traversing an envelope event?** This is an internal, AI-assisted research review and an author-review revision. It is not a journal referee report, an acceptance decision, or a claim that independent human peer review has occurred.

BHSM determines the research question and physical interpretation. The adopted envelope premise makes incoming-to-outgoing light transfer the central object. The already specified pure cubic R1 effective action supplies a **nonzero, normalized conditional source-to-optics calculation**. A positive conventional effective-mass representation reproduces the retained leading gravitational and optical response. This constructive equivalence advances the BHSM program without requiring a distinguishable departure from conventional gravity. The coupled effective source and its material flux have now been explicitly applied; the physical surface-to-R1 embedding remains unselected, with an explicit positive-lift nonuniqueness proof identifying what the retained coefficients do not fix.

The manuscript proceeds from envelope state and environmental preparation to action, normalization, source response and boundary data; then metric response, image directions and observer delay differences; finally conventional representation and checks. The adopted BHSM interpretation, derived effective equations, retained approximations and numerical results have separate status. The rolling scalar, target phase and physical metric coordinates retain their own independent variables and normalizations. An effective mass and its equivalent scalar contribution are counted once.

No submission target date is included. The [integration checkpoint](../STRONG_LENSING_BHSM_INTEGRATION_20261020.md) records current worktree validation and publication status; package evidence below remains historical to its original snapshot.

## 1. Evaluated envelope response and fixed-body reference

The local comparison keeps the finite matter profile, physical source displacement, measured redshifts, incident-light preparation, observer clock, and compatible background distance maps fixed. The lens has mass \(10^{11}M_\odot\), material radius \(3\) kpc, and redshift \(0.3\); the source redshift is \(1\). The R1 parameters and background normalization are reused without fitting a new local scalar amplitude.

| Observable | Fixed-body leading reference | Sourced R1 response | Offset from fixed-body reference |
|---|---:|---:|---:|
| Positive outer image | +0.835476034 arcsec | +0.835482157 arcsec | +6.1233 microarcsec |
| Central image | −0.082161940 arcsec | −0.082161152 arcsec | +0.7881 microarcsec |
| Negative outer image | −0.657495525 arcsec | −0.657499874 arcsec | −4.3486 microarcsec |
| Central-minus-positive delay | 21.716076566 days | 21.716377020 days | +25.9593 s |
| Negative-minus-positive delay | 14.271452387 days | 14.271552491 days | +8.6490 s |
| Ratio of these delays | 1.521644467 | 1.521654847 | +6.8212 ppm |

The force law also predicts a screening radius of **212.0292 kpc**, an asymptotic fractional force enhancement of **0.001507886**, and the fixed-epoch scaling \(r_V\propto M^{1/3}\). At the material surface, the local force correction is 5.0671 ppm and the projected deflection correction is 11.7467 ppm. The distinction matters because a light ray samples the extended field.

These are calculated values of a specified approximation. The printed digits support numerical reproduction; they are not statements of equally precise astrophysical accuracy.

The leading delay offsets are `25.95925815 s` and `8.64901938 s`. They compare the sourced response with the specified body alone. They are not an unavoidable difference between fully matched conventional and scalar representations: the positive effective mass described below reproduces the same retained leading response.

### Why the previous zero did not preclude this result

The earlier constant \(S^6\) target state has zero derivative current and stress in its source-free regular boundary problem. Its additional optical response remains zero.

The new calculation uses the separate rolling scalar \(T\) already present in the R1 covariant action. Its curvature coupling, combined with the Einstein constraint, supplies a matter source. No identification of \(T\) with the target phase, no duplicate metric scalar, and no new phenomenological scalar-matter amplitude is introduced.

The screened spherical branch is an established cubic-Galileon mechanism. The contribution here is the frozen R1 normalization, the finite-body optical application, and the explicit equivalence and identifiability analysis. References to the existing kinetic-braiding and Galileon-lensing literature have been added.

### Evaluated source and worked continuation equations

The [common-source derivation](BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
eliminates the same metric constraint to obtain
`J_eff = J_T - b_T(rho+rho_Sigma)/(2 Mpl²)` and the nonlinear canonical
scalar flux. On the already computed finite-body section it evaluates
`pi_T'(R) = -2.119825778845613e-4 Mpl H0` and the radial source integral
`integral_0^R r² J_eff dr = -1.216417408185084e-13` in the declared natural
units. The maximum scalar-flux relative residual is `4.15e-16`; the original
metric Gauss residual is zero at recorded precision; the material flux-jump
relative residual is `6.23e-16`. The exact-GR reference junction and the
separate leading-R1 junction have distinct extrinsic-curvature values and
zero jumps; the specified material seam charge is zero while the R1 volume
gravitational source remains nonzero. No exact modified-TOV matching is
asserted. These are explicit
applications of the source equations, not null fields inserted because a
microscopic source array was unavailable. Seven focused source tests passed.

Two positive same-domain quadratic lifts with the same frozen leading
`C,D` give different fourth-order coefficients and expansion scales.
They prove that the retained leading effective inputs do not uniquely
select the complement spectrum, physical length or embedding. Neither lift
is a selected galaxy, and neither adds an amplitude to the optical result.
The inherited nonzero `Theta/r²` seam source also cannot be identified with
the isolated R1 source by a constant field rescaling, because their flux
tails differ. This excludes that direct map; it leaves other action-owned
exports or nonlocal maps possible.

The [remaining-field derivation](R1_REMAINING_FIELD_DERIVATION_20261020.md)
now gives the exact cubic stress, scalar current, coupled spherical
Einstein/scalar rows, fluid continuity/Euler rows, environmental cross-Hessian,
overlap coefficients, matching conditions and null Hamiltonian. Twelve
independent symbolic identities pass. The rolling current shows why a frozen
exactly static isolated extension fails. The environmental equation consumes
the three-dimensional scalar Hessian; an optical convergence/shear does not
fix it. A compact-domain compatibility check also excludes transplanting the
isolated flat Poisson equation unchanged onto the closed slice.

These are derived and checked continuation equations. No new evolved
higher-order body, nonlinear environment or globally matched light-transfer
solution is claimed. That application consumes a material EOS/support relation,
environmental density and Cauchy data, and outer homogeneous/compensation and
matching prescription. The independently remaining microscopic relation is the
physical action embedding/source pullback into canonical `T`, its cubic vertex
and complement domain. The derivation has been advanced to these explicit
equations and nonuniqueness results; it is not left as an unworked source list.

## 2. Substantive corrections implemented

| Research concern | Correction in the revised manuscript |
|---|---|
| The old abstract centered on an unresolved local encounter and a null control. | Replaced with an approximately 200-word abstract centered on the computed local response and its observational-equivalence limit. |
| The relationship between the BHSM question and the R1 calculation could be misstated. | BHSM organizes the paper; the derived source and optics belong to the declared conditional R1 effective action. The source application and positive-lift proof now identify the remaining physical embedding and complement data. |
| Several different scalar or phase variables could be conflated. | Distinguished the independent rolling \(T\), target phase \(\eta\), and metric deformation \(\psi_{\rm geo}=-\Psi\). The detailed derivation also distinguishes \(\delta T\) from a cosmological time-displacement variable. |
| A matter source could be appended with a chosen normalization. | Derived the scalar source from the same action and metric constraint; retained the scalar-metric mixing contribution to the stiffness \(C\). |
| Exact GR fluid equations could be mistaken for exact modified-gravity equilibrium. | Retained the same prescribed density, calculated the leading pressure change, and stated that neither an exact modified TOV sequence nor a common fixed material equation of state has been constructed. |
| A constant background was described too broadly as having zero self response. | Distinguished zero background derivative stress from the positive stiffness of a nonconstant Dirichlet phase perturbation. |
| The nuisance equation omitted its baseline. | Defined residual data explicitly as \(\delta y=y-y_0\), specified a positive covariance for whitening, and explained what the projected rank means. |
| A one-column rank could be read as identification of a fundamental free parameter. | Identified the column as a fixed model-difference template and the surviving rank as a residual direction. |
| A delay ratio could appear to solve all degeneracies. | Explained that a common mass-sheet factor cancels, while source, center, radial-profile and environmental freedoms remain. |
| The conventional representation was absent from the central result. | Added the positive effective-mass construction and its proof to the main text as a constructive representation of the same retained response, including all source planes and leading matter dynamics. |
| Numerical residuals could be mistaken for physical accuracy. | Separated equation and quadrature errors from weak-field, finite-distance, environment and matching approximations; derived and applied the same finite-body ordinary-GR 2PM correction to both models, including its central-ray and delay-potential contributions. |
| Positive principal coefficients could be presented as a complete stability certificate. | Reported the positive leading coefficients together with the radial scalar characteristic speed squared of about 1.474 and the limited scope of this check. |
| A zero material seam charge was described as a zero exported source. | Corrected the technical source document to say that the specified material-surface seam contribution vanishes; the separate gravitational source is not thereby set to zero. |
| The draft contained internal project and checkpoint narrative. | Reorganized around the question, model, methods, results, degeneracies, discussion and conclusions. Detailed historical controls remain in the repository and appendices. |

## 3. What the observations constrain and what remains interchangeable

### Original single source

Allowing the source position, body mass and radius, a compensated mass-sheet scale, and the two center coordinates absorbs the linear screened-model difference: the projected rank is **zero**. Fixing the center leaves a nonzero residual, but that is a restriction on the nuisance model.

The forward calculation establishes delay offsets relative to the specified fixed-body reference with the same optical preparation. It does not establish a discrepancy between fully matched conventional and scalar representations. Their leading equivalence and the rank-zero result describe what observations constrain within the specified model; they do not constitute a failed BHSM prediction.

### Three-source finite-profile design

A bounded extension uses source redshifts \(0.8,1.2,1.8\), independent source coordinates and delay scales, shared body parameters and center, and physical convergence/shear nuisance columns. The fixed scalar response retains one projected direction even after independent source-plane convergences are added.

Only **0.1202%** of the raw whitened response remains in the expanded nuisance model. For the stated illustrative weighting, the residual coordinates are at the scale of **6 nanoarcseconds** and the residual delays at **0.086 seconds**. These are mathematical residuals, not an instrument forecast.

The tide columns are phenomenological optical controls. They do not constitute a solved nonlinear scalar field for an arbitrary external matter environment.

### Flexible conventional halo

Let \(f_M(r)=M(<r)/M\) and \(F(r)=\delta\Phi'/\Phi'_N\). Define

\[
\delta M_{\rm eff}(r)=M f_M(r)F(r),\qquad
\delta\rho_{\rm eff}(r)=\frac{M}{4\pi r^2}
\frac{d[f_M(r)F(r)]}{dr}\ge0.
\]

The revised article gives the positive-density proof. Because \(\Phi=\Psi\) at the retained order, this same conventional additional mass produces the complete scalar optical field, including relative delays and multiple source planes. Its asymptotic extra mass is

\[
\delta M_{\rm eff}(\infty)=1.50788624\times10^8M_\odot.
\]

This is broader than the mass-sheet degeneracy. Slow matter also responds to the same retained metric, so kinematics alone does not automatically distinguish an unrestricted effective halo from this scalar representation. Observations constrain the retained metric response. They constrain an underlying surface-state assignment only after its source map and any independent matter/environment information are supplied. A separately derived interface response could add observables, but no such response or distinctive signal is assumed here.

## 4. Optical accuracy: implemented correction and result

The earlier fluid calculation solved the GR matter equilibrium, but its optical encounter used a leading weak-field potential. At the two exterior image impacts, the ordinary Schwarzschild second-order term is 3.68 and 4.68 ppm of the leading conventional deflection, making it relevant at the scale of the reported scalar image shifts.

**The requested correction has now been calculated and applied to both models.** The finite-body second post-Minkowskian (2PM) deflection was derived from the exact Tolman–VII optical scattering integral. It includes the central ray through the fluid, the leading TOV pressure, and the conversion between areal radius and asymptotic impact. The exterior reduces to the standard expression

\[
\delta\hat\alpha_{\rm GR2}
=\frac{15\pi}{4}\left(\frac{m_g}{b_{\rm imp}}\right)^2.
\]

The analytic deflection is integrated to one consistent plane potential. The same potential, gradient and Hessian are then included in both the conventional and scalar-on image equations and arrival functional. The scalar parameters and leading radial profile remain fixed.

### Common GR2 control relative to the same fixed-body reference

| Sourced-minus-fixed-body observable | Leading models | Same GR2 term added to both | Change in their difference |
|---|---:|---:|---:|
| Central-minus-positive delay | 25.95925815 s | 25.95932714 s | +68.99 microseconds |
| Negative-minus-positive delay | 8.64901938 s | 8.64900924 s | −10.14 microseconds |
| Positive image shift | +6.12327625 microarcsec | +6.12326934 microarcsec | −0.00000691 microarcsec |
| Central image shift | +0.78810716 microarcsec | +0.78810711 microarcsec | −0.0000000558 microarcsec |
| Negative image shift | −4.34859153 microarcsec | −4.34858643 microarcsec | +0.00000509 microarcsec |

The conventional delays themselves move by **8.28473 s and 2.55245 s**, and the conventional coordinates by **+1.35456, −0.01519, −1.72123 microarcsec**. These absolute changes largely cancel between the matched models. The scalar delay differences change by only +2.66 and −1.17 ppm of their original values. The delay-ratio response remains 6.8212 ppm.

Direct exact static-GR null-geodesic integration checks the ordinary correction independently. At the three original impacts, the remainder after subtracting the 1PM and 2PM terms is below \(1.8\times10^{-16}\) rad for the actual compactness. An interior compactness extrapolation also recovers the analytic 2PM coefficient, and direct quadrature verifies the deflection-potential primitive.

This resolves the identified omission as a **common-conventional-term robustness check**. It does not complete scalar-dependent 2PM terms, scalar changes in the higher-order fluid geometry, or an exact finite-distance cosmological arrival calculation. The asymptotic scattering impact is identified with the physical thin-plane impact. The small numerical drift is not a bound on all omitted physical effects.

The positive effective-halo theorem and nuisance-rank statements are leading-order results. This common GR2 control neither extends the equivalence to all orders nor supplies a demonstrated way to break it.

The full derivation is in `FINITE_BODY_GR_2PM_DERIVATION_20261020.md`, the four-model calculation is in `gr2_control.json`, and the independent geodesic and composition checks are in the review outputs.

## 5. Verification already completed

The checks address distinct failure modes rather than simply counting tests:

- The integrated scalar flux equation has a maximum relative residual below \(8\times10^{-16}\).
- Independent direct line-of-sight and angular-substitution projections agree to \(4.5\times10^{-13}\) relative for the additional deflection.
- The re-solved image equation has a residual below \(8.5\times10^{-22}\) rad.
- Independently differentiating the potential and deflection verifies their gradient/Hessian consistency.
- The stationary delay response and finite re-solved response agree through the expected perturbative order; their differences are about 0.00016 and 0.000056 seconds.
- The common finite-body GR2 correction passes an exact static-geodesic check, an independent primitive/derivative check, and an independently implemented four-model comparison.
- The portable coefficient, projection and GR2 scripts use repository-relative paths.
- One portability replay of the three independent review scripts reproduced all 305 saved numerical values exactly in the working environment.

The supplied package analysis used Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. These identify that execution environment, not uniquely required or latest versions. Current worktree replay, focused-test and PDF evidence belongs to the integration checkpoint; the package evidence below remains attached to its original snapshot.

## 6. Galaxies preparation

The final destination is **Galaxies**, and the cover-letter draft names its “Strong Gravitational Lens Modeling” special issue. No submission target date is supplied.

The package records its check of the official instructions: free-format initial submission, a single-paragraph abstract of about 200 words maximum, 3–10 keywords, a reproducible source package and appropriate declarations. The source uses the official numeric-reference MDPI class with the Galaxies journal option and a review-draft footer. The current focused build and PDF inspection are recorded separately in the integration checkpoint.

Substantive AI assistance is disclosed in computational methods and acknowledgments. The disclosure includes mathematical and computational assistance, not merely language editing.

Official sources:

- [Galaxies author instructions](https://www.mdpi.com/journal/galaxies/instructions)
- [Strong Gravitational Lens Modeling](https://www.mdpi.com/journal/galaxies/special_issues/50J931CL3Q)
- [MDPI LaTeX instructions](https://www.mdpi.com/authors/latex)
- [Official numeric-reference template](https://res.mdpi.com/data/MDPI_template_ACS.zip?v=20260911)

## 7. Confirmed author declarations and publication identity

The author directed reuse of his other manuscripts' declarations and then
supplied the remaining correspondence and publication facts. The active
manuscript and unsent cover letter now identify Norman P. Carberry,
Independent Researcher, Oconomowoc, Wisconsin, USA; ORCID
`0009-0000-6650-3485`; correspondence `carberry8878130@hotmail.com`.
The declarations state no external funding, no conflicts and the confirmed
sole-author contribution roles, with substantive AI assistance disclosed
separately. The author confirmed that this manuscript is not published or
under consideration at another journal. Those facts are no longer placeholders.

The earlier public commit `49a1dab` supplies the prior reference calculation;
it does not contain the additions. The integrated calculation is pinned at
science commit `0a329fa13f6c90bb2e3b5999b0155aa1e81fe8bd`, containing 68 new
calculation, derivation, input, output and check files. The regenerated
manuscript cites that concrete revision in its bibliography and
data-availability statement. No archive DOI is claimed. This work updates
the existing draft PR; it does not merge it,
submit to the journal or send the cover letter. The physical source-map and
continuation choices above remain scientific scope, independently of these
completed author declarations.

## 8. Reproduction and integration

The companion package contains a recompilable manuscript source, scientific figures, new code and numerical outputs, review records, and a Codex integration prompt. The current [scope record](../STRONG_LENSING_SCOPE_20261020.md) supersedes earlier package framing wherever a distinctive departure was made the organizing objective. The retained technical derivations and numerical evidence keep their stated meanings.

The numerical additions are intended to be integrated into the existing Manuscript-Generation worktree based on commit 49a1dab or a compatible descendant. Existing source dependencies are identified by pinned revision and hashes. The package is not a duplicate repository and does not silently replace concurrent work.

From the repository root after integration:

~~~bash
python scripts/run_strong_lensing_r1_screened.py
python scripts/r1_screened_review/r1_qs_coefficients.py
python scripts/r1_screened_review/check_screened_projection.py
python scripts/r1_screened_review/multisource_crossratio.py
python scripts/r1_screened_review/shared_lens_rank.py
python scripts/run_strong_lensing_gr2_control.py
python scripts/r1_screened_review/check_gr_optical_control.py
python scripts/r1_screened_review/check_common_gr_composition.py
python scripts/plot_strong_lensing_r1_screened.py
latexmk -pdf -interaction=nonstopmode -halt-on-error -cd manuscript/strong_lensing_main.tex
~~~

At package creation, the previous CI job compiled the broader main.tex; its success did not independently establish a successful build of the focused strong-lensing manuscript. The package supplied a dedicated focused build definition and recorded a 13-page local build with no overfull boxes or undefined references. It also recorded rendered visual inspection of all pages and independent checking of equation and appendix pages. The package's focused GitHub workflow had not run at that historical point. The current integration checkpoint records the actual worktree PDF and CI outcomes after reconciliation.

The final current focused PDF has **18 pages**, built with the real Galaxies
MDPI class and dependencies. Its final log has no overfull/underfull boxes,
undefined references/citations or duplicate destinations. All 18 rendered
pages were inspected as a contact sheet, with full-page inspection of pages
1, 3, 5, 6, 7, 8, 10, 14, 15, 16, 17 and 18. Equations, result/source and
GR2 tables, figure, new source/continuation appendices, declarations,
references and journal metadata were checked. The correspondence email is
visible in the PDF; the science-commit citation is present. These current
checks supersede the package-era page count without rewriting its evidence.

At package creation, the connected GitHub integration rejected the attempted publication write with HTTP 403, “Resource not accessible by integration.” That historical attempt made no branch or PR update. The present worktree integration uses the existing draft PR and records its own resulting commit and remote outcome separately; the old package receipt is not current publication status.

## Primary scientific sources added or retained

1. Fleury, Larena and Uzan, *Gravitational lenses in arbitrary space-times*, Class. Quantum Grav. 38, 085002 (2021), [doi:10.1088/1361-6382/abea2d](https://doi.org/10.1088/1361-6382/abea2d).
2. Fleury, Larena and Uzan, *Line-of-sight effects in strong gravitational lensing*, JCAP 08, 024 (2021), [doi:10.1088/1475-7516/2021/08/024](https://doi.org/10.1088/1475-7516/2021/08/024).
3. Teodori et al., *Comments on the mass sheet degeneracy in cosmography analyses*, JCAP 07, 027 (2022), [doi:10.1088/1475-7516/2022/07/027](https://doi.org/10.1088/1475-7516/2022/07/027).
4. Deffayet et al., *Imperfect dark energy from kinetic gravity braiding*, JCAP 10, 026 (2010), [doi:10.1088/1475-7516/2010/10/026](https://doi.org/10.1088/1475-7516/2010/10/026).
5. Barreira et al., *Nonlinear structure formation in the Cubic Galileon gravity model*, JCAP 10, 027 (2013), [doi:10.1088/1475-7516/2013/10/027](https://doi.org/10.1088/1475-7516/2013/10/027).
6. Barreira et al., *Galaxy cluster lensing masses in modified lensing potentials*, MNRAS 454, 4085–4102 (2015), [doi:10.1093/mnras/stv2211](https://doi.org/10.1093/mnras/stv2211).
7. Iyer and Petters, *Light’s bending angle due to black holes: From the photon sphere to infinity*, Gen. Relativ. Gravit. 39, 1563–1582 (2007), [doi:10.1007/s10714-007-0481-8](https://doi.org/10.1007/s10714-007-0481-8).

The source action, frozen parameters and prior control results are pinned to the public [Manuscript-Generation revision 49a1dab](https://github.com/ncarberry64/Manuscript-Generation/tree/49a1dabbe148f8ba9da4ac079a8d93fe97c5de76). The full microscopic BHSM source is identified separately in the action-source manifest.
