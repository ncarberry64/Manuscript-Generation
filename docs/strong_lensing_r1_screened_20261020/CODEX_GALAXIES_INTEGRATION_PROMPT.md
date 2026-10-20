# Author handoff: BHSM strong lensing and the R1 effective response

## Governing scope

Working title: **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Effective Calculation**.

BHSM organizes this paper. The research question is **How do BHSM surface contributions determine gravitational lensing and time delays for light traversing an envelope event?** The central object is incoming-to-outgoing light transfer, conditional on envelope state and surrounding geometry. Follow the order: envelope and environmental preparation; action, normalization, source and boundary problem; metric response; image directions and observer delays; then conventional gravitational representation and checks.

The R1 effective source and its optical application have been derived in their declared approximation. Conventional matter and scalar descriptions are interpreted as representations of the same underlying surface contributions. The positive effective-mass construction establishes their equivalence for the retained leading gravitational and optical response. This is a constructive BHSM result. No distinguishable departure is required for success, and the fixed-body reference offsets are not disagreements between fully matched representations. Count the effective mass and its equivalent scalar response once.

The adopted BHSM surface interpretation and the evaluated R1 equations have different evidential status. The coupled effective source, nonlinear scalar flux and material junction have now been derived and applied to the existing finite-body section. Two positive same-domain lifts prove that the frozen leading coefficients do not uniquely choose the physical surface embedding, complement spectrum or length normalization. The remaining microscopic relation is that action-owned embedding and source pullback into canonical `T`, its cubic vertex and boundary domain; it does not identify `T` with the target phase or metric coordinate.

The [common-source derivation](BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
and [continuation equations](R1_REMAINING_FIELD_DERIVATION_20261020.md)
are now concrete work products. Exact scalar/metric/fluid applications,
environmental cross-Hessian and outer matching conditions are derived and
checked. Their physical application still consumes a material EOS/support,
environmental density/Cauchy state and outer homogeneous/compensation and
matching prescription. Derived equations are not an evolved physical solution.

This current scope supersedes the package's earlier objective wherever that objective centered on distinguishing BHSM from conventional gravity. Preserve the package's technical evidence and original numerical meanings. The [scope record](../STRONG_LENSING_SCOPE_20261020.md) and [integration checkpoint](../STRONG_LENSING_BHSM_INTEGRATION_20261020.md) identify the current manuscript; the instructions below retain the reproducible integration contract.

## Current completed checkpoint

Science commit `0a329fa13f6c90bb2e3b5999b0155aa1e81fe8bd` pins 68 new
calculation, derivation, input, output and check files. The current paper's
bibliography and data-availability statement cite that exact revision.
The final focused PDF has 18 pages and was built with the actual Galaxies
MDPI dependencies. All 18 rendered pages and full pages 1, 3, 5, 6, 7, 8,
10, 14, 15, 16, 17 and 18 were visually checked. No overfull/underfull boxes,
undefined references/citations or duplicate destinations remain. The visible
correspondence email and science-commit citation were verified. The final
manuscript/publication commit follows the science commit; no new build or
numerical rerun is required solely to repeat this already recorded evidence.

Continue the existing strong-lensing work in `ncarberry64/Manuscript-Generation`, draft PR #23. The reviewed source baseline is commit `49a1dabbe148f8ba9da4ac079a8d93fe97c5de76`, on `codex/strong-lensing-first-comparison-20261020`. The existing Windows worktree is:

`C:\Users\carbe\Manuscript-Generation\.worktrees\strong-lensing-first-comparison-20261020`

The attached `Galaxies_Lensing_Research_Package.zip` contains a reviewed overlay, complete focused manuscript sources, a patch for existing text changes, and new reproducible calculations. The journal is **Galaxies**, special issue **Strong Gravitational Lens Modeling**. Do not put a submission target date in any manuscript, cover letter, PR description or handoff.

## Objective

Integrate the research correction into the existing draft PR, verify the numerical and manuscript outputs, and provide a concrete reviewed diff. Preserve concurrent work. Keep the PR draft; do not merge or submit to the journal.

## First inspect the actual state

1. Read the applicable repository instructions and the existing PR/worktree state. Verify the current head and whether it is a descendant of the reviewed baseline. Inspect all local changes before editing.
2. Read the package `README.md`, `MANIFEST_SHA256.json`, `overlay/docs/strong_lensing_r1_screened_20261020/GALAXIES_RESEARCHER_REVIEW.md`, and the R1 and GR2 derivations. Do not infer that the new files were already present in commit 49a1dab.
3. Compare the overlay to the actual tree. Apply the existing-file patch only where it matches; otherwise merge the intended changes manually. Add the new files at their repository-relative paths. Never overwrite a newer scientific correction or a user edit without reconciling its meaning.

## Scientific content to preserve

- The fixed-body reference and sourced R1 calculation use the same incoming light preparation, density, physical source, redshifts, clock and compatible background maps. Fully matched conventional effective-mass and scalar representations reproduce the same retained leading transfer.
- The separate rolling scalar T comes from the already specified pure cubic R1 action and its frozen background. Do not fit a new local scalar amplitude or silently identify it with the previously constant target phase.
- The earlier source-free constant-target family remains a null control. Its zero material-surface seam contribution does not set every separate gravitational source to zero.
- The declared R1 example has a screening radius of 212.0292 kpc and relative-delay offsets of approximately 25.95925815 s and 8.64901938 s from the specified fixed-body reference. The image offsets are approximately +6.12327625, +0.78810716 and −4.34859153 microarcsec.
- The finite-body GR2 control is applied identically to both calculations, with a consistent potential, deflection and Hessian. It includes the central ray. Re-solving gives sourced-minus-fixed-body delays of 25.95932714 s and 8.64900924 s. Their changes are about +68.99 and −10.14 microseconds.
- This tests common ordinary-GR correction robustness. It is not a derivation of scalar-dependent 2PM field terms, an exact modified TOV family, or a complete finite-distance cosmological arrival function. Keep the asymptotic-impact/thin-plane identification explicit.
- The unrestricted single-source nuisance response has rank zero. The three-source finite-profile residual does not imply general identifiability. These state what the specified observations constrain and what remains interchangeable. The positive effective-mass equivalence remains central, including its consequence for leading slow-matter dynamics, and is asserted only at the retained order.
- Numerical convergence is not physical validation or an astrophysical error bar. AI-assisted independent code checks are not external human peer review. Keep the principal-coefficient, characteristic-speed, environment and global-matching qualifications.

## Verification and build

From the repository root, run the exact relevant reproduction commands in the package README. Compare generated numerical arrays with the supplied receipts using meaningful tolerances and account for ordinary floating-point variation; do not require PDF byte identity across TeX installations. In particular run:

```text
python scripts/run_strong_lensing_r1_screened.py
python scripts/run_strong_lensing_gr2_control.py
python scripts/r1_screened_review/r1_qs_coefficients.py
python scripts/r1_screened_review/check_screened_projection.py
python scripts/r1_screened_review/multisource_crossratio.py
python scripts/r1_screened_review/shared_lens_rank.py
python scripts/r1_screened_review/check_gr_optical_control.py
python scripts/r1_screened_review/check_common_gr_composition.py
python scripts/plot_strong_lensing_r1_screened.py
python scripts/build_strong_lensing_manuscript.py
```

Use a complete TeX Live/MiKTeX installation with the actual MDPI dependencies. The build script converts the official EPS logos with Ghostscript and builds `manuscript/strong_lensing_main.tex`; do not substitute a successful build of the broader `main.tex`. The included workflow builds the focused manuscript. Its unexecuted status in the original package is historical; report the actual worktree and remote CI outcomes at the integration commit.

Run existing focused tests that cover the touched interfaces, resolve any concrete regression, and distinguish an unrelated pre-existing broad CI import failure from the focused result. Render and inspect every PDF page for clipping, equation width, references, tables, scientific figure labels and journal metadata. Do not add broad redundant tests simply to increase the count.

## Finish the integration

The author has confirmed reuse of his prior no-external-funding, no-conflicts
and sole-author contribution declarations. His identification is Independent
Researcher, Oconomowoc, Wisconsin, USA; ORCID `0009-0000-6650-3485`;
correspondence `carberry8878130@hotmail.com`. He confirmed that this manuscript
is not published or under consideration at another journal. These declarations
are written in the active manuscript and cover letter; do not reopen the same
questions. The cover letter remains an unsent draft, and no journal submission
or PR merge is part of this integration. Historical AI labels are attributed
to their records without claiming an unsupported model for the current turn.

After validation, commit the integrated correction to the existing PR branch and update draft PR #23 if repository access is available. Explain the BHSM surface-response question, the retained R1 calculation, the constructive conventional representation, the common GR2 correction, what passed, and what remains conditional. Do not claim the public baseline already contained the new code. Identify the resulting immutable commit in the reproducibility record; if an archive DOI is needed, leave that factual step explicit rather than inventing an identifier.

Report the concrete science commit, changed files, scientific results,
source/continuation checks, focused build and PDF evidence, and precise
remaining physical source/domain prescriptions. Do not merge the PR,
submit the manuscript, or send the cover letter.
