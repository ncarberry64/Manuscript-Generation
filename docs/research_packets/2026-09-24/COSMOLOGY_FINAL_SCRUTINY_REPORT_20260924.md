# Final scientific scrutiny — Geometry Before Fields
Date: 2026-09-24

Public repository audited:
`ncarberry64/Manuscript-Generation`

Public baseline:
`main` at `083e14654d1b81b3853484b48e132bd5bc41567c`
(`Merge submission-ready geometry-to-prediction cosmology manuscript`)

This review treats the public source tree, not an older local PDF, as the
authoritative manuscript baseline.

## Executive assessment

The paper has a strong publishable scientific spine:

`S^3(R_H) -> physical n=2 sector -> canonical state -> constrained
metric/matter response -> luminosity-distance / transverse-BAO / growth /
Weyl / optical kernels -> prospective falsification`.

The public manuscript is already careful about conditional versus derived
claims and is stronger than the January preprint.  The remaining work is
primarily **consistency, presentation, and promotion of results that are
already present in the repository but underrepresented in the main text**.

The main scientific corrections should be made before submission.

## Priority 0 — corrections that should be made

### 1. Restore one canonical title everywhere

The public `manuscript/main.tex` currently uses:

`Constraint-Reduced Cosmology on a Closed Hypersphere: Berger--Hopf Optical
Transfer and Falsifiable Cross-Observable Predictions`

while the repository README and Universe cover letter use the stronger project
identity `Geometry Before Fields`.

Use one title everywhere:

**Geometry Before Fields: Constraint-Reduced Cosmology and Optical Predictions
in the Berger--Hopf Framework**

This is scientifically accurate, more memorable, and already matches the
cover letter.

### 2. The abstract is scientifically stale

The current abstract says that "complete matter--radiation stability" remains
open.  That is no longer the correct claim boundary after the committed
action-native ADM + Schutz--Sorkin matter/radiation work and the committed
principal gradient/dispersion audit.

The main text should instead say:

- the full frozen-R1 dust+radiation quadratic system has been action-reduced;
- the audited reduced kinetic matrix remains positive;
- no ghost or high-k scalar gradient instability is detected over the audited
  linear/quadratic R1 domain;
- nonlinear-fluid, ultraviolet and microscopic BHSM completion remain open.

Do **not** upgrade this to nonlinear stability or microscopic completion.

### 3. Promote the coupled environmental-state result into the cosmology paper

The public research provenance commit
`e32b13ad86d1b66889a4a0483fed9ef20ae84561`
establishes, within the unchanged frozen R1 six-state linear system,

\[
X(z)=U_{XX}(z)X_i+U_{XE}(z)E_i ,
\]

with `rank(U_XE)=2` at every sampled post-anchor epoch and a nonzero
matter-only `(delta_m,v_m)` determinant throughout the sampled interval.

The correct manuscript interpretation is:

- an isolated one-dimensional temporal attractor is **not required**;
- specified environmental density/velocity data can span both temporal
  topographic components;
- temporal realization and spatial profile/axis selection are separate;
- a single spatial profile is the special condition
  `rank_spatial(X)<=1`;
- the z=2.1 environmental initial state, full n=2 spatial tensors and the
  microscopic BHSM->R1 normalization remain open.

This is a real theoretical advance and belongs in the main paper.

### 4. Correct stale discussion language

`manuscript/sections/10_discussion.tex` still says:

> This is not a fully backreacting matter+radiation stability theorem.

and later lists full matter-radiation backreaction as an open scientific
input.

Replace that language with the current boundary:

- linear/quadratic action-native dust+radiation reduction and principal-symbol
  audit are closed at their stated level;
- environmental initial-state specification, early-time extension,
  spatial-profile reduction, nonlinear validity and UV completion remain open.

### 5. Reclassify the older matter-sourced appendix as historical compatibility work

`manuscript/sections/04b_matter_sourced_n2.tex` predates the action-native
six-state reduction.  It should remain for provenance, but its opening must
state that it is a historical curvature-aware compatibility derivation and is
not the current owner of R1 matter/radiation evolution.

The same clarification should be applied to the neighboring curved effective
dictionary section where necessary.

### 6. Scope the harmonic-filter claims

The current `Predictions and Falsification Criteria` section says that "the
bridge predicts selective n=2 softening" and states a large-n `C_n` law.
Those statements should be explicitly tied to the **retained illustrative
harmonic-filter realization**, not to the complete microscopic BHSM theory.
The physical n=2 harmonic result is stronger than the illustrative filter and
should not be conflated with it.

## Priority 1 — main-text scientific presentation

The paper should have a small number of information-dense, reproducible
scientific figures.  Decorative cosmology art, stock imagery, artistic black
holes, or pseudo-3D hypersphere renders should not be used.

Recommended six main publication figures:

1. **Geometry -> prediction chain** — retain the existing
   `constraint_prediction_chain`.
2. **Exact physical n=2 observer-shell geometry** — new composite:
   - exact monopole/dipole/quadrupole radial factors;
   - pure A4i sky dipole basis.
   The plotted axis is rotated to map center for visualization and is not an
   empirical axis.
3. **Coupled environmental realization** — new composite:
   - singular values of `U_XE`;
   - matter-only `det U_XE^(delta_m,v_m)`.
   This is the cleanest visual evidence for temporal rank two.
4. **R1 reference evolution** — promote the existing `r1_omegaT` and
   `r1_growth_ratio` from the appendix into one main-text composite figure.
5. **Principal-symbol stability** — new composite:
   - high-k scalar/radiation characteristic speeds versus redshift;
   - scaled positive-growth tail versus harmonic index at z=0 and z=2.1.
6. **Frozen SN comparison** — retain the existing observational summary figure.

Keep `r1_alphaB` and the illustrative `harmonic_filter` in the appendix.

A canonical-state covariance/calibration plot has also been prototyped during
this review, but it should not enter the final paper until its exact
high-z-state artifact is included in the public repository with the same
provenance standard as the other figures.

## Priority 2 — public audit/document consistency

Several public audit/readiness documents still repeat the now-stale statement
that full matter/radiation backreaction or stability is open.  Update the
current-facing documents:

- `docs/REVIEWER_RISK_AUDIT.md`
- `docs/HOSTILE_REFEREE_AUDIT.md`
- `docs/UNIVERSE_CLAIM_PROVENANCE.md`
- `submission/SUBMISSION_READINESS.md`
- `README.md` where relevant

Do not rewrite historical validation reports merely to make history look
cleaner.  Historical checkpoints should remain historical and be labeled as
such.

## References checked during scrutiny

The following high-risk records were checked against public sources and do
not require correction in this pass:

- Preprints.org 202601.1427 title and identifier;
- DESI DR2 Results II / Phys. Rev. D record and DOI;
- Specogna et al. 2026 closed-universe perturbation paper / Phys. Rev. D record
  and DOI;
- Gambino & Pace arXiv:2412.01781.

The bibliography should still receive the normal journal production pass.

## Submission blockers that require the author

The manuscript back matter still contains explicit placeholders for:

- CRediT / Author Contributions;
- Funding;
- Conflicts of Interest.

These must be confirmed by Norman P. Carberry before actual submission.
They should not be invented by an automated finalization pass.

A permanent public data/archive DOI is also still absent.  The GitHub
repository is useful and public, but a journal-quality archival deposit is
preferable if feasible.

## Final acceptance gates

Before merging the scrutiny branch:

1. Bring provenance commit `e32b13a...` onto the review branch without
   changing its scientific content.
2. Generate all new graphics from committed analytic formulas/artifacts.
3. Run focused cosmology tests.
4. Run the full repository test suite.
5. Build LaTeX from a clean worktree.
6. Require no undefined citations/references.
7. Require no overfull boxes.
8. Visually inspect every page of the rebuilt PDF.
9. Re-run the source ZIP/package builder.
10. Update submission-readiness page/figure counts only **after** the final
    build.
11. Confirm `RETUNING = FALSE`.
12. Do not merge automatically; present the final diff/PDF to the author.

## Bottom line

The paper does not need a new speculative mechanism before submission.
It needs the public manuscript brought into consistency with results the
repository already contains, especially the action-native matter/radiation
stability audit and the coupled environmental realization.

The strongest scientific story is:

**geometry selects the physical compact-space sector; the constrained R1
action transports a coupled state; matter and light read the same response;
and the resulting observable realization has explicit failure conditions.**
