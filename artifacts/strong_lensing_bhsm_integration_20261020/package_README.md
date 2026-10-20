# Galaxies strong-lensing research correction

## What this package contains

This package contains a revised 13-page author-review manuscript for **Galaxies**, special issue **Strong Gravitational Lens Modeling**, together with the new calculations, source derivations, numerical receipts, research review, cover-letter draft and a Codex integration prompt.

The conventional and scalar-on models use the same physical light preparation, source, density, measured redshifts, observer clock and compatible background maps. The separate rolling scalar follows the previously specified pure cubic R1 effective action with its frozen normalization. The earlier constant-target null result remains intact.

The leading R1 extra relative delays are 25.95925815 s and 8.64901938 s. Adding the same finite-body ordinary-GR second post-Minkowskian (2PM) correction to both models gives 25.95932714 s and 8.64900924 s. The changes in these model differences are approximately +68.99 and −10.14 microseconds. The common correction includes the central light ray and a consistent potential, deflection and Hessian.

The manuscript retains the physical scope: a conditional pure-cubic, leading quasi-static scalar model; missing scalar-dependent higher-order terms and global matching; leading effective-halo equivalence; and no claimed independent observational detection. The common GR2 control checks a specific omitted conventional term. Its small drift is not a bound on all uncalculated physics.

The public source baseline is:

- Repository: https://github.com/ncarberry64/Manuscript-Generation
- Existing draft PR: https://github.com/ncarberry64/Manuscript-Generation/pull/23
- Reviewed source commit: `49a1dabbe148f8ba9da4ac079a8d93fe97c5de76`
- Existing branch: `codex/strong-lensing-first-comparison-20261020`

The new files in this package are not asserted to be part of that earlier public commit. The connected GitHub integration rejected the attempted write with HTTP 403, “Resource not accessible by integration.” No branch or PR change, merge or journal submission was made. The complete local corrections are supplied here for integration through the existing Codex worktree.

## Layout

- `overlay/` mirrors repository-relative paths for the new and revised files. It is an integration overlay, not a replacement repository.
- `Galaxies_Strong_Lensing_Revised.pdf` is the compiled article.
- `existing_text_changes.patch` records the intended changes to the three existing text files, against the reviewed baseline.
- `MANIFEST_SHA256.json` contains hashes of the delivered files other than the manifest itself.
- `VALIDATION_RECEIPT.json` records the actual numerical/build checks and their scope.
- The main review and exact derivations are under `overlay/docs/strong_lensing_r1_screened_20261020/`.
- Full numerical values are under `overlay/artifacts/strong_lensing_r1_screened_20261020/`.

The three existing text files intentionally changed by the overlay are `README.md`, `manuscript/strong_lensing_main.tex` and `docs/STRONG_LENSING_LOCAL_FAMILY_20261020.md`. The tracked focused PDF is also replaced by the revised PDF. The main manuscript now compiles the new Galaxies section. The old historical strong-lensing sections remain in the repository; the source document correction limits its zero-source statement to the material seam contribution. Check the actual target tree before applying any file, since concurrent work may have added one of the new paths.

## Integration into the existing worktree

Inspect current Git status and head, read repository instructions, and preserve all existing edits. Compare this overlay with the current tree. The patch is a review aid for existing text changes; it is not permission to replace a newer manuscript without reconciling it. The new files can be added at their matching repository-relative locations.

Use `overlay/docs/strong_lensing_r1_screened_20261020/CODEX_GALAXIES_INTEGRATION_PROMPT.md` for the complete integration task. It targets the existing draft PR and asks Codex to validate and commit the integrated result while leaving the PR draft and the cover letter unsent.

For a clean numerical reproduction, obtain the public repository at the pinned commit and apply this overlay. The existing repository supplies the frozen background and optical code and the previously generated first/local comparison inputs. Source hashes and URLs are recorded in the numerical results and source manifest; this package does not duplicate the full source repositories.

## Numerical reproduction

The recorded environment was Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. `requirements-strong-lensing-review.txt` records these versions. Other compatible versions may work; exact floating-point byte identity is not required across platforms.

Run from the integrated repository root:

```bash
python scripts/run_strong_lensing_r1_screened.py
python scripts/run_strong_lensing_gr2_control.py
python scripts/r1_screened_review/r1_qs_coefficients.py
python scripts/r1_screened_review/check_screened_projection.py
python scripts/r1_screened_review/multisource_crossratio.py
python scripts/r1_screened_review/shared_lens_rank.py
python scripts/r1_screened_review/check_gr_optical_control.py
python scripts/r1_screened_review/check_common_gr_composition.py
python scripts/plot_strong_lensing_r1_screened.py
```

The first four-case control uses the common model implementation. The composition review uses independently implemented axial roots, potential quadrature and a mixed stationary-action derivative. The geodesic check evaluates exact static GR scattering through the solved finite body and verifies the interior coefficient by a compactness sequence. Those independent numerical methods test the retained equations; they do not provide independent empirical evidence.

The rank calculations are statements about their declared nuisance model. Their covariance entries are numerical weighting choices, not observed precision or a signal-to-noise forecast. The effective-halo theorem is stated at the retained leading scalar order. Do not present a surviving finite-profile template as a unique theory discriminator.

## Build the article independently

The focused LaTeX source and figure are complete inside `overlay/`; no numerical rerun is required to build the article. From `overlay/` or from the integrated repository root:

```bash
python scripts/build_strong_lensing_manuscript.py
```

This requires Python, Ghostscript, latexmk and a complete TeX installation with the official MDPI class dependencies. On an Ubuntu installation, the included CI workflow installs:

```text
ghostscript latexmk
texlive-latex-base texlive-latex-recommended texlive-latex-extra
texlive-fonts-recommended texlive-fonts-extra
texlive-science texlive-pictures texlive-humanities
```

The build script converts the original official EPS logos to PDF and compiles `manuscript/strong_lensing_main.tex`. The official `Definitions/mdpi.cls` and original logo assets are unmodified. The source uses the `galaxies` journal option; document-level settings identify this as an author-review draft and omit publication-date fields.

The class was obtained from the official numeric-reference MDPI template:
https://res.mdpi.com/data/MDPI_template_ACS.zip?v=20260911

The new `.github/workflows/strong-lensing-manuscript.yml` builds this focused manuscript. It has been created and checked locally, but has not yet executed on GitHub. Prior CI for the broader `main.tex` should not be treated as evidence that this new target was built. The separate broad numerical CI was not rerun or claimed repaired by this package.

## Author review and submission facts

The manuscript is a concrete review draft. The listed author still needs to confirm affiliation/correspondence, CRediT, funding, conflicts, historical Codex version details and final approval. The cover letter deliberately leaves publication/exclusivity declarations for factual confirmation. The new additions also need an immutable public release or inclusion as the actual submitted supplement; replace the draft data-availability statement with its real identifier.

The journal instructions and special-issue page were checked from official MDPI sources. The paper uses a single-paragraph abstract within the approximately 200-word limit, six keywords, numerical references and explicit substantive AI-use disclosure. No submission target date is included.

- https://www.mdpi.com/journal/galaxies/instructions
- https://www.mdpi.com/journal/galaxies/special_issues/50J931CL3Q
- https://www.mdpi.com/authors/latex
