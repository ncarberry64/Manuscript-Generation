# Final scrutiny validation — 2026-09-24

Branch: `review/cosmology-final-scrutiny-20260924`. Public baseline: `083e14654d1b81b3853484b48e132bd5bc41567c`.
The coupled provenance and corrected integration were carried forward as `e0794fa` and `850b1ab`. The final review commit is obtained with `git rev-parse HEAD` on this branch; this report does not embed its own self-referential hash.

Title: Geometry Before Fields: Constraint-Reduced Cosmology and Optical Predictions in the Berger--Hopf Framework

The review corrects stale matter/radiation claims to the audited frozen-R1 linear/quadratic result, separates temporal rank from spatial selection, identifies the retained GP equations as historical compatibility audits, and confines harmonic-filter extrapolations to the illustrative realization. Reference kernels and every frozen scientific input remain unchanged. No failed identity was patched and no observational fit was run.

The coupled table and figure use `(zeta,zeta_dot)` and a hatted transfer. Rank two is reproduced in both output bases at eight post-anchor epochs. Singular-value magnitudes are basis dependent; the anchor block is zero. `Pi2` is reference normalized and is not generally the coupled canonical momentum. Environmental alignment alone does not align an independent initial topographic state.

| Gate | Result |
|---|---|
| Title alignment / stale matter-radiation claims | PASS / CLEARED |
| Frozen scientific/test/protocol/sightline files | 191 byte-identical to the full-suite worktree |
| Focused action-native, environment and gradient checks | 14 passed |
| Full 102-test suite | 101 passed, one documented strict expected failure |
| Main figures | Six, plus two appendix figures |
| PDF / abstract / bibliography | 54 pages / 244 whitespace-delimited words / 22 references |
| LaTeX / isolated ZIP rebuild | PASS; no warnings, unresolved references/citations or overfull boxes |
| Visual page audit | PASS; every page rendered and visually checked |
| Retuning | FALSE |
| Ready for final author confirmation | TRUE |

The full suite is the run launched in `r1-coupled-environment-integration`, not a second fresh run on this review branch. All 102 test identities and 191 recorded scientific inputs match exactly, including the three retained JSON fixtures read by the environment and gradient tests. Its evidence is imported only after completion and checked for unexpected failures/skips. The new plotter was run separately. The focused tests were rerun here; the gradient tests validate the committed audit artifacts, not a newly rerun principal-symbol integration. Both figure generators ran successfully; the four regenerated legacy PNGs are unchanged.

The environmental replay recovers both topographic seeds to 6.01005e-14, inverts the seed basis to 2.22045e-16, and reproduces the saved z=1.5 propagator with maximum difference 0. It freshly integrates the full six-state fundamental matrix with the existing action-native DOP853 implementation and applies the constraint/basis maps. This is a repeat using the same numerical owner, not an independent integrator cross-check or reconstruction of physical environmental initial data.

The original plotted singular values/determinants and the current replay are not bitwise identical: maximum absolute difference 7.8509e-09, maximum relative difference 1.10053e-07. A trial 1e-12 absolute identity check failed. Both artifacts are preserved without modifying their values, the figures use the declared original input, and both audits retain the same sampled rank-two conclusion. This is not reported as exact numerical reproduction of every archived scalar.

`FINAL_SCIENTIFIC_FIGURES_RECEIPT.json` hashes the analytic/artifact inputs and all new figures. `FINAL_SCRUTINY_VISUAL_AUDIT.json` records PDF and per-page render hashes. All pages were inspected in 100-dpi contact sheets, with individual detail views for revised geometry/text. Crowded analytic longitude labels were corrected and rechecked. The isolated rebuild has identical extracted text and page count. The source ZIP verifies every internal hash.

Existing Python 3.14 warnings about legacy non-raw TeX string escapes were observed during collection. No scientific owner was changed to suppress them. The LaTeX builds have zero warnings.

The first GitHub numerical job failed during collection because the workflow omitted `PYTHONPATH=code` (`ModuleNotFoundError: background`). The workflow now supplies the same import path as the validated local command and runs the new figure generator as well. No test, solver, tolerance or scientific parameter was changed. Hosted CI status is separate from the local full-suite receipt; inspect the PR checks before merging.

Open science remains explicit: independently reconstructed environmental anchor data, spatial profile/axis selection, microscopic and absolute local-source normalization, nonlinear-fluid/UV completion, joint observational covariance, exact radial BAO and survey-ready RSD. No observational signal or prospective gate success is claimed.

Author-only CRediT, funding and conflict declarations remain unfilled. No archive DOI, journal submission, external message or merge was performed. Review and merge approval remain with the author, as explicitly required by the final scrutiny handoff.

Reproduce: run `python code/figures.py`, `python code/final_scientific_figures.py`, the focused pytest files listed in the JUnit receipt, `latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex` in `manuscript`, and `python scripts/package_submission.py`. The local receipt recorder additionally requires `pypdf` (`python -m pip install pypdf`) and the recorded sibling worktree, page-review evidence and isolated rebuild. With the sibling suite complete, run `python scripts/record_final_scrutiny_validation.py` to refresh the checked completion receipt. Full-suite reproduction: set `PYTHONPATH=code` and run `python -m pytest -q` (long integrations).
