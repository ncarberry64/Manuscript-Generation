# Blackhole Tracker: research paper and usable product

The current manuscript is **Blackhole Tracker: Calibrated Compact-Lens Inference and Admissible Boundary Responses in a Conditional BHSM Realization**. It accompanies a runnable specialist product for inference from resolved image pairs, flux ratios and delays. The earlier boundary-transfer and broad information experiment remain in the paper as distinct theoretical results; the preserved angular and Hopf papers remain supporting research.

The package establishes three connected results:

- A local GUI, CLI, Python API and two executed notebooks fit compact-lens mass, center and observer proper motion with the full supplied covariance and source positions free at every row. Negative or unresolved results are explicit outputs.
- A uniform sufficient bound controls all five normalized quadrupole coordinates throughout the retained R1 exterior, preserving its scalar principal signs within declared budgets. The exact Dirichlet/conormal pairing gives the corresponding source bound. This is a conditional local-action result, not full-system or horizon stability.
- The normalized boundary-to-metric-to-lensing transfer has explicit information losses: two line-of-sight modes vanish; the leading screened primary-pair delay response vanishes; and its first-order apparent-mass response is degenerate with external shear. The old joint-response experiment retains information discarded by the scalar statistic within its declared finite nuisance family.

## Use the product

Read `product/README.md` in the delivery bundle for the Windows double-click launcher, installation, local app, CSV/JSON import, exact covariance order, CLI fitting/profiles, Python API and notebooks. The canonical source-checkout package is the repository's top-level `product/`; the delivery includes a self-contained copy. It needs no private research archive, remote fitting service or API key.

The ordinary compact-lens mode is useful within its validated synthetic domain even when no BHSM-specific envelope can be measured. It can identify a supported mass/motion track, reject an inadequate model or report unresolved motion. It does not establish that a compact dark source has a horizon. All bundled data are synthetic, including the published-parameter stellar variation.

The frozen 64-draw full-covariance calibration completes without excluded fit failures. Nominal local 95% coverage is 61/64 for mass and both motion components, 59/64 for center x and 61/64 for center y. Exact binomial tests at the fixed alpha=0.01 reject none; finite-sample uncertainty and bias remain reported. The short cosmological control stays motion unresolved, and omitting conventional shear rejects the fit. The optional R1 envelope stays far below the nearby example's supplied noise at its admitted state cap.

## Reproduce the separate calculations

From the self-contained source product directory, install its development dependencies and run:

```text
python -m pip install -e ".[dev]"
python scripts/validate_product.py
python scripts/execute_notebooks.py
```

Product stage logs, hashes and machine receipts are under `product/validation/`; the paper's `review/` contains the protocol and delivery snapshots. The notebooks replay the workflow rather than replace calibration.

The original theory replay remains independent. Keep the sibling `bhsm_angular_lensing_20261021` directory containing its exact archived provider. From this paper directory, with `requirements.txt` installed, run:

```text
python scripts/reproduce.py
```

Its nine stages check pinned inputs, the source/radial/projection equations, estimator nulls, nonlinear image derivatives, nuisance projection and figures. Logs are under `results/` and `review/`. The old 54-component, 30-nuisance scan uses equal dimensionless errors, not a telescope covariance; its broad trajectory would require millions of years for the cosmological example. It is not the product's ten-year calibration fixture.

From `manuscript/`, compile with:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Read `RESULT_STATEMENT.md`, `LAY_REPORT.md`, `VALIDATION.md`, `CHAIN_AND_FREEDOMS.md`, the provenance map and internal reviews for scope and evidence. This is an author-review draft for *Galaxies*, not a submitted or accepted article. No journal submission or external contact has occurred.
