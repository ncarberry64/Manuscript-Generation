# Blackhole Tracker — central paper

The central manuscript is **Blackhole Tracker: Boundary-to-Lensing Transfer and Blind Directions in a Conditional BHSM Realization**. It supersedes the Hopf-charge paper as the priority of this branch, following Norman's explicit request. The two earlier paper directories remain preserved supporting research.

The checkable contribution is a normalized conditional boundary-to-metric-to-lensing transfer and its information losses. The adopted cubic R1 action fixes the response operator and metric coefficient. It does not select a microscopic BHSM identification, black-hole scalar flux or physical state. Deep screened quadrupoles can alter images and fluxes with a vanishing leading pair-delay contribution; their apparent-mass statistic is exactly shear-degenerate at first order. Joint data retain two transverse directions in a restricted synthetic information experiment. The tiny supplied physical state and unrealistic broad motion scan prevent a practical detection claim.

Read `RESULT_STATEMENT.md`, `CHAIN_AND_FREEDOMS.md`, `LAY_REPORT.md`, `VALIDATION.md`, the compact provenance map and the internal reviews before interpreting the curves.

## Reproduce

Keep the sibling `bhsm_angular_lensing_20261021` directory: its exact archived provider is the theoretical input, not an observational data set. From this directory, using Python with the requirements installed:

    python scripts/reproduce.py

The nine-stage replay verifies all pinned provider hashes and checks the source branch, radial equation, projection/derivative/normalization controls, estimator nulls, full primary-image responses, nuisance projection and figures. Logs and machine-readable receipts are under `results/` and `review/`. Random data are not needed. Source positions, matching fluxes and the state are supplied synthetic inputs.

To compile, from `manuscript/`:

    latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

No terminal TeX installation is needed to read the compiled PDF. The build uses standard LaTeX article packages; this is an author-review research draft for Galaxies, not a submitted or accepted journal article.

## Key numerical interpretation

The final experiment has 54 dimensionless observable components and 30 nuisance coordinates, including the open scalar flux and axial scalar quadrupole. Equal dimensionless errors are a design metric. The actual small-state example has optical amplitude about −8.89e−13. Large optical amplitudes used to test lens derivatives do not establish scalar admissibility at those enlarged states. The 70-digit projection refines arithmetic of fixed double-input response columns, not physical accuracy. Full line integrals approximate an isolated frozen local exterior; tail bounds cover tails only.

No observations, horizon evolution or empirical forecasts are included. Journal submission and external contact have not been performed.
