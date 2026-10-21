# Geometry Before Fields: Constraint-Reduced Cosmology and Optical Predictions in the Berger--Hopf Framework

Norman P. Carberry's conditional cosmology Article for **Universe (MDPI), Cosmology**.

The paper follows geometry -> action -> constraints/boundary conditions -> gauge quotient -> physical tangent space -> physical Hessian -> reduced response -> metric/matter -> optics -> observables. Field variables are response coordinates unless a stronger identification is independently established.

**Scope:** an explicit downstream response framework, not completed microscopic BHSM cosmology or a demonstrated solution of the Hubble tension. See the [submission checklist](docs/UNIVERSE_MDPI_SUBMISSION_CHECKLIST.md) and [hostile-referee audit](docs/HOSTILE_REFEREE_AUDIT.md).

## Scientific status

The separate **BHSM strong-lensing paper**, [focused manuscript](manuscript/strong_lensing_main.pdf), is prepared for **Galaxies**, special issue **Strong Gravitational Lens Modeling**. Its current empirical application is the [HE 0435-1223 observed-light baseline and withheld comparison](docs/STRONG_LENSING_OBSERVED_20261020.md): 902 historical COSMOGRAIL epochs train conventional and conditional R1 models; 100 later epochs test their frozen predictions. The primary test explicitly uses later image B as an observed driver to predict withheld A/C/D, so it is retrospective batch cross-image prediction.

The conventional/R1 primary joint-score difference is only **+0.00498** for R1, with **85.33% coverage of nominal 95% intervals** for both. This does not establish a BHSM preference, adequate predictive calibration or a unique envelope state. A source/optical forward audit also finds a fifth central image at about 6.3% of A's macro flux, without a fitted non-detection limit. The [current artifact index](artifacts/strong_lensing_observed_20261020/README.md) preserves the training, pre-exposure committed freeze, original residuals, reviews and one held-out application. These targets are consumed and cannot be retuned into a new untouched test.

The current manuscript is *Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Calculation and an Observed-Light Test*. Science checkpoint [`14ad777`](https://github.com/ncarberry64/Manuscript-Generation/tree/14ad777eaf8befda387b9c3cb4504cf499e76cfc) pins the observed code/data, source ancestors, scores and figure. The [publication manifest](artifacts/strong_lensing_observed_20261020/publication_manifest.json) distinguishes the 24-page real-MDPI build and rendered review from the historical 19-page control audit. The original pre-exposure freeze `e6e62bd` remains unchanged; the separate known-outcome replay is reproduction, not a second test.

The study connects the declared state/source through normalized metric and optical transfer to images and delays, then asks what observations identify under conventional lens freedoms. Source-derived positive effective mass reproduces the retained R1 response; it is counted once. The earlier [acceptance checkpoint](docs/STRONG_LENSING_ACCEPTANCE_20261020.md) records synthetic finite-refit, profile/coefficient and distributed-encounter controls. It remains historical evidence alongside the [scope record](docs/STRONG_LENSING_SCOPE_20261020.md). This observed pilot is conditional on a prescribed local source, remote optical screen and background, with physical admission false; it is not a selected microscopic BHSM state or complete galaxy reconstruction.

The inherited matched SIS transport and constant-target phase null remain
scientific controls. The screened example has radius `r_V = 212.0292 kpc` and
delay offsets `25.95925815 s` and `8.64901938 s` relative to its specified
fixed-body reference; fully matched conventional and scalar representations
share the demonstrated leading response. Reproduce the effective calculation
with `python scripts/run_strong_lensing_r1_screened.py` and its common finite-body
GR2 control with `python scripts/run_strong_lensing_gr2_control.py`. See the
[research review](docs/strong_lensing_r1_screened_20261020/GALAXIES_RESEARCHER_REVIEW.md)
and [historical acceptance records](artifacts/strong_lensing_acceptance_20261020/README.md)
for the derivations, limitations and focused
validation. Synthetic-control science checkpoint
[`87503ae`](https://github.com/ncarberry64/Manuscript-Generation/tree/87503ae9d7672405e07b573fd140c784b8ce6cbd)
pins the finite re-fit, profile/coefficient and distributed-encounter
applications and eighteen focused tests. This separate lensing work
preserves the Universe manuscript.

The [source application](docs/strong_lensing_r1_screened_20261020/BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
and [field continuation](docs/strong_lensing_r1_screened_20261020/R1_REMAINING_FIELD_DERIVATION_20261020.md)
state the original equations, computed flux/constraint residuals and precise
remaining physical data: an action-owned embedding, material EOS/support,
environmental state and outer Cauchy/matching prescription. No extra optical
amplitude or selected microscopic lift is introduced. These detailed source
and continuation records support the focused argument as supplementary
derivations; they are not a new evolved physical solution. The earlier
[integration checkpoint](docs/STRONG_LENSING_BHSM_INTEGRATION_20261020.md)
retains its original science and PDF evidence.

| Class | What the repository supports |
|---|---|
| DERIVED | Physical-domain Schur metric response; regular closed-FLRW harmonic exclusion; linear distance/optical identities; rank bound for a fixed finite latent-amplitude response. |
| CALIBRATED / CONDITIONAL | One selected physical n=2 profile, R1 parameters and reduced gravity-scalar coefficients, induced dust/growth/Weyl rows, SN amplitude/axis, coherent finite-channel stochastic realization. |
| NUMERICALLY AUDITED | Frozen R1 action-native quadratic dust/radiation backreaction; positive audited kinetic matrix and no detected high-k scalar gradient instability. Environmental transfer has temporal rank two at eight post-anchor epochs; no spatial profile/axis selection follows. |
| REJECTED / NOT ESTABLISHED | Reported homogeneous localization failure; frozen number-density, K-band luminosity and seam-charge foreground transfer failures. The original localization statistic is not archived; foreground numerical artifacts are retained. |
| OPEN | Normalized compact full-preimage background, local action matrix elements, seam export, source/rare-event statistics, absolute amplitudes, profile selection and survey-ready nuisance calibration. |

The full n=2 angular eigenspace is not a two-component state. The temporal pair applies after selecting one profile and axis. One SN amplitude leaves a state direction undetermined. The curved-EFT DAE's registered failure concerns that representation, not EFT/QFT in general.

## Frozen supernova sightline update

The [V1 report](analyses/sn_sightline_tomography_v1/FINAL_REPORT.md) adds epoch/band reconstruction, independent galaxy/PV information, full-covariance held-out tests and experimental foreground-cell tomography. With every scientific prediction frozen, the primary Pantheon held-out comparison gives Δχ² = −0.883, while the DES transfer gives +1.321 (negative favors the frozen template). A coherent residual is **not established**. The independent fixed-PV sensitivity and incomplete source/field error propagation are explicitly separated from the primary result.

The [Git snapshot guide](analyses/sn_sightline_tomography_v1/GIT_SNAPSHOT.md) identifies the complete local reproducible archive, its SHA-256, and the compact files retained here. Raw inputs and large intermediate matrices remain in that archive. This retrospective sightline test does not replace the prospective two-gate protocol or erase the historical negative results.

## Make-or-break prediction

**Gate A:** deterministic targets admit one common two-component state, `d=F X_2+epsilon`; reject that realization if the properly calibrated GLS orthogonal residual has **p<0.01**.

**Gate B:** the minimal coherent one-amplitude contribution has **rank <=1** after independently specified standard/noise treatment; reject that realization at **p<0.01** using a predeclared finite-sample procedure. One local stochastic field need not give rank one after path integration. No clipping, extra channels, per-observable amplitudes or post-result bin/kernel changes can rescue the test.

The gates have separate nulls and are not assumed statistically independent. The [structural protocol](preregistration/bhsm_geometry_first_make_or_break_v1.json) is frozen; [survey-design prerequisites](docs/BHSM_MAKE_OR_BREAK_PREDICTION_v1.md) are not. The Gaussian benchmark is synthetic and restricted to independently fixed covariance/loading.

## Reproduce

```powershell
python -m pip install -r requirements-dev.txt
$env:PYTHONPATH='code'
python -m pytest -q
python code/figures.py
python code/final_scientific_figures.py
cd manuscript
latexmk -C
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The suite includes long nested integrations and one documented strict expected failure. See [REPRODUCIBILITY](docs/REPRODUCIBILITY.md) for focused commands, tested versions, packaging and external-data limits. The PDF is [manuscript/main.pdf](manuscript/main.pdf).

## Repository map

| Path | Purpose |
|---|---|
| `code/` | Canonical numerical, transfer and statistical implementations |
| `tests/` | Numerical, algebraic, frozen-manifest and synthetic checks |
| `artifacts/` | Retained results and validation evidence; original public lens observations are in `data/` |
| `data/strong_lensing_observed_20261020/` | Original COSMOGRAIL bytes, provenance and published point astrometry for the observed-lens pilot |
| `preregistration/` | Historical freezes and new structural two-gate protocol |
| `manuscript/` | LaTeX, figures and PDF; detailed audits in appendices |
| `docs/` | [Reader map](docs/README.md), [claim ledger](docs/CLAIM_STATUS.md), reproducibility and journal documents |
| `external/` | [Input provenance](external/README.md); raw surveys not bundled |
| Root `apply_*` / `populate_*` scripts | Historical editing provenance; not current build entry points |

## Historical submission checkpoint

Upgrade branch: `theory/cosmology-universe-manuscript-final`, starting at **5dff33b**, which repaired and preserved the unfinished Universe work. Earlier lineage: **14fced6**, **22b58bd**, **129a88c**, **55e2134**, **203cd43**.

The intended review checkpoint is final HEAD on this branch and the draft PR. Run `git rev-parse HEAD`; the [validation report](docs/VALIDATION_REPORT.md) identifies the tested scientific snapshot. A permanent archive/DOI and author-approved submission commit remain TODO. No automatic merge.

That manuscript pass was merged as PR #21 into `theory/cosmology-universe-mdpi-submission-pass`. The subsequent sightline update is documented separately in [SN sightline validation](docs/SN_SIGHTLINE_VALIDATION_V1.md); the earlier full-suite evidence remains attached to its original scientific snapshot.

Related preprint: [Preprints.org 202601.1427](https://www.preprints.org/manuscript/202601.1427), under its verified published title. The supplied later Topographic Dark Energy title relationship needs author confirmation.

## Final Universe revision

See [submission readiness](submission/SUBMISSION_READINESS.md), [claim provenance](docs/UNIVERSE_CLAIM_PROVENANCE.md) and [reviewer risk audit](docs/REVIEWER_RISK_AUDIT.md). The final revision preserves the frozen science and adds no parameter rescue.

Current review branch: `review/cosmology-final-scrutiny-20260924`, based on public `083e14654d1b81b3853484b48e132bd5bc41567c`. See [final scrutiny validation](docs/FINAL_SCRUTINY_VALIDATION_20260924.md). No merge or submission without author approval.
