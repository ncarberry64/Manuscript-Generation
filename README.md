# Geometry Before Fields: Manuscript-Generation

Norman P. Carberry's conditional cosmology Article for **Universe (MDPI), Cosmology**.

The paper follows geometry -> action -> constraints/boundary conditions -> gauge quotient -> physical tangent space -> physical Hessian -> reduced response -> metric/matter -> optics -> observables. Field variables are response coordinates unless a stronger identification is independently established.

**Scope:** an explicit downstream response framework, not completed microscopic BHSM cosmology or a demonstrated solution of the Hubble tension. See the [submission checklist](docs/UNIVERSE_MDPI_SUBMISSION_CHECKLIST.md) and [hostile-referee audit](docs/HOSTILE_REFEREE_AUDIT.md).

## Scientific status

| Class | What the repository supports |
|---|---|
| DERIVED | Physical-domain Schur metric response; regular closed-FLRW harmonic exclusion; linear distance/optical identities; rank bound for a fixed finite latent-amplitude response. |
| CALIBRATED / CONDITIONAL | One selected physical n=2 profile, R1 parameters and reduced gravity-scalar coefficients, induced dust/growth/Weyl rows, SN amplitude/axis, coherent finite-channel stochastic realization. |
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
python -m pytest -q
python code/figures.py
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
| `artifacts/` | Retained results and validation evidence; not raw observations |
| `preregistration/` | Historical freezes and new structural two-gate protocol |
| `manuscript/` | LaTeX, figures and PDF; detailed audits in appendices |
| `docs/` | [Reader map](docs/README.md), [claim ledger](docs/CLAIM_STATUS.md), reproducibility and journal documents |
| `external/` | [Input provenance](external/README.md); raw surveys not bundled |
| Root `apply_*` / `populate_*` scripts | Historical editing provenance; not current build entry points |

## Submission checkpoint

Upgrade branch: `theory/cosmology-universe-manuscript-final`, starting at **5dff33b**, which repaired and preserved the unfinished Universe work. Earlier lineage: **14fced6**, **22b58bd**, **129a88c**, **55e2134**, **203cd43**.

The intended review checkpoint is final HEAD on this branch and the draft PR. Run `git rev-parse HEAD`; the [validation report](docs/VALIDATION_REPORT.md) identifies the tested scientific snapshot. A permanent archive/DOI and author-approved submission commit remain TODO. No automatic merge.

That manuscript pass was merged as PR #21 into `theory/cosmology-universe-mdpi-submission-pass`. The subsequent sightline update is documented separately in [SN sightline validation](docs/SN_SIGHTLINE_VALIDATION_V1.md); the earlier full-suite evidence remains attached to its original scientific snapshot.

Related preprint: [Preprints.org 202601.1427](https://www.preprints.org/manuscript/202601.1427), under its verified published title. The supplied later Topographic Dark Energy title relationship needs author confirmation.

## Final Universe revision

See [submission readiness](submission/SUBMISSION_READINESS.md), [claim provenance](docs/UNIVERSE_CLAIM_PROVENANCE.md) and [reviewer risk audit](docs/REVIEWER_RISK_AUDIT.md). The final revision preserves the frozen science and adds no parameter rescue.
