# Reader and manuscript map

The repository contains separate cosmology and strong-lensing manuscripts.
The existing cosmology material below remains the **Universe** paper. The
active lensing paper is for **Galaxies**, special issue **Strong Gravitational
Lens Modeling**, with no submission target date.

## BHSM strong lensing

Start with the [current scientific scope](STRONG_LENSING_SCOPE_20261020.md)
and the [focused manuscript](../manuscript/strong_lensing_main.pdf). The
[integration checkpoint](STRONG_LENSING_BHSM_INTEGRATION_20261020.md) records
the current reproducibility, focused PDF and remote outcomes. Its
organizing question is how BHSM surface contributions determine incoming-to-
outgoing light transfer through an envelope event. The order is preparation,
action/source/boundary problem, metric response, optical observables, and
conventional representation and checks.

The [R1 source derivation](strong_lensing_r1_screened_20261020/R1_LOCAL_COUPLING_DERIVATION_20261020.md)
gives the retained action, normalizations, metric constraint, regular screened
branch and approximation boundary. The
[research review](strong_lensing_r1_screened_20261020/GALAXIES_RESEARCHER_REVIEW.md)
records the finite-body images and delays, scoped rank results, and positive
effective-mass construction. That conventional construction reproduces the
same retained leading gravitational and optical response. The complete
normalized BHSM surface-to-effective-source map is the next derivation;
the leading equivalence does not assert that it is finished.

The subsequent [common-source application](strong_lensing_r1_screened_20261020/BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
evaluates the nonlinear flux and metric Gauss equation and proves a remaining
physical-map nonuniqueness through two positive common-domain lifts. The
[field continuation](strong_lensing_r1_screened_20261020/R1_REMAINING_FIELD_DERIVATION_20261020.md)
derives exact scalar/metric/fluid/environment/matching equations and checks
twelve component identities. It identifies the consumed EOS/support and outer
Cauchy data; these equations are not an evolved higher-order physical solution.

The [finite-body GR2 derivation](strong_lensing_r1_screened_20261020/FINITE_BODY_GR_2PM_DERIVATION_20261020.md)
and [composition review](strong_lensing_r1_screened_20261020/COMMON_GR_2PM_COMPOSITION_REVIEW_20261020.md)
identify the common conventional correction and its limited scope. The
[author handoff](strong_lensing_r1_screened_20261020/CODEX_GALAXIES_INTEGRATION_PROMPT.md)
and [unsent cover letter](strong_lensing_r1_screened_20261020/GALAXIES_COVER_LETTER_DRAFT.md)
carry the current scope and retain author-only declarations.

The earlier [first comparison](STRONG_LENSING_CHECKPOINT_20261020.md),
[constant-target family](STRONG_LENSING_LOCAL_FAMILY_20261020.md),
[local checkpoint](STRONG_LENSING_LOCAL_CHECKPOINT_20261020.md), and
[ADM/phase projection trail](STRONG_LENSING_ENVELOPE_DERIVATION_20261020.md)
remain historical scientific controls. Their numerical evidence is preserved.
In particular the constant-target phase null, rolling scalar and physical
metric coordinates are different objects. A scalar response and its equivalent
effective mass are counted once.

## Universe cosmology

Start with the [overview](../README.md), [current claims](CLAIM_STATUS.md), [two gates](BHSM_MAKE_OR_BREAK_PREDICTION_v1.md), and [skeptical audit](HOSTILE_REFEREE_AUDIT.md).

Main sequence: Introduction -> Methods -> claims -> closed-S3 parent -> Schur response -> conditional native n=2 transport -> propagator -> growth/Weyl -> distance -> optics -> coherent covariance -> representation audit -> R1 results -> negative results -> predictions/two gates -> Discussion -> one Conclusion.

Detailed Horndeski background, constraints, Bianchi dictionary, illustrative harmonic filter, conventional high-n projection, local-source candidates, constitutive/current projection contracts, numerical protocol and reference figures remain in the PDF as appendices.

| Object | Main code under `code/` | Artifact under `artifacts/` |
|---|---|---|
| R1 background/high-n growth | `background.py`, `growth.py` | Freeze in `preregistration/prediction_manifest.json` |
| Physical n=2 response | `r1_n2_reduced_propagator.py` | `R1_n2_reduced_propagator.json` |
| Matter/Weyl rows | `r1_n2_growth_lensing_transfer.py` | `R1_n2_growth_lensing_transfer.json` |
| Distance/calibration rank | `r1_n2_luminosity_distance_kernel.py` | `R1_n2_luminosity_distance_kernel.json` |
| Optical means/covariance | `bhsm_optical_transfer.py` and related modules | `BHSM_optical_*` contracts |
| Foreground failures | `r1_2mrs_*_transfer.py` | Three `R1_2MRS_*_v1.json` results |
| Prospective gates | `bhsm_make_or_break_prediction.py` | `BHSM_MAKE_OR_BREAK_PREDICTION_v1.json` (synthetic) |
| Retrospective SN sightlines | `analyses/sn_sightline_tomography_v1/` script snapshot | [Final report and compact results](../analyses/sn_sightline_tomography_v1/FINAL_REPORT.md); full local archive identified by hash |

Old numbered referee audits, accumulated derivation notes and unused `11_conclusions.tex` / `07a_geometry_first_implications.tex` record historical stages. They do not supersede main.tex or the current claim ledger. Generator `--tex` options also emit historical prose; regenerate numerical artifacts without overwriting editorial sections.

Operational documents: [reproduction](REPRODUCIBILITY.md), [validation](VALIDATION_REPORT.md), [checklist](UNIVERSE_MDPI_SUBMISSION_CHECKLIST.md), [cover letter](UNIVERSE_MDPI_COVER_LETTER_DRAFT.md), [external inputs](../external/README.md).
