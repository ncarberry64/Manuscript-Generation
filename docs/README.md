# Reader and manuscript map

The repository contains separate cosmology and strong-lensing manuscripts.
The existing cosmology material below remains the **Universe** paper. The
active lensing paper is for **Galaxies**, special issue **Strong Gravitational
Lens Modeling**, with no submission target date.

## BHSM strong lensing

Start with the [current observed-data checkpoint](STRONG_LENSING_OBSERVED_20261020.md), [observed artifact index](../artifacts/strong_lensing_observed_20261020/README.md), [data provenance](../data/strong_lensing_observed_20261020/README.md) and [focused manuscript](../manuscript/strong_lensing_main.pdf). The HE 0435-1223 pilot fits conventional and conditional R1 models to 902 historical epochs and tests 100 later epochs after a committed prediction freeze. The primary test conditions on all later B observations and withholds A/C/D; the secondary forecasts all four later images from historical training only.

The primary R1-minus-conventional score difference is **+0.00498**, with **85.33% nominal-95% marginal coverage** for both. The future-only test has **75.50%** coverage. The result provides no BHSM discrimination and exposes inadequate conditional predictive calibration. Both fitted models also predict a fifth central image, without an imposed detection limit. The frozen source/optical applications, original residuals and one score remain available; repairing and rescoring the same consumed targets cannot create a fresh held-out test.

The manuscript is **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Calculation and an Observed-Light Test**. Science checkpoint [`14ad777`](https://github.com/ncarberry64/Manuscript-Generation/tree/14ad777eaf8befda387b9c3cb4504cf499e76cfc) pins the observed application. The current real-MDPI build has 24 pages; [publication manifest](../artifacts/strong_lensing_observed_20261020/publication_manifest.json) and [contact sheet](../artifacts/strong_lensing_observed_20261020/pdf_contact_sheet.png) record the release/PDF review separately from the unchanged pre-exposure freeze `e6e62bd` and known-outcome numerical replay.

The argument remains specified state/source, normalized optical transfer, image/delay observables, then what those observables identify under conventional freedoms. [Scientific scope](STRONG_LENSING_SCOPE_20261020.md) and detailed derivations delimit the conditional physical assumptions. Physical admission remains false with the unchanged 2e-10 guard.

The earlier [acceptance checkpoint](STRONG_LENSING_ACCEPTANCE_20261020.md) and [artifact index](../artifacts/strong_lensing_acceptance_20261020/README.md) record the synthetic/control applications at science checkpoint [`87503ae`](https://github.com/ncarberry64/Manuscript-Generation/tree/87503ae9d7672405e07b573fd140c784b8ce6cbd). Their 19-page [publication manifest](../artifacts/strong_lensing_acceptance_20261020/publication_manifest.json) and rendered audit remain historical, distinct from both the older 18-page integration record and the new observed pilot.

The [R1 source derivation](strong_lensing_r1_screened_20261020/R1_LOCAL_COUPLING_DERIVATION_20261020.md)
gives the retained action, normalizations, metric constraint, regular screened
branch and approximation boundary. The
[research review](strong_lensing_r1_screened_20261020/GALAXIES_RESEARCHER_REVIEW.md)
records the finite-body images and delays, scoped rank results, and positive
effective-mass construction. That conventional construction reproduces the
same retained leading gravitational and optical response. The
[finite nonlinear re-fit](strong_lensing_r1_screened_20261020/NONLINEAR_LENS_REFIT_20261020.md)
executes conventional absorption for the declared axial example. The
[profile and coefficient controls](strong_lensing_r1_screened_20261020/LENS_MODEL_ROBUSTNESS_20261020.md)
apply a matched Plummer profile, one-input variations and state-coefficient
partials. The latter are not independently selected physical states: they
must be composed with an action-owned state-to-coefficient Jacobian. The
[literature and call assessment](strong_lensing_r1_screened_20261020/LITERATURE_AND_CALL_POSITIONING_20261020.md)
delimits the contribution and its journal fit.

The [distributed-encounter control](strong_lensing_r1_screened_20261020/ENCOUNTER_ACCURACY_CONTROL_20261020.md)
applies the same retained force and physical source between individually
derived closed-FLRW exterior wavefronts. Its common-window width/path result
is distinct from the common ordinary-GR2 control and from uncomputed scalar
higher-order or environmental effects.

The subsequent [common-source application](strong_lensing_r1_screened_20261020/BHSM_COMMON_SOURCE_DERIVATION_20261020.md)
evaluates the nonlinear flux and metric Gauss equation and proves a remaining
physical-map nonuniqueness through two positive common-domain lifts. The
[field continuation](strong_lensing_r1_screened_20261020/R1_REMAINING_FIELD_DERIVATION_20261020.md)
derives exact scalar/metric/fluid/environment/matching equations and checks
twelve component identities. It identifies the consumed EOS/support and outer
Cauchy data; these equations are not an evolved higher-order physical solution.
These detailed derivations and the long cap/source material are supplementary
to the main transfer/inference argument. The separate
[technical supplement source](../manuscript/sections/strong_lensing_technical_supplement_20261020.tex)
retains that material without adding it to the focused main document.

The [finite-body GR2 derivation](strong_lensing_r1_screened_20261020/FINITE_BODY_GR_2PM_DERIVATION_20261020.md)
and [composition review](strong_lensing_r1_screened_20261020/COMMON_GR_2PM_COMPOSITION_REVIEW_20261020.md)
identify the common conventional correction and its limited scope. The
[author handoff](strong_lensing_r1_screened_20261020/CODEX_GALAXIES_INTEGRATION_PROMPT.md)
and [unsent cover letter](strong_lensing_r1_screened_20261020/GALAXIES_COVER_LETTER_DRAFT.md)
carry the current scope and retain author-only declarations.

The earlier [package integration](STRONG_LENSING_BHSM_INTEGRATION_20261020.md),
[first comparison](STRONG_LENSING_CHECKPOINT_20261020.md),
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
