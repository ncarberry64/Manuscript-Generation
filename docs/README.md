# Reader and manuscript map

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
