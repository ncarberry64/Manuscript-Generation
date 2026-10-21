# HE 0435-1223 observed-light experiment

Current observed-data checkpoint: [study report](../../docs/STRONG_LENSING_OBSERVED_20261020.md). This supersedes the earlier synthetic acceptance checkpoint as the active empirical application; all earlier equations, supplied numerical outputs and PDF audits remain historical controls.

The predictions were frozen and committed at [`e6e62bd`](https://github.com/ncarberry64/Manuscript-Generation/commit/e6e62bddd793a6b5545040b2fc37b775876384da) before opening targets. The primary retrospective test conditions on the later B curve and predicts later A/C/D; the secondary uses training only for all later images. Original public light curves and attribution are in [data/strong_lensing_observed_20261020](../../data/strong_lensing_observed_20261020/README.md).

Science release [`14ad777`](https://github.com/ncarberry64/Manuscript-Generation/tree/14ad777eaf8befda387b9c3cb4504cf499e76cfc) contains the observed code/data, inherited source ancestors, all recorded scientific applications and the figure. It is distinct from protocol commit `e55725d` and prediction-freeze commit `e6e62bd`. The real-MDPI build is **24 pages**, **427799 bytes**, with no LaTeX warnings and PDF SHA-256 `5a031e394a1b9626f9eec4e5ba0cec6c20f97cc5320419771c7f8136e78c43f6`. Final release review/preservation status is recorded by the [publication manifest](publication_manifest.json), with [24-page rendered contact sheet](pdf_contact_sheet.png).

The manuscript is **Strong-Lensing Time Delays in the Berger–Hopf Framework: An R1 Calculation and an Observed-Light Test**. The [known-outcome replay](../../scripts/replay_he0435_scores.py) was executed once against these already exposed scores; it reproduced saved values within declared tolerances and writes separately under ignored build output by default. It is numerical reproduction, not another held-out experiment. Frozen module/parameters and the original score are unchanged.

| File | Meaning |
|---|---|
| `lightcurve_training.json`, `training_start_*.json` | Historical-only empirical OU/delay optimization, all four starts and selected result |
| `training_curvature.json`, `training_diagnostics.json` | Training-only local curvature and derivative diagnostics; not a calibrated parameter posterior |
| `code_normalization_audit.json` | Exact byte audit of the post-fit immutable-input copy fix; no numerical equation changed |
| `conventional_training.json`, `r1_training.json` | Historical astrometry plus original training-light likelihood macro fits, original source/field applications and residuals |
| `independent_data_likelihood_review.json` | Pre-exposure data, likelihood and guarded freeze/score review |
| `independent_macro_review.json` | Pre-exposure source, two-plane optics, clock, astrometry and fitting review |
| `macro_forward_audit.json` | Fixed-body R1 activation, distinct separately refitted changes, original source residuals and five-image census; no target values used |
| `frozen_predictions.json` | Exact fitted models, source/data/protocol hashes and numerical versus guard-edit provenance, committed before scoring |
| `scoring_started.json` | Persistent exclusive-create exposure marker; preparation stays locked after any scoring interruption |
| `withheld_scores.json` | One frozen application: all target values, predictions, residuals, marginal calibration and exact joint conditional scores |
| [publication_manifest.json](publication_manifest.json) | Final science/PDF provenance, rendered-review status and preservation checks |
| [pdf_contact_sheet.png](pdf_contact_sheet.png) | All 24 rendered pages of the observed-light revision |

The primary conventional/R1 joint log densities are **791.623845 / 791.628828**, a difference of only **+0.00498269** for R1. Both have **85.33% nominal-95% marginal coverage** and approximately **0.10133 mag RMSE**, with image A coverage only **62%**. The train-only secondary has **75.50%** overall nominal-95% coverage and about **0.19495 mag RMSE**. The shared OU transfer is not adequately calibrated by these data; correlated targets prevent treating coverage counts as independent-binomial tests. This is no evidence of a distinctive BHSM response.

Both fitted models additionally predict a fifth central stationary image at approximately **6.3% of A's macro flux**. No observational non-detection limit was fitted. The search uses an analytic outer bound and dense multiseed interior sampling, without a completeness proof for interior roots. Matching source-derived positive effective mass gives the same retained R1 transfer and reuses its score by structural identity; it is not an independent fitted confirmation or an extra halo.

The likelihood is classical empirical OU variability with diagonal measurement-noise approximation and correlated latent brightness covariance. It does not supply a BHSM quantum state kernel. Fixed local R1 state/background plus a conditional remote optical screen is not a solved global envelope, and physical admission remains **false** with the unchanged **2e-10 guard**.

Targets are consumed. Do not repair the model and relabel rescoring of these targets as a fresh held-out test. Preserve this result; use a new protocol and genuinely unused observations for a later confirmatory test. To reproduce the published statistics, inspect saved receipts or use an explicitly labeled numerical reproduction, preserving the original score and exposure marker.
