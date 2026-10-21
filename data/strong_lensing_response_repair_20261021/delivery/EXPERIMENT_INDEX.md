# Executed experiment index

All development calculations use the already consumed HE 0435 F115W data or
instrument calibration data. An optimizer success flag means numerical
convergence, not response adequacy. The absolute-error and physical-validation
flags remain false unless explicitly stated otherwise by a qualified receipt.

| Calculation | Data used for its parameters | Executed result | Integration status |
|---|---|---|---|
| Rate-mask adapter | Pipeline version, DQ and valid-rate provenance | 303 usable slopes restored; canonical G1 RMS 0.448457 mas | Apply the tested production patch; calibration remains open |
| Positive common field response | Four isolated-star exposures | Held-star reduced Q 5.73–13.09; transferred target reduced Q 6.224–18.314; raw G1 RMS 0.362383 mas | Partial mean-response improvement; development only |
| Along-track spatial variants | Isolated star only | Held Q improves only 0.40–1.16% | Insufficient; transverse extrapolation is unsupported |
| Finer field response | Three star exposures per fold | Held Q 24319.41; 1249 coefficients, 849 profiled data directions | Overfits; do not transfer |
| Finer local quasar response | Twelve consumed target stamps per fold | Profiled rank 1246/1249; held Q 39064.54 / 4475 nominal conditional dof; all fits fail | More calibration samples solve rank defect but not response; do not promote |
| Local-quasar geometry after calibration | Fixed completed calibration receipt, then fixed physical predictions | All-image RMS 0.47543 mas; cross-exposure RMS 0.42804 mas; exposure statistics 452.85 and 647.88 | Diagnostic only; lens scores never trained the response |
| Dated STPSF with equal 5× sampling | Measured WSS OPD; prescribed 4000 K / 6.2 mas diffusion | Field reduced Q 34.86–69.93 | Inadequate baseline |
| Shared stellar color and diffusion | Three star exposures per held fold | Fitted 4198 K / 8.05 mas; held reduced Q 32.13–66.32 | Tested physical family is insufficient |
| Eight bounded pupil-phase modes | Three star exposures per held fold | Rank 8, condition 3.12; exact nonlinear held reduced Q 20.50–35.23 | Candidate insufficient; not an exhaustive wavefront search |
| Limited raw-read diagnostic | First consumed F115W exposure | Wing count-rate trend under approximate reference correction | Historical; superseded for native-processing comparisons |
| Recorded detector/ramp replay | Exact recorded references, same consumed exposure | stcal 1.18 field rate/ERR differences 8.97e-8 / 1.65e-7; declared numerical gate passes | Reconciles the earlier failed 1.19.1 replay; archive version is inferred |
| Common-time first-two-read extraction | One consumed field star, 289 matched pixels | Normalized residual 7.6653% → 5.1787%; early Q/dof 1.20644, p=0.00984653 | Fails 0.01 adequacy; poorer precision, mixed wing result; development only |

Primary production receipts are in `mask_patch/artifacts/`.

Definitive empirical-response receipts are:

- `audit/positive_field_results.json`
- `audit/positive_target_transfer.json`
- `audit/spatial_field_results.json`
- `audit/fine_field_results.json`
- `audit/local_quasar_refined_results.json`
- `audit/local_quasar_geometry_results.json`

Definitive physical-response receipts are:

- `stpsf_field_calibration/results.json`
- `stpsf_field_calibration/phase_trial/results.json`

The initial unequal-sampling STPSF trial is superseded. The corrected physical
trials use equal FFT and detector oversampling and verify their accelerated
response against direct optical calculations. The phase candidate is evaluated
with exact nonlinear optics after a finite-difference proposal.

The empirical PSF covariance figures are conditional on the fitted response.
They are not a completed uncertainty propagation for shape calibration,
spatial extrapolation or shared instrument effects. Training residuals and
conditional nominal degrees of freedom are not independent validation.

Definitive detector-history receipt:

- `ramp_diagnostic/native_pipeline/common_time_field_stcal118_results.json`

The first `common_time_field_results.json` remains the archived unsuccessful
1.19.1 numerical replay. See the adjacent `NATIVE_READ_HISTORY_FINDINGS.md` and compact
corrected field-group arrays. The early/full images share reads, so this is
not an independent physical validation or a paired significance calculation.
