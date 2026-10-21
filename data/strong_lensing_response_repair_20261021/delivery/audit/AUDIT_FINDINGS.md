# HE 0435 F115W response audit

Source: Manuscript-Generation PR23, head `7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9`.

## Direct code findings

`native_pixel_mask` rejects `DO_NOT_USE | SATURATED` identically. This is the
priority implementation issue because SATURATED in the producing pipeline can
record saturation of any later ramp group, while the calibrated rate may still
contain a valid early-group slope. Root's independent pipeline-version and
native-SCI audit must decide eligibility; a valid rate must not be discarded
solely because later groups saturated. The old statement "the target cores have
no usable information" does not follow from the saved fit. The saved fit merely
discarded them. In exposure1 the old mask excludes 15,59,16,12 SAT pixels from
A,B,C,D respectively. The variable exclusion strongly changes the sampled PSF
wings and can magnify shape-dependent centroid bias.

The original sampler otherwise agrees with its declared convention: spatial
STDPSF, oversampling4, origin(50,50), pixel-integrated response sampled once,
correct x/y transpose, free flux normalization. The archived implementation
audit reproduces the official Photutils evaluator at floating precision and
the original chi-square. The code's model is a fixed PSF plus quadratic local
background; its covariance is conditional on an exact fixed PSF, and has no
instrument-response uncertainty or between-exposure covariance.

An image-dependent host-light error cannot by itself explain the separately
observed field-star failure. A distortion error likewise does not explain the
poor detector-space point-profile fit. Those effects may still be additional
limitations, but the instrumental response/mask must be repaired first.

## Shared mean-response experiment

`field_response_transfer.py` is a diagnostic, not production code. It learns one
smooth additive two-dimensional correction from the four archived field-star
wing datasets; physical lens predictions never enter fitting or selection.
Four leave-one-exposure-out checks are used to compare a finite nine-item
grid of basis spacing and regularization. The correction has exact zero added
integral and first moments on the saved101x101 PSF grid. It is then transferred
unchanged to all16 consumed quasar stamps, refitting only each point's flux,
centroid and existing diffuse background.

- Baseline field chi-square sum: 52,042.674.
- Best four-fold field held-exposure chi-square sum: 12,008.249.
- Held-exposure chi-square/dof: 10.906,15.095,6.268,13.080.
- Original16-target chi-square sum: 188,060.592.
- Shared field correction transferred to16 targets: 45,396.587.
- Corrected target chi-square/dof range: 6.201–18.191.
- Largest target centroid shift:0.08947nativepixel (~2.8mas at0.031arcsec/pixel).

This supports a coherent instrument-response component that transfers from an
isolated star to quasar wings. It does not identify its physical cause and does
not finish calibration. The 541-parameter exploratory additive correction has
negative low-level wings and must NOT be promoted to a physical PSF or a
validated measurement. The basis/regularization was selected on the same
four-fold development comparison, not an independent external test. Original
SAT masking was retained in this experiment because the archived field
residual arrays omit those pixels. Root's corrected rate-validity replay is
therefore the more direct and better-motivated repair.

Outputs are `field_response_results.json` and `field_response_correction.npy`.
All original inputs are unchanged. No remote write was performed. No unused
target band was opened.

## Positive response with valid early-group cores

`positive_field_response.py` reruns the field calibration using the complete
native CAL stamps and the corrected DNU-only validity mask. It restores18,16,
18,21 field-star pixels in the four exposures. It models one positive smooth
two-dimensional multiplicative response common to all exposures. A constant
pedestal makes the original signed library grid positive; the original library
itself has about4,950 negative samples and 2.3–2.4% integrated negative flux.
At zero correction the pedestal is exactly absorbed by the free point flux and
background, so it is an explicit background convention, not extra source data.
The shared log-response is regularized. Flux and first-moment directions are
removed to first order in fitting; normalization and exact finite-grid moment
offsets are separately reported for the fitted response.

With all valid cores retained, field baseline chi-square/dof is48.70–57.22.
The four leave-one-exposure-out corrected statistics are approximately10.98,
10.06,5.73,13.09. Thus the positive response removes most of the coherent
profile mismatch, but every held-field check still fails. It is inadequate.

`transfer_positive_field.py` applies this same learned response unchanged to
all16 consumed target stamps, using their original windows and unchanged ERR.
Only centroid, flux and the original quadratic diffuse background are fitted.
It uses the archived official GWCS. No physical lens prediction enters fitting.

- Pixel chi-square/dof becomes6.224–18.314; all still fail.
- Raw response-origin centroids give G1 five-coordinate RMS0.362383mas and
  exposure consistency chi-square145.866 (15degrees of freedom).
- The explicitly reported finite-grid moment convention gives RMS0.360240mas
  and exposure chi-square154.559.
- Formal G1 chi-square remains approximately829/5 in that convention.
- Centroid covariances are conditional on the learned response. Its calibration
  uncertainty and between-exposure instrument covariance remain unpropagated.

The correction improves both the profile residual and exposure consistency on
targets without using their lens residuals to learn it. It does not finish the
calibration. All validation flags remain false. No score was used to tune the
shape after transfer.

An important spatial limitation emerges from the full inputs: the four field
positions span approximately(900,876) to(1474,309), about192native pixels per
exposure, rather than only subpixel dithers. They approximately follow a line;
the targets are roughly400pixels from the star within each exposure. One shared
residual response cannot establish a two-dimensional detector-dependent
calibration. The existing spatial library describes its own calibration epoch;
its residual at this epoch may also vary spatially.

Outputs: `positive_field_results.json`, `positive_field_shape.npz`,
`positive_target_transfer.json`. `positive_response_synthetic_checks.json`
records exact fixed-response recovery for injected affine and quadratic
background cases: maximum centroid error2.73e-14pixel. This verifies the
implementation, not physical calibration.

## Remaining-response mechanism tests

`spatial_field_response.py` partitions residuals and tests two fixed, limited
first-order spatial alternatives using only the four field-star exposures.
The shared positive correction leaves6,746.567 total training chi-square. Of
that, only0.31% lies within radius1nativepixel and3.78% within radius2;
63.64% is at radius2–6. Restored SAT-flagged valid rates contribute10.11%.
Thus the remaining mismatch is chiefly a wing-response error, not a failure
dominated by the restored cores. Charge redistributed outward is not excluded
by this partition.

A single spatial amplitude gradient changes the four held-field chi-square
sum from11,279.766 to11,148.769 (1.16%); a smooth low-resolution spatial
log-response gradient changes it to11,235.202 (0.40%). All held-field fits
still fail. These limited gradients do not explain the remainder.

The calibrators span only about±0.165nativepixel perpendicular to their
dither track, while the quasar positions lie167.6–248.4pixels off that track.
The full2D affine design has condition number2,596.8, with target prediction
leverages290,715–638,839. A free transverse calibration inferred from these
four locations would be unsupported extrapolation.

`fine_field_response.py` tests exactly one finer common-response family,
spacing0.5 and width0.35pixel with unchanged priorSD0.2. It worsens held-field
chi-square to24,319.412. Each three-exposure training set has only849profiled
data directions for1,249coefficients, leaving400directions (including3gauges)
defined by the prior. The improved training fit is not a repair.

## Local quasar response development

`local_quasar_response.py` then addresses that specific data-rank limitation
using only the already-consumed16F115W quasar stamps. It uses the same finer
positive-response family, independently free image centroids/fluxes and the
original quadratic backgrounds. Four whole-exposure folds train on12stamps
and apply the frozen response to the other4. It never uses a physical lens
prediction or lens score for calibration.

`refine_local_quasar_response.py` finishes numerical convergence with the same
objective and each fold's own training-only iterate. All training and held
point fits converge. Every training set now has1,246profiled shape directions
for1,249coefficients: only the3explicit gauge directions remain. The star-only
rank deficiency has therefore been resolved, while the response still fails.

- Held-exposure chi-square sum39,064.541 over4,475nominal conditional degrees
  of freedom; individual held fit chi-square/dof4.524–16.004.
- The all16-image training fit has chi-square19,080.446. Its individual
  conditional statistics range2.722–6.522per degree of freedom. These are
  training quantities and omit shape-estimation degrees of freedom.
- `local_quasar_geometry.py` reads physical predictions only after the saved
  calibration receipt is finished. It performs no parameter fitting.
- All-image response: G1 raw-origin RMS0.475428mas, exposure statistic452.854.
- Cross-exposure response: G1 raw-origin RMS0.428035mas, exposure statistic647.878.
- Finite-grid moment versions are respectively0.473871mas/454.362 and
  0.429619mas/651.723.

The flexible local fit does not improve the earlier field-response geometry.
It must not be selected merely for its improved training profile statistic.
All errors remain native; calibration-shape/systematic uncertainties remain
unpropagated. All physical-validation flags are false. No unused band was read.

Definitive local outputs are `local_quasar_refined_results.json`,
`local_quasar_shape_refined.npz`, and `local_quasar_geometry_results.json`.
The latter binds the calibration receipt by SHA256. Provisional LSMR outputs
are retained separately and must not replace the converged results.
