# HE 0435 response repair: executed results and integration handoff

## Result

An actual pixel-selection error has been corrected and independently replayed.
The original code discarded calibrated slopes whenever a SATURATED flag was
present. In the producing JWST pipeline release, that flag can record a later
saturated group even when earlier unsaturated readings provide a valid slope.
The correction restores 303 usable pipeline-slope pixels across the sixteen F115W fits.
The same frozen G1 predictions then have a smaller raw geometry discrepancy.

A positive PSF correction learned solely from the isolated calibration star
also transfers to the quasar profiles. It reduces the profile and raw geometry
residuals, but does not yet satisfy the pixel or exposure-consistency checks.
The package therefore contains a production-ready mask repair and explicitly
labeled development response calculations. It does not contain a completed
precision calibration or an observational validation of BHSM.

The reviewed repository state is Manuscript-Generation commit
`7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9`, draft PR #23. The existing reference and
R1-inclusive columns remain controls within the retained BHSM optical response.
A fully matched effective-mass representation describes the same retained
surface contribution; this repair introduces no competition requirement or
double-counted gravitational term.

## 1. Confirmed mask error and implemented correction

The original `native_pixel_mask` rejected both DO_NOT_USE and SATURATED.
The four native files report CAL_VER 2.0.1, completed ramp fitting, SHALLOW2,
NFRAMES = 2 and nine groups. The producing release uses the newer any-group
saturation convention. Every restored target pixel has finite SCI, finite
positive ERR and no DO_NOT_USE flag.

The additive adapter checks native two-dimensional JWST/NIRCam CAL/rate
provenance and the audited release. It rejects raw/group-DQ inputs, mosaics and
unverified releases. The original fitter, science/error/DQ arrays, protocols
and numerical checkpoints are preserved. A private selection array connects
the corrected eligibility policy to the unchanged fitter.

The canonical replay uses full native arrays and global detector coordinates.

| Measurement | Published extraction | Corrected mask | Positive field-trained response |
|---|---:|---:|---:|
| G1 reference raw five-coordinate RMS (mas) | 0.779400 | 0.448457087207 | 0.362383367942 |
| Between-exposure statistic, 15 coordinates | 432.262011 | 253.027686060 | 145.866023 |
| Formal G1 geometry statistic, 5 coordinates | about 958.434 | 1033.113139179 | about 837.232 |
| Absolute pixel-response adequacy | Fails | Fails | Fails |

The last column is a separate development model, not part of the production
mask patch. Its centroid covariance is conditional on its learned response;
response uncertainty and shared instrument covariance are not yet propagated.
The quoted last-column RMS uses the raw response-origin convention. The
separately reported finite-grid moment convention gives 0.360240 mas.

Restoring informative pixels tightens the conditional centroid errors. Thus
raw RMS falls while the mask-only formal geometry statistic increases. Both
facts are retained. The smaller raw discrepancy cannot establish that the true
sky geometry has improved by the same amount.

All sixteen mask-corrected profile fits still fail; reduced chi-square spans
approximately 27.57–83.55. Native ERR values have not been multiplied or changed.
The R1 control is nearly identical: corrected-mask RMS 0.448562686643 mas.

The original centroids replay to 3.41e-12 pixel in canonical coordinates. A
portable small-stamp replay gives 0.448446682409 mas and exposure statistic
253.133300295; this small optimizer-coordinate difference is documented, not
called byte-identical or independent validation.

## 2. Field-trained positive response: a real but incomplete improvement

The positive response uses one shared, smooth multiplicative correction to the
library ePSF, constrained by the four isolated-star observations. Lens
predictions never enter this training. A documented background convention
handles negative samples in the empirical library while preserving a positive
final response and reporting its finite-grid normalization and position
convention. The fitted flux absorbs the normalization reparameterization.

Four leave-one-exposure-out field checks remain inadequate: reduced chi-square
is approximately 10.98, 10.06, 5.73 and 13.09. When transferred unchanged to all sixteen
already consumed quasar stamps, the corresponding values improve to 6.224–18.314.
Only point position, flux and the existing diffuse background are fitted during
that transfer. Every absolute response check still fails.

Residual localization narrows the remaining problem. After the shared positive
correction, about 63.64% of the field residual chi-square lies 2–6 native pixels
from the star center. Only 0.31% lies within 1 pixel, 3.78% within 2 pixels, and 10.11%
comes from the restored SAT-flagged pixels. The remaining failure is therefore
not dominated by the restored core pixels.

The four calibrator positions lie almost on a line: their maximum perpendicular
departure is about 0.165 pixel (full range about 0.323 pixel), whereas the quasar positions are roughly 168–248 pixels
off that track. A freely fitted two-dimensional correction cannot be reliably
extrapolated from this star alone. Bounded along-track variations reduce held
field chi-square by only about 0.4–1.2%.

A single finer response-basis trial was also checked. It improves training but
worsens held-field performance, with 1,249 coefficients and only 849 profiled data
directions in each three-exposure training set. That underconstrained trial is
rejected; it is not transferred to targets.

## 3. Physical telescope-response calculation

STPSF 2.2.0 was run with the measured WSS wavefront
`R2022101203-NRCA3_FP1-1.fits`, detector NRCB3, F115W, native detector position,
pixel integration, charge diffusion and IPC/PPC. The wavefront was measured
about 2 h 48 min after the first science exposure and is a pre-correction wavefront
measurement. The archived observation GWCS remains the geometry map.

An initial exploratory run used unequal FFT/detector oversampling. Inspection
found that this STPSF version uses the FFT header value for IPC tap spacing and
grid metadata. That exploratory calculation is superseded by equal 5× sampling;
the corrected physical baseline still fails. This sampling issue belongs to
the new exploratory STPSF calculation, not to the original archived library
PSF extraction.

A bounded physical family then fitted one common stellar blackbody temperature
and one common charge-diffusion width, with jitter fixed. The optimum is about
4198 K and 8.05 mas, inside the declared bounds. The four held-field reduced
chi-squares are 32.13, 38.95, 57.65 and 66.32. Simple color and isotropic diffusion do not
finish the repair.

An accelerated monochromatic construction agrees with direct broadband STPSF
to less than 6e-16 of the PSF peak on the fitting support. Thus that failure is
not explained by the tested spectral-mixture or convolution implementation.
The temperature model and parameter bounds are exploratory assumptions, not
measured stellar spectroscopy or externally calibrated detector coefficients.

A further physical trial added eight common Noll 4–11 pupil-phase modes to the
measured wavefront, with exploratory bounds of ±10 nm and a 5 nm regularizer.
The projected mode rank is 8 and condition number 3.12. Exact nonlinear optical
replays give held-field reduced chi-square 25.79, 20.50, 24.80 and 35.23; two
coefficients reach the bounds. A finite-difference proposal is scored by exact
nonlinear optics; this is not an exhaustive nonlinear wavefront search.
These are instrumental pupil modes, separate from the BHSM envelope/state
phase. The calculation improves some profiles but
does not calibrate the star. Its coefficient errors are conditional on an
inadequate response, not measurements of the telescope's mirror errors.

A local quasar-only calibration then trained the fixed finer positive family
on twelve already consumed stamps and assessed the other four, leaving out a
whole exposure at a time. The profiled shape rank is 1246 of 1249; the remaining
three directions are gauges. Thus the earlier field-only rank limitation is
removed. Nevertheless, the held profiles give total chi-square 39064.54 over 4475
nominal conditional degrees of freedom, with individual reduced values 4.52–16.00.
The all-image training fits also remain inadequate, 2.72–6.52 per conditional
degree of freedom.

Only a separate post-calibration calculation reads the fixed lens predictions.
The local all-image response gives RMS 0.47543 mas and exposure statistic 452.85;
the held-exposure response gives 0.42804 mas and 647.88. This more flexible family
is not promoted. It has not been selected by minimizing a lens residual.
All these trials preserve native errors and are development checks.

## 4. Executed detector-history test

The first already consumed F115W raw exposure was reprocessed with the recorded
JWST 2.0.1 configuration, full-frame SIRS reference-pixel correction, superbias,
linearity, saturation, exact dark reference and jump treatment. No reserved
band was opened. The complete official dark reference is provenance-bound;
exact first-42-frame and first-two-frame subsets retain every sample consumed
by the two calculations while avoiding a large memory allocation.

### Numerical provenance resolved

The native CAL header does not record its underlying stcal version. An initial
stcal 1.19.1 replay failed the declared numerical reproduction gate: the maximum
relative field rate difference was 1.86%, concentrated in pixels with jump-split
segments. A separate stcal 1.18.0 ramp-fit replay on the same corrected ramp
reproduces the archived field rate to a maximum relative difference of
8.97e-8 and its ERR to 1.65e-7. Both are inside the unchanged numerical tolerances
of 1e-4 and 1e-3. The earlier unsuccessful replay is retained.

This identifies a numerically matching ramp implementation; it does not prove
an unrecorded archive environment version. It is a check on these field pixels,
not a byte-identity claim for an entire exposure. The rate-mask repair does not
depend on the segment-weight difference between these versions.

### Common-time extraction

For this SHALLOW2 exposure, ZEROFRAME is the first individual read, not a
measurement at time zero. NFRAMES = 2 allows reconstruction, up to onboard
averaging/quantization, of the second raw read:

\[
R_1=Z,\qquad R_2\simeq 2G_0-Z.
\]

The reconstruction occurs before nonlinear calibration. The two reads receive
the recorded detector corrections separately, then are differenced over
10.73677 s and flat-fielded. A common reset pedestal cancels to leading order.
Noise is computed from the recorded gain, CDS read-noise, bias and flat maps;
no residual-based error multiplier is used.

The corrected native ramp directly confirms different exposure histories:
all three pixels within one native pixel of the star center have only two
usable unsaturated groups. All 21 pixels at radius 3–4 pixels and all 62 at
4–6 pixels retain nine groups. The early-read extraction instead measures all
289 retained pixels over the same interval. Different histories alone would
not change an ideal static linear response; their interaction with the real
detector and optics is the issue being tested.

Both images were fitted on exactly the same 289 pixels with the same dated,
fixed STPSF response, independently fitting only centroid, flux and the original
background plane.

| Field-star diagnostic | Archived full-rate image | First-two-read difference |
|---|---:|---:|
| Raw residual L2 / fitted source L2 | 7.66530% | 5.17873% |
| Raw residual RMS / fitted source peak | 0.544511% | 0.366093% |
| Raw residual RMS (MJy/sr) | 7.05958 | 4.59486 |
| Conditional chi-square / dof | 34.8634 | 1.20644 |
| Conditional centroid sigma, x (pixel) | 0.001677 | 0.005200 |
| Conditional centroid sigma, y (pixel) | 0.001623 | 0.005172 |

The normalized raw residual is 32.44% smaller, a descriptive mean-profile
change that cannot be summarized merely as a larger error bar. But the
chi-square reduction is also strongly affected by lower precision. Median
pixel ERR increases from 0.0928 to 2.810 MJy/sr; centroid standard errors are
about three times larger. The early fit has p = 0.00984653, still below the
predeclared 0.01 threshold. It does not pass that adequacy criterion.

The improvement is concentrated in the central profile; several wing annuli
worsen. This single comparison therefore does not establish a universal PSF
repair or a unique mechanism such as charge migration. Residual nonlinearity,
charge transport and time-dependent optics are still possible contributors.
The diagonal errors omit full SIRS-induced spatial covariance and complete
dark/reference covariance. The two images share reads, so a paired significance
claim would require their cross-covariance. All physical-validation flags remain
false. No quasar position was refitted in this raw-read test.

The definitive receipt is
`ramp_diagnostic/native_pipeline/common_time_field_stcal118_results.json`.
The earlier `common_time_field_results.json` explicitly preserves the failed
1.19.1 numerical gate. Compact corrected group stamps, fit receipts, software
versions, reference hashes and reproduction scripts accompany both results.
The earlier approximate reference-pixel diagnostic is retained only as
historical development evidence and is superseded for native-processing
comparisons by this full pipeline test.

All four original native configurations already specify OLS_C and
`suppress_one_group=False`. Repeating that setting is not a repair. The CAL SAT
bit alone does not identify which slopes used ZEROFRAME. No unreleased
ZEROFRAME variance patch or empirical error scaling was applied.

## 5. The functional calculation that needs completion

The observational prediction needs an explicit incoming-light-to-pixel
operator. Schematically,

\[
\mu_{ep}=\mathcal R_{ep}(\eta,\mathcal H_{ep})
 [I_{\rm sky}(\theta_{\rm BHSM},s)]+b_{ep}.
\]

Here the BHSM calculation and source state determine sky brightness;
the response parameters describe spectral weighting, telescope PSF, pixel
integration, detector response, and the coordinate map. The measured read history
$\mathcal H_{ep}$ records the usable groups, accumulated illumination and rate
estimation for each pixel; the response need not be linear in flux. A single
fixed spatial PSF does not explicitly represent those dependencies. The sky scene includes
the point source and extended source/galaxy light. The current local polynomial
background is a limited approximation to that extended scene.

The needed correction belongs in this measurement model. Its calibration must
come from field or development data and physical instrument information.
Lens-coordinate agreement must not be the fitting objective for the response.
Uncertainty in learned response parameters must enter the joint likelihood and
the correlated centroid covariance. Merely multiplying native errors does not
supply the missing mean response or constitute calibration.

The present evidence establishes a mask defect and a coherent response
mismatch. It does not establish that all remaining lens residuals are
instrumental. Once the measurement procedure is adequate, the fixed physical
prediction may pass or fail on genuinely unused observations.

## 6. Delivered implementation and verification

`mask_patch/jwst_rate_mask_repair.patch` adds ten files to the reviewed commit.
It includes the adapter, executable replay, four regression groups, immutable
input checks, output protection, numerical receipts, environment record and
Codex integration instructions. Its manifest records every file hash.

Executed checks include:

- Four focused regression groups, independently rerun by the root session.
- Rejection of altered baseline predictions and protected output paths.
- All original and corrected target fits in canonical and portable coordinates.
- Patch application checking, actual application and ten output-file hashes.
- Separate numerical and implementation reviews of the development trials.
- Native rate/ERR reconciliation with the archived field and an independently
  audited, matched-support first-two-read calculation.

The inherited 123-test suite and the journal PDF were not rerun for this local
repair package. Original physical admission remains false under the unchanged
2e-10 guard. The microscopic U,V preparation remains a separate unresolved
derivation; the present correction concerns the conditional classical optics
measurement and does not assert that microscopic construction.

## 7. Integration and publication status

The selected GitHub connection returned 403 for branch creation. No remote
branch, commit, PR comment, merge or journal submission was made. The patch and
the complete numerical handoff are ready for the user's external Codex checkout.
Follow `mask_patch/docs/CODEX_HANDOFF.md` and `NEXT_CODEX_TASK.md`.

Keep the original failed experiment intact. These are reanalyses of already
examined F115W photons. Freeze the completed measurement procedure before
opening any remaining reserved observations. The manuscript remains for
Galaxies, with no target submission date added.

## Primary references

- STScI stcal ramp fitting and data-quality semantics:
  https://stcal.readthedocs.io/en/1.17.x/stcal/ramp_fitting/description.html
- stcal release history, including the any-group SAT convention:
  https://stcal.readthedocs.io/en/latest/stcal/changes.html
- JWST 2.0.1 declared dependency range:
  https://github.com/spacetelescope/jwst/blob/2.0.1/pyproject.toml
- STPSF detector effects and adjustable diffusion approximation:
  https://stpsf.readthedocs.io/en/latest/jwst_detector_effects.html
- JWST polynomial linearity correction and reference assumptions:
  https://jwst-pipeline.readthedocs.io/en/latest/jwst/linearity/description.html
- JWST 2.0.1 detector processing and ZEROFRAME handling:
  https://jwst-pipeline.readthedocs.io/en/2.0.1/jwst/pipeline/calwebb_detector1.html
- NIRCam readout patterns and frame-zero limitations:
  https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-readout-patterns
- CDS read-noise reference definition:
  https://jwst-pipeline.readthedocs.io/en/2.0.1/jwst/references_general/readnoise_reffile.html
- Williams et al., TDCOSMO XXVII, Sections IV.5–IV.6: joint extended-light/PSF
  refinement and empirical PSF uncertainty for these systems:
  https://arxiv.org/html/2608.27566v1

The external lens study supplies methodological context, not validation of the
BHSM implementation or of any numerical result in this package.
