# Native read-history diagnostic: HE 0435 field star

## Finding and scope

The ordinary calibrated image combines different temporal support across the star: the three central pixels each retain only two usable groups, while the pixels at radii 3–8 pixels retain all nine. Reconstructing a common first-two-read interval, and applying the recorded full-detector calibrations to it, reduces the normalized field-star profile residual by **32.44%** on the same 289 pixels. The central residual improves, but several wing annuli worsen. The shorter exposure also gives roughly three times larger formal centroid uncertainties. This is evidence that the read-history-dependent response deserves explicit treatment; it does not identify one instrumental mechanism or complete the precision calibration.

The early-interval goodness-of-fit probability is **0.00984653077**, below the existing 0.01 threshold. It must not be rounded into a pass. All physical and precision-validation flags remain false.

Only the first already-consumed F115W exposure and its existing isolated field star enter this diagnostic. There is no quasar fit, lens-parameter change, new band, independently unused observation, residual variance inflation, or BHSM preference. The immutable repository baseline is `7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9` in `ncarberry64/Manuscript-Generation`. The original failures and frozen protocols remain preserved.

## What was calculated

The native SHALLOW2 exposure has nine groups, two reads per group, three skipped frames between groups, and a frame time of 10.73677 seconds. Its raw ZEROFRAME is the first read at time T, and its first averaged group is approximately `(read1 + read2) / 2`. We reconstructed `read2 = 2 * raw_group0 - raw_ZEROFRAME` **before** nonlinear correction. Onboard averaging and rounding prevent calling this a bit-exact reconstruction of the original second ADC sample.

A synthetic two-read RAPID model then passed through the public JWST 2.0.1 group scaling, DQ initialization, saturation, superbias, full-detector SIRS reference-pixel correction, linearity and matching rapid-dark subtraction. Steps recorded as skipped in the archived product remain skipped. The two reads each have their own validity/saturation checks, and the original first-group DO_NOT_USE/JUMP exclusions are retained. The estimator is the corrected read difference divided by T; this cancels a shared reset pedestal. A ZEROFRAME/T image would require additional reset and intercept assumptions and was not used for this comparison.

Both images were fit on exactly the same 17×17 support at global detector bounds `[892,868,909,885]` (zero-based, exclusive upper bounds). The same dated-WSS, fixed-temperature/diffusion/jitter STPSF model was used for both. Flux, centroid and a background plane are measurement nuisance parameters; the physical lens model remains untouched.

The full dark reference is 6,291,483,840 bytes. To stay within the 8 GiB memory limit, exact first-42-frame and first-two-frame subsets were made, retaining the original SCI, ERR, DQ and DQ_DEF arrays without approximating values. SHALLOW2 consumes rapid-frame pairs `[0,1]`, `[5,6]`, through `[40,41]`. The original full-file hash and derived subset hashes are in `dark_range_receipt.json` and `dark_subset_receipt.json`.

## Native numerical replay and the retained failed attempt

The archived CAL records JWST 2.0.1 and CRDS `jwst_1535.pmap`; it does not record its `stcal` version. An initial native replay using `stcal` 1.19.1 reproduced preprocessing but failed the declared rate comparison. The eight largest discrepancies carried JUMP_DET, localizing the disagreement to segmented slopes. The official 1.19.0 release changed segment-slope weights by removing Poisson variance from the weighting. Replaying only ramp fitting and gain scaling on the **same corrected ramp**, with compatible `stcal` 1.18.0, recovers the archived rate and ERR well inside the unchanged gates. This establishes compatible numerical behavior; it does not directly identify a version absent from the original metadata. [Official stcal change log](https://stcal.readthedocs.io/en/latest/stcal/changes.html)

| Replay on the field stamp | Initial stcal 1.19.1 | Compatible stcal 1.18.0 | Unchanged tolerance |
|---|---:|---:|---:|
| Maximum relative rate difference | 0.01859835 | 8.96568×10⁻⁸ | 1×10⁻⁴ |
| Maximum relative CAL ERR difference | 0.000474002 | 1.65304×10⁻⁷ | 1×10⁻³ |
| Rate/ERR gate | Failed | Passed | Both required |

The failed initial results are `common_time_field_results.json`; the final results are `common_time_field_stcal118_results.json`. No failed gate was overwritten or relaxed. `suppress_one_group=False` was already in the native settings. No unreleased ZEROFRAME variance patch was applied.

## Matched-support results

| Quantity | Archived full-slope CAL | Common first-two-read interval |
|---|---:|---:|
| Retained pixels | 289 | 289 |
| χ² / degrees of freedom | 34.86343742 | 1.206444889 |
| Fit probability | Numerically 0 | **0.00984653077** |
| Residual L2 / fitted source L2 | 0.07665300929 | 0.05178729343 |
| RMS residual / fitted source peak | 0.005445113916 | 0.003660932296 |
| Raw residual RMS, MJy/sr | 7.059584825 | 4.594862095 |
| Formal centroid σx, pixels | 0.001677229180 | 0.005199611246 |
| Formal centroid σy, pixels | 0.001622531240 | 0.005172133824 |

The normalized L2 reduction is **32.4393212%**. The fitted centroid changes by `(+0.01421214361, −0.01157624684)` detector pixels. Both optimizations succeed without active bounds. Formal centroid uncertainties increase by approximately 3.10× and 3.19×; they remain conditional estimates, not calibrated astrometric accuracy.

The mean-profile comparison must accompany the much smaller χ² because a 10.73677-second difference image has legitimately larger photon/read noise than the native slope image. A smaller χ² alone would not demonstrate repair. Here the raw normalized residual does improve, predominantly in the core, while the wings give a mixed result:

| Radius, pixels | Pixels | Full-slope summed residual / modeled source in annulus | Early-interval value |
|---|---:|---:|---:|
| 0–1 | 3 | −7.7531% | −3.4398% |
| 1–2 | 9 | +4.9833% | +7.4931% |
| 2–3 | 16 | +1.9502% | +3.6010% |
| 3–4 | 21 | +7.2004% | +23.0422% |
| 4–6 | 62 | +13.9113% | +19.1010% |
| 6–8 | 90 | −1.1452% | −37.6340% |
| 8–12 | 88 | −9.8023% | −27.3102% |

The outer annuli contain faint modeled source flux, so these fractions should be read together with the archived residuals normalized to the total fitted flux. No independent significance is assigned to the difference between these correlated images.

## Directly observed temporal support

Counting groups without GROUPDQ DO_NOT_USE or SATURATED gives the following. JUMP_DET starts a segment and is not itself removed by this support count.

| Radius, pixels | Distribution of usable groups |
|---|---|
| 0–1 | All 3 pixels: 2 groups |
| 1–2 | 6 pixels: 2; 2 pixels: 5; 1 pixel: 7 |
| 2–3 | 1 pixel: 5; 5 pixels: 7; 10 pixels: 9 |
| 3–4 | All 21 pixels: 9 groups |
| 4–6 | All 62 pixels: 9 groups |
| 6–8 | All 90 pixels: 9 groups |

Each of the three central pixels has native GROUPDQ sequence `0,0,1,2,2,2,2,2,2`. These counts come from corrected native group arrays, not from a hypothesized detector model. The associated functional missing from a single static image treatment is a response that accounts for each pixel's actual read history and calibration. This experiment does not yet determine whether charge migration, time-varying optical response, nonlinear group averaging, or another component dominates.

## Uncertainty model and limitations

The comparison propagates the actual gain, CDS read-noise, bias-error and flat-error maps. The field medians are approximately 2.05 electrons/DN, 9.338 DN CDS read noise and 0.706% fractional flat error. The READNOISE reference is a correlated-double-sample quantity, so each independent read contributes half its variance. [JWST 2.0.1 read-noise reference](https://jwst-pipeline.readthedocs.io/en/2.0.1/jwst/references_general/readnoise_reffile.html) and [gain reference](https://jwst-pipeline.readthedocs.io/en/2.0.1/jwst/references_general/gain_reffile.html)

For the rate difference `D=(L2−L1)/T`, the diagonal variance before flat correction contains photon variance `max(D+dark_rate,0)/(gain*T)`, read variance `R_CDS²*(L1prime²+L2prime²)/(2*T²)` and the common bias-reference term `ERR_bias²*(L2prime−L1prime)²/T²`. Flat propagation divides by `flat²` and adds `D²*ERR_flat²/flat⁴`; the photometric conversion multiplies the final variance by PHOTMJSR². There is no residual-derived rescaling.

These limitations remain explicit:

- Full SIRS/reference-induced spatial covariance and the paired covariance of the early/full estimators are unavailable here.
- Photon weights are estimated from observed early counts clipped at zero, so they are data dependent.
- Dark-reference error/covariance is not propagated. The first two dark SCI planes are zero across this stamp, while their ERR planes are nonzero and strongly alike. Adding the two ERR planes as independent variances would be unjustified.
- Reconstructing read2 retains onboard averaging/rounding error.
- The native dark step does not subtract dark from ZEROFRAME; the archived corrected ZEROFRAME must not be mistaken for a fully dark-subtracted image.
- A single already-consumed exposure cannot independently validate the proposed response or distinguish instrumental mechanisms.

The read timing and reset distinctions follow the [official NIRCam readout-pattern documentation](https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-detector-overview/nircam-detector-readout-patterns). No change to the BHSM microscopic preparation, physical admission guard, or lens geometry follows from this detector diagnostic.

## Files and reproduction

The compact files are sufficient to audit the saved support counts and arithmetic without downloading the large dark reference:

- `native_corrected_field_ramp.npz`: exact 66×66 native corrected groups, GROUPDQ, corrected ZEROFRAME, reconstructed corrected reads, both rate versions, archived CAL and reference-map stamps. Its receipt is `native_corrected_field_stamp_receipt.json`.
- `common_time_field_stcal118_results.json`: canonical final comparison.
- `common_time_field_stcal118_fit_receipts.json` and `common_time_field_stcal118_arrays.npz`: fitted pixels, model predictions and propagated variance components.
- `common_time_field_results.json` and associated arrays/fit receipts: preserved initial failed native replay gate.
- `first_two_reads_field_checkpoints.npz`: pre-linearity, post-linearity and post-dark field read arrays.
- `compact_field_audit.json`: independently recomputed saved arithmetic and group-support audit.
- `native_input_manifest.json`: exact inputs and hashes; `native_small_file_manifest.json`: small deliverable hashes.
- `environment.txt`: initial JWST 2.0.1/stcal 1.19.1 environment; `environment_stcal118.txt`: final compatible stcal 1.18.0 environment.

From this directory, a Python environment with NumPy can run:

```bash
python audit_compact_field_receipts.py
```

Full replay requires the original relative layout: this directory is `REPAIR/ramp_diagnostic/native_pipeline`, the repository snapshot is `REPAIR/repo`, and the fixed PSF is under `REPAIR/stpsf_equal_sampling`. Use the exact immutable repository head and input hashes. Download the three original bias/linearity/saturation references and raw input listed in `native_input_manifest.json` to `REPAIR/ramp_diagnostic`; download the six additional references to `native_pipeline/refs`. The CAL and fixed PSF locations are also recorded in that manifest.

To reproduce both the historical failed replay and its numerical reconciliation, start an isolated environment with the versions in `environment.txt` (JWST 2.0.1, stcal 1.19.1), then run:

```bash
python download_recorded_refs.py
python prepare_dark_subsets.py
python replay_detector1.py
python replay_first_two_reads.py
python compare_common_time_field.py
python -m pip install 'stcal==1.18.0'
python reramp_stcal118.py
python compare_common_time_field.py --rate-replay first_exposure_stcal118_rate.fits --output-prefix common_time_field_stcal118
python archive_corrected_field_stamp.py
python audit_compact_field_receipts.py
```

The installed detector pipeline may consult public CRDS server metadata despite explicit reference overrides; full replay therefore needs public network access or a populated CRDS cache. `resume_dark_ranges.py` records this run's interrupted-download recovery and is not required for a fresh complete download. Keep original receipts before rerunning in place: these scripts are replay utilities and write their named output files. Output FITS metadata may carry new runtime timestamps; numerical replay, original input hashes and the recorded tolerance gates are the reproducibility claims.

No full repository test suite was rerun for this detector experiment. The saved comparisons and compact audit were executed, and a separate read-only review independently confirmed the final replay errors, common support, normalized residual reduction, fit probability and uncertainty limitations.
