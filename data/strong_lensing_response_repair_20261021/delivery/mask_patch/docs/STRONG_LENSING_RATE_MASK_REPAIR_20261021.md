# HE 0435 native-rate saturation-mask repair

This is a development correction to the already consumed JWST F115W extraction.
It changes which existing calibrated slopes are retained. It does not refit the
lens, alter native SCI/ERR/DQ, supply a new PSF, or validate physical geometry.

## Identified cause

All four native files declare CAL_VER=2.0.1, DATAMODL=ImageModel, and
S_RAMP=COMPLETE. JWST 2.0.1 requires stcal>=1.18,<2. Since stcal 1.16.0,
the rate-product SATURATED bit can indicate that **any** group in the ramp
saturated. The slope may still be valid because the pipeline uses earlier
unsaturated segments. Invalid slopes remain excluded through DO_NOT_USE,
nonfinite SCI/ERR, and nonpositive ERR.

The archived point fitter excluded both DO_NOT_USE and SATURATED. This discarded
303 finite positive-error target slopes across all sixteen point stamps. The new
adapter retains these valid slopes under an explicit product/version check.
The original fitter and old protocol are unchanged.

Primary implementation sources:

- JWST 2.0.1 dependency:
  https://github.com/spacetelescope/jwst/blob/2.0.1/pyproject.toml
- Rate-product flag propagation:
  https://stcal.readthedocs.io/en/1.17.x/stcal/ramp_fitting/description.html
- stcal 1.16.0 flag change:
  https://stcal.readthedocs.io/en/latest/stcal/changes.html
- Upstream flag-only change:
  https://github.com/spacetelescope/stcal/pull/421

## Implementation

code/strong_lensing_jwst_rate_mask.py adds an explicit native 2D CAL/rate policy.
It checks the native-product header and the verified software convention. The
adapter sends a private boolean selection represented by 0/DO_NOT_USE to the
unchanged fitter. Original SCI/ERR/DQ arrays are not mutated, original DQ counts
are restored in the fit receipt, and the selected SAT-flagged count is recorded.

The legacy mask is still the default of the archived fitter. GROUPDQ/raw ramps,
resampled mosaics, missing provenance and unverified software conventions are
not silently interpreted as native 2D slope products. Extend the verified-release
table only after checking the release's actual flag semantics.

## Canonical native-coordinate replay

The replay consumes the four archived native CAL files and their exact saved
PSFs, GWCS and frozen origin. All 16 original fits are reproduced first.

| Quantity | Archived selection | Corrected selection |
|---|---:|---:|
| G1 conventional raw geometry RMS, mas | 0.779400 | 0.448457 |
| G1 R1 raw geometry RMS, mas | 0.779454 | 0.448563 |
| Between-exposure statistic, 15 coordinates | 432.262011 | 253.027686 |
| G1 conventional formal geometry statistic, 5 coordinates | 958.434 | 1033.113139 |
| G1 R1 formal geometry statistic, 5 coordinates | 958.428 | 1033.268517 |

The original global pixel centroids reproduce to 3.41e-12 pixel and original
five-coordinate geometry to 1.21e-10 arcsec. An independent complex-number
construction agrees with the geometry to 8.89e-16 arcsec; an independent finite
difference checks its Jacobian to 5.94e-10.

All sixteen corrected pixel fits still fail their known-ERR profile check:
chi2/dof ranges from 27.57 to 83.55. Between-exposure consistency still fails.
Restoring information reduces the raw geometry residual but also tightens
formal covariance, so the formal geometry statistic does not improve. These
scores are conditional diagnostics under the still inadequate PSF, not calibrated
physical rejection probabilities.

The portable small-stamp replay changes the fitter's numerical coordinate origin
and reproduces the corrected G1 RMS as 0.448447 mas and the exposure statistic as
253.133300. Its old-centroid reproduction differs by at most 1.23e-5 pixel and old
geometry by 0.000530 mas. This is a quantified numerical replay tolerance, not
byte-identical arithmetic. The native global-coordinate result is canonical.

## Existing early-read setting

Inspection of the native FITS ASDF pipeline configuration shows:

    algorithm: OLS_C
    suppress_one_group: False

The logs also record zero-frame processing. Therefore merely enabling
suppress_one_group=False again is not a correction for this dataset.

The current upstream development changelog separately mentions a correction for
read-noise variance of ZEROFRAME-derived OLS_C slopes (an NFRAMES factor).
NFRAMES=2 in these inputs, but the saved CAL DQ does not identify which individual
slopes used ZEROFRAME. No variance adjustment is applied here. That investigation
requires classifying the actual ramp/zero-frame usage and applying an appropriate
pipeline correction; it must not be replaced by an arbitrary error multiplier.

## Reproduction

Run the semantic, invalid-pixel, provenance and synthetic-point regressions:

    python tests/test_strong_lensing_jwst_rate_mask.py

Run the canonical replay with the four original CAL files:

    python scripts/replay_he0435_rate_mask.py --root . --native-dir data/strong_lensing_geometry_20261021/native --output artifacts/strong_lensing_geometry_20261021/mask_repair_native_replay.json

Run the portable replay using only the committed small stamps and GWCS:

    python scripts/replay_he0435_rate_mask.py --root . --output artifacts/strong_lensing_geometry_20261021/mask_repair_replay.json

The receipts record exact input SHA256 values, retained pixels, original-fit
reproduction tolerances, geometry crosschecks and unchanged scientific scope.
The baseline manifest locks the ten original JSON inputs and the normalized
legacy fitter text to the reviewed commit. The runner explicitly selects the
four consumed exposures and protects original inputs/code from output writes.
Four regression groups and two input/output integrity checks passed; the
inherited 123-test suite and journal build were not rerun for this package.
All physical-validation flags remain false. No unused band was opened.

