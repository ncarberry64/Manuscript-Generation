"""Eligibility regressions for partially saturated, valid JWST slope products."""
from pathlib import Path
import sys
import numpy as np
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "code"))
# Standalone patch test; a real checkout already has the original module.
ARCHIVE = HERE.parent / "repo"
if (ARCHIVE / "code").exists():
    sys.path.insert(1, str(ARCHIVE / "code"))

from strong_lensing_jwst_rate_mask import (
    jwst_native_rate_mask, fit_jwst_native_rate_point,
)
from strong_lensing_pixel_geometry import native_pixel_mask, fit_empirical_point

HEADER = {
    "TELESCOP": "JWST", "INSTRUME": "NIRCAM", "DATAMODL": "ImageModel",
    "S_RAMP": "COMPLETE", "CAL_VER": "2.0.1",
    "FILENAME": "jw01198007001_02101_00001_nrcb3_cal.fits",
}


def test_presaturation_slope_is_retained_but_unusable_slopes_are_rejected():
    sci = np.array([[10., 11., 12., np.nan, 14., 15., 16.]])
    err = np.array([[1., 2., 3., 1., 0., np.inf, 1.]])
    dq = np.array([[0, 2, 3, 2, 2, 2, 4]], dtype=np.uint32)
    copies = [x.copy() for x in (sci, err, dq)]
    good, receipt = jwst_native_rate_mask(sci, err, dq, primary_header=HEADER)
    assert good.tolist() == [[True, True, False, False, False, False, True]]
    assert receipt["retained_sat_flagged_slopes"] == 1
    for x, before in zip((sci, err, dq), copies):
        np.testing.assert_array_equal(x, before)
    assert native_pixel_mask(sci, err, dq).tolist() == [
        [True, False, False, False, False, False, True]
    ]


def test_raw_mosaic_or_unverified_semantics_are_not_silently_accepted():
    changes = [
        {"CAL_VER": "1.20.2"}, {"CAL_VER": "99.0.0"}, {"CAL_VER": ""},
        {"DATAMODL": "RampModel"}, {"FILENAME": "jw01198_i2d.fits"},
        {"FILENAME": "jw01198_uncal.fits"}, {"S_RAMP": "SKIPPED"},
        {"INSTRUME": "NIRISS"},
    ]
    for change in changes:
        with unittest.TestCase().assertRaises(ValueError):
            jwst_native_rate_mask(
                np.ones((3,3)), np.ones((3,3)), np.zeros((3,3), np.uint32),
                primary_header=HEADER | change,
            )


def test_group_or_invalid_flag_arrays_are_rejected():
    cases = [np.zeros((1,2,3),np.uint32), np.full((2,3),-1,np.int32),
             np.zeros((2,3),float)]
    for dq in cases:
        with unittest.TestCase().assertRaises(ValueError):
            jwst_native_rate_mask(
                np.ones(dq.shape),np.ones(dq.shape),dq,primary_header=HEADER
            )


def test_adapter_preserves_original_dq_and_legacy_behavior_on_synthetic_point():
    from scipy.special import erf
    yy, xx = np.mgrid[:21, :21]
    # Exact pixel-integrated Gaussian grid sampled at quarter pixels.
    gridcoord = (np.arange(101)-50)/4
    z1 = erf((gridcoord+.5)/(np.sqrt(2)*.95))
    z0 = erf((gridcoord-.5)/(np.sqrt(2)*.95))
    grid = np.outer((z1-z0)/2, (z1-z0)/2)
    from strong_lensing_pixel_geometry import empirical_psf_values
    truecenter = (10.12, 9.83)
    data = 2000*empirical_psf_values(grid,4,(50.,50.),xx,yy,truecenter) + 2.
    error = np.ones_like(data)*.5
    dq = np.zeros(data.shape, np.uint32)
    dq[9:12, 9:12] = 2
    dq[4, 4] = 1
    before = dq.copy()
    old = fit_empirical_point(data,error,dq,(10,10),grid,4,(50,50))
    new = fit_jwst_native_rate_point(
        data,error,dq,(10,10),grid,4,(50,50),primary_header=HEADER
    )
    np.testing.assert_array_equal(dq,before)
    assert new["usable_pixels"] == old["usable_pixels"]+9
    assert new["pixel_mask_receipt"]["fitted_stamp_retained_sat_flagged_slopes"] == 9
    assert new["dq_flag_counts_in_stamp"]["2"] == 9
    np.testing.assert_allclose(new["centroid_xy_pixel"],truecenter,atol=2e-6)
    # More correctly retained information must sharpen this correctly specified test.
    assert np.trace(new["centroid_covariance_pixel2"]) < np.trace(old["centroid_covariance_pixel2"])



if __name__ == "__main__":
    tests = [value for name,value in globals().copy().items()
             if name.startswith("test_") and callable(value)]
    suite = unittest.TestSuite(unittest.FunctionTestCase(test) for test in tests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    raise SystemExit(not result.wasSuccessful())
