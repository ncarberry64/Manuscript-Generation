"""Explicit masking of usable JWST rate/CAL slopes with any-group SAT flags.

This adapter leaves the archived fitter and its default mask byte-for-byte
unchanged. It selects finite, positive-error slopes whose original DQ lacks
DO_NOT_USE. SATURATED alone is retained only for a verified native product and
calibration-software convention. This is not a GROUPDQ/raw-ramp policy.
"""
from __future__ import annotations

import numpy as np

from strong_lensing_pixel_geometry import fit_empirical_point

DO_NOT_USE = 1
SATURATED = 2

# Audited release metadata: JWST 2.0.1 requires stcal>=1.18,<2.
# stcal 1.16.0 changed rate-product SAT to mean ANY ramp group saturated.
# Add releases only after checking their actual stcal dependency/semantics.
_VERIFIED_CAL_VERSIONS = {
    "2.0.1": {
        "stcal_requirement": ">=1.18,<2",
        "dependency_source": (
            "https://github.com/spacetelescope/jwst/blob/2.0.1/pyproject.toml"
        ),
    }
}
_SATURATION_SEMANTICS_SOURCE = (
    "https://stcal.readthedocs.io/en/1.17.x/stcal/ramp_fitting/description.html"
)


def _native_slope_provenance(primary_header):
    """Check native 2D slope-product provenance, not a filename guess alone."""
    header = dict(primary_header)
    for name, expected in (
        ("TELESCOP", "JWST"),
        ("INSTRUME", "NIRCAM"),
        ("DATAMODL", "ImageModel"),
        ("S_RAMP", "COMPLETE"),
    ):
        if str(header.get(name, "")).strip() != expected:
            raise ValueError(f"Expected native JWST/NIRCam slope {name}={expected}")
    filename = str(header.get("FILENAME", "")).strip()
    if filename.endswith("_cal.fits"):
        product = "cal"
    elif filename.endswith("_rate.fits"):
        product = "rate"
    else:
        raise ValueError("Only native _cal.fits and _rate.fits products are supported")
    version = str(header.get("CAL_VER", "")).strip()
    if version not in _VERIFIED_CAL_VERSIONS:
        raise ValueError(
            f"Unverified CAL_VER={version!r}; check its SAT flag semantics first"
        )
    return {
        "policy": "jwst_native_rate_any_group_saturation_v1",
        "native_filename": filename,
        "native_product": product,
        "cal_ver": version,
        "s_ramp": header["S_RAMP"],
        "sat_flag_semantics": "Any ramp group saturated; usable early-group slope may remain",
        "bad_pixel_rule": "DO_NOT_USE, nonfinite SCI/ERR, or nonpositive ERR excluded",
        "saturation_semantics_source": _SATURATION_SEMANTICS_SOURCE,
        **_VERIFIED_CAL_VERSIONS[version],
    }


def jwst_native_rate_mask(data, error, dq, *, primary_header):
    """Return a boolean usable-slope mask and its verified provenance receipt.

    Pixel DQ is informational about the ramp history. The pipeline has already
    fitted unsaturated segments. This function never turns a missing, nonfinite,
    zero-error or DO_NOT_USE slope into a measurement and never alters its ERR.
    It is deliberately inapplicable to raw/group arrays and resampled mosaics.
    """
    receipt = _native_slope_provenance(primary_header)
    data, error, dq = np.asarray(data), np.asarray(error), np.asarray(dq)
    if data.ndim != 2 or data.shape != error.shape or data.shape != dq.shape:
        raise ValueError("Native SCI, ERR and pixel DQ must have the same 2D shape")
    if not np.issubdtype(dq.dtype, np.integer):
        raise ValueError("Pixel DQ must be an integer flag array")
    if np.issubdtype(dq.dtype, np.signedinteger) and np.any(dq < 0):
        raise ValueError("Pixel DQ flags cannot be negative")
    flags = dq.astype(np.uint64, copy=False)
    good = (
        np.isfinite(data)
        & np.isfinite(error)
        & (error > 0)
        & ((flags & DO_NOT_USE) == 0)
    )
    receipt.update({
        "input_pixels": int(data.size),
        "selected_pixels": int(good.sum()),
        "retained_sat_flagged_slopes": int((good & ((flags & SATURATED) != 0)).sum()),
        "original_dq_modified": False,
        "science_or_error_modified": False,
        "validation_status": "MASK_SEMANTICS_CHECK_ONLY",
    })
    return good, receipt


def fit_jwst_native_rate_point(data, error, dq, guess, psf_grid, oversampling,
                              psf_origin, *, primary_header, stamp_radius=8,
                              centroid_bound=5., background_degree=1):
    """Apply the unchanged point fitter with the verified native-slope selection.

    A private 0/DO_NOT_USE *selection* array bridges the archived fitter's legacy
    interface. It is not a replacement for the input DQ and is never persisted
    as an observational product. All reported DQ counts use the original flags.
    This fixes eligibility, not PSF adequacy, centroid calibration or lens physics.
    """
    good, receipt = jwst_native_rate_mask(
        data, error, dq, primary_header=primary_header
    )
    selection_dq = np.where(good, 0, DO_NOT_USE).astype(np.uint32)
    result = fit_empirical_point(
        data, error, selection_dq, guess, psf_grid, oversampling, psf_origin,
        stamp_radius=stamp_radius, centroid_bound=centroid_bound,
        background_degree=background_degree,
    )
    x0, y0, x1, y1 = result["stamp_bounds_xyxy"]
    original_flags = np.asarray(dq)[y0:y1, x0:x1].astype(np.uint64)
    selected = good[y0:y1, x0:x1]
    result["dq_flag_counts_in_stamp"] = {
        str(1 << bit): int(np.count_nonzero(original_flags & (np.uint64(1) << np.uint64(bit))))
        for bit in range(64)
        if np.any(original_flags & (np.uint64(1) << np.uint64(bit)))
    }
    receipt["fitted_stamp_selected_pixels"] = int(selected.sum())
    receipt["fitted_stamp_retained_sat_flagged_slopes"] = int(
        (selected & ((original_flags & SATURATED) != 0)).sum()
    )
    result["pixel_mask_receipt"] = receipt
    return result

