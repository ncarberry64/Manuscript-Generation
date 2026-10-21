"""Compare fixed field-star PSF on matched support after a common-time extraction.

Both images are already consumed development data. No lens parameters change.
Noise is calculated from count differences and recorded detector references;
no residual rescaling is used. Absolute validation stays false because the
full spatial/reference covariance is not calibrated by this one exposure.
"""
from pathlib import Path
import os
import sys
import json
import argparse
import numpy as np
from astropy.io import fits

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
sys.path.insert(0, str(ROOT / "repo/code"))
from strong_lensing_pixel_geometry import fit_empirical_point

CAL = ROOT / "repo/data/strong_lensing_geometry_20261021/native/jw01198007001_02101_00001_nrcb3_cal.fits"
PSF = ROOT / "stpsf_equal_sampling/jw01198007001_02101_00001_nrcb3_FIELD_stpsf_epsf.fits"
SL = np.s_[868:885, 892:909]
OFFSET = np.array([892., 868.])
GUESS = np.array([899.7549609755204, 875.9549994961739]) - OFFSET


def image(path, name):
    with fits.open(path, memmap=False) as hd:
        return np.array(hd[name].section[SL], dtype=float)


def summarize(fit):
    raw = np.array(fit["pixel_data"])
    model = np.array(fit["pixel_model"])
    x, y = np.array(fit["pixel_x"]), np.array(fit["pixel_y"])
    b = np.array(fit["background_plane_SCI"])
    background = b[0] + b[1] * (x - GUESS[0]) + b[2] * (y - GUESS[1])
    source = model - background
    residual = raw - model
    answer = {key: value for key, value in fit.items()
              if not key.startswith("pixel_") and key != "standardized_residual"}
    answer["centroid_xy_pixel"] = (np.array(answer["centroid_xy_pixel"]) + OFFSET).tolist()
    answer["stamp_bounds_xyxy"] = (np.array(answer["stamp_bounds_xyxy"]) + np.tile(OFFSET, 2)).astype(int).tolist()
    answer["coordinate_frame"] = "Global native detector pixels, zero-based"
    answer["background_reference_xy_pixel"] = (GUESS + OFFSET).tolist()
    answer["chi2_per_dof"] = fit["chi2"] / fit["dof"]
    answer["residual_RMS_MJy_sr"] = float(np.sqrt(np.mean(residual**2)))
    answer["residual_RMS_over_fitted_source_peak"] = float(np.sqrt(np.mean(residual**2)) / source.max())
    answer["residual_L2_over_fitted_source_L2"] = float(np.linalg.norm(residual) / np.linalg.norm(source))
    answer["source_peak_MJy_sr"] = float(source.max())
    answer["centroid_sigma_pixel"] = np.sqrt(np.diag(fit["centroid_covariance_pixel2"])).tolist()
    radius = np.hypot(x - GUESS[0], y - GUESS[1])
    radial = []
    edges = [0, 1, 2, 3, 4, 6, 8, 12]
    for low, high in zip(edges[:-1], edges[1:]):
        choose = (radius >= low) & (radius < high)
        if choose.any():
            radial.append({"radius_pixel": [low, high], "pixels": int(choose.sum()),
                           "summed_residual_over_bin_source": float(residual[choose].sum() / source[choose].sum()),
                           "summed_residual_over_total_fitted_flux": float(residual[choose].sum() / fit["flux_in_SCI_pixel_units"])})
    answer["coherent_radial_residuals"] = radial
    answer["physical_response_validated"] = False
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rate-replay", default="first_exposure_native_rate.fits")
    parser.add_argument("--output-prefix", default="common_time_field")
    args = parser.parse_args()
    with fits.open(CAL, memmap=False) as hd:
        cal = np.array(hd["SCI"].section[SL], float)
        cal_error = np.array(hd["ERR"].section[SL], float)
        cal_dq = np.array(hd["DQ"].section[SL], np.uint32)
        photmjsr = float(hd["SCI"].header["PHOTMJSR"])
    with fits.open(OUT / "first_two_reads_corrected_ramp.fits", memmap=False) as hd:
        reads = np.array(hd["SCI"].section[0, :, SL[0], SL[1]], float)
        groupdq = np.array(hd["GROUPDQ"].section[0, :, SL[0], SL[1]], np.uint32)
        pixeldq = np.array(hd["PIXELDQ"].section[SL], np.uint32)
        interval = float(hd[0].header["TFRAME"])
    refs = OUT / "refs"
    gain = image(refs / "jwst_nircam_gain_0094.fits", "SCI")
    readnoise_cds = image(refs / "jwst_nircam_readnoise_0282.fits", "SCI")
    flat = image(refs / "jwst_nircam_flat_0667.fits", "SCI")
    flat_error = image(refs / "jwst_nircam_flat_0667.fits", "ERR")
    flat_dq = image(refs / "jwst_nircam_flat_0667.fits", "DQ").astype(np.uint32)
    bias_error = image(OUT.parent / "jwst_nircam_superbias_0180.fits", "ERR")
    with fits.open(OUT.parent / "jwst_nircam_linearity_0054.fits", memmap=False) as hd:
        coefficients = np.array(hd["COEFFS"].section[:, SL[0], SL[1]], float)
    checkpoint = np.load(OUT / "first_two_reads_field_checkpoints.npz")
    before_linearity = checkpoint["refpix_field_reads"]
    derivative = np.zeros_like(before_linearity, dtype=float)
    for power in range(1, len(coefficients)):
        derivative += power * coefficients[power] * before_linearity ** (power - 1)
    dark_removed = checkpoint["linearity_field_reads"] - checkpoint["dark_current_field_reads"]
    dark_rate = (dark_removed[1] - dark_removed[0]) / interval
    rate = (reads[1] - reads[0]) / interval
    # Photon noise uses the total pre-dark photon rate. Readnoise reference
    # is CDS DN: each single read contributes R_CDS^2/2 before linearity.
    variance_photon = np.maximum(rate + dark_rate, 0.) / (gain * interval)
    variance_read = readnoise_cds**2 * np.sum(derivative**2, axis=0) / (2 * interval**2)
    variance_bias = bias_error**2 * (derivative[1] - derivative[0])**2 / interval**2
    variance_unflat = variance_photon + variance_read + variance_bias
    variance_flat = rate**2 * flat_error**2 / flat**4
    variance = variance_unflat / flat**2 + variance_flat
    early = rate / flat * photmjsr
    early_error = np.sqrt(variance) * photmjsr
    common = (np.all((groupdq & 3) == 0, axis=0) & ((pixeldq & 1) == 0)
              & ((cal_dq & 1) == 0) & ((flat_dq & 1) == 0)
              & np.isfinite(cal) & np.isfinite(cal_error) & (cal_error > 0)
              & np.isfinite(early) & np.isfinite(early_error) & (early_error > 0)
              & np.isfinite(gain) & (gain > 0) & (flat > 0))
    selection_dq = np.where(common, 0, 1).astype(np.uint32)
    psf = fits.getdata(PSF)
    origin = ((psf.shape[1] - 1) / 2, (psf.shape[0] - 1) / 2)
    fullfit = fit_empirical_point(cal, cal_error, selection_dq, GUESS, psf, 5, origin)
    earlyfit = fit_empirical_point(early, early_error, selection_dq, GUESS, psf, 5, origin)
    # The original saved CAL provides the independent numerical-replay target.
    with fits.open(OUT / args.rate_replay, memmap=False) as hd:
        native_rate = np.array(hd["SCI"].section[SL], float)
        native_var_p = np.array(hd["VAR_POISSON"].section[SL], float)
        native_var_r = np.array(hd["VAR_RNOISE"].section[SL], float)
    equivalent_rate = cal * flat / photmjsr
    comparable = np.isfinite(native_rate) & np.isfinite(equivalent_rate) & ((cal_dq & 1) == 0)
    difference = native_rate[comparable] - equivalent_rate[comparable]
    reconstructed_cal_variance = photmjsr**2 * ((native_var_p + native_var_r) / flat**2
                              + native_rate**2 * flat_error**2 / flat**4)
    relative_error_diff = np.sqrt(reconstructed_cal_variance[comparable]) / cal_error[comparable] - 1
    rate_relative_difference = np.abs(difference) / np.maximum(np.abs(equivalent_rate[comparable]), 1e-10)
    replay_pass = bool(np.max(rate_relative_difference) < 1e-4 and np.max(np.abs(relative_error_diff)) < 1e-3)
    with fits.open(refs / "derived_jwst_nircam_dark_0337_first2.fits", memmap=False) as hd:
        dark_error = np.array(hd["ERR"].section[:, SL[0], SL[1]], float)
    report = {
        "status": "COMMON_TIME_FIELD_RESPONSE_DEVELOPMENT_DIAGNOSTIC",
        "source_observation": CAL.name, "fixed_psf": str(PSF),
        "native_pipeline_replay": {
            "rate_file": args.rate_replay,
            "field_rate_RMS_difference_DN_s": float(np.sqrt(np.mean(difference**2))),
            "field_rate_max_abs_difference_DN_s": float(np.max(np.abs(difference))),
            "field_rate_max_relative_difference": float(np.max(rate_relative_difference)),
            "CAL_ERR_max_relative_replay_difference": float(np.max(np.abs(relative_error_diff))),
            "numerical_replay_gate": {"rate_max_relative_tolerance": 1e-4, "ERR_max_relative_tolerance": 1e-3, "pass": replay_pass},
        },
        "interval_s": interval,
        "matched_pixels": int(common.sum()),
        "same_pixel_support": fullfit["pixel_x"] == earlyfit["pixel_x"] and fullfit["pixel_y"] == earlyfit["pixel_y"],
        "original_native_CAL_same_support": summarize(fullfit),
        "first_two_reads_same_support": summarize(earlyfit),
        "noise_model": {
            "gain_e_per_DN_median": float(np.median(gain[common])),
            "readnoise_CDS_DN_median": float(np.median(readnoise_cds[common])),
            "flat_fractional_error_median": float(np.median((flat_error / flat)[common])),
            "dark_reference_ERR_max_DN_on_retained_pixels": float(np.max(np.abs(dark_error[:, common]))),
            "formula": "Var(delta_L/T)/flat^2 + rate^2*ERR_flat^2/flat^4; Poisson=max(rate+dark_rate,0)/(gain*T); read=R_CDS^2*(L1prime^2+L2prime^2)/(2*T^2); bias=ERR_bias^2*(L2prime-L1prime)^2/T^2",
            "residual_rescaling": False,
            "limits": ["Conditional diagonal approximation: full SIRS/reference-induced spatial covariance is unavailable", "Photon variance is estimated from observed early counts clipped at zero, so weights are data-dependent", "Dark-reference covariance/error is not propagated into this diagonal calculation", "Onboard averaging reconstruction rounding is not independent photon noise", "This fit does not certify precision astrometry or an adequate PSF", "Lower chi-square at10.73677seconds alone does not demonstrate repaired mean response"],
        },
        "physical_validation": False,
        "new_observations": False,
        "lens_parameter_refitting": False,
    }
    a = report["original_native_CAL_same_support"]
    b = report["first_two_reads_same_support"]
    report["raw_shape_change"] = {
        "residual_fraction_L2_ratio_early_over_full": b["residual_L2_over_fitted_source_L2"] / a["residual_L2_over_fitted_source_L2"],
        "peak_normalized_RMS_ratio_early_over_full": b["residual_RMS_over_fitted_source_peak"] / a["residual_RMS_over_fitted_source_peak"],
        "centroid_shift_pixel": (np.array(b["centroid_xy_pixel"]) - a["centroid_xy_pixel"]).tolist(),
        "interpretation": "Mean-profile diagnostics on identical support; a larger p value alone is not evidence of a repaired response",
    }
    np.savez_compressed(OUT / (args.output_prefix + "_arrays.npz"), full=cal, full_error=cal_error,
        early=early, early_error=early_error, common=common, reads=reads, groupdq=groupdq,
        gain=gain, readnoise_cds=readnoise_cds, flat=flat, flat_error=flat_error,
        variance_photon=variance_photon, variance_read=variance_read, variance_bias=variance_bias,
        linearity_derivative=derivative, x_offset=OFFSET[0], y_offset=OFFSET[1])
    (OUT / (args.output_prefix + "_results.json")).write_text(json.dumps(report, indent=2) + "\n")
    (OUT / (args.output_prefix + "_fit_receipts.json")).write_text(json.dumps({"full": fullfit, "early": earlyfit}, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
