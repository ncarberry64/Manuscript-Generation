"""Measure raw science/reference covariance from independent detector darks.

Source support is the actual first18 RAPID frames selected as nine SHALLOW2
pairs. SIRS is applied to each TWO-READ AVERAGE before measuring covariance.
No corrected-CDS variance is added again, no noise is fitted to lens residuals.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/"data/strong_lensing_response_completion_20261021"
DATA = BASE/"dark_calibration"
RAW = BASE/"native_cache/dark_calibration"
OUT = ROOT/"artifacts/strong_lensing_response_completion_20261021/dark_calibration"
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_reference_covariance import conditional_sirs


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8*1024**2), b""):
            h.update(b)
    return h.hexdigest()


def dump(path, value):
    def default(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, default=default, allow_nan=False)+"\n",
                         encoding="utf-8", newline="\n")


def paired_integration_difference(first, second, *, allow_masked_nonfinite=False):
    """Cov[(integration1-integration0)/sqrt2]=C if independent equal-C ramps.

    Same-pixel fixed bias/dark mean cancels within an epoch. Independence and
    equality of the two integrations' noise are declared calibration assumptions.
    """
    a, b = np.asarray(first, float), np.asarray(second, float)
    if (a.shape != b.shape or a.ndim < 2
            or (not allow_masked_nonfinite and (not np.isfinite(a).all() or not np.isfinite(b).all()))):
        raise ValueError("Matching finite calibration integrations required")
    return (b-a)/np.sqrt(2.)


def pooled_time_covariance(values, valid):
    """Spatially pooled second moment under explicit zero-mean/stationarity.

    One epoch supplies ONE independent integration-difference realization.
    Averaging different detector pixels is an additional stationary ensemble
    assumption. This is not an empirically known full pixelwise covariance.
    """
    a = np.asarray(values, float)
    good = np.asarray(valid, bool)
    if a.ndim != 3 or good.shape != a.shape[1:] or not good.any():
        raise ValueError("Group,row,column field and nonempty calibration mask required")
    selected = a[:, good]
    if not np.isfinite(selected).all():
        raise ValueError("Finite selected calibration pixels required")
    second = selected@selected.T/selected.shape[1]
    means = selected.mean(axis=1)
    centered = (selected-means[:, None])@(selected-means[:, None]).T/selected.shape[1]
    return {"zero_mean_second_moment": second, "spatially_centered_covariance": centered,
            "spatial_mean": means, "number_of_spatial_samples": selected.shape[1]}


def spatial_lag_second_moment(image, valid, dy, dx):
    """Actual finite-lag scalar moment, with both endpoints valid; no wrap."""
    a = np.asarray(image, float)
    m = np.asarray(valid, bool)
    if a.ndim != 2 or a.shape != m.shape or dy < 0 or dx < 0 or (dy == 0 and dx == 0):
        raise ValueError("A finite image/mask and positive nonwrapped lag required")
    yend, xend = a.shape[0]-dy, a.shape[1]-dx
    left, right = a[:yend, :xend], a[dy:, dx:]
    good = m[:yend, :xend] & m[dy:, dx:]
    if not good.any():
        raise ValueError("No valid finite-lag pairs")
    return {"dy": dy, "dx": dx, "pairs": int(good.sum()),
        "uncentered_cross_DN2": float(np.mean(left[good]*right[good])),
        "adjacent_difference_filtered_variance_DN2": float(np.mean((left[good]-right[good])**2)/2)}


def native_group_sirs(raw_reads, bias, control_xy, control_dq, kernels):
    """Same group-averaging/bias/SIRS channel, preserving native buffer stages."""
    ny, nx = bias.shape
    if raw_reads.shape != (18, ny, nx) or (ny, nx) != (2048, 2048):
        raise ValueError("Exact18 full-frame individual raw ADC reads required")
    groups = raw_reads.reshape(9, 2, ny, nx).astype(np.float32).mean(axis=1)
    pre = (groups-bias).astype(np.float32)
    post = np.empty_like(pre)
    controls = np.empty((9, len(control_xy)))
    clipping = []
    xx = np.arange(nx)
    amp, parity = xx//512, xx % 2
    for group in range(9):
        dtype = np.float32 if group == 0 else np.float64
        science = np.asarray(pre[group], dtype=dtype)
        r = science[control_xy[:, 1], control_xy[:, 0]]
        # A single dummy science feature is sufficient for the forward offsets
        # and whole-detector templates; no16k-feature covariance/eigenfloor.
        op = conditional_sirs(r, control_xy, control_dq, np.array([[4, 4]]), kernels,
            np.zeros(len(r)), native_work_dtype=dtype)
        first = np.asarray(science-op["native_offsets"][amp, parity][None], dtype=dtype)
        template = np.repeat(op["side_template"][:, ::-1].T, 512, axis=1)
        post[group] = np.asarray(first-template, np.float32)
        controls[group] = r
        clipping.append({"group": group, "native_work_dtype": np.dtype(dtype).name,
            "top_bottom_clip_counts": op["top_bottom_clip_counts"],
            "left_side_kept": int(op["side_good_left"].sum()),
            "right_side_kept": int(op["side_good_right"].sum())})
    return pre, post, controls, clipping


def measure():
    from astropy.io import fits
    from stdatamodels.jwst import datamodels
    from jwst.refpix.optimized_convolution import make_kernels
    acquisition = json.loads((OUT/"raw_dark_acquisition.json").read_text())
    cache = BASE/"native_cache"
    with datamodels.SuperBiasModel(cache/"jwst_nircam_superbias_0180.fits") as m:
        bias = np.asarray(m.data, np.float32).copy()
        bdq = np.asarray(m.dq, np.uint32).copy()
    with datamodels.MaskModel(cache/"jwst_nircam_mask_0075.fits") as m:
        dq = np.asarray(m.dq, np.uint32).copy() | bdq
    with datamodels.SaturationModel(cache/"jwst_nircam_saturation_0100.fits") as m:
        dq |= np.asarray(m.dq, np.uint32)
        threshold = np.asarray(m.data).copy()
    with datamodels.SIRSKernelModel(cache/"jwst_nircam_sirskernel_0002.asdf") as m:
        kernels = make_kernels(m, "NRCB3", 1., 30)
    # Metadata only: coordinate extraction never consults any lens residual.
    native = np.load(BASE/"histories/jw01198007001_02101_00001_nrcb3_histories.npz")
    control_xy = native["reference_control_native_xy"]
    control_dq = dq[control_xy[:, 1], control_xy[:, 0]]
    xy = []
    for e in range(1, 5):
        with np.load(BASE/f"histories/jw01198007001_02101_0000{e}_nrcb3_histories.npz") as h:
            xy.append(np.stack((h["x_pixel"][:, 24:41, 24:41], h["y_pixel"][:, 24:41, 24:41]), axis=-1))
    xy = np.asarray(xy, int)
    flatxy = xy.reshape(-1, 2)
    report = []
    epoch_arrays = []
    started = time.perf_counter()
    for epoch, record in enumerate(acquisition["records"]):
        pre_ramps, post_ramps, control_ramps, clipping = [], [], [], []
        saturation_counts = []
        for subset in record["subsets"]:
            rawpath = ROOT/subset["file"]
            if sha(rawpath) != subset["sha256"]:
                raise RuntimeError("Exact calibration raw subset changed")
            raw = np.load(rawpath, mmap_mode="r")
            saturation_counts.append(np.any(raw >= threshold[None], axis=0))
            pre, post, controls, masks = native_group_sirs(raw, bias, control_xy, control_dq, kernels)
            pre_ramps.append(pre); post_ramps.append(post); control_ramps.append(controls)
            clipping.append(masks)
            print("DARK_NATIVE_GROUP_SIRS", epoch, subset["integration_zero_based"], flush=True)
        pre_delta = paired_integration_difference(*pre_ramps, allow_masked_nonfinite=True)
        post_delta = paired_integration_difference(*post_ramps, allow_masked_nonfinite=True)
        correction_delta = pre_delta-post_delta
        control_delta = paired_integration_difference(*control_ramps)
        good = (((dq & 1) == 0) & ~np.logical_or.reduce(saturation_counts)
                & np.all(np.isfinite(pre_delta) & np.isfinite(post_delta), axis=0))
        # This fixed particle/hot-pixel screen is independent of lens residuals.
        # Both unscreened and screened moments are retained; no mask is chosen
        # by covariance adequacy. Truncation is a declared conditional limit.
        stable = np.max(abs(pre_delta-pre_delta[0]), axis=0) < 500.
        clean = good & stable
        arrays = {"matched_patch_native_xy": xy, "control_native_xy": control_xy,
            "pre_SIRS_reference_difference_DN": control_delta,
            "pre_SIRS_matched_patch_difference_DN": pre_delta[:, flatxy[:, 1], flatxy[:, 0]].reshape(9, *xy.shape[:-1]),
            "post_SIRS_matched_patch_difference_DN": post_delta[:, flatxy[:, 1], flatxy[:, 0]].reshape(9, *xy.shape[:-1]),
            "reference_correction_matched_patch_difference_DN": correction_delta[:, flatxy[:, 1], flatxy[:, 0]].reshape(9, *xy.shape[:-1]),
            "matched_patch_valid_DQ": clean[flatxy[:, 1], flatxy[:, 0]].reshape(xy.shape[:-1]),
            "native_frame_indices": np.array([[5*j, 5*j+1] for j in range(9)]),
            "SIRS_kernel_left": np.asarray(kernels[0]), "SIRS_kernel_right": np.asarray(kernels[1])}
        amp_reports = []
        for amp in range(4):
            sl = np.s_[4:2044, 512*amp+4:512*(amp+1)-4]
            mask = clean[sl]
            pre = pre_delta[:, sl[0], sl[1]]
            post = post_delta[:, sl[0], sl[1]]
            correction = correction_delta[:, sl[0], sl[1]]
            ps = pooled_time_covariance(pre, mask)
            qs = pooled_time_covariance(post, mask)
            cs = pooled_time_covariance(correction, mask)
            cross = pre[:, mask]@correction[:, mask].T/mask.sum()
            identity = ps["zero_mean_second_moment"]+cs["zero_mean_second_moment"]-cross-cross.T
            zpre, zpost = pre[1]-pre[0], post[1]-post[0]
            early_variance = lambda c: float(c[1, 1]+c[0, 0]-2*c[0, 1])
            zero_moment = qs["zero_mean_second_moment"]
            # Profiling a first-group intercept leaves8 group contrasts. These
            # are measured source matrices, no raw READNOISE identification.
            contrast = np.c_[-np.ones(8), np.eye(8)]
            contrast_cov = contrast@zero_moment@contrast.T
            arrays[f"amp{amp}_pre_time_second_moment_DN2"] = ps["zero_mean_second_moment"]
            arrays[f"amp{amp}_post_time_second_moment_DN2"] = zero_moment
            arrays[f"amp{amp}_correction_time_second_moment_DN2"] = cs["zero_mean_second_moment"]
            arrays[f"amp{amp}_pre_correction_cross_DN2"] = cross
            arrays[f"amp{amp}_post_group_contrast_covariance_DN2"] = contrast_cov
            arrays[f"amp{amp}_post_spatially_centered_time_covariance_DN2"] = qs["spatially_centered_covariance"]
            # Long-range coherent row behavior stays visible; not silently
            # absorbed in a white per-pixel variance estimate.
            row_count = mask.sum(axis=1)
            row_mean = np.where(mask[None], post, 0.).sum(axis=2)/row_count[None]
            arrays[f"amp{amp}_post_row_means_DN"] = row_mean
            row_cov = row_mean@row_mean.T/row_mean.shape[1]
            arrays[f"amp{amp}_post_row_mean_second_moment_DN2"] = row_cov
            lags = [spatial_lag_second_moment(zpost, mask, dy, dx)
                    for dy, dx in ((0, 1), (1, 0), (0, 2), (2, 0), (0, 8), (8, 0),
                                   (0, 16), (16, 0), (0, 32), (32, 0), (0, 64), (64, 0))]
            unscreened = pooled_time_covariance(post, good[sl])
            amp_reports.append({"amp_zero_based": amp, "detector_region_xyxy": [512*amp+4, 4, 512*(amp+1)-4, 2044],
                "valid_native_DQ_unsaturated_pixels": int(good[sl].sum()),
                "fixed500DN_particle_screen_excluded": int((good[sl] & ~stable[sl]).sum()),
                "screened_spatial_samples": int(mask.sum()),
                "pre_time_second_moment_DN2": ps["zero_mean_second_moment"],
                "post_time_second_moment_DN2": zero_moment,
                "post_spatially_centered_time_covariance_DN2": qs["spatially_centered_covariance"],
                "post_spatial_mean_DN": qs["spatial_mean"],
                "pre_correction_cross_DN2": cross,
                "correction_time_second_moment_DN2": cs["zero_mean_second_moment"],
                "full_cross_identity_max_abs_DN2": float(abs(identity-zero_moment).max()),
                "pre_early_group_contrast_RMS_DN": np.sqrt(early_variance(ps["zero_mean_second_moment"])),
                "post_early_group_contrast_RMS_DN": np.sqrt(early_variance(zero_moment)),
                "post_early_group_contrast_RMS_unscreened_DN": np.sqrt(early_variance(unscreened["zero_mean_second_moment"])),
                "post_group_contrast_covariance_min_eigenvalue_DN2": float(np.linalg.eigvalsh(contrast_cov).min()),
                "post_row_mean_early_contrast_RMS_DN": np.sqrt(early_variance(row_cov)),
                "post_early_spatial_lags": lags})
        array_path = DATA/f"epoch{epoch}_measured_covariance.npz"
        np.savez_compressed(array_path, **arrays)
        # Preserve full fields for independent PSD/kernel estimation without
        # any new source read or target-data use. No generated covariance floor.
        for kind, value in (("pre", pre_delta), ("post", post_delta)):
            np.save(RAW/f"epoch{epoch}_{kind}_integration_difference_DN.npy", value.astype(np.float32))
        np.save(DATA/f"epoch{epoch}_calibration_valid.npy", clean)
        row = {"source": record["source"], "source_EXPSTART": record["primary_header"]["EXPSTART"],
            "array_file": str(array_path.relative_to(ROOT)).replace('\\', '/'), "array_sha256": sha(array_path),
            "amplifiers": amp_reports, "native_SIRS_clipping": clipping,
            "independent_integration_difference_realizations": 1,
            "pixelwise_full_covariance_empirical_rank_at_most": 1,
            "source_complete_instrument_covariance": False}
        report.append(row); epoch_arrays.append(arrays)
        dump(OUT/"measured_covariance_progress.json", report)
    held = []
    for amp in range(4):
        ctrain = epoch_arrays[0][f"amp{amp}_post_group_contrast_covariance_DN2"]
        ctest = epoch_arrays[1][f"amp{amp}_post_group_contrast_covariance_DN2"]
        expectation = float(np.trace(linalg.solve(ctrain, ctest, assume_a="pos")))
        held.append({"amp_zero_based": amp, "earliest_training_latest_held": True,
            "held_expected_Q_per_8group_contrast_pixel": expectation,
            "expected_if_same_stationary_covariance": 8.,
            "held_to_train_early_contrast_variance_ratio": float(ctest[0, 0]/ctrain[0, 0]),
            "not_a_physical_or_formal_independent_pixel_chi2_p_value": True})
    final = {"status": "actual public calibration raw ADC/reference/post-SIRS covariance measured",
        "records": report, "held_epoch_diagnostic": held, "seconds": time.perf_counter()-started,
        "read_support": "nine2-read averages ofRAPID[0,1],[5,6],...,[40,41], before SAME bias/SIRS channel; no per-RAPID SIRS averaging substitution",
        "units": "ADC DN and DN^2 at pre-linearity/post-SIRS stage; no IPC correction/deconvolution applied",
        "covariance_identity": "Cpost=Css+Ccc-Csc-Ccs from the SAME calibration realization; no correctedCDS+newreference double count",
        "ensemble_assumptions": ["within-epoch integrations independent and same mean/noise law; difference/sqrt2 removes fixed bias/dark mean",
            "per-amplifier spatial pooling is an explicit stationary noise approximation; pixels are not declared independent for significance",
            "earliest epoch training, latest epoch held for same-noise-state diagnostics; no lens-residual tuning",
            "500DN fixed particle/hot-pixel screen conditional; unscreened variance also preserved",
            "onboard integer averaging/ADC quantization cannot be exactly reconstructed from ground RAPID pair averages; difference documented"],
        "remaining_relation": ["one difference realization per epoch cannot identify unrestricted pixelwise covariance; raw spatial/time fields and cross terms are now available for a separately specified stationary response/noise model",
            "exposure-state dependence and finite calibration sampling uncertainty remain; stationary pooled moment matrix is not itself full cross-stamp covariance",
            "photon charge-to-voltage IPC/NLIPC covariance kernel and collection/fullwell law remain distinct from raw read/dark covariance",
            "this detector counting calibration is not a BHSM dressed optical-state kernel or a preparation contraction"],
        "input_acquisition_sha256": sha(OUT/"raw_dark_acquisition.json"),
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (
            Path(__file__), ROOT/"code/strong_lensing_reference_covariance.py")},
        "source_complete_instrument_covariance": False, "noise_rescaling": False,
        "no_reserved_science_observations": True, "physical_guard": 2e-10, "physical_validation": False}
    dump(OUT/"measured_covariance.json", final)
    print("DARK_MEASURED", json.dumps(held, default=lambda a: a.item()), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute-independent-dark-calibration", action="store_true")
    args = parser.parse_args()
    if not args.execute_independent_dark_calibration:
        parser.error("Explicit independent calibration-data execution flag required")
    measure()


if __name__ == "__main__":
    main()
