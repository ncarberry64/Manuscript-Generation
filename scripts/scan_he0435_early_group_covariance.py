"""Actual two-group A/FIELD conditional mean/covariance scan, no reserved data.

This is not a response calibration: the fixed native rate pilot has consumed
later groups.  It checks a fully stated small-signal counting law on the first
two unsaturated group averages, preserving raw residuals and original κ=0 noise.
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
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_group_calibration import averaged_read_observation, second_order_mean_correction
from strong_lensing_stochastic_response import (
    BoundaryMomentModel, nearest_neighbor_edges, pixel_read_projection,
    project_collection_covariance, independent_adc_read_covariance,
    shared_calibration_covariance, fixed_residual_covariance_application,
    harmonic_rate_initializer,
    within_pixel_read_covariance_first_order,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def application(native, patch, kappa):
    # A 49x49 finite charge context excludes the one negative FIELD native-rate
    # pixel far outside the 17x17 output. No negative rate is clipped or floored.
    sl = np.s_[8:57, 8:57]
    shape = (49, 49)
    rate_DN_s = np.asarray(native["native_rate_pre_gain_scale_DN_s"][patch][sl], float)
    gain = np.asarray(native["gain_e_per_DN"][patch][sl], float)
    known_rate = np.isfinite(rate_DN_s) & (rate_DN_s >= 0)
    rate, continuation = harmonic_rate_initializer(gain*rate_DN_s, known_rate)
    rate = rate.ravel()
    read_times = np.asarray(native["native_read_times_s"][:2], float)
    times = read_times.ravel()
    # Readout calibration acts only on the observed 17x17 pixels. Context
    # carriers still propagate through the full 49x49 graph. No invalid native
    # polynomial at a distant DNU context pixel is replaced by an identity map.
    coeff = np.asarray(native["linearity_coefficients"][patch][:, 24:41, 24:41], float)
    anchor = np.asarray(native["pre_linearity_group_DN"][patch, :2, 24:41, 24:41], float)
    data = np.asarray(native["corrected_group_DN"][patch, :2, 24:41, 24:41], float)
    dark_group = np.asarray(native["native_subtracted_dark_group_DN"][patch, :2, 24:41, 24:41], float)
    dark_indices = np.asarray(native["native_dark_RAPID_frame_indices"][:2], int)
    dark_read = np.asarray(native["dark_RAPID42_DN"][patch, dark_indices, 24:41, 24:41], float)
    raw = np.asarray(native["raw_group_DN"][patch, :2, 24:41, 24:41], float)
    sat = np.asarray(native["saturation_threshold_DN"][patch, 24:41, 24:41], float)
    # Saturation is recorded in raw group-DN coordinates. The native additive
    # superbias/refpix corrections provide the matching pre-linearity threshold.
    upper = np.min(sat[None]-(raw-anchor), axis=0)
    model = BoundaryMomentModel(rate, kappa, nearest_neighbor_edges(shape),
        max_time_s=float(times[-1]),
        domain_source="first two native group supports only; no SAT flags; individual-read ADC ceiling additionally checked")
    source_context = np.stack([model.original_leading_mean(t).reshape(shape)/gain for t in times]).reshape(2, 2, *shape)
    if np.any(source_context < 0):
        raise ValueError("Original leading source leaves the nonnegative charge domain")
    source = source_context[:, :, 16:33, 16:33]
    latent_without_intercept = source+dark_read
    initial_ADC = np.repeat(anchor[:, None], 2, axis=1)
    intercept = data[0]-np.mean(latent_without_intercept[0], axis=0)
    for _ in range(12):
        observed = averaged_read_observation(latent_without_intercept+intercept,
            coeff, initial_ADC, dark_group_dn=dark_group, adc_upper=upper)
        delta = (observed["mean_DN"][0]-data[0])/observed["read_to_group_jacobian"][0].sum(axis=0)
        intercept -= delta
        if np.max(abs(delta)) < 1e-8:
            break
    else:
        raise RuntimeError("First-group electronic-intercept solve failed")
    observed = averaged_read_observation(latent_without_intercept+intercept,
        coeff, initial_ADC, dark_group_dn=dark_group, adc_upper=upper)
    # Output pixels must be unflagged in BOTH groups, including observed jumps.
    dq = native["groupdq"][patch, :2, 24:41, 24:41]
    pdq = native["pixeldq"][patch, 24:41, 24:41]
    good = np.all((dq & 7) == 0, axis=0) & ((pdq & 1) == 0)
    full_indices = np.arange(49*49).reshape(shape)[16:33, 16:33][good]
    native_xy = np.column_stack((native["x_pixel"][patch, 24:41, 24:41][good],
                                native["y_pixel"][patch, 24:41, 24:41][good]))
    output_indices = np.flatnonzero(good.ravel())
    jac = observed["read_to_group_jacobian"].reshape(2, 2, -1)[:, :, output_indices]
    derivative_ratio = jac[1].sum(axis=0)/jac[0].sum(axis=0)
    photon_weights = np.c_[-derivative_ratio[:, None]*jac[0].T, jac[1].T]/gain.ravel()[full_indices, None]
    result = project_collection_covariance(model, times,
        pixel_read_projection(full_indices, photon_weights, len(rate)),
        provenance="actual native polynomial read derivatives/gain; group1 minus group0 intercept-propagation contrast; frozen native adaptive rate pilot")
    group_d = observed["native_group_linearity_derivative"].reshape(2, -1)[:, output_indices]
    adc_weights = np.c_[np.repeat((-derivative_ratio*group_d[0]/2)[:, None], 2, axis=1),
                       np.repeat((group_d[1]/2)[:, None], 2, axis=1)]
    cds = np.asarray(native["readnoise_CDS_DN"][patch][sl], float).ravel()
    adc_stage = ("explicit incomplete approximation: identify reference/IPC-corrected native0282 CDS marginals "
                 "with independent raw ADC reads; source does not establish this identification")
    adc = independent_adc_read_covariance(pixel_read_projection(full_indices, adc_weights, len(rate)),
                                         cds, stage_assumption=adc_stage)
    # Keep repeated SAME physical bias error shared between both groups.
    # Its cross-pixel covariance is not supplied by the native ERR map. This
    # diagonal-in-reference-pixel component is therefore explicitly conditional.
    bias_derivative = -group_d[1]+derivative_ratio*group_d[0]
    bias_sigma = native["superbias_ERR_DN"][patch][sl].ravel()[full_indices]
    bias_j = np.diag(bias_derivative)
    bias = shared_calibration_covariance(bias_j, np.diag(np.asarray(bias_sigma, float)**2),
        source="actual native bias ERR marginals; diagonal spatial prior is an explicit incomplete approximation; physical pixel shared across groups")
    covariance = result["covariance"]+adc+bias
    residual = (data[1]-observed["mean_DN"][1]).ravel()[output_indices]
    # Preserve a κ=0 native Poisson/CDS comparator even for nonzero scan points.
    zero = BoundaryMomentModel(rate, 0., model.edges, max_time_s=float(times[-1]),
                               domain_source=model.domain_source)
    original_poisson = project_collection_covariance(zero, times,
        pixel_read_projection(full_indices, photon_weights, len(rate)),
        provenance="same native contrast/J/g map; independent κ=0 photon comparator")["covariance"]
    frozen_error = np.sqrt(np.diag(original_poisson+adc))
    score = fixed_residual_covariance_application(residual, covariance, frozen_error,
        covariance_scope="conditional local boundary counting + native CDS single-read + spatial-diagonal shared bias approximation; SIRS/IPC/refcov missing")
    # The comparator is derived group-domain Poisson/CDS noise, not an ERR
    # array supplied by native RampModel (which has none at this stage).
    score["same_mapping_kappa_zero_diagonal_Poisson_ADC_chi2"] = score.pop("native_ERR_chi2")
    chol = linalg.cho_factor(covariance, lower=True)
    logdet = float(2*np.log(np.diag(chol[0])).sum())
    score["conditional_gaussian_Q_plus_logdet"] = score["conditional_stochastic_chi2"]+logdet
    margin = upper[None, None]-observed["inverse_ADC_reads_DN"]
    correction = result["resummed_minus_original_mean"]
    # Execute the defining Hessian contraction as a separate, UNAPPLIED
    # diagnostic. This does not silently change either frozen source mean.
    read_cov = within_pixel_read_covariance_first_order(model, times, full_indices)
    read_cov_DN = read_cov/gain.ravel()[full_indices][None, None]**2
    within_group = np.stack([read_cov_DN[2*j:2*j+2, 2*j:2*j+2] for j in range(2)])
    hessian = observed["read_to_group_hessian"].reshape(2, 2, 2, -1)[:, :, :, output_indices]
    photon_second_mean = second_order_mean_correction(hessian, within_group)
    # ADC noise is injected after inverse-L; its independent average variance
    # is CDS^2/(2*N), with N=2 here.
    group_dd = observed["native_group_linearity_second_derivative"].reshape(2, -1)[:, output_indices]
    adc_second_mean = .5*group_dd*cds[full_indices][None]**2/4
    second_contrast = (photon_second_mean[1]+adc_second_mean[1]
                      -derivative_ratio*(photon_second_mean[0]+adc_second_mean[0]))
    report = {"label": str(native["labels"][patch]), "kappa_e_inverse": kappa,
        "number_of_output_contrasts": len(residual), "number_of_excluded_17x17_pixels": int(289-good.sum()),
        "finite_charge_context_shape": list(shape), "output_native_xy": native_xy.tolist(),
        "undefined_context_rate_continuation": continuation,
        "minimum_inferred_individual_read_ADC_saturation_margin_DN": float(np.min(margin)),
        "first_group_intercept_max_abs_residual_DN": float(np.max(abs(observed["mean_DN"][0]-data[0]))),
        "electronic_intercept_range_linear_DN": [float(intercept.min()), float(intercept.max())],
        "photon_source_rate_e_s_range": [float(rate.min()), float(rate.max())],
        "gain_e_per_DN_range": [float(gain.min()), float(gain.max())],
        "ADC_noise_stage_assumption": adc_stage,
        "source_complete_instrument_covariance": False,
        "dark_read_DN_max_abs": float(np.max(abs(dark_read))),
        "group_dark_DN_max_abs": float(np.max(abs(dark_group))),
        "dark_read_ERR_DN_max_abs": float(np.max(abs(native["dark_RAPID42_ERR_DN"][patch, dark_indices, 24:41, 24:41]))),
        "full_context_native_SAT_flags_first_two_groups": int(np.sum((native["groupdq"][patch, :2, 8:57, 8:57]&2)!=0)),
        "max_resummed_minus_original_contrast_mean_DN": float(np.max(abs(correction))),
        "executed_but_unapplied_native_mean_curvature": {
            "physical_count_covariance_order": "O(kappa); distinct from resummed covariance score",
            "photon_group_correction_DN_range": [float(photon_second_mean.min()), float(photon_second_mean.max())],
            "ADC_group_correction_DN_range": [float(adc_second_mean.min()), float(adc_second_mean.max())],
            "max_abs_contrast_correction_DN": float(abs(second_contrast).max()),
            "mean_was_changed": False},
        "components_variance_DN2_ranges": {"photon_boundary": [float(np.diag(result["covariance"]).min()), float(np.diag(result["covariance"]).max())],
            "independent_ADC": [float(np.diag(adc).min()), float(np.diag(adc).max())],
            "shared_bias_conditional_diagonal_spatial_prior": [float(np.diag(bias).min()), float(np.diag(bias).max())]},
        "projection_diagnostics": result["diagnostics"], "residual_application": score,
        "first_group_was_used_for_intercept": True,
        "native_rate_pilot_consumed_later_groups": True,
        "not_an_independent_calibration_or_validated_test": True,
        "physical_validation": False}
    arrays = {"residual_DN": residual, "native_xy": native_xy,
        "photon_read_weights_DN_per_e": photon_weights, "ADC_read_weights": adc_weights,
        "bias_derivative": bias_derivative, "photon_boundary_covariance_DN2": result["covariance"],
        "independent_ADC_covariance_DN2": adc, "shared_bias_conditional_covariance_DN2": bias,
        "conditional_covariance_DN2": covariance, "frozen_poisson_ADC_error_DN": frozen_error,
        "first_group_mean_DN": observed["mean_DN"][0], "second_group_mean_DN": observed["mean_DN"][1],
        "source_original_leading_read_DN": source, "source_original_leading_context_read_DN": source_context,
        "electronic_intercept_linear_DN": intercept,
        "inferred_ADC_reads_DN": observed["inverse_ADC_reads_DN"], "native_read_times_s": read_times,
        "context_rate_e_s": rate.reshape(shape), "context_gain_e_per_DN": gain,
        "native_group_jacobian": observed["read_to_group_jacobian"],
        "native_group_hessian": observed["read_to_group_hessian"]}
    arrays.update(photon_native_second_order_group_mean_DN=photon_second_mean,
        ADC_native_second_order_group_mean_DN=adc_second_mean,
        unapplied_second_order_contrast_mean_DN=second_contrast,
        physical_within_group_read_covariance_DN2=within_group)
    return report, arrays


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute-native-small-signal", action="store_true")
    args = parser.parse_args()
    if not args.execute_native_small_signal:
        parser.error("Actual consumed-native application requires --execute-native-small-signal")
    source = ROOT/"data/strong_lensing_response_completion_20261021/histories/jw01198007001_02101_00001_nrcb3_histories.npz"
    output = ROOT/"artifacts/strong_lensing_response_completion_20261021/stochastic/native_early_groups"
    output.mkdir(parents=True, exist_ok=True)
    records = []
    start = time.perf_counter()
    with np.load(source) as native:
        for patch in (0, 4):
            for kappa in (0., 1e-7, 2e-7, 4e-7):
                record, arrays = application(native, patch, kappa)
                name = f"{record['label']}_kappa_{kappa:.0e}".replace("+", "")
                array_path = output/(name+".npz")
                np.savez_compressed(array_path, **arrays)
                record["array_file"] = str(array_path.relative_to(ROOT)).replace('\\', '/')
                record["array_sha256"] = sha(array_path)
                records.append(record)
                (output/"progress.json").write_text(json.dumps(records, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
                print(record["label"], kappa, record["residual_application"], flush=True)
    report = {"status": "actual first-exposure target A and matching FIELD early-group conditional scan executed",
        "records": records, "seconds": time.perf_counter()-start,
        "input": str(source.relative_to(ROOT)).replace('\\', '/'), "input_sha256": sha(source),
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (
            Path(__file__), ROOT/"code/strong_lensing_stochastic_response.py", ROOT/"code/strong_lensing_group_calibration.py")},
        "individually_unresolved": [
            "NRCB3 collection-boundary law and kappa not identified by a source-mean fit; scan does not identify them",
            "actual SIRS/control covariance and its shared projection, IPC/NLIPC, complete spatial reference covariance",
            "actual0282 CDS reference was already reference/bias/IPC corrected; raw ADC covariance is not uniquely recoverable from those marginal CDS values",
            "finite-context incident rates at native DQ-marked context pixels; pilot is conditional, not a new physical state solve",
            "post-saturation collection/overflow: groups after first two explicitly not extrapolated",
            "native nonlinearity Hessian mean contraction is executed separately at O(kappa) and not silently added to the deterministic map; exact higher moments require a fully selected counting law",
            "gain and native coefficient reference covariance not supplied",
            "BHSM U,V source contractions and statistical-state kernel are distinct unresolved preparation inputs"],
        "physical_guard": 2e-10, "physical_validation": False,
        "no_reserved_observations_consumed": True,
        "no_error_rescaling": True, "no_physical_lens_parameters_changed": True}
    (output/"results.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
