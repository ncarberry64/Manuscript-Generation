"""Project independently measured READ+DARK spectra onto actual native outputs.

The calibration spectrum is fixed from the earliest dark epoch. The later dark
epoch and saved HE0435 means are scored without any noise/physical mean refit.
Original ERR residual statistics remain alongside every conditional noise score.
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
from strong_lensing_dark_covariance import (
    dark_field_to_clock_templates, template_spectra, pool_spectral_matrices,
    BinnedSpectralModel, native_to_clock_coordinates,
)
from strong_lensing_native_ols import piecewise_ols_arrays


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8*1024**2), b""):
            h.update(b)
    return h.hexdigest()


def dump(path, value):
    def default(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, default=default, allow_nan=False)+"\n",
                         encoding="utf-8", newline="\n")


def spectral_model():
    protocol = {"status": "fixed before spectral projection and target residual whitening",
        "training": "earliest selected1453 epoch0 integration difference/sqrt2 only",
        "calibration_holdout": "latest selected1453 epoch1; already summarized separately; not historically unused science",
        "clock_cell": [2040, 512], "NRCB3_clock": "native rows reversed; odd amplifiers columns reversed; INTERNAL amplifier science boundaries preserved",
        "channels": "fouramplifiers×ninegroups, same innovation field; cross-amplifier spectra retained",
        "law": "zero-mean stationary periodic Gaussian electronic/read+dark field with PSD matrices measured from independent dark Fourier outer products",
        "frequency_bins": "signed powers-of-two Fourier-index bins, DC/Nyquist separate; piecewise constant perbin; no fitted smoothing choice",
        "mask": "previously fixed mappedDQ/unsaturated/500DN calibration-only screen; masked periodogram normalized by coverage; no window deconvolution",
        "variance": "spectral sums preserved by bin averaging; no noise multiplier, eigenfloor or additional correctedCDS/reference noise",
        "limits": "one independent integration-difference realization per epoch; frequency pooling/stationarity/periodic boundary and finite sampling not identified instrument truth",
        "HE0435_residuals_used_to_calibrate_noise": False, "physical_validation": False}
    protocolpath = OUT/"stationary_spectral_protocol_512.json"
    if not protocolpath.exists(): dump(protocolpath, protocol)
    path = DATA/"earliest_epoch_binned_cross_spectrum_512.npz"
    if not path.exists():
        field = np.load(RAW/"epoch0_post_integration_difference_DN.npy", mmap_mode="r")
        good = np.load(DATA/"epoch0_calibration_valid.npy")
        template, mask = dark_field_to_clock_templates(field, good)
        f, coverage = template_spectra(template, mask)
        pooled = pool_spectral_matrices(f)
        np.savez_compressed(path, **pooled, calibration_window_coverage=coverage)
        print("FIXED_CROSS_SPECTRUM_READY", path.stat().st_size, flush=True)
    with np.load(path) as f:
        model = BinnedSpectralModel(f["spectral_matrices"], f["y_frequency_bin"],
            f["x_frequency_bin"], f["spatial_shape"], int(f["namp"]), int(f["ngroup"]))
        metadata = {"file": str(path.relative_to(ROOT)).replace('\\', '/'), "sha256": sha(path),
            "bytes": path.stat().st_size, "spectral_matrix_shape": list(f["spectral_matrices"].shape),
            "calibration_window_coverage": f["calibration_window_coverage"],
            "minimum_Fourier_modes_per_bin": int(f["bin_frequency_counts"].min()),
            "maximum_Fourier_modes_per_bin": int(f["bin_frequency_counts"].max()),
            "protocol_sha256": sha(protocolpath)}
    return model, metadata


def score(residual, covariance):
    r = np.asarray(residual, float)
    chol = linalg.cho_factor(covariance, lower=True)  # No arbitrary jitter.
    white = linalg.solve_triangular(chol[0], r, lower=True)
    q = float(white@white)
    logdet = float(2*np.log(np.diag(chol[0])).sum())
    return {"n": len(r), "raw_RMS": float(np.sqrt(np.mean(r*r))), "Q": q,
        "Q_per_observation": q/len(r), "logdet": logdet,
        "Gaussian_NLL": .5*(q+logdet+len(r)*np.log(2*np.pi)), "physical_validation": False}, white


def requests_and_saved_residuals(derivative_map="analytic"):
    static = ROOT/"artifacts/strong_lensing_response_completion_20261021/static_diagnostics"
    target = np.load(static/"all_target.npz")
    field = np.load(static/"all_field.npz")
    heldtarget = np.load(static/"held_target_strict.npz")
    trainingfield = np.load(static/"training_field.npz")
    late = np.load(DATA/"epoch1_measured_covariance.npz")
    requests, records = [], []
    for exposure in range(1, 5):
        stem = f"jw01198007001_02101_0000{exposure}_nrcb3"
        historypath = BASE/"histories"/(stem+"_histories.npz")
        dqpath = BASE/"histories"/(stem+"_native_diagnostics_v2.npz")
        h = np.load(historypath); d = np.load(dqpath)
        t = np.flatnonzero(np.char.startswith(target["name"].astype(str), stem))
        f = np.flatnonzero(np.char.startswith(field["name"].astype(str), stem))
        if len(t) != 4 or len(f) != 1:
            raise ValueError("Saved static residual order does not match one native exposure")
        masks = np.concatenate((target["mask"][t], field["mask"][f]), axis=0)
        xy = np.stack((h["x_pixel"][:, 24:41, 24:41], h["y_pixel"][:, 24:41, 24:41]), axis=-1)
        native_xy = xy[masks]
        amp, clock = native_to_clock_coordinates(native_xy)
        finite_jac = np.moveaxis(h["native_rate_jacobian_group_and_zeroframe_h1_per_s"][:, :9, 24:41, 24:41], 1, -1)[masks]
        native_match = None
        branches = []
        if derivative_map == "analytic":
            values, analytic, branches = piecewise_ols_arrays(
                np.moveaxis(h["corrected_group_DN"][:, :, 24:41, 24:41], 1, -1)[masks].T,
                np.moveaxis(h["groupdq"][:, :, 24:41, 24:41], 1, -1)[masks].T,
                gain=h["gain_e_per_DN"][:, 24:41, 24:41][masks],
                readnoise_cds_dn=h["readnoise_CDS_DN"][:, 24:41, 24:41][masks],
                group_time=float(h["group_time_header_s"]), frame_time=float(h["frame_time_s"]), nframes=2,
                dark_dn_s=h["average_dark_current_DN_s"][:, 24:41, 24:41][masks],
                zeroframe=h["native_corrected_zeroframe_DN"][:, 24:41, 24:41][masks],
                pixeldq=h["pixeldq"][:, 24:41, 24:41][masks], cast_input=True)
            actual = h["native_rate_pre_gain_scale_DN_s"][:, 24:41, 24:41][masks]
            native_match = bool(np.array_equal(values.astype(np.float32), actual.astype(np.float32)))
            if not native_match:
                raise ValueError("Source-derived native OLS values do not reproduce originalC float32 on consumed pixels")
            jac, zerojac = analytic[:9].T, analytic[9]
        elif derivative_map == "h1":
            jac = finite_jac
            zerojac = h["native_rate_jacobian_group_and_zeroframe_h1_per_s"][:, 9, 24:41, 24:41][masks]
        else:
            raise ValueError("Explicit native derivative source required")
        if np.any(zerojac != 0):
            raise ValueError("Actual selected native branch consumes ZEROFRAME; unsupported by9group calibration projection")
        group_Lprime = np.moveaxis(h["linearity_polynomial_first_derivative"][:, :, 24:41, 24:41], 1, -1).copy()
        no_lin = (d["linearity_DQ_mapped"][:, 24:41, 24:41] & 1048576) != 0
        group_Lprime[no_lin] = 1.
        flat = h["flat"][:, 24:41, 24:41][masks]
        phot = float(h["photmjsr"])
        gain_scale = (h["native_rate_DN_s"][:, 24:41, 24:41]/h["native_rate_pre_gain_scale_DN_s"][:, 24:41, 24:41])[masks]
        if not np.all(np.isfinite(gain_scale)):
            raise ValueError("Native gain_scale derivative is undefined")
        weights = jac*group_Lprime[masks]*(phot*gain_scale/flat)[:, None]
        if not np.isfinite(weights).all():
            raise ValueError("Nonfinite actual native/linearity projection derivative")
        requests.append({"amp": amp, "clock_xy": clock, "group_weights": weights})
        residual = np.concatenate((target["raw_residual"][t], field["raw_residual"][f]), axis=0)[masks]
        native_error = np.concatenate((target["ERR"][t], field["ERR"][f]), axis=0)[masks]
        # Existing native Poisson branch is retained as a declared diagonal
        # comparator until matching IPC/charge counting propagation is supplied.
        poisson = h["native_rate_VAR_POISSON"][:, 24:41, 24:41][masks]*(phot*gain_scale/flat)**2
        rate = h["native_rate_DN_s"][:, 24:41, 24:41][masks]
        flat_error = h["flat_ERR"][:, 24:41, 24:41][masks]
        flat_marginal = (rate*phot*flat_error/flat**2)**2
        if np.any(poisson < 0) or np.any(flat_marginal < 0) or not np.isfinite(poisson+flat_marginal).all():
            raise ValueError("Native photon/reference variance branch is invalid")
        latest = late["post_SIRS_matched_patch_difference_DN"][:, exposure-1]
        latest_weighted = np.sum(np.moveaxis(latest, 0, -1)[masks]*weights, axis=1)
        count_target = int(masks[:4].sum())
        strict = None
        if exposure == 4:
            ti = np.flatnonzero(np.char.startswith(heldtarget["name"].astype(str), stem))
            fi = np.flatnonzero(np.char.startswith(trainingfield["name"].astype(str), stem))
            if not np.array_equal(heldtarget["mask"][ti], masks[:4]) or not np.array_equal(trainingfield["mask"][fi], masks[4:]):
                raise ValueError("Held/training fixed-mean mask differs from calibrated projection order")
            strict = np.concatenate((heldtarget["raw_residual"][ti], trainingfield["raw_residual"][fi]), axis=0)[masks]
        records.append({"stem": stem, "residual": residual, "native_error": native_error,
            "native_xy": native_xy, "masks": masks, "poisson": poisson, "flat_marginal": flat_marginal,
            "late_dark_draw": latest_weighted, "strict_held_residual": strict,
            "count_target": count_target, "native_input_sha256": sha(historypath),
            "native_DQ_sidecar_sha256": sha(dqpath),
            "derivative_map": derivative_map, "analytic_native_value_matches_originalC_float32": native_match,
            "analytic_branch_nondifferentiable_pixels": sum(not b["locally_unique_derivative"] for b in branches),
            "analytic_vs_h1_relative_Frobenius_difference": float(np.linalg.norm(jac-finite_jac)/np.linalg.norm(jac)),
            "native_J_h1_h2_max_difference_per_s": float(np.max(abs(np.moveaxis(
                h["native_rate_jacobian_group_and_zeroframe_h1_per_s"][:, :9, 24:41, 24:41]-
                h["native_rate_jacobian_group_and_zeroframe_h2_per_s"][:, :9, 24:41, 24:41], 1, -1)[masks])))})
    return requests, records, {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (
        static/"all_target.npz", static/"all_field.npz", static/"held_target_strict.npz", static/"training_field.npz")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute-fixed-calibration-projection", action="store_true")
    parser.add_argument("--native-derivative-map", choices=("analytic", "h1"), default="analytic")
    args = parser.parse_args()
    if not args.execute_fixed_calibration_projection:
        parser.error("Explicit source-fixed calibration projection flag required")
    started = time.perf_counter()
    model, spectral = spectral_model()
    requests, records, means = requests_and_saved_residuals(args.native_derivative_map)
    data_output = DATA/"analytic_native_projection" if args.native_derivative_map == "analytic" else DATA
    receipt_output = OUT/"analytic_native_projection" if args.native_derivative_map == "analytic" else OUT
    data_output.mkdir(exist_ok=True); receipt_output.mkdir(exist_ok=True)
    matrices = model.project(requests, progress=lambda n: print("SOURCE_CROSS_KERNELS", n, flush=True))
    outputs = []
    for request, record, electronic in zip(requests, records, matrices):
        covariance = electronic+np.diag(record["poisson"]+record["flat_marginal"])
        result, white = score(record["residual"], covariance)
        original, _ = score(record["residual"], np.diag(record["native_error"]**2))
        calibration_held, _ = score(record["late_dark_draw"], electronic)
        n = record["count_target"]
        target_result, _ = score(record["residual"][:n], covariance[:n, :n])
        field_result, _ = score(record["residual"][n:], covariance[n:, n:])
        strict_result = None
        if record["strict_held_residual"] is not None:
            strict_result, _ = score(record["strict_held_residual"][:n], covariance[:n, :n])
        path = data_output/(record["stem"]+"_spectral_rate_covariance.npz")
        np.savez_compressed(path, native_xy=record["native_xy"], valid_mask=record["masks"],
            output_order=np.array(["A", "B", "C", "D", "FIELD"]), count_target=n,
            native_ADC_group_to_SCI_rate_weights=request["group_weights"],
            source_READ_DARK_covariance_SCI2=electronic,
            retained_native_diagonal_Poisson_SCI2=record["poisson"],
            conditional_diagonal_flat_reference_SCI2=record["flat_marginal"],
            conditional_total_covariance_SCI2=covariance,
            unchanged_static_raw_residual_SCI=record["residual"], unchanged_native_ERR=record["native_error"],
            conditional_whitened_static_residual=white,
            latest_calibration_epoch_projected_draw_SCI=record["late_dark_draw"])
        outputs.append({"exposure": record["stem"], "file": str(path.relative_to(ROOT)).replace('\\', '/'),
            "sha256": sha(path), "n": len(electronic), "count_target": n,
            "source_native_history_sha256": record["native_input_sha256"],
            "native_DQ_sidecar_sha256": record["native_DQ_sidecar_sha256"],
            "native_J_h1_h2_max_difference_per_s": record["native_J_h1_h2_max_difference_per_s"],
            "native_derivative_map": record["derivative_map"],
            "analytic_native_value_matches_originalC_float32": record["analytic_native_value_matches_originalC_float32"],
            "analytic_branch_nondifferentiable_pixels": record["analytic_branch_nondifferentiable_pixels"],
            "analytic_vs_h1_relative_Frobenius_difference": record["analytic_vs_h1_relative_Frobenius_difference"],
            "raw_original_native_ERR": original, "conditional_full_covariance": result,
            "conditional_target_only": target_result, "conditional_field_only": field_result,
            "strict_held_exposure_fixed_mean_target": strict_result,
            "latest_dark_epoch_fixed_earliest_spectrum": calibration_held,
            "max_abs_target_field_electronic_cross_covariance_SCI2": float(abs(electronic[:n, n:]).max()),
            "source_READ_DARK_variance_SCI2_range": [float(np.diag(electronic).min()), float(np.diag(electronic).max())],
            "source_complete_instrument_covariance": False, "physical_validation": False})
        print("FIXED_COVARIANCE_APPLICATION", record["stem"], result["Q"], original["Q"], flush=True)
    final = {"status": "actual source-calibrated full cross-stamp/amp READ+DARK covariance projected and saved fixed means whitened; no mean or noise refit",
        "spectral_model": spectral, "exposures": outputs,
        "source_mean_SHA256": means, "seconds": time.perf_counter()-started,
        "projection": f"native OLS_C {args.native_derivative_map} Jacobian × native group Lprime × gain_scale ×PHOTMJSR/flat; ZEROFRAME derivative proven zero on consumed masks",
        "components": "source-calibrated postSIRS READ+DARK REPLACES native diagonal VAR_RNOISE; retained native diagonal VAR_POISSON plus explicit diagonal flat-ERR prior; no original ERR added again",
        "cross_correlations": "within-exposure target A/B/C/D/FIELD plus all four read-clock amplifier channels share source spectra; not five independent noise copies",
        "conditional_limits": ["analytic map is the continuous extension on actual fixed DQ/SNR/median branches, with originalC float32 values proved; boundary uniqueness reported. h1 option is a finite-step diagnostic that may cross branches, not a demonstrated local Jacobian",
            "stationary periodic log-frequency-bin calibration model; masked spectral window not deconvolved; one difference realization per epoch; calibration uncertainty/state dependence not known exactly",
            "latest dark epoch is an internal calibration holdout already summarized; no new science validation is claimed",
            "native diagonal Poisson is a producing-branch comparator, not the complete collection/IPC/fullwell statistical law; matching photon IPC projection is a separate incoming application",
            "measured READ+DARK includes physical dark shot noise; retaining nativeVAR_POISSON may duplicate its small dark component unless photon-only/dark counting rates are separately identified",
            "flat ERR supplies marginal uncertainty only; full shared spatial flat/kernel calibration covariance absent",
            "saved static mean is unchanged; lower Q under calibrated noise is not physical validation or a correction to failed mean residuals",
            "these scores are conditional fixed-mean Gaussian noise diagnostics; no BHSM optical-state Wick/Gaussian closure or prepared U,V solve is implied"],
        "source_sha256": {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in (
            Path(__file__), ROOT/"code/strong_lensing_dark_covariance.py", ROOT/"code/strong_lensing_native_ols.py")},
        "physical_guard": 2e-10, "physical_validation": False,
        "reserved_science_observations_unused": True, "noise_rescaling": False}
    dump(receipt_output/"fixed_spectral_projection.json", final)


if __name__ == "__main__":
    main()
