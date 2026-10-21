"""Apply the native averaged-read map to actual consumed target-history anchors.

This is a local, data-anchored calibration diagnostic, not a fit of the physical
illumination or of segment intercepts. It does not turn archived rate estimates
into independent training data and does not apply a saturated polynomial branch.
"""
from pathlib import Path
import sys, json, hashlib
import numpy as np
from stdatamodels.jwst.datamodels import dqflags

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_group_calibration import averaged_read_observation, polynomial_derivatives


def main():
    source = ROOT/"data/strong_lensing_response_completion_20261021/histories"
    out = ROOT/"artifacts/strong_lensing_response_completion_20261021/group_average"
    out.mkdir(exist_ok=True)
    reports = []
    for path in sorted(source.glob("*_histories.npz")):
        a = np.load(path, allow_pickle=False)
        time = a["native_read_times_s"]
        dark_indices = a["native_dark_RAPID_frame_indices"]
        rows = []; maps = np.full_like(a["corrected_group_DN"], np.nan, dtype=float)
        for k, label in enumerate(a["labels"]):
            coeff = a["linearity_coefficients"][k].astype(float)
            native_adc = a["pre_linearity_group_DN"][k].astype(float)
            rate = a["native_rate_DN_s"][k].astype(float)
            yy, xx = np.indices(rate.shape)
            gx, gy = a["prior_guess_native_xy"][k]-a["bounds_xyxy"][k, :2]
            roi = ((np.abs(xx-np.rint(gx)) <= 8) & (np.abs(yy-np.rint(gy)) <= 8))
            valid = (((a["groupdq"][k] & 7) == 0)
                     & ((a["pixeldq"][k] & (1 | dqflags.pixel["NO_LIN_CORR"])) == 0)[None]
                     & np.isfinite(rate)[None] & roi[None])
            # A native skipped/bad polynomial is not replaced by an extrapolation.
            finite_coeff = np.all(np.isfinite(coeff), axis=0)
            valid &= finite_coeff[None]
            group, y, x = np.nonzero(valid)
            c = coeff[:, y, x]
            anchor = native_adc[group, y, x]
            center_q = polynomial_derivatives(c, anchor)[0]
            dark = a["dark_RAPID42_DN"][k][dark_indices[group], y[:, None], x[:, None]]
            delta_dark = dark-dark.mean(axis=1, keepdims=True)
            reads = (center_q[:, None]+rate[y, x, None]*(time[group]-time[group].mean(axis=1, keepdims=True))
                     + delta_dark)
            answer = averaged_read_observation(reads.T[None], c,
                                                np.broadcast_to(anchor[None, None], (1, 2, len(anchor))))
            correction = answer["mean_DN"][0]-center_q
            adc_read = answer["inverse_ADC_reads_DN"][0].T
            # Native saturation is evaluated before bias/refpix subtraction.
            ref_correction = a["pre_refpix_group_DN"][k][group, y, x]-native_adc[group, y, x]
            adc_ceiling = (a["saturation_threshold_DN"][k][y, x]
                           -a["superbias_DN"][k][y, x]-ref_correction)
            in_domain = np.all(adc_read < adc_ceiling[:, None], axis=1)
            maps[k, group[in_domain], y[in_domain], x[in_domain]] = correction[in_domain]
            safe = correction[in_domain]
            illumination = rate[y[in_domain], x[in_domain]]*time[group[in_domain]].mean(axis=1)
            radius = np.hypot(x[in_domain]-gx, y[in_domain]-gy)
            annuli = []
            for lo, hi in [(0, 2), (2, 3), (3, 6), (6, 12)]:
                choose = (radius >= lo) & (radius < hi)
                annuli.append({"radius_pixel": [lo, hi], "group_samples": int(choose.sum()),
                               "max_abs_group_bias_DN": float(np.max(np.abs(safe[choose]))) if choose.any() else None,
                               "bias_L2_over_illumination_L2": float(np.linalg.norm(safe[choose])/max(np.linalg.norm(illumination[choose]), 1e-30))})
            rows.append({"label": str(label), "used_group_samples": int(in_domain.sum()),
                         "individual_read_calibration_limit_exceeded": int((~in_domain).sum()),
                         "excluded_jump_groups": int(np.sum(((a["groupdq"][k] & 4) != 0) & roi[None])),
                         "max_abs_group_bias_DN": float(np.max(np.abs(safe))),
                         "median_abs_group_bias_DN": float(np.median(np.abs(safe))),
                         "annuli": annuli})
        array_path = out/(path.stem+"_group_bias.npz")
        np.savez_compressed(array_path, correction_DN=maps, labels=a["labels"], bounds_xyxy=a["bounds_xyxy"])
        reports.append({"source": str(path.relative_to(ROOT)),
                        "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                        "array": str(array_path.relative_to(ROOT)), "patches": rows})
    receipt = {"status": "EXECUTED_NATIVE_CALIBRATION_MAP_DIAGNOSTIC_ONLY",
               "equation": "q_jr=L(native_ADC_group_j)+f_native*(t_jr-mean(t_j))+dark_jr-mean(dark_jr); delta_j=L(mean_r L^-1(q_jr))-L(native_ADC_group_j)",
               "meaning": "Actual native pre-linearity anchors and consumed rate pilots quantify the local read-averaging mean effect. Group intercept/state is anchored separately at each observation, not solved as a physical preparation or a joint likelihood. This is not an independently calibrated rate correction.",
               "selection": "Original DNU/SAT/JUMP groups and science-pixel DNU/NO_LIN_CORR flags excluded; raw reference DQ encoding is not used as global flags. Only original17x17 ROI. Individual inverse read exceeding the native saturation calibration threshold is excluded from this diagnostic rather than extrapolated.",
               "unresolved": "Actual readwise charge and jump timing, segment intercept uncertainty and charge-dependent collection law still belong in the coupled fit. No optical response, physical lens prediction or ERR is changed here.",
               "applications": reports, "physical_admission": False, "physical_guard": 2e-10,
               "no_reserved_observations": True,
               "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__), ROOT/"code/strong_lensing_group_calibration.py"]}}
    (out/"diagnostic.json").write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"exposures": len(reports), "patches": sum(len(r["patches"]) for r in reports),
                      "status": receipt["status"]}))


if __name__ == "__main__":
    main()
