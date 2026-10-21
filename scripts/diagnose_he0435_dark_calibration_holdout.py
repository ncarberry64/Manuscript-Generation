"""Apply the SAME predeclared calibration screen to the later dark epoch.

All-pixel held-dark scores are preserved in fixed_spectral_projection.json.
This separate conditional bulk-noise diagnostic removes only the held dark's
previously frozen DQ/unsaturated/500DN screen; no spectrum or mean is refitted.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from scipy import linalg

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/"data/strong_lensing_response_completion_20261021/dark_calibration"
OUT = ROOT/"artifacts/strong_lensing_response_completion_20261021/dark_calibration"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--native-derivative-map", choices=("analytic", "h1"), default="analytic")
    args = parser.parse_args()
    matrix_directory = DATA/"analytic_native_projection" if args.native_derivative_map == "analytic" else DATA
    receipt_directory = OUT/"analytic_native_projection" if args.native_derivative_map == "analytic" else OUT
    late = np.load(DATA/"epoch1_measured_covariance.npz")
    rows = []
    for exposure in range(1, 5):
        path = matrix_directory/f"jw01198007001_02101_0000{exposure}_nrcb3_spectral_rate_covariance.npz"
        with np.load(path) as a:
            good = late["matched_patch_valid_DQ"][exposure-1][a["valid_mask"]]
            c = a["source_READ_DARK_covariance_SCI2"][np.ix_(good, good)]
            r = a["latest_calibration_epoch_projected_draw_SCI"][good]
            z = linalg.solve_triangular(linalg.cholesky(c, lower=True), r, lower=True)
            rows.append({"exposure": exposure, "original_native_outputs": len(good),
                "same_fixed_calibration_screen_retained": int(good.sum()),
                "excluded_dark_particle_or_DQ": int((~good).sum()),
                "conditional_calibration_held_Q": float(z@z),
                "Q_per_observation": float(z@z/len(z)),
                "projected_covariance_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    report = {"status": "same previously fixed500DN/DQ calibration screen applied to held dark; no spectrum adjustment",
        "records": rows, "allnative_unfiltered_held_scores_preserved": "fixed_spectral_projection.json",
        "native_derivative_map": args.native_derivative_map,
        "interpretation": "Unscreened dark particle events are outside the calibration bulk-noise conditioning; same-screen transfer still varies by detector location/epoch and is not assumed universally calibrated.",
        "mask_not_selected_by_PSFor_lens_agreement": True, "physical_validation": False,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (receipt_directory/"latest_dark_same_screen_diagnostic.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n",
                                                             encoding="utf-8", newline="\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
