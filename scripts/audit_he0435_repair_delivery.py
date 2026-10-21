"""Read-only verification of the delivered repair and compact read-history data."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(delivery):
    delivery = Path(delivery)
    manifest = json.loads((delivery / "DELIVERY_MANIFEST.json").read_text(encoding="utf-8"))
    verified = []
    for item in manifest["files"]:
        path = delivery / item["path"]
        if path.stat().st_size != item["size_bytes"] or digest(path) != item["sha256"]:
            raise ValueError("Delivery bytes changed: " + item["path"])
        verified.append(item["path"])
    native = delivery / "ramp_diagnostic/native_pipeline"
    with np.load(native / "native_corrected_field_ramp.npz", allow_pickle=False) as ramp:
        x, y = ramp["x_pixel"], ramp["y_pixel"]
        selected = (x >= 892) & (x < 909) & (y >= 868) & (y < 885)
        if selected.sum() != 289:
            raise ValueError("Original matched support is incomplete")
        flags = ramp["native_groupdq"][:, selected]
        support = ((flags & 3) == 0).sum(axis=0)
        radius = np.hypot(x[selected] - 899.7549609755204, y[selected] - 875.9549994961739)
        support_rows = []
        for lower, upper in zip([0, 1, 2, 3, 4, 6, 8], [1, 2, 3, 4, 6, 8, 12]):
            use = (radius >= lower) & (radius < upper)
            values, counts = np.unique(support[use], return_counts=True)
            support_rows.append(dict(radius_pixel=[lower, upper], pixels=int(use.sum()),
                                     usable_groups={str(v): int(n) for v, n in zip(values, counts)}))
        if support_rows[0]["usable_groups"] != {"2": 3}:
            raise ValueError("Delivered central group support differs")
        if any(row["usable_groups"] != {"9": row["pixels"]} for row in support_rows[3:6]):
            raise ValueError("Delivered wing group support differs")
        times = ramp["native_read_times_s"].tolist()
    result = json.loads((native / "common_time_field_stcal118_results.json").read_text(encoding="utf-8"))
    fits = json.loads((native / "common_time_field_stcal118_fit_receipts.json").read_text(encoding="utf-8"))
    diagnostics = {}
    for label, key in [("full", "original_native_CAL_same_support"),
                       ("early", "first_two_reads_same_support")]:
        fit = fits[label]
        xx, yy = np.asarray(fit["pixel_x"]), np.asarray(fit["pixel_y"])
        bg = np.asarray(fit["background_plane_SCI"])
        plane = bg[0] + bg[1] * (xx - (899.7549609755204 - 892)) + bg[2] * (yy - (875.9549994961739 - 868))
        model, data = np.asarray(fit["pixel_model"]), np.asarray(fit["pixel_data"])
        normalized = float(np.linalg.norm(data - model) / np.linalg.norm(model - plane))
        np.testing.assert_allclose(normalized, result[key]["residual_L2_over_fitted_source_L2"], rtol=1e-13)
        diagnostics[label] = dict(normalized_residual_L2=normalized, chi2=float(fit["chi2"]),
                                  dof=int(fit["dof"]), p=float(fit["p_value"]))
    if not result["native_pipeline_replay"]["numerical_replay_gate"]["pass"]:
        raise ValueError("Delivered matching-version numerical gate did not pass")
    if diagnostics["early"]["p"] >= .01 or result["physical_validation"]:
        raise ValueError("Original failed adequacy/admission state changed")
    return dict(status="PASS_DELIVERY_AND_COMPACT_ARITHMETIC_ONLY", files_verified=len(verified),
                delivery_manifest_sha256=digest(delivery / "DELIVERY_MANIFEST.json"),
                matched_pixels=289, read_times_s=times, group_support=support_rows, diagnostics=diagnostics,
                normalized_L2_reduction_fraction=1 - diagnostics["early"]["normalized_residual_L2"] / diagnostics["full"]["normalized_residual_L2"],
                detector_pipeline_rerun=False, new_independent_observations=False,
                physical_admission=False, physical_guard=2e-10,
                limitations=["Arithmetic verification does not repeat or independently validate native detector processing.",
                             "The matching stcal replay establishes numerical compatibility, not an unrecorded archive installation version.",
                             "Shared SIRS, dark/reference and cross-read covariance remain incomplete."])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.resolve().is_relative_to(args.delivery.resolve()):
        raise ValueError("Do not write inside the immutable delivered package")
    result = audit(args.delivery)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ["status", "files_verified", "normalized_L2_reduction_fraction"]}))


if __name__ == "__main__":
    main()
