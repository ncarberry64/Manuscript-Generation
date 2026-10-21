"""Verify arithmetic and read support from the small, saved field-only archive.

This does not refit any PSF, target, lens parameter, or noise scale. It checks
the decisive quantities without requiring the multi-gigabyte references.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

OUT = Path(__file__).resolve().parent


def main():
    stamp_path = OUT / "native_corrected_field_ramp.npz"
    stamp_receipt = json.loads((OUT / "native_corrected_field_stamp_receipt.json").read_text())
    assert hashlib.sha256(stamp_path.read_bytes()).hexdigest() == stamp_receipt["sha256"]
    stamp = np.load(stamp_path)
    results = json.loads((OUT / "common_time_field_stcal118_results.json").read_text())
    receipts = json.loads((OUT / "common_time_field_stcal118_fit_receipts.json").read_text())
    arrays = np.load(OUT / "common_time_field_stcal118_arrays.npz")
    assert arrays["common"].shape == (17, 17)
    assert arrays["common"].all()
    x = stamp["x_pixel"]
    y = stamp["y_pixel"]
    choose_stamp = (x >= 892) & (x < 909) & (y >= 868) & (y < 885)
    assert choose_stamp.sum() == 289
    groupdq = stamp["native_groupdq"][:, choose_stamp]
    # DNU and SAT exclude a group. JUMP begins a new ramp segment and is
    # not automatically excluded from this support count.
    usable_groups = ((groupdq & (1 | 2)) == 0).sum(axis=0)
    centre = np.array([899.7549609755204, 875.9549994961739])
    radius = np.hypot(x[choose_stamp] - centre[0], y[choose_stamp] - centre[1])
    support = []
    edges = [0, 1, 2, 3, 4, 6, 8, 12]
    for low, high in zip(edges[:-1], edges[1:]):
        selection = (radius >= low) & (radius < high)
        values, counts = np.unique(usable_groups[selection], return_counts=True)
        support.append({"radius_pixel": [low, high], "pixels": int(selection.sum()),
                        "count_of_pixels_by_usable_groups": {
                            str(int(v)): int(n) for v, n in zip(values, counts)}})
    assert support[0]["count_of_pixels_by_usable_groups"] == {"2": 3}
    assert all(row["count_of_pixels_by_usable_groups"] == {"9": row["pixels"]}
               for row in support[3:6])

    diagnostics = {}
    for name, key in [("full", "original_native_CAL_same_support"),
                      ("early", "first_two_reads_same_support")]:
        fit = receipts[name]
        xx, yy = np.asarray(fit["pixel_x"]), np.asarray(fit["pixel_y"])
        data, model = np.asarray(fit["pixel_data"]), np.asarray(fit["pixel_model"])
        background = np.asarray(fit["background_plane_SCI"])
        local_centre = centre - np.array([892., 868.])
        plane = background[0] + background[1] * (xx - local_centre[0]) + background[2] * (yy - local_centre[1])
        source = model - plane
        value = float(np.linalg.norm(data - model) / np.linalg.norm(source))
        np.testing.assert_allclose(value, results[key]["residual_L2_over_fitted_source_L2"], rtol=1e-13)
        assert fit["valid_measurement_fit"] and fit["optimizer_success"]
        assert fit["usable_pixels"] == 289 and not fit["active_bound_parameter_indices"]
        diagnostics[name] = {"normalized_residual_L2": value,
                             "chi2_per_dof": float(fit["chi2"] / fit["dof"]),
                             "p_value": fit["p_value"]}
    assert results["native_pipeline_replay"]["numerical_replay_gate"]["pass"]
    assert not results["physical_validation"]
    assert diagnostics["early"]["p_value"] < 0.01
    np.testing.assert_array_equal(receipts["full"]["pixel_x"], receipts["early"]["pixel_x"])
    np.testing.assert_array_equal(receipts["full"]["pixel_y"], receipts["early"]["pixel_y"])
    summary = {"status": "SAVED_FIELD_ARITHMETIC_AND_SUPPORT_AUDIT_PASS",
               "not_an_independent_observation": True,
               "field_stamp_sha256": stamp_receipt["sha256"],
               "group_exclusions": "GROUPDQ DO_NOT_USE|SATURATED; JUMP starts a segment",
               "group_support_by_radius": support,
               "core_groupdq_sequences": groupdq[:, radius < 1].T.astype(int).tolist(),
               "matched_pixels": 289, "diagnostics": diagnostics,
               "normalized_L2_reduction_fraction": 1 - diagnostics["early"]["normalized_residual_L2"] / diagnostics["full"]["normalized_residual_L2"],
               "early_p_passes_0_01": False,
               "precision_validation": False}
    (OUT / "compact_field_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
