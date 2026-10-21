"""Plot the once-executed JWST diagnostic from saved JSON only; never fit pixels.

The failed independent PSF and exposure-consistency gates are retained. The
reported q offsets are descriptive, conditional on the inadequate fixed PSF.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ART = ROOT / "artifacts/strong_lensing_geometry_20261021"
FIG = ROOT / "manuscript/figures/he0435_jwst_geometry_diagnostic"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def residual_summary(point: dict) -> dict:
    """Descriptive summaries of already-saved, known-ERR residuals."""
    z = np.asarray(point["standardized_residual"], dtype=float)
    residual = np.asarray(point["pixel_residual_data_minus_model"], dtype=float)
    error = np.asarray(point["pixel_error"], dtype=float)
    x, y = np.asarray(point["pixel_x"]), np.asarray(point["pixel_y"])
    cx, cy = point["centroid_xy_pixel"]
    radius = np.hypot(x - cx, y - cy)
    # Reproduce the saved original residual statistic; no new likelihood fit.
    np.testing.assert_allclose(residual / error, z, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(np.dot(z, z), point["chi2"], rtol=1e-10)
    annuli = []
    for lo, hi in [(0, 2), (2, 4), (4, 6), (6, 12)]:
        use = (radius >= lo) & (radius < hi)
        annuli.append({
            "radius_native_pixel": [lo, hi],
            "usable_pixels": int(use.sum()),
            "mean_standardized_residual": float(z[use].mean()) if use.any() else None,
            "mean_SCI_residual": float(residual[use].mean()) if use.any() else None,
            "chi2_fraction": float(np.dot(z[use], z[use]) / np.dot(z, z)),
        })
    return {
        "label": point["label"], "chi2": point["chi2"], "dof": point["dof"],
        "chi2_per_dof": point["chi2"] / point["dof"],
        "p_value": point["p_value"], "usable_pixels": len(z),
        "standardized_residual_rms": float(np.sqrt(np.mean(z * z))),
        "fraction_abs_standardized_residual_gt_3": float(np.mean(abs(z) > 3)),
        "fraction_abs_standardized_residual_gt_5": float(np.mean(abs(z) > 5)),
        "maximum_abs_standardized_residual": float(max(abs(z))),
        "annuli": annuli,
        "pixel_adequacy_pass": point["pixel_model_adequate_at_frozen_threshold"],
    }


def main() -> None:
    aggregate_path = ART / "jwst_target_diagnostic/pixel_geometry_results.json"
    field_path = ART / "standard_psf/field_psf_calibration.json"
    aggregate, field = read_json(aggregate_path), read_json(field_path)
    exposure_paths = sorted((ART / "jwst_target_diagnostic/exposures").glob("*.json"))
    exposures = [read_json(p) for p in exposure_paths]
    assert len(exposures) == 4 and len(field["field_fits"]) == 4
    assert aggregate["diagnostic_mode"] and not aggregate["calibration_adequacy_pass"]
    assert not aggregate["all_pixel_models_adequate"]
    assert not aggregate["exposures_consistent"]
    assert not aggregate["physical_admission"]
    assert aggregate["unchanged_physical_guard"] == 2e-10
    assert all(not m["conditional_geometry_validation_pass"] for m in aggregate["models"].values())
    assert all(not m["conditional_geometry_validation_pass"]
               for e in exposures for m in e["model_scores"].values())
    summaries = []
    for path, exposure in zip(exposure_paths, exposures):
        assert [p["label"] for p in exposure["points"]] == list("ABCD")
        summaries.append({"exposure": path.stem, "points": [residual_summary(p) for p in exposure["points"]]})
    assert all(not p["pixel_adequacy_pass"] for e in summaries for p in e["points"])

    q_mas = np.asarray([e["intrinsic_q_arcsec"] for e in exposures]) * 1000
    models = aggregate["models"]
    model_differences = {}
    for old, new, label in [("old_frozen_conventional", "old_frozen_R1", "old"),
                            ("revised_G1_conventional", "revised_G1_R1", "G1")]:
        a = np.asarray(models[old]["predicted_intrinsic_q_arcsec"])
        b = np.asarray(models[new]["predicted_intrinsic_q_arcsec"])
        model_differences[label] = float(np.max(abs(a - b)) * 1e6)
    source_paths = [aggregate_path, field_path, *exposure_paths]
    findings = {
        "status": "DIAGNOSTIC_ONLY_NO_VALIDATION",
        "inputs_sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): digest(p) for p in source_paths},
        "analysis": "Saved-result summaries only; no FITS read, parameter estimation, response fitting or covariance rescaling.",
        "freeze_commit": aggregate["freeze_commit"],
        "q_definition": "u=(B-A)/|B-A|, v=J u; q=(|B-A|, (C-A).u, (C-A).v, (D-A).u, (D-A).v); translation/proper rotation removed, scale and chirality retained.",
        "model_scores": models,
        "maximum_conventional_R1_prediction_difference_microarcsec": model_differences,
        "q_exposure_span_mas": np.ptp(q_mas, axis=0).tolist(),
        "q_exposure_sample_sd_mas": np.std(q_mas, axis=0, ddof=1).tolist(),
        "q_exposure_consistency_chi2": aggregate["combined"]["exposure_consistency_chi2"],
        "q_exposure_consistency_dof": aggregate["combined"]["exposure_consistency_dof"],
        "q_exposure_consistency_p_value": aggregate["combined"]["exposure_consistency_p_value"],
        "target_residuals": summaries,
        "field_control_fit_statistics": [{k: f[k] for k in ["exposure", "chi2", "dof", "p_value", "absolute_adequacy_pass"]} for f in field["field_fits"]],
        "interpretation": "G1 reduces the raw five-coordinate offset descriptively. Coherent positive wing residuals and exposure scatter demonstrate inadequate measurement response. Neither model is validated or physically rejected by these conditional formal scores; no BHSM-specific signal is established.",
        "remaining_dependency": "A demonstrably adequate instrument/scene mean response, independently calibrated and with its uncertainty characterized, is required before calibrated geometry inference from this dataset. Finite PSF uncertainty is not supplied by native ERR. These target photons have now been consumed for diagnosis and are not unused future validation data.",
        "physical_admission": False, "unchanged_physical_guard": 2e-10,
    }
    (ART / "target_response_findings.json").write_text(json.dumps(findings, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    plt.rcParams.update({"font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
                         "legend.fontsize": 8, "pdf.fonttype": 42, "ps.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.35))
    fig.subplots_adjust(left=.075, right=.985, bottom=.23, top=.77, wspace=.29)
    fig.suptitle("HE 0435−1223: JWST F115W diagnostic — no validation", fontsize=13, y=.97)
    ax = axes[0]
    x = np.arange(5)
    specs = [("old_frozen_conventional", "#286a9d", "o", "-", "Old conventional"),
             ("old_frozen_R1", "#286a9d", "s", "--", "Old R1"),
             ("revised_G1_conventional", "#b44c19", "o", "-", "G1 conventional"),
             ("revised_G1_R1", "#b44c19", "s", "--", "G1 R1")]
    for name, color, marker, line, label in specs:
        ax.plot(x, models[name]["intrinsic_residual_mas"], color=color, marker=marker,
                linestyle=line, linewidth=1.25, markersize=4 if marker == "o" else 6,
                markerfacecolor=color if marker == "o" else "none", label=label)
    ax.axhline(0, color=".55", linewidth=.8)
    ax.set_xticks(x, [r"$|B-A|$", r"$C_u$", r"$C_v$", r"$D_u$", r"$D_v$"])
    ax.set_ylabel("Observed − predicted q (mas)")
    ax.set_xlabel("Intrinsic image geometry")
    ax.set_title("(a) Raw offsets: old vs two-plane G1", loc="left", pad=44)
    ax.legend(ncol=2, loc="lower left", bbox_to_anchor=(-.015, 1.005), borderaxespad=0,
              frameon=False, handlelength=2.1, columnspacing=1.1)
    ax.text(.02, .96, "Coordinate RMS: 3.31 → 0.78 mas\nConventional and R1 nearly coincide", transform=ax.transAxes,
            fontsize=8, va="top", bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    ax.set_ylim(-5.45, 6.25)
    ax = axes[1]
    colors = ["#286a9d", "#b44c19", "#22845b", "#9163a4", "#333333"]
    offsets = np.linspace(-.23, .23, 5)
    for j, label in enumerate(["A", "B", "C", "D", "Field star"]):
        values = [e["points"][j]["chi2"] / e["points"][j]["dof"] for e in exposures] if j < 4 else [f["chi2"] / f["dof"] for f in field["field_fits"]]
        ax.scatter(np.arange(1, 5) + offsets[j], values, color=colors[j], s=25,
                   marker="D" if j == 4 else "o", label=label, zorder=3)
    ax.axhline(1, color=".4", linestyle=":", linewidth=1, label=r"$E[\chi^2/\nu]=1$")
    ax.set_xticks(range(1, 5), ["1", "2", "3", "4"])
    ax.set_xlabel("Exposure")
    ax.set_ylabel(r"Native-pixel $\chi^2/\nu$")
    ax.set_ylim(-2, 98)
    ax.set_title("(b) Independent PSF fails field and target", loc="left", pad=44)
    ax.legend(ncol=3, loc="lower left", bbox_to_anchor=(-.02, 1.005), borderaxespad=0,
              frameon=False, handlelength=1.3, columnspacing=.65, handletextpad=.4)
    ax.text(.02, .96, "All 20 fits fail fixed-ERR adequacy", transform=ax.transAxes, va="top", fontsize=8)
    fig.text(.075, .095, "Independent PSF calibration and exposure-consistency gates both fail; every validation flag is false.", color="#a02e2e", fontsize=9, fontweight="bold")
    fig.text(.075, .045, r"Conditional known-ERR statistics; no covariance inflation. Exposure consistency: $Q=432.3$ for 15 dof. No response refit.", fontsize=8)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG.with_suffix(".pdf"), metadata={"Title": "HE 0435 JWST geometry diagnostic: no validation", "Creator": "plot_he0435_jwst_geometry.py; saved JSON only"})
    fig.savefig(FIG.with_suffix(".png"), dpi=220)
    plt.close(fig)
    for path in [FIG.with_suffix(".pdf"), FIG.with_suffix(".png"), ART / "target_response_findings.json"]:
        print(path.relative_to(ROOT), digest(path))


if __name__ == "__main__":
    main()
