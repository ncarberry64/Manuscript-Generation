"""Replay the mask correction with an explicit, repository-verified baseline.

This additive integration runner preserves the delivered runner and manifest.
Its replay operations are inherited unchanged; only provenance validation is
made explicit to correct the delivered legacy-source hash inconsistency.

No download, physical refit, uncertainty rescaling, or new-observation validation.
This runner requires the archived native stamps, exact GWCS and original headers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess

import asdf
from astropy.io import fits
import numpy as np
from scipy.linalg import block_diag

CODE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CODE_ROOT / "code"))


SUPPLIED_BASELINE_SHA256 = "08ff73d45de3e27746648e58c83cbfe620d39d9d7e0b36c8428242255a9b32d3"
REVIEWED_COMMIT = "7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9"
REVIEWED_FITTER = "code/strong_lensing_pixel_geometry.py"


def verified_baseline(manifest_path, legacy_path, repository_root):
    """Verify the unchanged delivered baseline or an explicit Git-bound correction.

    Alternate manifests cannot change predictions or observation identities.
    The corrected source must equal the file at the immutable reviewed commit.
    """
    supplied_path = CODE_ROOT / "artifacts/strong_lensing_geometry_20261021/mask_baseline_manifest.json"
    supplied_bytes = supplied_path.read_bytes()
    supplied_sha = hashlib.sha256(supplied_bytes).hexdigest()
    if supplied_sha != SUPPLIED_BASELINE_SHA256:
        raise RuntimeError("Delivered baseline manifest changed")
    supplied = json.loads(supplied_bytes)
    manifest_path, legacy_path = Path(manifest_path).resolve(), Path(legacy_path).resolve()
    value = json.loads(manifest_path.read_text(encoding="utf-8"))
    if value["source_commit"] != REVIEWED_COMMIT:
        raise RuntimeError("Unexpected baseline manifest source commit")
    if value["json_semantic_sha256"] != supplied["json_semantic_sha256"]:
        raise RuntimeError("Integration baseline cannot change original JSON identities")
    if value["canonical_json_rule"] != supplied["canonical_json_rule"]:
        raise RuntimeError("Integration baseline cannot change canonical JSON rule")
    legacy_bytes = legacy_path.read_bytes()
    legacy_hash = hashlib.sha256(legacy_path.read_text(encoding="utf-8").encode("utf-8")).hexdigest()
    if manifest_path == supplied_path.resolve():
        if legacy_hash != value["legacy_fitter_normalized_text_sha256"]:
            raise RuntimeError("Legacy fitter no longer matches the immutable baseline")
        return value, {"mode": "SUPPLIED_MANIFEST", "supplied_manifest_sha256": supplied_sha}
    correction = value.get("integration_provenance_correction", {})
    required = {
        "kind": "IMMUTABLE_REPOSITORY_SOURCE_CORRECTION_NOT_FITTER_CHANGE",
        "supplied_manifest_path": "artifacts/strong_lensing_geometry_20261021/mask_baseline_manifest.json",
        "supplied_manifest_sha256": supplied_sha,
        "supplied_legacy_fitter_normalized_text_sha256": supplied["legacy_fitter_normalized_text_sha256"],
        "verified_git_commit": REVIEWED_COMMIT,
        "verified_git_path": REVIEWED_FITTER,
    }
    if any(correction.get(k) != v for k, v in required.items()):
        raise RuntimeError("Invalid explicit integration provenance correction")
    git_bytes = subprocess.check_output(
        ["git", "show", REVIEWED_COMMIT + ":" + REVIEWED_FITTER],
        cwd=repository_root,
    )
    git_sha = hashlib.sha256(git_bytes).hexdigest()
    git_normalized = hashlib.sha256(git_bytes.decode("utf-8").replace("\r\n", "\n").encode("utf-8")).hexdigest()
    if legacy_bytes != git_bytes:
        raise RuntimeError("Current legacy fitter differs from immutable reviewed Git object")
    if correction.get("verified_git_object_raw_sha256") != git_sha:
        raise RuntimeError("Incorrect declared immutable Git object hash")
    if correction.get("verified_git_object_normalized_text_sha256") != git_normalized:
        raise RuntimeError("Incorrect declared normalized Git source hash")
    if value["legacy_fitter_normalized_text_sha256"] != git_normalized or legacy_hash != git_normalized:
        raise RuntimeError("Integration legacy-source hash differs from verified Git source")
    return value, {
        "mode": "EXPLICIT_REPOSITORY_VERIFIED_CORRECTION",
        "supplied_manifest_sha256": supplied_sha,
        "supplied_legacy_fitter_hash": supplied["legacy_fitter_normalized_text_sha256"],
        "actual_verified_legacy_fitter_hash": legacy_hash,
        "verified_git_commit": REVIEWED_COMMIT,
        "supplied_ten_JSON_identities_preserved": True,
        "original_fitter_changed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=CODE_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--native-dir", type=Path,
                        help="Optional full CAL directory for canonical global-coordinate fitting")
    parser.add_argument("--baseline-manifest", type=Path,
                        help="Explicit integration baseline; default preserves the supplied strict guard")
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path.append(str(root / "code"))
    from strong_lensing_jwst_rate_mask import fit_jwst_native_rate_point
    import strong_lensing_pixel_geometry as pg

    base = root / "artifacts/strong_lensing_geometry_20261021"
    manifest_path = (args.baseline_manifest or
                     CODE_ROOT/"artifacts/strong_lensing_geometry_20261021/mask_baseline_manifest.json")
    baseline_manifest, baseline_audit = verified_baseline(manifest_path, pg.__file__, root)
    output = args.output.resolve()
    protected = [root/"data",root/"preregistration",base/"jwst_target_diagnostic",
                 root/"code",root/"scripts",CODE_ROOT/"code",CODE_ROOT/"scripts"]
    if any(output.is_relative_to(path.resolve()) for path in protected):
        raise RuntimeError("Output would overwrite protected original input or code")
    inputs = {}
    def track(path, expected=None):
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected is not None and actual != expected:
            raise RuntimeError(f"Archived input hash mismatch: {path}")
        try:
            relative = path.resolve().relative_to(root).as_posix()
        except ValueError:
            relative = "external_native/"+path.name
        inputs[relative] = actual
        return path

    def trusted_json(path):
        relative = path.resolve().relative_to(root).as_posix()
        expected = baseline_manifest["json_semantic_sha256"].get(relative)
        if expected is None:
            raise RuntimeError(f"Unlisted immutable JSON input: {relative}")
        value = json.loads(track(path).read_text(encoding="utf-8"))
        canonical = json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False,allow_nan=False)
        if hashlib.sha256(canonical.encode("utf-8")).hexdigest() != expected:
            raise RuntimeError(f"JSON input differs from immutable baseline: {relative}")
        return value

    original = trusted_json(base/"jwst_target_diagnostic/pixel_geometry_results.json")
    protocol = trusted_json(root/"preregistration/he0435_geometry_jwst_diagnostic_20261021.json")
    origin = protocol["tangent_origin_radec_deg"]
    rows = []
    baseline_values, baseline_covariances = [], []
    largest_original_centroid_error = 0.
    largest_original_q_error = 0.
    largest_geometry_crosscheck_error = 0.
    largest_jacobian_crosscheck_error = 0.
    total_restored = 0

    def complex_geometry(points):
        z = np.asarray(points) @ np.array([1,1j])
        b = z[1]-z[0]
        rotation = b.conjugate()/abs(b)
        c,d = (z[2:]-z[0])*rotation
        return np.array([abs(b),c.real,c.imag,d.real,d.imag])

    exposure_paths = [base/"jwst_target_diagnostic/exposures"/
                      f"jw01198007001_02101_{i:05d}_nrcb3_cal.json" for i in range(1,5)]
    for filename in exposure_paths:
        old = trusted_json(filename)
        basename = filename.stem
        header_path = root/"data/strong_lensing_geometry_20261021/metadata"/f"{basename}_headers_only.json"
        header = trusted_json(header_path)["primary_header"]
        if header["FILTER"] != "F115W" or header["DETECTOR"] != "NRCB3":
            raise RuntimeError("Only the already-consumed F115W NRCB3 exposure is allowed")
        if [point["label"] for point in old["points"]] != list("ABCD"):
            raise RuntimeError("Expected exactly ordered A/B/C/D point records")
        if header["FILENAME"] != basename+".fits":
            raise RuntimeError("Header/native exposure mismatch")
        native_arrays = None
        if args.native_dir is not None:
            native_path = track(args.native_dir/(basename+".fits"),old["sha256"])
            with fits.open(native_path,memmap=False) as h:
                if h[0].header["CAL_VER"] != header["CAL_VER"]:
                    raise RuntimeError("Native/header calibration provenance mismatch")
                native_arrays = [np.array(h[name].data) for name in ("SCI","ERR","DQ")]
        positions, covariances, old_positions, old_covariances = [],[],[],[]
        points = []
        wcs_archive = old["official_GWCS_archive"]
        wcs_path = track(base/"jwst_target_diagnostic"/wcs_archive["path"],wcs_archive["sha256"])
        with asdf.open(wcs_path) as af:
            wcs = af["wcs"]
            for point in old["points"]:
                label = point["label"]
                archive = old["native_stamp_archives"][label]
                path = track(base/"jwst_target_diagnostic"/archive["path"],archive["sha256"])
                with fits.open(path) as h:
                    if h[0].header["CAL_SHA"] != old["sha256"]:
                        raise RuntimeError("Native stamp/source provenance mismatch")
                    data,error,dq = [np.array(h[name].data) for name in ("SCI","ERR","DQ")]
                    offset = np.array([h[0].header["X0"],h[0].header["Y0"]])
                path = track(root/point["field_psf_path"],point["field_psf_sha256"])
                with fits.open(path) as h:
                    grid = np.array(h[0].data)
                    oversampling = h[0].header["OVERSMP"]
                    psf_origin = (h[0].header["ORIGINX"],h[0].header["ORIGINY"])
                if native_arrays is not None:
                    data,error,dq = native_arrays
                    offset = np.zeros(2)
                guess = np.array(point["initial_guess_xy_pixel"])-offset
                oldfit = pg.fit_empirical_point(
                    data,error,dq,guess,grid,oversampling,psf_origin,background_degree=2
                )
                fit = fit_jwst_native_rate_point(
                    data,error,dq,guess,grid,oversampling,psf_origin,
                    primary_header=header,background_degree=2,
                )
                oldxy = np.array(oldfit["centroid_xy_pixel"])+offset
                xy = np.array(fit["centroid_xy_pixel"])+offset
                centroid_error = float(np.max(np.abs(oldxy-point["centroid_xy_pixel"])))
                if centroid_error > 1e-4:
                    raise RuntimeError("Original centroid replay failed")
                if abs(oldfit["chi2"]-point["chi2"]) > 1e-4*point["chi2"]:
                    raise RuntimeError("Original pixel score replay failed")
                largest_original_centroid_error = max(largest_original_centroid_error,centroid_error)
                p,c,_ = pg.gwcs_centroid(wcs,xy,fit["centroid_covariance_pixel2"],origin)
                po,co,_ = pg.gwcs_centroid(wcs,oldxy,oldfit["centroid_covariance_pixel2"],origin)
                positions.append(p);covariances.append(c)
                old_positions.append(po);old_covariances.append(co)
                restored = fit["usable_pixels"]-oldfit["usable_pixels"]
                total_restored += restored
                points.append({
                    "label":label,"restored_pixels":restored,
                    "centroid_xy_pixel":xy.tolist(),
                    "centroid_covariance_pixel2":fit["centroid_covariance_pixel2"],
                    "tangent_east_north_arcsec":p.tolist(),
                    "tangent_covariance_arcsec2":c.tolist(),
                    "old_chi2":oldfit["chi2"],"old_dof":oldfit["dof"],
                    "chi2":fit["chi2"],"dof":fit["dof"],
                    "p_value":fit["p_value"],
                    "pixel_model_adequate":bool(fit["p_value"]>=.001),
                    "pixel_mask_receipt":fit["pixel_mask_receipt"],
                    "dq_flag_counts_in_stamp":fit["dq_flag_counts_in_stamp"],
                })
        q,cov = pg.intrinsic_with_covariance(positions,block_diag(*covariances))
        oldq,oldcov = pg.intrinsic_with_covariance(old_positions,block_diag(*old_covariances))
        q_error = float(np.max(np.abs(oldq-old["intrinsic_q_arcsec"])))
        if q_error > 1e-6:
            raise RuntimeError("Original GWCS/geometry replay failed")
        largest_original_q_error = max(largest_original_q_error,q_error)
        largest_geometry_crosscheck_error = max(
            largest_geometry_crosscheck_error,float(np.max(np.abs(q-complex_geometry(positions)))))
        step = 1e-6
        x = np.array(positions).reshape(-1)
        jac = np.column_stack([
            (complex_geometry((x+step*np.eye(8)[i]).reshape(4,2))-
             complex_geometry((x-step*np.eye(8)[i]).reshape(4,2)))/(2*step)
            for i in range(8)])
        largest_jacobian_crosscheck_error = max(
            largest_jacobian_crosscheck_error,
            float(np.max(np.abs(jac-pg.intrinsic_jacobian(positions)))))
        rows.append({"exposure":basename,"points":points,
                     "intrinsic_q_arcsec":q.tolist(),"intrinsic_covariance_arcsec2":cov.tolist()})
        baseline_values.append(oldq);baseline_covariances.append(oldcov)

    if len(rows)!=4 or sum(len(row["points"]) for row in rows)!=16:
        raise RuntimeError("Expected all four exposures and sixteen points")
    combined = pg.combine_exposures(
        [row["intrinsic_q_arcsec"] for row in rows],
        [row["intrinsic_covariance_arcsec2"] for row in rows])
    baseline_combined = pg.combine_exposures(baseline_values,baseline_covariances)
    scores = {}
    for name,model in original["models"].items():
        prediction = model["predicted_intrinsic_q_arcsec"]
        residual = np.array(combined["mean"])-prediction
        score = pg.gaussian_score(residual,combined["covariance"])
        score["intrinsic_coordinate_rms_mas"] = float(np.sqrt(np.mean(residual**2))*1000)
        score["conditional_geometry_validation_pass"] = False
        scores[name] = score
    result = {
        "status":"DEVELOPMENT_MASK_REPAIR_REPLAY_NOT_NEW_VALIDATION",
        "source_commit":"7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9",
        "coordinate_mode":"native_global" if args.native_dir is not None else "portable_stamp_local",
        "replay_tolerance":{"centroid_pixel":1e-4,"relative_pixel_chi2":1e-4,"geometry_arcsec":1e-6},
        "retained_extra_pixels":total_restored,"exposures":rows,
        "combined":combined,"models":scores,
        "baseline_combined_replay":baseline_combined,
        "audit":{
            "maximum_original_centroid_replay_error_pixel":largest_original_centroid_error,
            "maximum_original_geometry_replay_error_arcsec":largest_original_q_error,
            "independent_complex_geometry_error_arcsec":largest_geometry_crosscheck_error,
            "independent_finite_difference_jacobian_error":largest_jacobian_crosscheck_error,
            "original_fitter_modified":False,"original_dq_modified":False,
            "native_errors_modified":False,"lens_predictions_refitted":False,
        },
        "input_sha256":inputs,
        "code_sha256":{
            "code/strong_lensing_pixel_geometry.py":hashlib.sha256(Path(pg.__file__).read_bytes()).hexdigest(),
            "code/strong_lensing_jwst_rate_mask.py":hashlib.sha256((CODE_ROOT/"code/strong_lensing_jwst_rate_mask.py").read_bytes()).hexdigest(),
            "scripts/replay_he0435_verified_rate_mask.py":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        "physical_admission":False,"unchanged_physical_guard":2e-10,
        "interpretation":"Mask eligibility fixed; PSF adequacy and absolute centroid calibration remain open.",
    }
    result["baseline_manifest_provenance"] = baseline_audit
    result["baseline_manifest_sha256"] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if largest_geometry_crosscheck_error > 1e-8 or largest_jacobian_crosscheck_error > 1e-6:
        raise RuntimeError("Independent geometry/Jacobian crosscheck failed")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"retained_extra_pixels":total_restored,
        "exposure_consistency_chi2":combined["exposure_consistency_chi2"],
        "models":{name:score["intrinsic_coordinate_rms_mas"] for name,score in scores.items()},
        "audit":result["audit"]},indent=2))


if __name__ == "__main__":
    main()
