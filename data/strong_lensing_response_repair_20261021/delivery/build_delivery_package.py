"""Assemble the reviewed repair and definitive development receipts.

Large public observations, detector references, environments, provisional
optimizers and reproducible caches are deliberately excluded. Run after the
common-time diagnostic and root report have been finalized.
"""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

ROOT = Path(__file__).resolve().parent
BASE_COMMIT = "7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    entries = {}

    def add(relative, destination=None):
        path = ROOT / relative
        if not path.is_file() or path.is_symlink():
            raise RuntimeError(f"Required regular file unavailable: {relative}")
        if path.stat().st_size > 20 * 1024**2:
            raise RuntimeError(f"Unexpected large file: {relative}")
        entries[destination or relative] = path

    for name in ["HE0435_RESPONSE_REPAIR_REPORT.md", "NEXT_CODEX_TASK.md",
                 "EXPERIMENT_INDEX.md", "build_delivery_package.py"]:
        add(name)
    readme = "README_DELIVERY.md" if (ROOT / "README_DELIVERY.md").exists() else "README.md"
    add(readme, "README.md")
    for path in sorted((ROOT / "mask_patch").rglob("*")):
        if path.is_file() and "__pycache__" not in path.parts:
            add(path.relative_to(ROOT).as_posix())

    for name in [
        "field_response_transfer.py", "positive_field_response.py",
        "transfer_positive_field.py", "spatial_field_response.py",
        "fine_field_response.py", "local_quasar_response.py",
        "refine_local_quasar_response.py", "local_quasar_geometry.py",
        "positive_field_results.json", "positive_target_transfer.json",
        "spatial_field_results.json", "fine_field_results.json",
        "local_quasar_refined_results.json", "local_quasar_geometry_results.json",
        "positive_field_shape.npz", "local_quasar_shape_refined.npz",
        "local_quasar_fold1_refined.npz", "local_quasar_fold2_refined.npz",
        "local_quasar_fold3_refined.npz", "local_quasar_fold4_refined.npz",
        "fine_field_shape.npz", "AUDIT_FINDINGS.md",
        "positive_response_synthetic_checks.json",
    ]:
        add("audit/" + name)

    for name in ["download_used_cal.py", "download_used_ramp.py",
                 "download_stpsf_reference.py", "diagnose_used_ramp.py",
                 "run_stpsf_field_trial.py", "run_stpsf_equal_sampling_trial.py"]:
        add(name)
    for relative in [
        "stpsf_reference/download_receipt.json",
        "stpsf_trial/field_trial.json",
        "stpsf_trial/wss_opds/R2022101203-NRCA3_FP1-1.fits",
        "stpsf_equal_sampling/field_trial.json",
        "stpsf_equal_sampling/jw01198007001_02101_00001_nrcb3_FIELD_stpsf_epsf.fits",
        "stpsf_field_calibration/run_bounded_field_calibration.py",
        "stpsf_field_calibration/run_bounded_phase_trial.py",
        "stpsf_field_calibration/protocol.json",
        "stpsf_field_calibration/results.json",
        "stpsf_field_calibration/coarse_grid.json",
        "stpsf_field_calibration/spectral_weights.npz",
        "stpsf_field_calibration/final_field_pixel_receipts.json",
        "stpsf_field_calibration/phase_trial/protocol.json",
        "stpsf_field_calibration/phase_trial/results.json",
        "ramp_diagnostic/download_receipt.json",
        "ramp_diagnostic/ramp_diagnostic.json",
        "ramp_diagnostic/field_ramp_stamp.npz",
    ]:
        add(relative)
    for index in range(1, 5):
        add(f"stpsf_field_calibration/field_{index}_selected_epsf.fits")
        add(f"stpsf_field_calibration/phase_trial/field_{index}_global_exact_epsf.fits")

    native = ROOT / "ramp_diagnostic/native_pipeline"
    # A finished result is required; never issue a package claiming an unrun test.
    final_receipt = native / "common_time_field_stcal118_results.json"
    if not final_receipt.is_file():
        raise RuntimeError("Common-time diagnostic has not finished")
    replay_gate = json.loads(final_receipt.read_text())["native_pipeline_replay"]["numerical_replay_gate"]
    if not replay_gate["pass"]:
        raise RuntimeError("Final native field replay remains unresolved")
    for pattern in ["*.py", "*.json", "*.npz", "*.md", "*.txt", "*.log", "refs/*.json"]:
        for path in sorted(native.glob(pattern)):
            add(path.relative_to(ROOT).as_posix())

    manifest = {
        "package": "HE0435_RESPONSE_REPAIR_20261021",
        "reviewed_repository": "ncarberry64/Manuscript-Generation",
        "reviewed_commit": BASE_COMMIT,
        "original_pr": "https://github.com/ncarberry64/Manuscript-Generation/pull/23",
        "remote_write_performed": False,
        "precision_validation_claimed": False,
        "excluded": ["Full native and raw observations", "Large official detector references",
                     "Python environments", "Regenerable PSF caches", "Provisional optimizers"],
        "files": [],
    }
    payloads = {}
    for name, path in sorted(entries.items()):
        data = path.read_bytes()
        payloads[name] = data
        manifest["files"].append({"path": name, "size_bytes": len(data),
                                  "sha256": hashlib.sha256(data).hexdigest()})
    payloads["DELIVERY_MANIFEST.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix + ".partial")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=6) as archive:
        for name, data in sorted(payloads.items()):
            info = zipfile.ZipInfo(name, date_time=(2026, 10, 21, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, data)
    with zipfile.ZipFile(temporary) as archive:
        bad = archive.testzip()
        if bad:
            raise RuntimeError(f"Archive integrity failed: {bad}")
        for record in manifest["files"]:
            if hashlib.sha256(archive.read(record["path"])).hexdigest() != record["sha256"]:
                raise RuntimeError(f"Archive payload hash mismatch: {record['path']}")
    temporary.replace(args.output)
    print(json.dumps({"path": str(args.output.resolve()), "files": len(payloads),
                      "bytes": args.output.stat().st_size,
                      "sha256": hashlib.sha256(args.output.read_bytes()).hexdigest()}, indent=2))


if __name__ == "__main__":
    main()
