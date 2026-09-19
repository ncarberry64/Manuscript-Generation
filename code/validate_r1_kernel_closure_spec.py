from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = ROOT / "docs" / "r1_kernel_closure_spec.json"
MANIFEST = ROOT / "preregistration" / "prediction_manifest.json"

EXPECTED_MANIFEST_SHA256 = (
    "0a995d15171f3b133edaf203869329e9e"
    "ec7913a9ce8e23e8469175bc7538bfa"
)

REQUIRED_KERNELS = {
    "F_DL",
    "F_H",
    "F_DA",
    "F_fsigma8",
    "F_lens",
}


def canonical_json_digest(path: Path) -> str:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    payload = json.dumps(
        obj,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def validate() -> dict:
    spec = json.loads(SPEC.read_text(encoding="utf-8"))

    kernels = set(spec["required_kernels"])
    if kernels != REQUIRED_KERNELS:
        raise AssertionError(
            f"kernel set mismatch: {sorted(kernels)}"
        )

    forbidden = set(spec["forbidden_operations"])
    required_forbidden = {
        "independent_directional_amplitude_refit",
        "post_hoc_kernel_retuning_to_prospective_data",
        "duplicate_physical_scalar_mode",
        "changing_frozen_R1_manifest",
    }
    if not required_forbidden.issubset(forbidden):
        raise AssertionError("missing no-retuning protections")

    digest = canonical_json_digest(MANIFEST)
    if digest != EXPECTED_MANIFEST_SHA256:
        raise AssertionError(
            f"frozen manifest changed: {digest}"
        )

    calibration = spec["calibration"]
    if calibration["observable"] != "SN_distance_modulus_dipole":
        raise AssertionError("unexpected calibration observable")

    return {
        "status": "R1_KERNEL_CLOSURE_PROTOCOL_VALID",
        "manifest_sha256": digest,
        "required_kernels": sorted(kernels),
        "calibration": calibration,
    }


if __name__ == "__main__":
    print(json.dumps(validate(), indent=2, sort_keys=True))
