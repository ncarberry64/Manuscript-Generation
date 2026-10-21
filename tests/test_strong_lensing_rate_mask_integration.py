"""Provenance guards for the additive repository-verified replay runner."""
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "verified_mask_replay", ROOT / "scripts/replay_he0435_verified_rate_mask.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
SUPPLIED = ROOT / "artifacts/strong_lensing_geometry_20261021/mask_baseline_manifest.json"
EXPLICIT = ROOT / "artifacts/strong_lensing_response_repair_20261021/mask_verified_baseline_manifest.json"
FITTER = ROOT / "code/strong_lensing_pixel_geometry.py"


class IntegrationGuardTests(unittest.TestCase):
    def setUp(self):
        # Unit guards must run in a shallow CI checkout. The native integration
        # receipt separately exercises the real immutable Git object locally.
        frozen = FITTER.read_bytes()
        self.assertEqual(hashlib.sha256(frozen).hexdigest(),
                         "41fa74373d8e7caf658af50fa7fd55d714f0c2d36af3e9e6660a75a6ded15f0a")
        self.git_patch = patch.object(MODULE.subprocess, "check_output", return_value=frozen)
        self.git_show = self.git_patch.start()
        self.addCleanup(self.git_patch.stop)

    def test_delivered_default_still_rejects_source_hash_discrepancy(self):
        with self.assertRaisesRegex(RuntimeError, "Legacy fitter no longer matches"):
            MODULE.verified_baseline(SUPPLIED, FITTER, ROOT)

    def test_explicit_git_bound_correction_preserves_all_ten_json_identities(self):
        value, audit = MODULE.verified_baseline(EXPLICIT, FITTER, ROOT)
        supplied = json.loads(SUPPLIED.read_text())
        self.assertEqual(value["json_semantic_sha256"], supplied["json_semantic_sha256"])
        self.assertEqual(len(value["json_semantic_sha256"]), 10)
        self.assertEqual(audit["mode"], "EXPLICIT_REPOSITORY_VERIFIED_CORRECTION")
        self.assertFalse(audit["original_fitter_changed"])
        self.git_show.assert_called_once_with(
            ["git", "show", MODULE.REVIEWED_COMMIT + ":" + MODULE.REVIEWED_FITTER], cwd=ROOT
        )

    def test_current_source_tamper_rejected_against_immutable_git_object(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tampered_fitter.py"
            path.write_bytes(FITTER.read_bytes() + b"\n# altered source\n")
            with self.assertRaisesRegex(RuntimeError, "differs from immutable reviewed Git object"):
                MODULE.verified_baseline(EXPLICIT, path, ROOT)

    def test_modified_observation_or_prediction_identity_rejected(self):
        value = json.loads(EXPLICIT.read_text())
        first = next(iter(value["json_semantic_sha256"]))
        value["json_semantic_sha256"][first] = "0" * 64
        self.check_bad_manifest(value, "cannot change original JSON identities")

    def test_false_supplied_or_source_hash_claim_rejected(self):
        value = json.loads(EXPLICIT.read_text())
        wrong_supplied = deepcopy(value)
        wrong_supplied["integration_provenance_correction"]["supplied_manifest_sha256"] = "0" * 64
        self.check_bad_manifest(wrong_supplied, "Invalid explicit integration provenance")
        wrong_source = deepcopy(value)
        wrong_source["integration_provenance_correction"]["verified_git_object_raw_sha256"] = "0" * 64
        self.check_bad_manifest(wrong_source, "Incorrect declared immutable Git object hash")
        wrong_legacy = deepcopy(value)
        wrong_legacy["legacy_fitter_normalized_text_sha256"] = "0" * 64
        self.check_bad_manifest(wrong_legacy, "differs from verified Git source")

    def test_numerical_replay_body_is_identical_to_delivered_runner(self):
        delivered = (ROOT / "scripts/replay_he0435_rate_mask.py").read_text()
        integrated = (ROOT / "scripts/replay_he0435_verified_rate_mask.py").read_text()
        start, stop = '    original = trusted_json(', '    result = {'
        self.assertEqual(delivered[delivered.index(start):delivered.index(stop)],
                         integrated[integrated.index(start):integrated.index(stop)])

    def check_bad_manifest(self, value, message):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tampered_manifest.json"
            path.write_text(json.dumps(value))
            with self.assertRaisesRegex(RuntimeError, message):
                MODULE.verified_baseline(path, FITTER, ROOT)


if __name__ == "__main__":
    unittest.main(verbosity=2)
