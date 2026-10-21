# Verified integration of the HE 0435 rate-mask package

The supplied ten-file patch is preserved byte for byte, including its original
manifest, replay receipts and environment record. All 109 entries in the ZIP's
delivery manifest were verified before additive application. The package SHA256
is `28ee0e7d343a428fc81a0ee91ab907111d02170372caa2e9db358ef70a07dad3`.

The first replay correctly stopped at its source-identity guard. The supplied
baseline declares legacy fitter SHA256 `e7333e45…`; the unchanged fitter at the
declared Git commit `7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9` has SHA256
`41fa74373d8e7caf658af50fa7fd55d714f0c2d36af3e9e6660a75a6ded15f0a`.
Its worktree bytes equal the immutable Git object, and all ten supplied JSON
semantic hashes match. The discrepancy is recorded in
`artifacts/strong_lensing_response_repair_20261021/mask_first_replay_guard_failure.json`.
The delivered source hash is not represented as reproduced.

The additive `replay_he0435_verified_rate_mask.py` runner accepts an explicit
integration manifest. It verifies the original manifest's exact hash, retains
all ten JSON identities, and checks the actual fitter bytes against the file at
the immutable reviewed Git commit. The delivered default guard still rejects
the mismatch. The numerical replay body is identical to the supplied runner;
the original fitter, selection policy, observations, predictions and protocol
are unchanged. Six additional regressions exercise these provenance guards.
Those unit tests use an exact SHA-bound Git-response fixture so they also run
in shallow CI checkouts; the actual native replay uses the real immutable Git
object and records that separate integration check.

Executed in the existing Windows geometry environment:

```powershell
.venv-geometry/Scripts/python.exe tests/test_strong_lensing_jwst_rate_mask.py
.venv-geometry/Scripts/python.exe tests/test_strong_lensing_rate_mask_integration.py
.venv-geometry/Scripts/python.exe scripts/replay_he0435_verified_rate_mask.py --root . --baseline-manifest artifacts/strong_lensing_response_repair_20261021/mask_verified_baseline_manifest.json --native-dir data/strong_lensing_geometry_20261021/native --output artifacts/strong_lensing_response_repair_20261021/mask_verified_native_replay.json
```

The four supplied regression groups and six integration regressions pass. The
canonical native replay restores 303 valid slopes and exactly reproduces all
16 original centroids and their original intrinsic geometry in this runtime.
The corrected G1 reference/R1 raw RMS values are 0.448457055291 and
0.448562654687 mas. Exposure consistency is 253.027694606 for 15 degrees of
freedom. All 16 absolute profile checks fail; their reduced chi-square spans
27.56851–83.55374. The reference formal geometry statistic is 1033.11355549
for five degrees of freedom. Native ERR has not been rescaled.

The actual runtime is Python 3.14.3, NumPy 2.3.5, SciPy 1.17.0, Astropy 7.2.0,
ASDF 6.0.0 and GWCS 1.0.3. It is distinct from the delivered Linux environment.
The reference RMS differs from the supplied canonical replay by approximately
3.2e-8 mas; the exposure statistic differs by 8.5e-6. The actual native hashes
and producing CAL_VER=2.0.1 headers were checked. The producing stcal version
is not recorded in those headers. This execution performs point fitting, not
JWST/stcal detector-pipeline reprocessing.

`mask_integration_review.json` and `mask_patch_integration_environment.json`
record the executed checks, package versions, native provenance and result
hashes. The inherited full suite and manuscript build are not claimed by this
integration audit.

This is a correction to already consumed F115W development data. PSF/scene mean
response and calibrated correlated uncertainty remain inadequate. There is no
new observational validation, physical admission remains false, and the 2e-10
physical guard is unchanged. No reserved band was opened.
