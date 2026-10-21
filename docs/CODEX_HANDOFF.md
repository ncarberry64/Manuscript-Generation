# Codex handoff: integrate the verified HE 0435 rate-mask repair

## Task

Integrate and publish the supplied, tested JWST native-rate mask correction in
the manuscript Git workflow, then continue the instrumental-response repair
before making a precision-validation claim. The user has authorized Git work.
This package was prepared locally after the selected GitHub connector returned
403 for branch creation; no remote branch, commit, PR comment or merge was made.

Repository: ncarberry64/Manuscript-Generation
Existing draft PR: https://github.com/ncarberry64/Manuscript-Generation/pull/23
Existing branch: codex/strong-lensing-first-comparison-20261020
Immutable reviewed starting point: 7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9

Read the repository's current instructions and inspect its current state first.
If the branch has advanced, inspect and preserve those changes. Do not reset,
force-push, or overwrite another session's work. Use an isolated worktree for
this repair.

The BHSM framing remains primary: the reference and R1-inclusive calculations
are controls within the specified response calculation. A fully matched
conventional effective-mass representation describes the same retained optics.
No requirement that BHSM must beat conventional gravity is introduced.

## Confirmed issue and correction

The old native_pixel_mask excluded DO_NOT_USE | SATURATED. In these native
CAL_VER=2.0.1 files, SATURATED can mean that a later ramp group saturated while
earlier unsaturated groups still supplied a valid slope. All 303 restored
target pixels have finite SCI, finite positive ERR and no DO_NOT_USE flag.

The additive adapter checks native 2D JWST/NIRCam CAL/rate provenance, completed
ramp fitting and the audited software convention. It retains valid SAT-only
slopes through a private selection array passed to the unchanged fitter.
The original SCI/ERR/DQ, original fitter and frozen protocol are preserved.

The strict release guard currently accepts only CAL_VER=2.0.1. Do not remove
that guard to accommodate another product. Check the other release's actual
stcal convention first and document the evidence before extending the table.
This policy is not applicable to GROUPDQ/raw ramps or resampled mosaics.

## Apply the package

The package contains jwst_rate_mask_repair.patch and mask_patch_manifest.json.
The patch only adds files. Confirm its hashes against the manifest.

A possible worktree sequence, executed from an existing repository checkout:

    git status --short
    git fetch origin
    git worktree add -b codex/he0435-jwst-rate-mask-repair ../he0435-jwst-rate-mask-repair 7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9

Inside that worktree, replace the patch path with its actual extracted location:

    git apply --check /path/to/jwst_rate_mask_repair.patch
    git apply /path/to/jwst_rate_mask_repair.patch

If a branch or worktree already exists, inspect and resume it instead of
recreating it. These commands do not merge the repair into main.

The original 123-test statement belongs to the published checkpoint. This
repair package passed four regression groups plus the specified numerical
replays. It has not rerun the entire inherited 123-test suite or journal build.

## Reproduce using the existing geometry Python environment

The executable below is written as python. In the Windows worktree, use the
existing .venv-geometry\Scripts\python.exe if that is the established runtime.
Record the actual package versions; do not silently relabel an old test run.

Run the focused mask regressions:

    python tests/test_strong_lensing_jwst_rate_mask.py

Run the canonical replay against the four original CAL files:

    python scripts/replay_he0435_rate_mask.py --root . --native-dir data/strong_lensing_geometry_20261021/native --output artifacts/strong_lensing_geometry_20261021/mask_repair_native_replay.json

The native files are the same already consumed F115W NRCB3 exposures. Their
hashes are checked against the archived receipts. Do not substitute a later
MAST reprocessing with a changed hash, and do not open an unused band.

If the full native inputs are unavailable locally, run the portable replay:

    python scripts/replay_he0435_rate_mask.py --root . --output artifacts/strong_lensing_geometry_20261021/mask_repair_replay.json

The runner checks an immutable baseline manifest for the original predictions,
exposure results, headers and protocol. It explicitly selects the four consumed
F115W exposures and refuses outputs inside protected original-input/code paths.
It cannot relabel a new prediction file as the old baseline.

## Expected result

Canonical global-detector-coordinate calculation:

- Restored slopes: 303 across all sixteen target stamps.
- G1 reference raw geometry RMS: 0.448457087207 mas.
- G1 R1 raw geometry RMS: 0.448562686643 mas.
- Between-exposure statistic: 253.027686060 for 15 coordinates.
- G1 reference formal geometry statistic: 1033.113139179 for 5 coordinates.
- All sixteen pixel-profile adequacy checks still fail.
- Corrected chi2/dof range: approximately 27.57 to 83.55.

The original 0.779400/0.779454 mas extraction and exposure statistic 432.262011
remain historical results. Report the mask-only correction alongside them.
Restoring information tightens conditional errors, so the lower raw RMS does
not imply a lower formal geometry statistic or validated astrometry.

The native-coordinate replay reproduces original centroids within 3.41e-12
pixel. The small-stamp portable replay gives 0.448446682409 mas and exposure
statistic 253.133300295 because optimizer coordinates change. Its quantified
difference is numerical conditioning, not independent validation. Retain the
native-coordinate values as canonical.

Keep all original failures and original validation flags. This is reanalysis of
used F115W observations; it does not create another untouched test. The original
2e-10 guard and physical_admission=false remain unchanged.

## Do not repeat an already active pipeline setting

All four actual native ASDF configurations already have:

    algorithm: OLS_C
    suppress_one_group: False

Their logs also record zero-frame processing. Merely enabling
suppress_one_group=False again will repeat the existing setup.

The current upstream development changelog mentions an NFRAMES correction for
read-noise variance of ZEROFRAME-derived OLS_C slopes. These inputs have
NFRAMES=2, but which individual slopes used ZEROFRAME has not been established
from CAL DQ alone. No unreleased variance patch or error multiplier is applied
in this package. Investigate actual read provenance before any such correction.

## Continue the response repair

Use the corrected eligibility policy consistently for field-star calibration
and quasar extraction. The remaining task is an adequate instrument/extended-light
mean and calibrated astrometric error propagation. No additional gravitational
term is inferred from a PSF-shape residual.

Any accompanying positivity-preserving empirical PSF experiment remains
development work. A better field or target residual alone does not establish
transferable calibration. Check that the correction preserves positivity,
normalization and the chosen position convention; separate its training from
its assessment; and propagate its uncertainty and shared correlations into
centroid comparisons. Do not fit PSF/host freedom using agreement with the
frozen lens coordinates as an objective.

Keep the remaining untouched bands untouched until a concrete measurement
procedure and its assessment are fixed. The existing F115W observations can
be used openly for development.

## Git and manuscript integration

Inspect the added code and receipts, rerun the focused checks, and commit the
repair with the original failed experiment intact. Publish the repair branch
and a draft PR targeting the existing manuscript branch, or integrate into
draft PR #23 according to the current authorized repository workflow. Do not
merge main or submit the manuscript automatically.

Update the manuscript's current-results discussion to distinguish the corrected
pixel-selection policy from the still unresolved PSF/host/distortion response.
Do not replace old numerical tables as though the correction existed at the
original freeze. If the manuscript is rebuilt, report the build and visual
inspection actually performed on the new commit. Keep the journal Galaxies;
do not add a target submission date.

Return the published commit/PR, exact reproduction commands, what was corrected,
and the remaining specific instrumental dependency.

