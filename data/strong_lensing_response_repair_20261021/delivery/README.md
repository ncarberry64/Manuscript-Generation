# HE 0435 response repair package

This package contains an executed JWST rate-mask correction and the subsequent
instrument-response development calculations for Manuscript-Generation PR #23.
Start with **HE0435_RESPONSE_REPAIR_REPORT.md**, then **NEXT_CODEX_TASK.md**.

## Production patch

- `mask_patch/jwst_rate_mask_repair.patch`: ten additive repository files.
- `mask_patch/mask_patch_manifest.json`: immutable source commit and file hashes.
- `mask_patch/docs/CODEX_HANDOFF.md`: application, reproduction and Git workflow.

The reviewed source commit is
`7ec9cea7b4d1c642e04407a6ef2f0a4eaef1baf9` in
https://github.com/ncarberry64/Manuscript-Generation/pull/23.

The patch restores 303 usable pipeline slopes, gives canonical mask-only
G1-reference RMS 0.448457087207 mas, and preserves all original checkpoints.
It does not claim completed PSF calibration or physical validation.

## Experiment layout

The scripts were run in a working directory with these relative locations:

- `repo/`: a checkout or read-only snapshot of the reviewed manuscript commit,
  including its original numerical artifacts and the four original native
  F115W CAL files.
- `mask_patch/`: the adapter, replay and regression tests delivered here.
- `audit/`: positive response, spatial, finer-basis and local-quasar experiments.
- `stpsf_field_calibration/`: bounded color/diffusion and pupil-phase experiments.
- `stpsf_trial/wss_opds/`: the actual WSS wavefront used in physical PSF runs.
- `ramp_diagnostic/`: historical limited diagnostic and public input receipts.
- `ramp_diagnostic/native_pipeline/`: completed recorded-pipeline and common-time
  checks, compact corrected group stamps, software/reference records and scripts.

The large native observations, detector references and Python environments are
not duplicated in this archive. Their exact filenames and hashes are recorded.
The public acquisition scripts and the existing repository receipts provide
the inputs. Do not substitute a later archive reprocessing with different
bytes. No reserved target band was opened for this work.

For development replay, make a separate working copy of this package and put
the reviewed repository at `repo/` in that copy. The development scripts write
new result files beside themselves. Preserve the delivered result receipts
and `DELIVERY_MANIFEST.json` as the record of this executed pass.

For the refined local-quasar calculation, run `audit/local_quasar_response.py`
before `audit/refine_local_quasar_response.py`: the first script regenerates
the optimizer starting points. The definitive refined receipts and response
arrays are included; superseded optimizer outputs are not. Field leave-one-out
shape arrays were not separately saved, so reproducing those folds requires
rerunning the field scripts rather than loading the all-field shape.

The empirical reproduction order is:

```bash
python audit/positive_field_response.py
python audit/transfer_positive_field.py
python audit/spatial_field_response.py
python audit/fine_field_response.py
python audit/local_quasar_response.py
python audit/refine_local_quasar_response.py
python audit/local_quasar_geometry.py
```

`audit/field_response_transfer.py` is retained for an imported FITS helper;
its superseded additive trial does not need to run. The supplementary synthetic
check receipt records an inline implementation check; its original inline
generator is not part of the standalone reproduction commands.

The STPSF calculations use STPSF 2.2.0 and its official 2.2.0 reference data,
with equal FFT/detector oversampling. `download_stpsf_reference.py` records the
reference archive URL and digest. The initial unequal-sampling STPSF experiment
is explicitly superseded; its receipt is retained for audit purposes.

The production mask runtime is recorded in
`mask_patch/artifacts/strong_lensing_geometry_20261021/mask_repair_environment.json`.
The completed full-ramp diagnostic records its separate software/reference
configuration and reproduction sequence in
`ramp_diagnostic/native_pipeline/NATIVE_READ_HISTORY_FINDINGS.md`.
The definitive receipt is `common_time_field_stcal118_results.json`; the earlier
`common_time_field_results.json` retains the unsuccessful 1.19.1 numerical gate.
Do not assume a new rerun is identical merely because the primary JWST package
version agrees.

## Status and interpretation

The mask correction is tested and ready to integrate. The response trials are
development evidence with their limitations and failed adequacy flags intact.
The positive field-trained mean response reduces the raw geometry discrepancy
to about 0.36 mas, while formal measurement adequacy remains unresolved.
A separate one-star common-time read extraction reduces raw profile residual
by 32.44%, but still fails its declared adequacy threshold and loses precision.
It is a measured lead on the remaining response, not a completed repair.

No GitHub branch or commit was written by this session: branch creation through
the selected connection was denied with HTTP 403. The external Codex checkout
can apply and publish the supplied patch using the already authorized workflow.
No manuscript submission or main-branch merge is part of this package.
