# BHSM SN sightline tomography V1

Read **FINAL_REPORT.md** first. The analysis preserves the previous scientific prediction. It does not establish a robust coherent topographic residual.

Install the pinned requirements in Python 3.14, then run:

```powershell
python -m pip install -r requirements.txt
python run_all.py
python package_outputs.py
```

Numerical stages use only bundled data. Prior result packages are checked when present; their original manifests are bundled for portable reproduction. Do not use files in `audit/SUPERSEDED_free_PV_feature.zip` as final results.

Scientific outputs are deterministic at numerical tolerance; gzip/ZIP timestamps and some acquisition metadata need not be byte-identical after a rerun. `RESULT_MANIFEST.json` hashes the delivered snapshot. The archive excludes itself and bytecode caches.
