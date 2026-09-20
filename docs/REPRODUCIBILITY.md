# Reproducing the manuscript and calculations

Run commands from the repository root unless stated otherwise. The branch is `theory/cosmology-universe-manuscript-final`. The release report records the tested scientific fingerprint and checkpoint; the intended submission is the final commit on this branch, pending author approval and archive. Retrieve it with `git rev-parse HEAD` rather than inferring it from timestamps.

## Environment

```powershell
python -m venv .venv
.venv/Scripts/python -m pip install -r requirements-dev.txt
.venv/Scripts/python -m pytest -q
```

On POSIX use `.venv/bin/python`. The validation environment's direct versions are in `requirements-validation.txt`; Python 3.14 on Windows was used locally. This is not a complete transitive lock. The CI workflow independently selects Python 3.12 on Linux. CI status must be read from the PR; local success does not imply remote success. Pandas, previously missing from requirements, is needed by foreground-test imports; matplotlib is needed for figures.

## Focused checks and full suite

```powershell
python -m pytest -q tests/test_bhsm_make_or_break_prediction.py tests/test_bhsm_optical_low_rank_covariance.py tests/test_bhsm_optical_transfer.py
python -m pytest -q tests/test_r1_n2_reduced_propagator.py tests/test_r1_n2_growth_lensing_transfer.py tests/test_r1_n2_luminosity_distance_kernel.py
python -m pytest -q --junitxml=artifacts/submission_pytest.xml -o cache_dir=manuscript/build/pytest-cache
```

The full suite can take a long time, particularly the nested distance integrals. Do not repeatedly restart it. Exactly one existing **strict xfail** is registered: `test_dae_preserves_00_0i_constraints_from_index_consistent_seed`, for the curved-EFT compatibility defect. An XPASS is a failure under strict mode. No unexpected failure should be relabeled xfail. The JUnit record describes software checks, not a new validation of all observational artifacts.

## Figures, numerical outputs, and safe regeneration

```powershell
python code/figures.py
python code/bhsm_make_or_break_prediction.py --output manuscript/build/make-or-break-synthetic.json
python code/r1_n2_reduced_propagator.py --artifact manuscript/build/r1-propagator.json
python code/r1_n2_growth_lensing_transfer.py --artifact manuscript/build/r1-growth-weyl.json
python code/r1_n2_luminosity_distance_kernel.py --artifact manuscript/build/r1-distance.json
```

Compare numerical JSON values with tolerances from the focused tests, not raw byte equality across numerical libraries. Existing `--tex` generator options emit historical section prose and can overwrite the editorial pass: use artifact-only output for reproduction. The generator flux escape bug was repaired without changing equations/numerical parameters. Historical populate/apply PowerShell scripts are provenance, not current build entry points.

The old R1 manifest is checked by canonical-JSON SHA-256 `0a995d15171f3b133edaf203869329e9eec7913a9ce8e23e8469175bc7538bfa`. Its numerical parameters and old forecasts were not retuned. The new two-gate manifest has a separate raw-byte SHA-256 and LF normalization. It freezes structural choices, not a fully specified survey analysis. Never run a prediction-freezing command against an existing manifest simply to refresh it.

## LaTeX

Install a working TeX Live or MiKTeX distribution with latexmk, BibTeX, Latin Modern, amsmath, natbib, graphicx, hyperref and booktabs. From `manuscript`:

```powershell
latexmk -C
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The source uses the generic article class, lowercase title and maketitle. The PDF is `manuscript/main.pdf`. No force flag is used. Inspect `main.log` for undefined commands/references/citations, natbib/fatal warnings and overfull boxes. The generic numbered bibliography is acceptable for the journal's free-format route; final MDPI ACS styling is a separate formatting step.

Render with `pdftoppm -r 90 -png manuscript/main.pdf manuscript/build/page` for visual QA. Do not treat successful compilation as proof of readable equations and tables.

## Submission source archive

After the final PDF build, run `python scripts/package_submission.py` from the repository root. It produces `submission/universe-source.zip`, containing only the active LaTeX input graph, bibliography, figures, PDF, build instructions and a per-file SHA-256 manifest. The ZIP is deterministic for fixed input bytes and is checked after writing. Extract into a fresh directory and run latexmk there to verify the archive is self-contained.

## External observations and missing provenance

Raw observations are **not bundled**. See [external-input map](../external/README.md). Foreground scripts accept explicit `--sn` and `--gal` paths; their tests use controlled synthetic inputs. The old result JSONs do not contain full raw/derived input checksums. The historical temporary photometric CSV and original cleaned SN/galaxy files must be recovered and hashed before claiming byte-reproducible observational reruns. Do not silently substitute a later catalog release.

Example, only after the exact frozen inputs are recovered:

```powershell
python code/r1_2mrs_foreground_transfer.py --sn external/PantheonPlus_frozen.csv --gal external/2MRS_frozen.csv --permutations 999 --artifact manuscript/build/number-density.json --report manuscript/build/number-density.md
python code/r1_2mrs_luminosity_transfer.py --sn external/PantheonPlus_frozen.csv --gal external/2MRS_photometry_frozen.csv --permutations 999 --artifact manuscript/build/luminosity.json
python code/r1_2mrs_bhsm_seam_transfer.py --sn external/PantheonPlus_frozen.csv --gal external/2MRS_photometry_frozen.csv --permutations 999 --artifact manuscript/build/seam.json
```

These filenames are placeholders, not supplied files. The carried-forward SN amplitude/axis likelihood is not implemented as a fresh fit here. The prior homogeneous localization failure is reported in documentation; its original numerical statistic/threshold are absent. These limitations are explicitly retained in the hostile-referee audit.
