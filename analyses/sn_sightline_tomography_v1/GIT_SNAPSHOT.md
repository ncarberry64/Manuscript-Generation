# Git snapshot and full reproduction package

This directory preserves the exact analysis scripts, final report, compact result/sightline tables, figures, and provenance. It is a review snapshot, not the complete data-bearing working directory. The full archive is identified by `PACKAGE_REFERENCE.json` and `ARCHIVE_SHA256.txt`.

The 1.16 GB full reproducible package is stored locally at `C:\Users\carbe\Downloads\BHSM_SN_SIGHTLINE_TOMOGRAPHY_V1`. It includes all downloaded public inputs, epoch records, dense covariance/GLS operators, sparse tomography matrices and mock/permutation arrays. Extract that archive and run `python run_all.py` **inside the extracted package**, after installing its pinned requirements. Do not run the numerical stages in this lightweight Git snapshot without the complete inputs.

`RESULT_MANIFEST.json` describes the full delivered package; it is not a list of files claimed to be committed here. `GIT_SNAPSHOT_MANIFEST.json` records the subset committed here. No public permanent archive is claimed. The superseded development ZIP is retained only in the complete local archive and is not scientific evidence.
