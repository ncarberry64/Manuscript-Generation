# Frozen SN sightline V1 update

The user requested source reconstruction, independent foreground/PV interpretation, covariance/rank diagnostics and an update of the manuscript repository without retuning the scientific prediction. The [final report](../analyses/sn_sightline_tomography_v1/FINAL_REPORT.md) is authoritative for the new observational results. The older full numerical validation remains historical evidence for its unchanged canonical scientific snapshot.

The primary held-out Pantheon result is Δχ²=−0.882838 (conditional Gaussian-mock p=0.07246); the overlap-clean DES transfer is +1.321176. The fixed independent-PV sensitivity is −4.908836, with incomplete external/source joint covariance. Neither aggregate foreground-cell localization nor a coherent cross-survey residual is established. The fitted-only scalar amplitude is −0.0328705 ± 0.0123979 mag; it never replaces the frozen predictive amplitude −0.04321623435997208 mag.

The package checks 2,619 input hashes, unchanged original full-covariance/PV manifests, exact original likelihood reproduction, covariance row-order invariance, training/test SN separation, source parameter injection, orthogonal epoch decomposition, prediction-error covariance, mock moments and path-length conservation. It includes 132,494 individual-MJD epoch deletions (132,491 identifiable), checked against 25 direct retained-epoch linear refits to 7.6e−15 maximum parameter discrepancy. All 3,521 nominal rows have a local source reconstruction; 41 non-converged fits remain flagged. No free PV normalization or raw common-amplitude nuisance predictor remains.

The full archive passed ZIP CRC verification; its SHA-256 and exact byte size are in [PACKAGE_REFERENCE.json](../analyses/sn_sightline_tomography_v1/PACKAGE_REFERENCE.json). The Git copy is deliberately a compact review snapshot. It is not a claim that raw photometry or the large covariance/mock arrays are committed here.

Repository integration validation:

- **25 focused tests passed** for the prospective gates, optical covariance and optical transfer. One non-fatal pytest cache-permission warning did not affect results.
- The canonical `code/`, `tests/` and `preregistration/` hashes exactly match the prior tested scientific snapshot. No previous 87-pass full-suite result is relabeled as a new full run.
- The **41-page manuscript builds successfully**, with no undefined references/citations, overfull boxes or LaTeX warnings. The added section, results table, calibration distinction, conclusion and data statement were visually reviewed.
- The updated generic source ZIP contains **44 verified entries** and rebuilds after cleaning in an isolated extraction directory.
- **132 Git-snapshot file hashes** match their export manifest. The 1.16 GB full archive remains outside Git and passed its separate CRC check.
- The manuscript explicitly distinguishes the historical −0.0412 mag calibration from the later frozen full-covariance amplitude; neither is retuned.

The added interpretation does not test or alter the prospective multi-observable rank-one rejection gate. A public permanent data archive, author declarations and the existing submission checklist remain open.
