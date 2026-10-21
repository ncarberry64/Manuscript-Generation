"""Execute the archived IPC inverse and conditional low-charge photon map.

No fit, native processing, uncertainty inflation, mean reconvolution, or new
science-data acquisition. The actual already-computed pilot means are inputs.
"""
from pathlib import Path
import hashlib
import json
import sys
import argparse
import numpy as np
from scipy.signal import convolve2d
from scipy import sparse
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/"data/strong_lensing_response_completion_20261021"
OUT = ROOT/"artifacts/strong_lensing_response_completion_20261021/ipc"
sys.path.insert(0, str(ROOT/"code"))
from strong_lensing_ipc_response import (
    inverse_deconvolution_kernel, convolution_matrix, photon_to_raw_dn_jacobian)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(match_native_early=False):
    OUT.mkdir(exist_ok=True, parents=True)
    ref = BASE/"native_cache/jwst_nircam_ipc_0022.fits"
    g = np.asarray(fits.getdata(ref), float)
    inversion = inverse_deconvolution_kernel(g)
    h = inversion.pop("kernel")
    h33 = convolution_matrix(h, (33, 33))
    roi_indices = np.arange(33**2).reshape(33, 33)[8:25, 8:25].ravel()
    hroi = h33[roi_indices]
    replay = convolve2d(g, h)
    impulse = np.zeros_like(replay); impulse[9, 9] = 1.
    rows = []
    for e in range(1, 5):
        stem = f"jw01198007001_02101_0000{e}_nrcb3"
        hp = BASE/f"histories/{stem}_histories.npz"
        pp = ROOT/f"artifacts/strong_lensing_response_completion_20261021/common_early/pilot_means/training_plus_held_{stem}_pilot.npz"
        with np.load(hp) as native, np.load(pp) as pilot:
            times = native["native_read_times_s"][:2].ravel()
            interval = float(np.diff(native["native_read_times_s"][:2].mean(axis=1))[0])
            denominator = float(native["group_time_header_s"]) if match_native_early else interval
            physical_rate_conversion = denominator/interval if match_native_early else 1.
            weights = np.array([-.5, -.5, .5, .5])/denominator
            time_coefficient = float(weights@np.minimum.outer(times, times)@weights)
            saved = {"reference_deconvolution_G": g, "conditional_forward_H_order8": h,
                "individual_read_times_s": times, "early_rate_read_weights_per_s": weights,
                "early_poisson_time_coefficient_per_s": time_coefficient,
                "pilot_reported_rate_to_physical_rate": physical_rate_conversion,
                "early_channel_denominator_s": denominator}
            summaries = []
            for label_index, label in enumerate(native["labels"]):
                mu = pilot["point_plus_sky_DN_s"][label_index]
                # This is the low-charge effective-DN plane. G is not applied
                # to an already IPC-free optical PSF, nor is H convolved twice.
                before = convolve2d(mu, g, mode="same", boundary="fill")
                gain = native["gain_e_per_DN"][label_index]
                lam = before*gain*physical_rate_conversion
                pre33 = lam[16:49, 16:49]
                d = native["dark_RAPID42_DN"][label_index]
                stats = {"label": str(label), "pre_IPC_rate_min_inner17_e_s": float(lam[24:41,24:41].min()),
                    "pre_IPC_negative_inner17_pixels": int((lam[24:41,24:41] < 0).sum()),
                    "pre_IPC_negative_inner33_pixels": int((pre33 < 0).sum()),
                    "pre_IPC_rate_min_inner33_e_s": float(pre33.min()),
                    "effective_rate_negative_inner33_pixels": int((mu[16:49,16:49] < 0).sum()),
                    "pre_IPC_negative_inner63_pixels": int((lam[1:-1,1:-1] < 0).sum()),
                    "native_dark_nonzero_inner17_pixels": int(np.any(d[:,24:41,24:41] != 0, axis=0).sum())}
                saved[f"{label}_pre_IPC_pilot_e_s_inner63"] = lam[1:-1,1:-1]
                if np.min(pre33) >= 0:
                    # Low-charge Linv'=1. H acts in raw DN AFTER source gain;
                    # native early linear read difference follows on the left.
                    j = photon_to_raw_dn_jacobian(hroi, gain[16:49,16:49], np.ones((33,33)))
                    cov_dn = (j@sparse.diags(pre33.ravel())@j.T).toarray()*time_coefficient
                    sci_factor = (native["photmjsr"]/native["flat"][label_index,24:41,24:41]).ravel()
                    cov = sci_factor[:,None]*cov_dn*sci_factor[None,:]
                    saved[f"{label}_low_charge_early_photon_covariance_SCI2"] = cov
                    eig = np.linalg.eigvalsh(cov)
                    diag_independent = (mu[24:41,24:41]/gain[24:41,24:41]).ravel()*physical_rate_conversion*time_coefficient*sci_factor**2
                    ratio = np.diag(cov)/diag_independent
                    correlation = cov/np.sqrt(np.diag(cov)[:,None]*np.diag(cov)[None,:])
                    np.fill_diagonal(correlation, 0.)
                    stats["conditional_photon_application"] = {"status": "computed low-charge constant-reference Poisson photon projection",
                        "dimension": 289, "minimum_eigenvalue_SCI2": float(eig.min()),
                        "maximum_offdiagonal_correlation": float(np.max(abs(correlation))),
                        "variance_ratio_to_independent_effective_rate_min_max": [float(ratio.min()),float(ratio.max())],
                        "source_complete_or_validated_instrument_covariance": False}
                else:
                    stats["conditional_photon_application"] = {"status": "pre_IPC_context_not_a_positive_photon_state",
                        "reason": "G applied to the delivered effective field template gives negative rates on104 inner33 pixels; no clipping or Gaussian replacement used",
                        "required_relation": "positive pre_IPC optical/context mean whose H-image matches the effective field response, with boundary/support prescribed"}
                summaries.append(stats)
            suffix = "conditional_ipc_native_early_projection" if match_native_early else "conditional_ipc_projection"
            destination = OUT/f"{stem}_{suffix}.npz"
            np.savez_compressed(destination, **saved)
            rows.append({"exposure": stem, "pilot_file": str(pp.relative_to(ROOT)).replace("\\","/"),
                "pilot_sha256": sha(pp), "history_file": str(hp.relative_to(ROOT)).replace("\\","/"),
                "history_sha256": sha(hp), "array_file": str(destination.relative_to(ROOT)).replace("\\","/"),
                "array_sha256": sha(destination), "individual_read_times_s": times.tolist(),
                "early_channel_denominator_s": denominator,
                "pilot_reported_rate_to_physical_rate": physical_rate_conversion,
                "early_poisson_time_coefficient_per_s": time_coefficient,"patches": summaries})
    headers = {}
    for name in ("jwst_nircam_ipc_0022.fits", "jwst_nircam_linearity_0054.fits",
                 "derived_jwst_nircam_dark_0337_first42.fits"):
        header = fits.getheader(BASE/"native_cache"/name)
        headers[name] = {"cards": [str(c) for c in header.cards],
            "science_array_not_opened_for_header_audit": True}
    result = {"status": "archived reference inverse and actual pilot photon projection executed",
        "reference_file": ref.name, "reference_sha256": sha(ref), "reference_bytes": ref.stat().st_size,
        "inverse": inversion,"reference_G": g.tolist(), "inverse_H_center": float(h[8,8]),
        "G_H_identity_max_abs": float(abs(replay-impulse).max()),
        "reference_sum_differs_from_one_by": float(g.sum()-1.),
        "physical_guard": 2e-10, "reference_renormalized": False,
        "calibration_source_relation": "linearity0054 S_IPC=COMPLETE,S_REFPIX=COMPLETE,S_BIAS=COMPLETE; HISTORY explicitly applies them before fitting native L",
        "physical_order": "Q_e/g -> L_cal^-1 -> H(raw-DN IPC) -> individual-read/group-average -> actual bias/reference channel -> L_native -> native OLS -> /flat*PHOTMJSR",
        "low_charge_optical_state_relation": "lambda_pre_e_s=diag(gain)*G*mu_effective_DN_s; H must not reconvolve the effective PSF mean",
        "temporal_channel": "matching native early two-group product with rounded headerTGROUP denominator and physical-time pilot conversion" if match_native_early else "explicit unrounded physical-interval diagnostic; not the rounded native early-rate denominator",
        "temporal_conversion_relation": "for matching native early product: lambda_physical=(headerTGROUP/delta_mean_readtime)*diag(gain)*G*mu_reported_DN_s, while read weights divide by headerTGROUP; recorded individual times remain unrounded",
        "photon_covariance_relation": "Jraw=H*diag((L_cal^-1)'/gain); Cphoton=J_native*Jraw*Ccollection*Jraw.T*J_native.T. Here only the stated low-charge independent-Poisson/first-two-group channel is applied",
        "dark_source_relation": "dark0337 contains thresholded post-linearity slope-times-time at hot pixels and zeros at nominal pixels; header has no explicit S_IPC/S_REFPIX provenance. Its pre_IPC high-charge relation is therefore conditional, not derived from absent keywords",
        "electronic_noise_accounting": "empirical post-SIRS calibration dark covariance already includes actual electronics/IPC; no second H or independent CDS+reference variance is added",
        "not_identified": ["flight signal-dependent/nonlinear IPC and its uncertainty", "pre_IPC field-template support/positivity", "hot-dark pre_IPC/readwise relation", "full nonlinear/correlated photon-count law from a static ground correction kernel"],
        "headers": headers, "records": rows,
        "sources": ["https://jwst-pipeline.readthedocs.io/en/stable/jwst/ipc/description.html",
            "https://jwst-pipeline.readthedocs.io/en/stable/jwst/ipc/reference_files.html",
            "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_nircam_ipc_0022.fits"],
        "source_sha256": {str(p.relative_to(ROOT)).replace("\\","/"): sha(p) for p in (
            Path(__file__),ROOT/"code/strong_lensing_ipc_response.py",ROOT/"tests/test_strong_lensing_ipc_response.py")},
        "native_pipeline_rerun": False, "new_optical_fit": False, "noise_rescaled": False,
        "physical_validation": False, "no_reserved_science_observations": True}
    name = "ipc_native_early_projection_relation.json" if match_native_early else "ipc_projection_relation.json"
    (OUT/name).write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("IPC_PROJECTION_DONE",json.dumps(inversion),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--match-native-early-channel", action="store_true")
    main(parser.parse_args().match_native_early_channel)
