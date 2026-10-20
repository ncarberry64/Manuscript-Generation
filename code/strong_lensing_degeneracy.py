"""Analytic scalar-map nuisance factorization; audit saved comparison arrays.

This module does not rerun cosmological propagation, change event roots, or
apply the proposition to an additional intrinsic local encounter. Observable
coordinates and nuisance normalization are inherited from the saved comparison.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


MINIMAL_NUISANCES = (
    "source_transverse_0_per_arcsec",
    "log_galaxy_deflection_strength",
    "log_galaxy_core_radius",
    "independent_mass_sheet_log_scale",
)


def scalar_map_coefficients(distances, map_response, beta_arcsec, *, core):
    """Return R in B=N_min R for the axial softened-isothermal family.

    Distances and derivatives use the same length units. The source column is
    per arcsecond; all other selected nuisance coordinates are dimensionless.
    The state changes scalar maps only, at fixed measured lens redshift,
    dimensional transverse source, physical deflection strength and core.
    """
    distances = np.asarray(distances, dtype=float)
    response = np.asarray(map_response, dtype=float)
    if distances.shape != (3,) or response.ndim != 2 or response.shape[0] != 3:
        raise ValueError("three scalar distances and a three-row response required")
    if not np.isfinite(distances).all() or np.any(distances <= 0):
        raise ValueError("positive finite scalar distances required")
    if not np.isfinite(response).all() or not np.isfinite(beta_arcsec):
        raise ValueError("finite response and source coordinate required")
    f_ol, f_os, f_ls = response / distances[:, None]
    rows = [-float(beta_arcsec)*f_os, f_ls-f_os]
    if core:
        rows.append(-f_ol)
    rows.append(f_ol+f_os-f_ls)
    return np.vstack(rows)


def factorization_check(B, N_min, R, covariance):
    """Check factorization and orthogonal projection for an SPD covariance."""
    B, N_min, R, covariance = (
        np.asarray(value, dtype=float) for value in (B, N_min, R, covariance)
    )
    if B.ndim != 2 or N_min.ndim != 2 or R.ndim != 2:
        raise ValueError("matrix responses required")
    if N_min.shape[0] != B.shape[0] or R.shape != (N_min.shape[1], B.shape[1]):
        raise ValueError("matching observable, nuisance and state dimensions required")
    if covariance.shape != (B.shape[0], B.shape[0]):
        raise ValueError("matching observational covariance required")
    if not all(np.isfinite(value).all() for value in (B, N_min, R)):
        raise ValueError("finite response matrices required")
    if not np.isfinite(covariance).all() or not np.allclose(covariance, covariance.T, atol=1e-14, rtol=0):
        raise ValueError("symmetric finite observational covariance required")
    L = np.linalg.cholesky(covariance)
    WB, WN = np.linalg.solve(L, B), np.linalg.solve(L, N_min)
    difference = np.linalg.solve(L, B-N_min@R)
    scale = max(float(np.linalg.norm(WB)), np.finfo(float).tiny)
    norms = np.linalg.norm(WN, axis=0)
    keep = norms > 0
    U, singular_values, _ = np.linalg.svd(WN[:, keep]/norms[keep], full_matrices=False)
    rank = int(np.sum(singular_values > singular_values[0]*1e-10)) if keep.any() else 0
    projected = WB-U[:, :rank]@(U[:, :rank].T@WB)
    return dict(
        relative_whitened_factorization_discrepancy=float(np.linalg.norm(difference)/scale),
        maximum_absolute_unwhitened_discrepancy=float(np.max(np.abs(B-N_min@R))),
        relative_whitened_projected_response=float(np.linalg.norm(projected)/scale),
        minimal_nuisance_rank=rank,
        coefficient_matrix=R.tolist(),
        maximum_projected_whitened_singular_value=float(np.linalg.svd(projected, compute_uv=False)[0]),
    )


def saved_array_receipt(input_dir):
    """Evaluate the identity using only first-checkpoint JSON and NPZ inputs."""
    input_dir = Path(input_dir)
    result_path, array_path = (input_dir/name for name in ("results.json", "response_arrays.npz"))
    report = json.loads(result_path.read_text(encoding="utf-8"))
    background = report["transport"]["background"]
    distances = []
    for name in ("D_ol", "D_os", "D_ls"):
        matrix = np.asarray(background[name], dtype=float)
        if matrix.shape != (2, 2) or not np.allclose(matrix, np.eye(2)*matrix[0, 0], atol=1e-15):
            raise ValueError("saved proposition requires scalar screen maps")
        distances.append(matrix[0, 0])
    beta = report["experiment"]["reference_beta_arcsec"]
    if beta[1] != 0:
        raise ValueError("saved proposition uses the declared axial source")
    families = {}
    with np.load(array_path, allow_pickle=False) as arrays:
        for family, core, prefix, container in (
            ("SIS", False, "", report["identifiability"]),
            ("smooth_core", True, "smooth_", report["smooth_check"]),
        ):
            names = tuple(name for name in MINIMAL_NUISANCES if core or name != MINIMAL_NUISANCES[2])
            indices = [container["nuisance_names"].index(name) for name in names]
            B, N_min = arrays[prefix+"B"], arrays[prefix+"N"][:, indices]
            R = scalar_map_coefficients(distances, arrays["d_maps_d_anchor"], beta[0], core=core)
            if core:
                n_images = len(container["solution"]["images_rad"])
                # Same illustrative precision as the saved runner; its central
                # image is an idealized control, not an observational forecast.
                sigma = np.r_[np.full(2*n_images, .001), np.full(n_images-1, .1)]
                covariance = np.diag(sigma**2)
            else:
                covariance = arrays["covariance"]
            check = factorization_check(B, N_min, R, covariance)
            check.update(nuisance_names=names, nuisance_indices=indices,
                         observable_rows=B.shape[0], state_columns=B.shape[1])
            families[family] = check
    trial = report["transport"]["finite_trials"][0]
    affine_shift = (trial["plus_endpoints"]["source"]["affine"]
                    - trial["minus_endpoints"]["source"]["affine"])/(2*trial["step"])
    radial_shift = report["transport"]["d_events_d_anchor"][1][0]
    chi_prime = 1/background["source"]["a"]**2
    return dict(
        schema="strong-lensing-scalar-map-nuisance-factorization-v1",
        provenance_anchor="97bf5571335cf60e15e2e486e1197ca226198ed2",
        method="saved arrays only; no transport or numerical event-root change",
        scope="scalar-map SIS and softened-isothermal controls; no additional intrinsic encounter",
        input_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in (result_path, array_path)},
        scalar_distances_in_c_over_H0=distances,
        beta_arcsec=beta[0],
        families=families,
        source_endpoint_unit_zeta=dict(
            delta_chi_event=radial_shift,
            delta_lambda_event=affine_shift,
            background_chi_prime_lambda=chi_prime,
            converted_affine_displacement=chi_prime*affine_shift,
            inferred_delta_chi_ray_fixed_lambda=radial_shift-chi_prime*affine_shift,
            interpretation="last value follows from the chain rule and saved moving-event responses",
        ),
        physical_equivalence_claim=False,
        covariance_claim="P W B=0 for every positive-definite observational covariance in the exact model",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        default=Path("artifacts/strong_lensing_first_comparison_20261020"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.input_dir/"nuisance_factorization.json"
    receipt = saved_array_receipt(args.input_dir)
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    for family, check in receipt["families"].items():
        print(f"{family}: relative whitened factorization discrepancy "
              f"{check['relative_whitened_factorization_discrepancy']:.12g}")


if __name__ == "__main__":
    main()
