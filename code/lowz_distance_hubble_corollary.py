from __future__ import annotations

import json
import math


A_MU = -0.0412
SIGMA_A_MU = 0.009
H0_REF = 67.4


def fractional_distance_from_mag(A_mu: float) -> float:
    return (math.log(10.0) / 5.0) * A_mu


def fractional_distance_sigma(sigma_A_mu: float) -> float:
    return (math.log(10.0) / 5.0) * sigma_A_mu


def diagnostics() -> dict:
    dL = fractional_distance_from_mag(A_MU)
    sdL = fractional_distance_sigma(SIGMA_A_MU)

    # At fixed sufficiently low z, D_L ~= cz/H0, so
    # delta D_L/D_L = -delta H0/H0 + O(z).
    dH = -dL
    sdH = sdL

    return {
        "status": "LOW_Z_DISTANCE_HUBBLE_COROLLARY_EVALUATED",
        "A_mu": A_MU,
        "A_mu_sigma": SIGMA_A_MU,
        "delta_DL_over_DL": dL,
        "delta_DL_over_DL_sigma": sdL,
        "delta_DA_over_DA": dL,
        "delta_DA_over_DA_sigma": sdL,
        "delta_H0_over_H0": dH,
        "delta_H0_over_H0_sigma": sdH,
        "H0_reference_km_s_Mpc": H0_REF,
        "delta_H0_km_s_Mpc": dH * H0_REF,
        "delta_H0_km_s_Mpc_sigma": sdH * H0_REF,
        "claim_boundary": (
            "D_A equality follows from Etherington reciprocity at fixed "
            "observed redshift; H0 relation is the z->0 cosmographic limit, "
            "not the full-redshift H(z) kernel"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(diagnostics(), indent=2, sort_keys=True))
