from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import r1_2mrs_foreground_transfer as base
import r1_2mrs_luminosity_transfer as lum
import bhsm_los_topographic_transfer as los


def seam_charge_weight(mag: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Relative V_Sigma^2 proxy from the BHSM spherical seam scaling.

    Frozen chain:
      M_b âˆ L_K,
      V_Sigma^2 âˆ sqrt(M_b),
      L_K âˆ 10^(-0.4 m_K) d_L^2,

    therefore
      V_Sigma^2 âˆ 10^(-0.2 m_K) d_L.
    """
    mag = np.asarray(mag, dtype=float)
    z = np.asarray(z, dtype=float)
    dl = np.array(
        [(1.0 + zz) * los.chi_of_z(float(zz)) for zz in z],
        dtype=float,
    )
    return (10.0 ** (-0.2 * mag)) * np.maximum(dl, 1.0e-9)


def weighted_gaussian_sky(
    sn_u: np.ndarray,
    gal_u: np.ndarray,
    weights: np.ndarray,
    sigma_deg: float = 10.0,
    chunk: int = 3000,
) -> np.ndarray:
    sigma = math.radians(float(sigma_deg))
    weights = np.asarray(weights, dtype=float)
    out = np.zeros(sn_u.shape[0], dtype=float)
    for j0 in range(0, gal_u.shape[0], chunk):
        g = gal_u[j0:j0+chunk]
        w = weights[j0:j0+chunk]
        dot = np.clip(sn_u @ g.T, -1.0, 1.0)
        theta = np.arccos(dot)
        out += (
            np.exp(-0.5 * (theta / sigma) ** 2)
            * w[None, :]
        ).sum(axis=1)
    return out


def shell_seam_contrasts(sn_ra, sn_dec, gal, kmag_col):
    sn_u = base.radec_to_unit(sn_ra, sn_dec)
    fields = []
    meta = []

    for zlo, zhi in zip(
        base.SHELL_EDGES[:-1],
        base.SHELL_EDGES[1:],
    ):
        g = gal[
            (gal["_z"] >= zlo)
            & (gal["_z"] < zhi)
        ].copy()

        g["_kmag"] = pd.to_numeric(
            g[kmag_col], errors="coerce"
        )
        g = g.dropna(subset=["_kmag"])

        if len(g) < 20:
            raise ValueError(
                f"Too few galaxies in shell {zlo}-{zhi}"
            )

        wt = seam_charge_weight(
            g["_kmag"].to_numpy(),
            g["_z"].to_numpy(),
        )

        med = float(np.nanmedian(wt))
        if not np.isfinite(med) or med <= 0:
            raise ValueError(
                f"Invalid seam-weight median in {zlo}-{zhi}"
            )

        # Only numerical normalization; no residual information.
        wt = wt / med

        gu = base.radec_to_unit(
            g["_ra"].to_numpy(),
            g["_dec"].to_numpy(),
        )
        rho = weighted_gaussian_sky(
            sn_u, gu, wt, sigma_deg=10.0
        )

        mu = float(np.mean(rho))
        sd = float(np.std(rho, ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(
                f"Degenerate seam proxy in {zlo}-{zhi}"
            )

        fields.append((rho - mu) / sd)
        meta.append(
            {
                "zlo": float(zlo),
                "zhi": float(zhi),
                "n_galaxies": int(len(g)),
                "median_relative_seam_weight": 1.0,
            }
        )

    return np.column_stack(fields), meta


def analyze(sn_path: str, gal_path: str, n_perm: int) -> dict:
    sn = pd.read_csv(sn_path)
    gal = pd.read_csv(gal_path)

    zc = base.find_column(
        sn, ["zCMB","zcmb","zHD","z","redshift"]
    )
    rac = base.find_column(sn, ["RA","ra","RA_DEG"])
    decc = base.find_column(
        sn, ["DEC","dec","DEC_DEG","DECL"]
    )
    muc = base.find_column(
        sn,
        ["MU_SH0ES","MU","mu","distance_modulus"],
    )
    ec = base.find_column(
        sn,
        [
            "MU_SH0ES_ERR_DIAG",
            "MUERR",
            "muerr",
            "MU_ERR",
            "distance_modulus_err",
        ],
    )

    gz = base.find_column(
        gal, ["z","zcmb","redshift","Z"]
    )
    gra = base.find_column(gal, ["RA","ra","RA_DEG"])
    gdec = base.find_column(
        gal, ["DEC","dec","DEC_DEG","DECL"]
    )
    kmag = lum.find_kmag(gal)

    sn = sn.copy()
    sn["_z"] = pd.to_numeric(sn[zc], errors="coerce")
    sn["_ra"] = pd.to_numeric(sn[rac], errors="coerce")
    sn["_dec"] = pd.to_numeric(sn[decc], errors="coerce")
    sn["_mu"] = pd.to_numeric(sn[muc], errors="coerce")
    sn["_err"] = pd.to_numeric(sn[ec], errors="coerce")
    sn = (
        sn.replace([np.inf,-np.inf], np.nan)
        .dropna(
            subset=["_z","_ra","_dec","_mu","_err"]
        )
    )
    sn = sn[
        (sn["_z"] >= base.SN_Z_MIN)
        & (sn["_z"] <= base.SN_Z_MAX)
        & (sn["_err"] > 0)
    ].reset_index(drop=True)

    gal = gal.copy()
    gal["_z"] = pd.to_numeric(gal[gz], errors="coerce")
    gal["_ra"] = pd.to_numeric(gal[gra], errors="coerce")
    gal["_dec"] = pd.to_numeric(gal[gdec], errors="coerce")
    gal = (
        gal.replace([np.inf,-np.inf], np.nan)
        .dropna(subset=["_z","_ra","_dec"])
    )
    gal = gal[
        (gal["_z"] >= base.SHELL_EDGES[0])
        & (gal["_z"] < base.SHELL_EDGES[-1])
    ].reset_index(drop=True)

    baseline = base.fit_global_baseline(
        sn["_z"].to_numpy(),
        sn["_mu"].to_numpy(),
        sn["_err"].to_numpy(),
    )
    mu0 = (
        base.mu_flat(
            sn["_z"].to_numpy(),
            70.0,
            baseline["omega_m"],
        )
        + baseline["deltaM"]
    )
    residual = sn["_mu"].to_numpy() - mu0

    shell_delta, shell_meta = shell_seam_contrasts(
        sn["_ra"].to_numpy(),
        sn["_dec"].to_numpy(),
        gal,
        kmag,
    )

    raw = base.foreground_raw_score(
        sn["_z"].to_numpy(),
        shell_delta,
    )
    score, score_mean, score_sd = base.standardize(raw)

    controls = base.build_controls(sn)

    null, full, dchi2, pval, _ = base.blocked_permutation(
        sn["_z"].to_numpy(),
        score,
        residual,
        sn["_err"].to_numpy(),
        controls,
        n_perm,
    )

    beta = full["coeff"]["foreground"]
    beta_err = full["stderr"]["foreground"]
    tval = beta / beta_err if beta_err > 0 else float("nan")

    jackknife = []
    for leave in range(shell_delta.shape[1]):
        dtmp = shell_delta.copy()
        dtmp[:, leave] = 0.0
        raw_j = base.foreground_raw_score(
            sn["_z"].to_numpy(),
            dtmp,
        )
        score_j, _, _ = base.standardize(raw_j)
        Xj, namesj = base.design_matrix(
            controls, score_j
        )
        fitj = base.wls(
            residual,
            sn["_err"].to_numpy(),
            Xj,
            namesj,
        )
        bj = fitj["coeff"]["foreground"]
        jackknife.append(
            {
                "left_out_shell": [
                    float(base.SHELL_EDGES[leave]),
                    float(base.SHELL_EDGES[leave+1]),
                ],
                "beta_foreground": float(bj),
                "sign_matches_primary": bool(
                    np.sign(bj) == np.sign(beta)
                    or bj == 0
                    or beta == 0
                ),
            }
        )

    stable = all(
        row["sign_matches_primary"]
        for row in jackknife
    )

    if pval < 0.01 and stable:
        classification = "strong"
    elif pval < 0.05:
        classification = "provisional"
    else:
        classification = "not_established"

    transfer = base.low_to_high_transfer(
        sn["_z"].to_numpy(),
        score,
        residual,
        sn["_err"].to_numpy(),
        controls,
    )

    return {
        "schema": "r1-2mrs-bhsm-seam-transfer-result-v1",
        "classification": classification,
        "data": {
            "n_sn": int(len(sn)),
            "n_galaxies": int(len(gal)),
            "k_magnitude_column": kmag,
        },
        "source_chain": {
            "mass_proxy": "M_b proportional to L_K",
            "seam_scaling": "V_Sigma^2 proportional to sqrt(M_b)",
            "object_weight": "10^(-0.2 m_K) d_L",
            "weight_exponent_fitted": False,
        },
        "foreground": {
            "score_mean_raw": float(score_mean),
            "score_sd_raw": float(score_sd),
            "shells": shell_meta,
        },
        "primary": {
            "beta_foreground_mag_per_score_sigma": float(beta),
            "beta_foreground_err": float(beta_err),
            "t_foreground": float(tval),
            "delta_chi2": float(dchi2),
            "blocked_permutation_p": float(pval),
            "sign_stable_all_shell_dropouts": bool(stable),
            "leave_one_shell_out": jackknife,
        },
        "low_to_high_transfer": transfer,
        "claim_boundary": [
            "This is a BHSM seam-charge proxy, not an action-normalized de-encapsulation map.",
            "K-band luminosity is only a baryonic-mass proxy.",
            "No SMBH mass or accretion information is used.",
        ],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sn", required=True)
    ap.add_argument("--gal", required=True)
    ap.add_argument(
        "--permutations",
        type=int,
        default=999,
    )
    ap.add_argument("--artifact", required=True)
    args = ap.parse_args()

    result = analyze(
        args.sn,
        args.gal,
        args.permutations,
    )
    Path(args.artifact).write_text(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        ) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
