from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import r1_2mrs_foreground_transfer as base
import bhsm_los_topographic_transfer as los


K_CANDIDATES = [
    "k_m_k20fe",
    "K_m_k20fe",
    "k_m_ext",
    "K_m_ext",
    "k_m",
    "K_m",
    "Kmag",
    "kmag",
    "K_MAG",
    "K",
    "Ks",
    "K_s",
]


def find_kmag(df: pd.DataFrame) -> str:
    lower = {str(c).lower(): c for c in df.columns}
    for c in K_CANDIDATES:
        if c in df.columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    likely = [c for c in df.columns if ("k" in str(c).lower() and "mag" in str(c).lower())]
    raise KeyError(
        "No supported K-band magnitude column found. "
        f"Likely columns: {likely}. All columns: {list(df.columns)}"
    )


def luminosity_weight(mag: np.ndarray, z: np.ndarray) -> np.ndarray:
    # Relative luminosity is enough because each shell map is standardized.
    # d_L in arbitrary H0 units; overall constants cancel.
    z = np.asarray(z,float)
    mag = np.asarray(mag,float)
    dl = np.array([
        (1.0+zz) * los.chi_of_z(float(zz))
        for zz in z
    ],dtype=float)
    flux = 10.0 ** (-0.4 * mag)
    return flux * np.maximum(dl,1e-9)**2


def weighted_gaussian_sky(sn_u, gal_u, weights, sigma_deg=10.0, chunk=3000):
    sigma = math.radians(float(sigma_deg))
    out = np.zeros(sn_u.shape[0],dtype=float)
    weights=np.asarray(weights,float)
    for j0 in range(0,gal_u.shape[0],chunk):
        g=gal_u[j0:j0+chunk]
        w=weights[j0:j0+chunk]
        dot=np.clip(sn_u @ g.T,-1.0,1.0)
        theta=np.arccos(dot)
        out += (np.exp(-0.5*(theta/sigma)**2) * w[None,:]).sum(axis=1)
    return out


def shell_luminosity_contrasts(sn_ra,sn_dec,gal,kmag_col):
    sn_u=base.radec_to_unit(sn_ra,sn_dec)
    arr=[]
    meta=[]
    for zlo,zhi in zip(base.SHELL_EDGES[:-1],base.SHELL_EDGES[1:]):
        g=gal[(gal["_z"]>=zlo)&(gal["_z"]<zhi)].copy()
        g["_kmag"]=pd.to_numeric(g[kmag_col],errors="coerce")
        g=g.dropna(subset=["_kmag"])
        if len(g)<20:
            raise ValueError(f"Too few galaxies with K magnitude in shell {zlo}-{zhi}: {len(g)}")
        wt=luminosity_weight(g["_kmag"].to_numpy(),g["_z"].to_numpy())
        med=np.nanmedian(wt)
        if not np.isfinite(med) or med<=0:
            raise ValueError(f"Invalid shell luminosity normalization {zlo}-{zhi}")
        # Normalize weights by shell median to avoid arbitrary numerical scaling.
        wt=wt/med
        gu=base.radec_to_unit(g["_ra"].to_numpy(),g["_dec"].to_numpy())
        rho=weighted_gaussian_sky(sn_u,gu,wt)
        mu=float(np.mean(rho)); sd=float(np.std(rho,ddof=1))
        if sd<=0 or not np.isfinite(sd):
            raise ValueError("Degenerate luminosity map")
        arr.append((rho-mu)/sd)
        meta.append({
            "zlo":float(zlo),
            "zhi":float(zhi),
            "n_galaxies_with_kmag":int(len(g)),
            "median_relative_luminosity_weight":1.0,
        })
    return np.column_stack(arr),meta


def analyze(sn_path,gal_path,n_perm):
    sn=pd.read_csv(sn_path)
    gal=pd.read_csv(gal_path)

    zc=base.find_column(sn,["zCMB","zcmb","zHD","z","redshift"])
    rac=base.find_column(sn,["RA","ra","RA_DEG"])
    decc=base.find_column(sn,["DEC","dec","DEC_DEG","DECL"])
    muc=base.find_column(sn,["MU_SH0ES","MU","mu","distance_modulus"])
    ec=base.find_column(sn,["MU_SH0ES_ERR_DIAG","MUERR","muerr","MU_ERR","distance_modulus_err"])

    gz=base.find_column(gal,["z","zcmb","redshift","Z"])
    gra=base.find_column(gal,["RA","ra","RA_DEG"])
    gdec=base.find_column(gal,["DEC","dec","DEC_DEG","DECL"])
    kmag=find_kmag(gal)

    sn=sn.copy()
    sn["_z"]=pd.to_numeric(sn[zc],errors="coerce")
    sn["_ra"]=pd.to_numeric(sn[rac],errors="coerce")
    sn["_dec"]=pd.to_numeric(sn[decc],errors="coerce")
    sn["_mu"]=pd.to_numeric(sn[muc],errors="coerce")
    sn["_err"]=pd.to_numeric(sn[ec],errors="coerce")
    sn=sn.replace([np.inf,-np.inf],np.nan).dropna(subset=["_z","_ra","_dec","_mu","_err"])
    sn=sn[(sn["_z"]>=base.SN_Z_MIN)&(sn["_z"]<=base.SN_Z_MAX)&(sn["_err"]>0)].reset_index(drop=True)

    gal=gal.copy()
    gal["_z"]=pd.to_numeric(gal[gz],errors="coerce")
    gal["_ra"]=pd.to_numeric(gal[gra],errors="coerce")
    gal["_dec"]=pd.to_numeric(gal[gdec],errors="coerce")
    gal=gal.replace([np.inf,-np.inf],np.nan).dropna(subset=["_z","_ra","_dec"])
    gal=gal[(gal["_z"]>=base.SHELL_EDGES[0])&(gal["_z"]<base.SHELL_EDGES[-1])].reset_index(drop=True)

    bfit=base.fit_global_baseline(
        sn["_z"].to_numpy(),sn["_mu"].to_numpy(),sn["_err"].to_numpy()
    )
    mu0=base.mu_flat(sn["_z"].to_numpy(),70.0,bfit["omega_m"])+bfit["deltaM"]
    resid=sn["_mu"].to_numpy()-mu0

    delta,shell_meta=shell_luminosity_contrasts(
        sn["_ra"].to_numpy(),sn["_dec"].to_numpy(),gal,kmag
    )
    raw=base.foreground_raw_score(sn["_z"].to_numpy(),delta)
    score,mean,sd=base.standardize(raw)
    controls=base.build_controls(sn)

    null,full,dchi2,pval,perms=base.blocked_permutation(
        sn["_z"].to_numpy(),score,resid,sn["_err"].to_numpy(),controls,n_perm
    )
    beta=full["coeff"]["foreground"]
    berr=full["stderr"]["foreground"]
    t=beta/berr if berr>0 else float("nan")

    jack=[]
    for leave in range(delta.shape[1]):
        dtmp=delta.copy()
        dtmp[:,leave]=0.0
        rawj=base.foreground_raw_score(sn["_z"].to_numpy(),dtmp)
        sj,_,_=base.standardize(rawj)
        Xj,nj=base.design_matrix(controls,sj)
        fj=base.wls(resid,sn["_err"].to_numpy(),Xj,nj)
        bj=fj["coeff"]["foreground"]
        jack.append({
            "left_out_shell":[float(base.SHELL_EDGES[leave]),float(base.SHELL_EDGES[leave+1])],
            "beta_foreground":float(bj),
            "sign_matches_primary":bool(np.sign(bj)==np.sign(beta) or bj==0 or beta==0)
        })
    stable=all(x["sign_matches_primary"] for x in jack)

    if pval<0.01 and stable:
        classification="strong"
    elif pval<0.05:
        classification="provisional"
    else:
        classification="not_established"

    transfer=base.low_to_high_transfer(
        sn["_z"].to_numpy(),score,resid,sn["_err"].to_numpy(),controls
    )

    return {
        "schema":"r1-2mrs-luminosity-transfer-result-v1",
        "classification":classification,
        "k_magnitude_column":kmag,
        "data":{"n_sn":int(len(sn)),"n_galaxies":int(len(gal))},
        "baseline":bfit,
        "foreground":{
            "score_mean_raw":float(mean),
            "score_sd_raw":float(sd),
            "shells":shell_meta,
            "weight_law":"10^(-0.4 m_K) d_L^2, normalized by shell median"
        },
        "primary":{
            "beta_foreground_mag_per_score_sigma":float(beta),
            "beta_foreground_err":float(berr),
            "t_foreground":float(t),
            "delta_chi2":float(dchi2),
            "blocked_permutation_p":float(pval),
            "sign_stable_all_shell_dropouts":bool(stable),
            "leave_one_shell_out":jack,
        },
        "low_to_high_transfer":transfer,
        "claim_boundary":[
            "K-band luminosity is a stellar-mass/gravitational-well proxy, not an SMBH/de-encapsulation measurement.",
            "No K-band weighting exponent was fit.",
            "A positive result still requires separation from conventional local-structure effects."
        ]
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sn",required=True)
    ap.add_argument("--gal",required=True)
    ap.add_argument("--permutations",type=int,default=999)
    ap.add_argument("--artifact",required=True)
    args=ap.parse_args()
    r=analyze(args.sn,args.gal,args.permutations)
    Path(args.artifact).write_text(json.dumps(r,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(r,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
