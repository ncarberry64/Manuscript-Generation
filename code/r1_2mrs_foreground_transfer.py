from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import bhsm_los_topographic_transfer as los


SHELL_EDGES = np.array([0.0, 0.01, 0.02, 0.03, 0.04, 0.05], dtype=float)
SMOOTH_SIGMA_DEG = 10.0
AXIS_RA_DEG = 211.48
AXIS_DEC_DEG = -12.81
SN_Z_MIN = 0.01
SN_Z_MAX = 0.15
PERM_Z_EDGES = np.array([0.01,0.02,0.03,0.04,0.06,0.10,0.15], dtype=float)
RNG_SEED = 20260919


def find_column(df: pd.DataFrame, candidates: Iterable[str], required: bool = True):
    lower = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c in df.columns:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    if required:
        raise KeyError(f"Could not find any of {list(candidates)} in columns {list(df.columns)}")
    return None


def radec_to_unit(ra_deg, dec_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    cd = np.cos(dec)
    return np.column_stack((cd*np.cos(ra), cd*np.sin(ra), np.sin(dec)))


def axis_cosine(ra_deg, dec_deg):
    u = radec_to_unit(ra_deg, dec_deg)
    a = radec_to_unit([AXIS_RA_DEG], [AXIS_DEC_DEG])[0]
    return u @ a


def gaussian_sky_counts(sn_u, gal_u, sigma_deg=SMOOTH_SIGMA_DEG, chunk=4000):
    sigma = math.radians(float(sigma_deg))
    out = np.zeros(sn_u.shape[0], dtype=float)
    for j0 in range(0, gal_u.shape[0], chunk):
        g = gal_u[j0:j0+chunk]
        dot = np.clip(sn_u @ g.T, -1.0, 1.0)
        theta = np.arccos(dot)
        out += np.exp(-0.5 * (theta / sigma) ** 2).sum(axis=1)
    return out


def shell_density_contrasts(sn_ra, sn_dec, gal):
    sn_u = radec_to_unit(sn_ra, sn_dec)
    contrasts = []
    meta = []
    for zlo, zhi in zip(SHELL_EDGES[:-1], SHELL_EDGES[1:]):
        g = gal[(gal["_z"] >= zlo) & (gal["_z"] < zhi)]
        if len(g) < 20:
            raise ValueError(f"Too few galaxies in shell {zlo}-{zhi}: {len(g)}")
        gu = radec_to_unit(g["_ra"].to_numpy(), g["_dec"].to_numpy())
        rho = gaussian_sky_counts(sn_u, gu)
        mean = float(np.mean(rho))
        sd = float(np.std(rho, ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"Degenerate shell density map {zlo}-{zhi}")
        delta = (rho - mean) / sd
        contrasts.append(delta)
        meta.append({
            "zlo": float(zlo),
            "zhi": float(zhi),
            "n_galaxies": int(len(g)),
            "density_mean_at_sn_positions": mean,
            "density_sd_at_sn_positions": sd,
        })
    return np.column_stack(contrasts), meta


def chi_array(z):
    return np.array([los.chi_of_z(float(x)) for x in np.asarray(z, float)], dtype=float)


def foreground_raw_score(sn_z, shell_delta):
    sn_z = np.asarray(sn_z, dtype=float)
    chi_sn = chi_array(sn_z)
    chi_edges = chi_array(SHELL_EDGES)
    score = np.zeros_like(sn_z)

    for s in range(len(SHELL_EDGES)-1):
        lo = chi_edges[s]
        hi = chi_edges[s+1]
        visible = np.clip(np.minimum(chi_sn, hi) - lo, 0.0, hi-lo)
        score += shell_delta[:, s] * visible

    radial = np.array([los.C_over_S(float(c)) for c in chi_sn])
    return radial * score


def standardize(x):
    x = np.asarray(x, dtype=float)
    m = float(np.mean(x))
    s = float(np.std(x, ddof=1))
    if not np.isfinite(s) or s <= 0:
        raise ValueError("Cannot standardize degenerate predictor")
    return (x-m)/s, m, s


def mu_flat(z, h0=70.0, om=0.3):
    z = np.asarray(z, dtype=float)
    grid = np.linspace(0.0, max(0.16, float(np.max(z))), 4096)
    E = np.sqrt(om*(1.0+grid)**3 + 1.0-om)
    dz = grid[1]-grid[0]
    integ = np.zeros_like(grid)
    integ[1:] = np.cumsum(0.5*(1/E[1:]+1/E[:-1])*dz)
    chi = np.interp(z, grid, integ)
    dl = (1+z) * 299792.458/h0 * chi
    return 5*np.log10(dl)+25


def fit_global_baseline(z, mu, err):
    w = 1/np.maximum(err, 1e-6)**2
    best = None
    for om in np.linspace(0.05,0.6,221):
        base = mu_flat(z, 70.0, om)
        dm = np.sum(w*(mu-base))/np.sum(w)
        r = mu-(base+dm)
        chi2 = float(np.sum(w*r*r))
        if best is None or chi2 < best[0]:
            best = (chi2,float(om),float(dm))
    return {"chi2":best[0],"omega_m":best[1],"deltaM":best[2]}


def build_controls(df):
    controls = {"axis_cosine": axis_cosine(df["_ra"], df["_dec"])}
    host = find_column(df, ["HOST_LOGMASS","HOSTMASS","host_logmass","host_mass"], required=False)
    if host is not None:
        x = pd.to_numeric(df[host], errors="coerce").to_numpy(float)
        if np.isfinite(x).sum() > 20 and np.nanstd(x) > 0:
            fill = np.nanmedian(x)
            controls["host_mass"] = np.where(np.isfinite(x), x, fill)

    survey = find_column(df, ["IDSURVEY","SURVEY","survey","survey_id"], required=False)
    if survey is not None:
        vals = df[survey].astype(str)
        dummies = pd.get_dummies(vals, prefix="survey", drop_first=True)
        for c in dummies.columns:
            x = dummies[c].to_numpy(float)
            if x.sum() >= 3:
                controls[str(c)] = x
    return controls


def design_matrix(controls, foreground=None, mask=None):
    keys = list(controls.keys())
    n = len(next(iter(controls.values())))
    if mask is None:
        mask = np.ones(n, dtype=bool)
    cols = [np.ones(mask.sum())]
    names = ["Intercept"]
    if foreground is not None:
        cols.append(np.asarray(foreground)[mask])
        names.append("foreground")
    for k in keys:
        cols.append(np.asarray(controls[k])[mask])
        names.append(k)
    return np.column_stack(cols), names


def wls(y, err, X, names, fixed_foreground_beta=None):
    y = np.asarray(y,float)
    err = np.asarray(err,float)
    w = 1/np.maximum(err,1e-6)**2

    if fixed_foreground_beta is not None:
        j = names.index("foreground")
        yfit = y - fixed_foreground_beta*X[:,j]
        keep = [k for k in range(X.shape[1]) if k != j]
        Xfit = X[:,keep]
        namesfit = [names[k] for k in keep]
    else:
        yfit = y
        Xfit = X
        namesfit = names

    sw = np.sqrt(w)
    A = Xfit*sw[:,None]
    b = yfit*sw
    beta, *_ = np.linalg.lstsq(A,b,rcond=None)
    pred = Xfit@beta
    if fixed_foreground_beta is not None:
        pred = pred + fixed_foreground_beta*X[:,j]

    resid = y-pred
    chi2 = float(np.sum(w*resid*resid))
    dof = max(1,len(y)-Xfit.shape[1])
    sigma2 = chi2/dof
    cov = np.linalg.pinv(A.T@A)*sigma2
    se = np.sqrt(np.maximum(np.diag(cov),0))
    coeff = {n:float(v) for n,v in zip(namesfit,beta)}
    stderr = {n:float(v) for n,v in zip(namesfit,se)}
    if fixed_foreground_beta is not None:
        coeff["foreground"] = float(fixed_foreground_beta)
        stderr["foreground"] = float("nan")
    return {
        "chi2":chi2,
        "dof":dof,
        "coeff":coeff,
        "stderr":stderr,
        "prediction":pred,
        "residual":resid,
    }


def blocked_permutation(z, score, resid_mu, err, controls, n_perm):
    X0,n0 = design_matrix(controls, None)
    X1,n1 = design_matrix(controls, score)
    null = wls(resid_mu,err,X0,n0)
    full = wls(resid_mu,err,X1,n1)
    observed = null["chi2"]-full["chi2"]

    blocks = np.digitize(z, PERM_Z_EDGES[1:-1], right=False)
    rng = np.random.default_rng(RNG_SEED)
    vals = []
    for _ in range(int(n_perm)):
        p = np.asarray(score).copy()
        for b in np.unique(blocks):
            idx = np.flatnonzero(blocks==b)
            if len(idx) >= 2:
                p[idx] = p[rng.permutation(idx)]
        Xp,npn = design_matrix(controls,p)
        fp = wls(resid_mu,err,Xp,npn)
        vals.append(null["chi2"]-fp["chi2"])
    vals=np.asarray(vals)
    pval=(1+np.sum(vals>=observed))/(1+len(vals))
    return null,full,float(observed),float(pval),vals


def low_to_high_transfer(z, score, resid_mu, err, controls):
    cal = (z>=0.01)&(z<=0.03)
    val = (z>0.03)&(z<=0.15)
    Xc,nc = design_matrix(controls,score,cal)
    fitc = wls(resid_mu[cal],err[cal],Xc,nc)
    beta = fitc["coeff"]["foreground"]

    # Validation: allow nuisance/intercept refit, but foreground amplitude is frozen.
    controls_v={k:np.asarray(v)[val] for k,v in controls.items()}
    Xv,nv=design_matrix(controls_v,score[val])
    fullv=wls(resid_mu[val],err[val],Xv,nv,fixed_foreground_beta=beta)
    X0,n0=design_matrix(controls_v,None)
    nullv=wls(resid_mu[val],err[val],X0,n0)
    return {
        "n_calibration":int(cal.sum()),
        "n_validation":int(val.sum()),
        "beta_foreground_calibration":float(beta),
        "calibration_chi2":float(fitc["chi2"]),
        "validation_chi2_null":float(nullv["chi2"]),
        "validation_chi2_fixed_foreground":float(fullv["chi2"]),
        "validation_delta_chi2":float(nullv["chi2"]-fullv["chi2"]),
    }


def analyze(sn_path, gal_path, n_perm):
    sn=pd.read_csv(sn_path)
    gal=pd.read_csv(gal_path)

    zc=find_column(sn,["zCMB","zcmb","zHD","z","redshift"])
    rac=find_column(sn,["RA","ra","RA_DEG"])
    decc=find_column(sn,["DEC","dec","DEC_DEG","DECL"])
    muc=find_column(sn,["MU_SH0ES","MU","mu","distance_modulus"])
    ec=find_column(sn,["MU_SH0ES_ERR_DIAG","MUERR","muerr","MU_ERR","distance_modulus_err"])

    gz=find_column(gal,["z","zcmb","redshift","Z"])
    gra=find_column(gal,["RA","ra","RA_DEG"])
    gdec=find_column(gal,["DEC","dec","DEC_DEG","DECL"])

    sn=sn.copy()
    sn["_z"]=pd.to_numeric(sn[zc],errors="coerce")
    sn["_ra"]=pd.to_numeric(sn[rac],errors="coerce")
    sn["_dec"]=pd.to_numeric(sn[decc],errors="coerce")
    sn["_mu"]=pd.to_numeric(sn[muc],errors="coerce")
    sn["_err"]=pd.to_numeric(sn[ec],errors="coerce")
    sn=sn.replace([np.inf,-np.inf],np.nan).dropna(subset=["_z","_ra","_dec","_mu","_err"])
    sn=sn[(sn["_z"]>=SN_Z_MIN)&(sn["_z"]<=SN_Z_MAX)&(sn["_err"]>0)].reset_index(drop=True)

    gal=gal.copy()
    gal["_z"]=pd.to_numeric(gal[gz],errors="coerce")
    gal["_ra"]=pd.to_numeric(gal[gra],errors="coerce")
    gal["_dec"]=pd.to_numeric(gal[gdec],errors="coerce")
    gal=gal.replace([np.inf,-np.inf],np.nan).dropna(subset=["_z","_ra","_dec"])
    gal=gal[(gal["_z"]>=SHELL_EDGES[0])&(gal["_z"]<SHELL_EDGES[-1])].reset_index(drop=True)

    base=fit_global_baseline(sn["_z"].to_numpy(),sn["_mu"].to_numpy(),sn["_err"].to_numpy())
    mu0=mu_flat(sn["_z"].to_numpy(),70.0,base["omega_m"])+base["deltaM"]
    resid=sn["_mu"].to_numpy()-mu0

    delta,shell_meta=shell_density_contrasts(
        sn["_ra"].to_numpy(),sn["_dec"].to_numpy(),gal
    )
    raw=foreground_raw_score(sn["_z"].to_numpy(),delta)
    score,score_mean,score_sd=standardize(raw)
    controls=build_controls(sn)

    null,full,dchi2,pval,perms=blocked_permutation(
        sn["_z"].to_numpy(),score,resid,sn["_err"].to_numpy(),controls,n_perm
    )
    beta=full["coeff"]["foreground"]
    beta_err=full["stderr"]["foreground"]
    t=beta/beta_err if beta_err>0 else float("nan")

    jack=[]
    for leave in range(delta.shape[1]):
        d2=np.delete(delta,leave,axis=1)
        edges2=[]
        for s in range(delta.shape[1]):
            if s!=leave:
                edges2.append(s)
        # Recompute using original shell geometry but zero the omitted shell.
        dtmp=delta.copy()
        dtmp[:,leave]=0.0
        rawj=foreground_raw_score(sn["_z"].to_numpy(),dtmp)
        sj,_,_=standardize(rawj)
        Xj,nj=design_matrix(controls,sj)
        fj=wls(resid,sn["_err"].to_numpy(),Xj,nj)
        bj=fj["coeff"]["foreground"]
        jack.append({
            "left_out_shell":[float(SHELL_EDGES[leave]),float(SHELL_EDGES[leave+1])],
            "beta_foreground":float(bj),
            "sign_matches_primary":bool(np.sign(bj)==np.sign(beta) or bj==0 or beta==0)
        })

    sign_stable=all(x["sign_matches_primary"] for x in jack)
    if pval < 0.01 and sign_stable:
        classification="strong"
    elif pval < 0.05:
        classification="provisional"
    else:
        classification="not_established"

    transfer=low_to_high_transfer(
        sn["_z"].to_numpy(),score,resid,sn["_err"].to_numpy(),controls
    )

    return {
        "schema":"r1-2mrs-foreground-transfer-result-v1",
        "classification":classification,
        "data":{
            "sn_path":str(sn_path),
            "galaxy_path":str(gal_path),
            "n_sn":int(len(sn)),
            "n_galaxies_zlt0p05":int(len(gal)),
        },
        "baseline":base,
        "foreground_definition":{
            "shell_edges":SHELL_EDGES.tolist(),
            "angular_smoothing_sigma_deg":SMOOTH_SIGMA_DEG,
            "fixed_axis_ra_deg":AXIS_RA_DEG,
            "fixed_axis_dec_deg":AXIS_DEC_DEG,
            "score_mean_raw":score_mean,
            "score_sd_raw":score_sd,
            "shells":shell_meta,
        },
        "primary":{
            "beta_foreground_mag_per_score_sigma":float(beta),
            "beta_foreground_err":float(beta_err),
            "t_foreground":float(t),
            "chi2_null":float(null["chi2"]),
            "chi2_full":float(full["chi2"]),
            "delta_chi2":float(dchi2),
            "blocked_permutation_p":float(pval),
            "n_permutations":int(n_perm),
            "leave_one_shell_out":jack,
            "sign_stable_all_shell_dropouts":bool(sign_stable),
        },
        "low_to_high_transfer":transfer,
        "claim_boundary":[
            "This tests a 2MRS density proxy for LOS topography, not a BHSM-specific seam/de-encapsulation catalog.",
            "A positive result can also arise from conventional local-structure effects unless separately controlled.",
            "The 10-degree smoothing scale is carried forward from earlier exploratory work and was not rescanned here.",
        ]
    }


def write_report(result, path):
    p=result["primary"]
    t=result["low_to_high_transfer"]
    lines=[
        "# Frozen 2MRS foreground-transfer test",
        "",
        f"Classification: **{result['classification']}**",
        "",
        f"- N_SN: {result['data']['n_sn']}",
        f"- N_2MRS(z<0.05): {result['data']['n_galaxies_zlt0p05']}",
        f"- foreground coefficient: {p['beta_foreground_mag_per_score_sigma']:.6g} Â± {p['beta_foreground_err']:.6g} mag / score Ïƒ",
        f"- t: {p['t_foreground']:.3f}",
        f"- Î”Ï‡Â² beyond fixed-axis + controls: {p['delta_chi2']:.6g}",
        f"- blocked permutation p: {p['blocked_permutation_p']:.6g}",
        f"- sign stable under all leave-one-shell-out checks: {p['sign_stable_all_shell_dropouts']}",
        "",
        "## Low-z calibration â†’ higher-z transfer",
        f"- calibration N: {t['n_calibration']}",
        f"- validation N: {t['n_validation']}",
        f"- frozen low-z foreground Î²: {t['beta_foreground_calibration']:.6g}",
        f"- held-out Î”Ï‡Â² (null - fixed foreground): {t['validation_delta_chi2']:.6g}",
        "",
        "## Claim boundary",
        "This is a foreground-density transfer test. It is not yet a BHSM-specific black-hole/de-encapsulation source test, and conventional peculiar-velocity/lensing/selection explanations must still be separated.",
    ]
    path.write_text("\n".join(lines)+"\n",encoding="utf-8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--sn",required=True)
    ap.add_argument("--gal",required=True)
    ap.add_argument("--permutations",type=int,default=999)
    ap.add_argument("--artifact",required=True)
    ap.add_argument("--report",required=True)
    args=ap.parse_args()
    result=analyze(args.sn,args.gal,args.permutations)
    Path(args.artifact).parent.mkdir(parents=True,exist_ok=True)
    Path(args.artifact).write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    write_report(result,Path(args.report))
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()
