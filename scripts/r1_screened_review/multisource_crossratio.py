"""Conditional multiple-source distance invariant, using the frozen R1 background.

No observed data or fitted perturbation amplitude is used. The direct angular
formula assumes one spherical finite lens, two exterior images per source,
one common physical sheet at the lens, and a determined angular lens centre.
It is NOT invariant to arbitrary independent sheets for each source plane.
Shared shear / uncertain-centre robustness is a separate inference question.
"""
from pathlib import Path
import json
import sys
import numpy as np
from scipy.integrate import quad

REPO = Path(__file__).resolve().parents[2]
REVIEW = REPO / 'artifacts/strong_lensing_r1_screened_20261020/review'
sys.path.insert(0, str(REPO / 'code'))
import action_native_matter_n2 as an
import background as bg

H0 = 67.4
C_KMS = 299792.458
G = 4.300917270036279e-9
ARCSEC = np.pi / (180 * 3600)


def invariant(q):
    inv = 1 / np.asarray(q)
    return float((inv[0] - inv[1]) / (inv[0] - inv[2]))


def distances(E, zl, source_redshifts):
    chi = lambda z: quad(lambda t: 1 / E(t), 0, z, epsabs=2e-12, epsrel=2e-12)[0]
    sk = lambda x: np.sin(np.sqrt(bg.KCURV) * x) / np.sqrt(bg.KCURV)
    cl = chi(zl)
    Dol = C_KMS / H0 * sk(cl) / (1 + zl)
    rows = []
    for zs in source_redshifts:
        cs = chi(zs)
        Dos = C_KMS / H0 * sk(cs) / (1 + zs)
        Dls = C_KMS / H0 * sk(cs - cl) / (1 + zs)
        rows.append(dict(z_s=zs, D_os_Mpc=Dos, D_ls_Mpc=Dls, q=Dls / Dos))
    return Dol, rows


def exterior_pairs(q, Dol, mass_Msun=1e11, R_Mpc=.003, sheet_per_Mpc=0.):
    m = G * mass_Msun / C_KMS**2
    q = np.asarray(q)
    E2 = 4 * m * q / Dol
    beta = .1 * np.sqrt(E2)
    a = 1 - q * Dol * sheet_per_Mpc
    disc = np.sqrt(beta**2 + 4 * a * E2)
    plus = (beta + disc) / (2 * a)
    minus = (beta - disc) / (2 * a)
    impacts = np.c_[plus, abs(minus)] * Dol / R_Mpc
    if not np.all(impacts > 1):
        raise ValueError('at least one nominal exterior ray intersects the finite body')
    product = -plus * minus
    inv = 1 / product
    return dict(images_arcsec=np.c_[plus, minus].tolist(),
        images_unit='radians; converted to arcsec below', impact_over_R=impacts.tolist(),
        P_arcsec_squared=(product / ARCSEC**2).tolist(),
        angular_invariant=float((inv[0]-inv[1])/(inv[0]-inv[2])))


def main(output=None):
    zl = .3
    zs = [.8, 1.2, 1.8]
    sol = an.frozen_background()
    Er1 = lambda z: float(sol.sol(-np.log1p(z))[2])
    olambda = 1 - bg.OMEGA_M0 - bg.OMEGA_R0 - bg.OMEGA_K0
    Elcdm = lambda z: np.sqrt(bg.OMEGA_M0*(1+z)**3 + bg.OMEGA_R0*(1+z)**4
                            + bg.OMEGA_K0*(1+z)**2 + olambda)
    Dol, r1 = distances(Er1, zl, zs)
    Dllcdm, lcdm = distances(Elcdm, zl, zs)
    q = np.array([x['q'] for x in r1])
    c1 = invariant(q)
    c0 = invariant([x['q'] for x in lcdm])
    examples = []
    for mass, kappa_mid in [(1e11, 0.), (1e11, .08), (1.4e11, -.04)]:
        S = kappa_mid / (Dol*q[1])
        row = exterior_pairs(q, Dol, mass, sheet_per_Mpc=S)
        row['images_arcsec'] = (np.asarray(row['images_arcsec']) / ARCSEC).tolist()
        row['images_unit'] = 'arcsec'
        row.update(mass_Msun=mass, sheet_per_Mpc=S, convergence_at_middle_source=kappa_mid)
        examples.append(row)
    # Arbitrary independent source sheets have three independent nuisance rows.
    # This computes their local effect on the very same angular statistic.
    h = 1e-6
    def with_independent_kappa(kappa):
        P = (4*G*1e11/C_KMS**2/Dol)*q/(1-np.asarray(kappa))
        ip = 1/P
        return float((ip[0]-ip[1])/(ip[0]-ip[2]))
    dk = []
    for j in range(3):
        d = np.eye(3)[j]*h
        dk.append((with_independent_kappa(d)-with_independent_kappa(-d))/(2*h))
    result = dict(status='conditional distance-ratio prediction; not a detection or novelty claim',
        lens_redshift=zl, source_redshifts=zs, H0_km_s_Mpc=H0,
        frozen_R1=dict(D_ol_Mpc=Dol, sources=r1, crossratio=c1),
        comparison_closed_LCDM=dict(D_ol_Mpc=Dllcdm, sources=lcdm, crossratio=c0,
            Omega_m0=bg.OMEGA_M0, Omega_r0=bg.OMEGA_R0, Omega_k0=bg.OMEGA_K0,
            Omega_Lambda0=olambda,
            comparison_status='declared same-H0 same-density/curvature cosmological comparator; not observationally refitted'),
        difference=c1-c0, fractional_difference=c1/c0-1,
        common_mass_and_sheet_controls=examples,
        derivative_wrt_independent_source_convergences=dk,
        assumptions=['one common spherical finite mass outside all measured rays',
            'known or independently determined lens centre',
            'one shared physical sheet at the primary lens plane',
            'no independent arbitrary per-source-plane convergence',
            'shared shear and profile uncertainty require additional observables and a joint physical model'],
        invariant='(1/P1-1/P2)/(1/P1-1/P3)=(1/q1-1/q2)/(1/q1-1/q3), P=-theta_plus*theta_minus, q=Dls/Dos')
    path = output or REVIEW / 'multisource_crossratio_results.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, help='Write a separate replay receipt; default preserves the original command behavior.')
    args = parser.parse_args()
    main(args.output)
