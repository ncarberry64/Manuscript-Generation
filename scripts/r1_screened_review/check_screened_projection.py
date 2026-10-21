"""Independent direct-z optical projection of the supplied R1 radial force.

This does not audit the action reduction: it independently integrates its
radial force along z, rather than reusing the theta-substituted projection.
It re-solves the inherited image branches and stationary delay differences.
"""
from pathlib import Path
import json,sys
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

REPO=Path(__file__).resolve().parents[2]
REVIEW=REPO/'artifacts/strong_lensing_r1_screened_20261020/review'
sys.path.insert(0,str(REPO/'code'))
from strong_lensing_r1_screened import R1ScreenedLens
from strong_lensing_comparison import (ARCSEC,MPC_OVER_C_DAYS,complete_map_response,
    response_vector,nuisance_projected_information)
from strong_lensing_local_family import TolmanVIILens,local_nuisance_columns

REF=REPO/'artifacts/strong_lensing_local_family_20261020/results.json'


def main(output=None):
    old=json.loads(REF.read_text());inputs=old['physical_inputs'];sol=old['optical_encounter']['solution']
    lens=R1ScreenedLens();c=lens.coefficients;m=lens.mass_length_Mpc;R=lens.radius_Mpc
    Dol=inputs['physical_maps_Mpc']['D_ol'][0][0];Dos=inputs['physical_maps_Mpc']['D_os'][0][0];Dls=inputs['physical_maps_Mpc']['D_ls'][0][0]
    q=Dls/Dos; beta=inputs['source_transverse_Mpc'][0]/Dos
    T=(1+inputs['z_l'])*Dol/q
    oldtheta=np.array(sol['images_rad'])[:,0]
    def mu(r):
        x=r/R
        return .5*(5*x**3-3*x**5) if x<1 else 1.
    def delta_phi_prime(r):
        enclosed=mu(r)
        disc=1+16*c.b*m*c.hubble_length_Mpc**2*enclosed/(c.lambda_over_MplH02*c.C**2*r**3)
        F=c.b**2/(c.C*(1+np.sqrt(disc)))
        return m*enclosed*F/r**2
    @lru_cache(maxsize=4096)
    def alpha_z(b):
        # Finite z intervals resolve the body and screening transition; the
        # final interval integrates the remaining tail to infinity directly.
        cuts=[0.]
        if b<R:cuts.append(np.sqrt(R*R-b*b))
        cuts+=sorted(set([R,lens.rV_Mpc,20.]))
        cuts=sorted(set(cuts))
        total=0.;error=0.
        integrand=lambda z:delta_phi_prime(np.hypot(b,z))/np.hypot(b,z)
        for lo,hi in zip(cuts,cuts[1:]+[np.inf]):
            value,err=quad(integrand,lo,hi,epsabs=1e-21,epsrel=2e-10,limit=200)
            total+=value;error+=err
        return 4*b*total
    @lru_cache(maxsize=1024)
    def psi_z(b):
        return quad(alpha_z,R,b,epsabs=1e-23,epsrel=2e-10,limit=100)[0]
    def F2(y):
        return 1. if y>=1 else -np.expm1(2.5*np.log1p(-y*y))
    def base_alpha(theta):
        impact=Dol*theta
        return 4*m*F2(abs(impact)/R)/impact
    def residual(theta):
        return theta-beta-q*(base_alpha(theta)+np.sign(theta)*alpha_z(abs(Dol*theta)))
    roots=np.array([brentq(residual,t-.001*ARCSEC,t+.001*ARCSEC,xtol=2e-20,rtol=1e-14) for t in oldtheta])
    def arrival(t):
        return (.5*T*(t-beta)**2-(1+inputs['z_l'])*(lens.base.application([Dol*t,0.])[0]+psi_z(abs(Dol*t))))*MPC_OVER_C_DAYS
    arrival_new=np.array([arrival(t) for t in roots]);delays=arrival_new[1:]-arrival_new[0]
    first_time=-(1+inputs['z_l'])*np.array([psi_z(abs(Dol*t)) for t in oldtheta])*MPC_OVER_C_DAYS
    first_delays=first_time[1:]-first_time[0]
    first_images=[]
    for t in oldtheta:
        J=1-q*Dol*lens.base.application([Dol*t,0.])[2][0,0]
        first_images.append(q*np.sign(t)*alpha_z(abs(Dol*t))/J/ARCSEC)
    impacts=abs(Dol*oldtheta)
    az=np.array([alpha_z(b) for b in impacts]);am=np.array([lens.extra_alpha(b) for b in impacts])
    pz=np.array([psi_z(b) for b in impacts]);pm=np.array([lens.extra_potential(b) for b in impacts])
    result=dict(status='independent direct-z projection; conditional same radial equation and body',
        scalar_coefficients={'b':c.b,'C':c.C,'F_infinity':c.unscreened_extra,'rV_Mpc':lens.rV_Mpc},
        impact_Mpc=impacts.tolist(),direct_z_extra_alpha_rad=az.tolist(),module_extra_alpha_rad=am.tolist(),
        max_relative_alpha_discrepancy=float(np.max(abs(az-am)/az)),
        direct_z_extra_potential_Mpc=pz.tolist(),module_extra_potential_Mpc=pm.tolist(),
        max_absolute_potential_discrepancy_Mpc=float(np.max(abs(pz-pm))),
        new_images_arcsec=(roots/ARCSEC).tolist(),image_changes_arcsec=((roots-oldtheta)/ARCSEC).tolist(),
        first_order_image_changes_arcsec=first_images,
        new_relative_delays_days=delays.tolist(),
        exact_model_delay_changes_seconds=((delays-np.array(sol['relative_delays_days']))*86400).tolist(),
        first_order_delay_changes_seconds=(first_delays*86400).tolist(),
        max_image_equation_residual_rad=float(max(abs(residual(t)) for t in roots)),
        precision_scope='numerical agreement of local QS weak-field model, not accuracy of its truncations')
    output = output or REVIEW/'screened_projection_check_results.json'
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n');print(json.dumps(result,indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, help='Write a separate replay receipt; default preserves the original command behavior.')
    args = parser.parse_args()
    main(args.output)
