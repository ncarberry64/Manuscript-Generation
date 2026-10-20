"""Bounded shared-lens nuisance check, with three finite-body image systems.

The source redshifts and offset sizes are fixed before the rank calculation.
This uses the frozen unperturbed R1 distances. Generic fractional
q_s=D_ls/D_os columns test propagation-map identifiability. A separate
fixed screened-R1 model-difference column tests the derived local response;
neither is a newly fitted fundamental scalar amplitude.
One shared physical convergence/shear is retained, and the centre is free.
Independent source-plane convergence is then added to expose the limitation.
"""
from pathlib import Path
import json
import numpy as np
from multisource_crossratio import ARCSEC, G, C_KMS, H0
REPO = Path(__file__).resolve().parents[2]
REVIEW = REPO / 'artifacts/strong_lensing_r1_screened_20261020/review'
from strong_lensing_local_family import TolmanVIILens
from strong_lensing_r1_screened import R1ScreenedLens
from strong_lensing_comparison import (LensGeometry, solve_images,
    complete_map_response, response_vector, nuisance_projected_information)


class SharedLens:
    def __init__(self, body, sheet):
        self.body=body; self.sheet=sheet; self.b=body.b
    def application(self, x):
        psi,alpha,hess=self.body.application(x)
        v=np.asarray(x)-np.asarray(self.body.centre_Mpc)
        return psi+.5*v@self.sheet@v,alpha+self.sheet@v,hess+self.sheet


def main(output=None):
    base=json.loads((REVIEW/'multisource_crossratio_results.json').read_text())
    rows=base['frozen_R1']['sources']; Dol=base['frozen_R1']['D_ol_Mpc']
    q=np.array([r['q'] for r in rows]); ns=len(q)
    centre=np.array([.02,-.01])*ARCSEC*Dol
    body=TolmanVIILens(centre_Mpc=tuple(centre))
    screened=R1ScreenedLens(body)
    mats=[np.eye(2),np.diag([1.,-1.]),np.array([[0.,1.],[1.,0.]])]
    sheet_scale=1/(Dol*q[1])
    sheet=sheet_scale*(.03*mats[0]+.01*mats[1]+.007*mats[2])
    lens=SharedLens(body,sheet)
    geometries=[]; sources=[]; solutions=[]; offsets=[0]
    for row,phi in zip(rows,[.15,1.1,2.2]):
        g=LensGeometry(Dol*np.eye(2),row['D_os_Mpc']*np.eye(2),row['D_ls_Mpc']*np.eye(2),.3)
        E=np.sqrt(4*body.mass_length_Mpc*row['q']/Dol)
        beta=centre/Dol+.1*E*np.array([np.cos(phi),np.sin(phi)])
        source=g.D_os@beta
        solution=solve_images(g,lens,source)
        geometries.append(g);sources.append(source);solutions.append(solution)
        offsets.append(offsets[-1]+3*len(solution['images_rad'])-1)
    size=offsets[-1]; zero=np.zeros((3,2,2))
    names=[f'source_{i}_{axis}' for i in range(ns) for axis in ['x','y']]
    names+=['shared_log_M','shared_log_R','shared_center_x','shared_center_y',
        'shared_physical_convergence','shared_physical_shear1','shared_physical_shear2']
    names+=[f'free_delay_scale_{i}' for i in range(ns)]
    N=np.zeros((size,len(names)));B=np.zeros((size,ns));Nlos=np.zeros((size,ns));sig=[]
    Bscreen=np.zeros((size,1))
    for j,(g,source,sol) in enumerate(zip(geometries,sources,solutions)):
        sl=slice(offsets[j],offsets[j+1]);images=sol['images_rad'];ni=len(images)
        sig.extend([.001]*(2*ni)+[.1]*(ni-1))
        def response(partial=None,ds=None,maps=zero):
            return response_vector(complete_map_response(images,g,lens,source,maps,
                d_source_Mpc=ds,intrinsic_encounter=partial))
        for k in range(2):N[sl,2*j+k]=response(ds=g.D_os[:,k]*ARCSEC)
        for k,key in enumerate(['log_mass','log_radius']):
            N[sl,2*ns+k]=response(partial=body.intrinsic(key))
        for k in range(2):
            dc=Dol*ARCSEC*np.eye(2)[k]
            def centre_partial(x,dc=dc):
                psi,alpha,hess=lens.application(x)
                return -alpha@dc,-hess@dc
            N[sl,2*ns+2+k]=response(partial=centre_partial)
        for k,mat in enumerate(mats):
            def tidal_partial(x,mat=mat):
                v=np.asarray(x)-centre
                return .5*sheet_scale*v@mat@v,sheet_scale*mat@v
            N[sl,2*ns+4+k]=response(partial=tidal_partial)
        delay=np.zeros(3*ni-1);delay[2*ni:]=sol['relative_delays_days']
        N[sl,2*ns+7+j]=delay
        B[sl,j]=response(maps=np.array([np.zeros((2,2)),np.zeros((2,2)),g.D_ls]))
        Bscreen[sl,0]=response(partial=screened.intrinsic_extra)
        def independent_sheet(x):
            v=np.asarray(x)-centre
            scale=1/(Dol*q[j])
            return .5*scale*(v@v),scale*v
        Nlos[sl,j]=response(partial=independent_sheet)
    C=np.diag(np.array(sig)**2)
    fixed=nuisance_projected_information(B,np.empty((size,0)),C)
    shared=nuisance_projected_information(B,N,C)
    arbitrary=nuisance_projected_information(B,np.c_[N,Nlos],C)
    screened_fixed=nuisance_projected_information(Bscreen,np.empty((size,0)),C)
    screened_shared=nuisance_projected_information(Bscreen,N,C)
    screened_independent=nuisance_projected_information(Bscreen,np.c_[N,Nlos],C)
    fields=['nuisance_rank','signal_rank','raw_signal_singular_values','projected_signal_singular_values']
    pack=lambda x:{k:(x[k].tolist() if isinstance(x[k],np.ndarray) else x[k]) for k in fields}
    whitened_residual=screened_independent['projected_response'][:,0]
    residual_norm=np.linalg.norm(whitened_residual)
    statistic_weights=whitened_residual/(np.asarray(sig)*residual_norm)
    result=dict(status='conditional identifiability design and fixed R1 local response; no data detection or forecast',
        source_redshifts=base['source_redshifts'], source_offset='0.1 theta_E at angles (0.15,1.1,2.2) rad relative to shared center',
        image_counts=[len(s['images_rad']) for s in solutions],
        images_arcsec=[s['images_arcsec'].tolist() for s in solutions],
        impact_over_R=[(np.linalg.norm(Dol*s['images_rad']-centre,axis=1)/body.radius_Mpc).tolist() for s in solutions],
        nuisance_names=names,signal_names=[f'log_q_{j}' for j in range(ns)],
        fixed_nuisance=pack(fixed),shared_lens_and_tides_free_delay_scales=pack(shared),
        additionally_independent_source_convergence=pack(arbitrary),
        local_screened_signal=dict(
            source='same R1ScreenedLens intrinsic_extra; no geometry or nuisance changes',
            fixed_nuisance=pack(screened_fixed),
            shared_lens_and_tides_free_delay_scales=pack(screened_shared),
            additionally_independent_source_convergence=pack(screened_independent),
            shared_retained_norm_fraction=float(np.linalg.norm(screened_shared['projected_response'])/np.linalg.norm(Bscreen/np.asarray(sig)[:,None])),
            independent_source_convergence_retained_norm_fraction=float(np.linalg.norm(screened_independent['projected_response'])/np.linalg.norm(Bscreen/np.asarray(sig)[:,None])),
            projected_response_units='whitened by 0.001 arcsec / 0.1 day; not original physical coordinates',
            nuisance_orthogonal_statistic=dict(
                definition='S = weights dot (observables - conventional baseline); weights annihilate every admitted nuisance column to first order',
                observable_order='for each source in increasing redshift: image0_x,image0_y,image1_x,image1_y,image2_x,image2_y,delay1,delay2; angles arcsec, delays days',
                weights=statistic_weights.tolist(),
                expected_screened_response=float(residual_norm),
                units='dimensionless numerical normalization; illustrative covariance is not a measured likelihood'),
            response_vector=Bscreen[:,0].tolist(),
            projected_shared_response=screened_shared['projected_response'][:,0].tolist(),
            projected_independent_source_convergence_response=screened_independent['projected_response'][:,0].tolist()),
        covariance='illustrative 1 mas coordinate and 0.1 day delay; not forecast',
        effective_halo_degeneracy=dict(
            mass_profile='delta M_eff(r) = M mu(r) F(r); delta rho_eff = M/(4 pi r^2) d[mu F]/dr >= 0',
            total_extra_effective_mass_Msun=float(body.mass_Msun*screened.coefficients.unscreened_extra),
            scope='because Phi=Psi at this order, a sufficiently flexible conventional spherical halo reproduces the complete optical modification; the retained rank is only relative to the stated finite profile family'),
        qualification='same finite spherical radial profile, common physical sheet/shear/center; source-dependent convergence removes the cosmographic q response but does not remove the computed local radial-shape response')
    path = output or REVIEW/'shared_lens_rank_results.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n', encoding='utf-8', newline='\n');print(json.dumps(result,indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, help='Write a separate replay receipt; default preserves the original command behavior.')
    args = parser.parse_args()
    main(args.output)
