"""Replay the conditional frozen-R1 local screened lens prediction."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))

import numpy as np
from scipy.integrate import quad
import action_native_matter_n2 as an
import background as bg
from strong_lensing_comparison import (ARCSEC, solve_images, complete_map_response,
    observable_vector, response_vector, nuisance_projected_information)
from strong_lensing_local_family import (TolmanVIILens, saved_geometry,
    local_nuisance_columns)
from strong_lensing_r1_screened import R1ScreenedLens, R1Coefficients, solve_axis_images


def encode(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    raise TypeError(type(x).__name__)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--reference',type=Path,default=ROOT/'artifacts/strong_lensing_first_comparison_20261020/results.json')
    parser.add_argument('--output',type=Path,default=ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json')
    args=parser.parse_args()
    old,geometry,source=saved_geometry(args.reference)
    base=TolmanVIILens()
    coefficients=R1Coefficients.frozen(geometry.z_l,old['experiment']['H0_km_s_Mpc'])
    lens=R1ScreenedLens(base,coefficients)
    conventional=solve_images(geometry,base,source)
    screened=solve_axis_images(geometry,lens,source,conventional['images_rad'])
    zero=np.zeros((3,2,2))
    linear=complete_map_response(conventional['images_rad'],geometry,base,source,zero,
                               intrinsic_encounter=lens.intrinsic_extra)
    impacts=np.linalg.norm(conventional['images_rad']@geometry.D_ol.T,axis=1)
    changes=dict(images_microarcsec=(screened['images_arcsec']-conventional['images_arcsec'])*1e6,
        relative_delays_seconds=(screened['relative_delays_days']-conventional['relative_delays_days'])*86400,
        relative_delays_fractional=screened['relative_delays_days']/conventional['relative_delays_days']-1,
        linear_relative_delays_seconds=linear['relative_delays_days']*86400,
        linear_images_microarcsec=linear['images_arcsec']*1e6)
    ratio0=conventional['relative_delays_days'][0]/conventional['relative_delays_days'][1]
    ratio1=screened['relative_delays_days'][0]/screened['relative_delays_days'][1]
    N,names=local_nuisance_columns(geometry,base,source,conventional)
    B=response_vector(linear)[:,None]
    # This covariance is an explicit numerical scale for rank calculations;
    # it is not an observational forecast or an actual data likelihood.
    scales=np.r_[np.full(6,.001),np.full(2,.1)]
    nuisance=nuisance_projected_information(B,N,np.diag(scales**2))
    step=1e-4
    centers=[]
    from dataclasses import replace
    for direction in np.eye(2):
        shift=geometry.D_ol@(direction*ARCSEC*step)
        p=solve_images(geometry,replace(base,centre_Mpc=tuple(shift)),source,seeds=conventional['images_rad'])
        m=solve_images(geometry,replace(base,centre_Mpc=tuple(-shift)),source,seeds=conventional['images_rad'])
        centers.append((observable_vector(p)-observable_vector(m))/(2*step))
    expanded=nuisance_projected_information(B,np.column_stack([N,*centers]),np.diag(scales**2))
    radii_kpc=[.1,.3,1.,3.,10.,30.,100.,212.0292249771296,300.,1000.]
    radial=[]
    for rk in radii_kpc:
        r=rk/1000
        _,alpha,_=base.application([r,0.])
        radial.append(dict(radius_kpc=rk,local_force_extra_fraction=lens.extra_force(r),
            projected_deflection_extra_fraction=lens.extra_alpha(r)/alpha[0]))
    # Leading hydrostatic response for the same density, rather than reusing
    # the exact GR TOV pressure as an exact modified-gravity equilibrium.
    pressure=[]
    for x in [0.,.25,.5,.75,1.]:
        delta=quad(lambda t:15/(8*np.pi)*(1-t*t)*t*(2.5-1.5*t*t)*
                   lens.extra_force(base.radius_Mpc*t),x,1.,epsabs=1e-17,epsrel=1e-12)[0]
        pN=15/(32*np.pi)*(2-5*x*x+4*x**4-x**6)
        pressure.append(dict(r_over_R=x,newtonian_pbar=pN,extra_pbar=delta,
                             extra_fraction=None if pN==0 else delta/pN))
    history=an.frozen_background()
    redshift_law=[]
    for z in [0.,.1,.3,.5,.8,1.,1.5]:
        n=-np.log1p(z); y=history.sol(n); v,H=y[1:]; d=H*bg.rhs(n,y)[1]
        mixing=2*v*v/bg.LAMBDA
        stiffness=1+4*(d+2*H*v)/bg.LAMBDA-2*v**4/bg.LAMBDA**2
        rv=(16*mixing*base.mass_length_Mpc*coefficients.hubble_length_Mpc**2/
            (bg.LAMBDA*stiffness**2))**(1/3)
        redshift_law.append(dict(z=z,unscreened_extra_fraction=mixing**2/(2*stiffness),
                                rV_kpc=rv*1000))
    provenance={}
    for path in [ROOT/'preregistration/prediction_manifest.json',ROOT/'code/background.py',
                 ROOT/'manuscript/sections/04_scalar_metric_bridge.tex',args.reference]:
        provenance[str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    experiment={key:old['experiment'][key] for key in (
        'H0_km_s_Mpc','light','observer_wavelength_nm','source_rest_wavelength_nm',
        'source_transverse_Mpc','reference_beta_arcsec','z_l','z_s')}
    experiment.update(
        fixed='same prescribed finite density excess, physical source displacement, observed redshifts, observer clock and compatible background distance maps',
        varied='activate the independently specified frozen-R1 cubic scalar response of the same body; re-solve the same image branches',
        envelope_role='conditional R1 metric response sourced by the finite body; distinct from the previously null constant target-phase family',
        double_counting='smooth R1 background enters the distance maps once; body potential plus its derived scalar response enters the encounter once',
        incoming_light_distinct_from=old['experiment']['incoming_light_distinct_from'])
    result=dict(status='CONDITIONAL_R1_CUBIC_EFT_LOCAL_QUASISTATIC_PREDICTION',
        source_commit='49a1dabbe148f8ba9da4ac079a8d93fe97c5de76',
        experiment=experiment,physical_maps_Mpc=old['physical_maps_Mpc'],
        coefficients=asdict(coefficients),body=dict(mass_Msun=base.mass_Msun,radius_Mpc=base.radius_Mpc,
            density='15 M/(8 pi R^3) (1-r^2/R^2) inside R; prescribed excess over smooth background'),
        screening=dict(rV_Mpc=lens.rV_Mpc,unscreened_extra_fraction=coefficients.unscreened_extra,
            H_l_rV_over_c=coefficients.H_over_H0*lens.rV_Mpc/coefficients.hubble_length_Mpc,
            effective_extra_mass_infinity_Msun=base.mass_Msun*coefficients.unscreened_extra,
            mass_scaling_at_fixed_z='rV proportional to M^(1/3)',
            effective_halo_degeneracy='delta M_eff(r)=M_enclosed(r)*extra_force(r), with nonnegative density; flexible ordinary halos reproduce this retained Phi=Psi field'),
        conventional=conventional,screened=screened,changes=changes,
        delay_ratio=dict(conventional=ratio0,screened=ratio1,change=ratio1-ratio0,
                         fractional_change=ratio1/ratio0-1,
                         qualification='common MST scale cancels; source/profile/centering freedoms remain'),
        image_impact_radii_Mpc=impacts,radial_profile=radial,
        frozen_redshift_law=redshift_law,
        leading_hydrostatic_response=dict(pbar_units='p=(G M^2/R^4)*pbar',values=pressure,
            equation='delta p(r)=integral_r^R rho(s) GM_enc(s)/s^2 F(s) ds',
            qualification='leading weak-field pressure correction; not an exact modified TOV solution'),
        checks=lens.checks(impacts),
        identifiability=dict(illustrative_image_sigma_arcsec=.001,illustrative_delay_sigma_days=.1,
            fixed_center_nuisance_names=names,fixed_center=nuisance,
            free_center_nuisance_names=names+['lens_center_x_per_arcsec','lens_center_y_per_arcsec'],
            free_center=expanded,
            qualification='rank concerns this finite spherical profile only; covariance is a numerical scaling, not observed precision'),
        assumptions=[
            'Pure cubic R1 compatibility EFT with the previously frozen background; no retuned parameters.',
            'The independent rolling scalar T is not the S6 target phase and not a duplicate metric scalar.',
            'Leading weak-field, local quasi-static approximation with nonlinear second spatial derivatives retained.',
            'Regular scalar branch, pi prime(0)=0, decaying isolated excess at large radius; no extra interface charge.',
            'Same prescribed density, physical source, measured redshifts and three compatible background maps.',
            'No optional high-spatial-derivative topographic kernel is included in this named EFT specialization.',
            'Line projection treats the local coefficients at lens time; global closed-cosmology matching is not completed.',
            'Numerical precision is not a bound on omitted post-Newtonian, time-dependent, or additional-sector physics.',
            'The existing zero additional constant-eta response is preserved; this is a separate previously specified scalar sector.'
        ],
        provenance_sha256=provenance)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,default=encode)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({k:result[k] for k in ['status','coefficients','screening','changes','delay_ratio','checks']},indent=2,default=encode))
    print('fixed-center projected rank',nuisance['signal_rank'],'free-center projected rank',expanded['signal_rank'])
    print('output',args.output)


if __name__=='__main__': main()
