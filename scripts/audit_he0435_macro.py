"""Audit frozen macro predictions without fitting or consuming withheld light."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import root

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import strong_lensing_observed_macro as macro
from strong_lensing_comparison import ARCSEC, lens_application
from strong_lensing_robustness import ProfileR1Lens

OUT=ROOT/'artifacts/strong_lensing_observed_20261020'


def encode(value):
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,np.generic): return value.item()
    raise TypeError(type(value).__name__)


def census(parameters, model):
    p=np.asarray(parameters); d=macro.frozen_distances()
    geometry,U=macro.screen_geometry(*p[4:7]); source=d['os']*ARCSEC*p[2:4]
    body=macro.body_from_parameters(p)
    exact=body if model=='conventional' else ProfileR1Lens(body,macro.frozen_coefficients())
    temporary=exact if model=='conventional' else macro.CompositeContinuationLens(body,macro.FrozenExtraProjection(exact,count=240))
    beta=np.linalg.solve(geometry.D_os,source)
    center=np.asarray(body.centre_Mpc)/d['ol']
    # Positive circular projected mass is bounded by its finite total mass.
    # This excludes roots outside this circle, but a grid search inside is
    # still not an analytic proof that every stationary point was located.
    epsilon=0. if model=='conventional' else exact.coefficients.unscreened_extra
    strength=np.linalg.norm(geometry.G,2)*4*body.mass_length_Mpc*(1+epsilon)/d['ol']
    bn,cn=np.linalg.norm(beta),np.linalg.norm(center)
    bound=.5*(bn+cn+np.sqrt((bn-cn)**2+4*strength))/ARCSEC
    axis=np.linspace(-1.03*bound,1.03*bound,23)
    seeds=[np.array([x,y]) for x in axis for y in axis]
    angles=np.linspace(0.,2*np.pi,48,endpoint=False)
    for radius in (.001,.01,.03,.07,.15,.3,.6,.9,1.2,1.6,2.2):
        seeds.extend(center/ARCSEC+radius*np.c_[np.cos(angles),np.sin(angles)])
    seeds.extend([center/ARCSEC,beta/ARCSEC,np.zeros(2)])
    def kernel(q, lens):
        point=d['ol']*q*ARCSEC
        _,alpha,hessian=lens.application(point)
        residual=q-beta/ARCSEC-geometry.G@alpha/ARCSEC
        jacobian=np.eye(2)-geometry.G@hessian*d['ol']
        return residual,jacobian
    roots=[]; rejected=0
    for seed in seeds:
        solution=root(lambda q:kernel(q,temporary)[0],seed,jac=lambda q:kernel(q,temporary)[1],tol=1e-10)
        if (np.max(abs(kernel(solution.x,temporary)[0]))>2e-7 or
            np.linalg.norm(solution.x)>bound*1.001 or not np.isfinite(solution.x).all()):
            rejected+=1; continue
        if not any(np.linalg.norm(solution.x-item)<1e-5 for item in roots): roots.append(solution.x)
    # Every newly located section is refined using the original, noninterpolated
    # source projection before magnifications/residuals are admitted to this audit.
    refined=[]
    for seed in roots:
        solution=root(lambda q:kernel(q,exact)[0],seed,jac=lambda q:kernel(q,exact)[1],tol=2e-11)
        if np.max(abs(kernel(solution.x,exact)[0]))>2e-7:
            raise RuntimeError('original projection did not verify a searched image')
        if not any(np.linalg.norm(solution.x-item)<1e-5 for item in refined): refined.append(solution.x)
    observed=macro.macro_application(p,model=model)
    A=geometry.D_os/d['os']
    reference=observed['images_light_centered_arcsec'][0]
    reference_jac=kernel(reference,exact)[1]
    reference_mu=1/np.linalg.det(A@reference_jac)
    images=[]
    for position in sorted(refined,key=np.linalg.norm):
        residual,jacobian=kernel(position,exact)
        signed_mu=1/np.linalg.det(A@jacobian)
        distances=np.linalg.norm(observed['images_light_centered_arcsec']-position,axis=1)
        nearest=int(np.argmin(distances))
        label='ABCD'[nearest] if distances[nearest]<.05 else 'additional_stationary_image'
        original=lens_application(position*ARCSEC,geometry,exact,source)
        images.append(dict(label=label,position_light_centered_arcsec=position,
                           position_relative_A_arcsec=position-reference,
                           signed_magnification_relative_unperturbed_FLRW=float(signed_mu),
                           absolute_macro_flux_relative_A=float(abs(signed_mu/reference_mu)),
                           image_parity=int(np.sign(signed_mu)),
                           normalized_lens_jacobian_singular_values=np.linalg.svd(jacobian,compute_uv=False),
                           original_lens_residual_rad=float(max(abs(original['residual']))),
                           original_arrival_stationarity_Mpc=float(max(abs(original['stationarity'])))))
    return dict(searched_stationary_image_count=len(images),images=images,
                seed_count=len(seeds),rejected_seed_count=rejected,
                rigorous_outer_image_radius_bound_arcsec=float(bound),
                bound='For R>|beta_eff|,|center|: (R-|beta_eff|)(R-|center|)<=||G||4m_total/D_l; m_total includes the derived positive R1 effective mass exactly once.',
                search_scope='Dense Cartesian and logarithmic polar multiseed search inside the analytic outer radial bound; no proof that every interior root basin was sampled.',
                magnification_convention='Actual source angular Jacobian=(D_os_eff/D_s) times normalized main-lens Jacobian. Ratios exclude microlensing, extinction, finite source and galaxy contamination.',
                observational_scope='Only published four point-image positions trained the model. No central-image detection limit or non-detection likelihood has been imposed.')


def main():
    conventional_path=OUT/'conventional_training.json'; r1_path=OUT/'r1_training.json'
    conventional=json.loads(conventional_path.read_text()); fitted=json.loads(r1_path.read_text())
    p=np.asarray(conventional['parameter_vector'])
    reference=macro.macro_application(p,model='conventional')
    activated=macro.macro_application(p,model='r1')
    impacts=np.linalg.norm(activated['images_light_centered_arcsec']*ARCSEC*macro.frozen_distances()['ol']-
                           np.asarray(activated['lens'].centre_Mpc),axis=1)
    fixed= dict(parameters=p,images_light_centered_arcsec=activated['images_light_centered_arcsec'],
                relative_delays_days=activated['delays_days'],
                image_change_microarcsec=(activated['images_light_centered_arcsec']-reference['images_light_centered_arcsec'])*1e6,
                delay_change_seconds=(activated['delays_days']-reference['delays_days'])*86400,
                original_source_checks=activated['lens'].checks(impacts),
                max_lens_residual_rad=activated['max_lens_residual_rad'],
                max_arrival_stationarity_Mpc=activated['max_arrival_stationarity_Mpc'],
                rV_kpc=activated['lens'].rV_Mpc*1000,
                scope='Activate the original fixed-state R1 source at identical fitted ordinary body, source, screen, redshifts, maps and observer clock; no nuisance refit.')
    refitted=(np.asarray(fitted['predicted_delays_days'])-np.asarray(conventional['predicted_delays_days']))*86400
    result=dict(status='FROZEN_TRAINING_PREDICTIONS_FORWARD_AUDIT_NO_TARGET_CONSUMPTION',
                fixed_body_R1_activation=fixed,
                separately_refitted_R1_minus_conventional_delay_seconds=refitted,
                refit_minus_fixed_body_delay_change_seconds=refitted-fixed['delay_change_seconds'],
                interpretation='The separately refitted difference includes parameter reallocation and finite optimizer termination. It is not the fixed-body physical R1 signal. No further fit or model selection was performed in this audit.',
                conventional_census=census(conventional['parameter_vector'],'conventional'),
                conditional_R1_census=census(fitted['parameter_vector'],'r1'),
                input_sha256={path.relative_to(ROOT).as_posix():hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in (conventional_path,r1_path,Path(__file__),ROOT/'code/strong_lensing_observed_macro.py')},
                target_values_consumed=False,physical_admission=False)
    path=OUT/'macro_forward_audit.json'
    if path.exists(): raise RuntimeError('Preserve existing forward audit rather than silently overwrite it')
    path.write_text(json.dumps(result,indent=2,default=encode,allow_nan=False)+'\n',encoding='utf8',newline='\n')
    print(json.dumps(dict(path=str(path),fixed_body_delay_change_seconds=fixed['delay_change_seconds'],
                          separately_refitted_delay_change_seconds=refitted,
                          conventional_image_count=result['conventional_census']['searched_stationary_image_count'],
                          R1_image_count=result['conditional_R1_census']['searched_stationary_image_count']),default=encode),flush=True)


if __name__=='__main__': main()
