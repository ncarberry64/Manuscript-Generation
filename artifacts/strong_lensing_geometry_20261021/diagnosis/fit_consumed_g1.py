"""New geometry development fit. Every consumed observation is already exposed."""
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import least_squares,minimize

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts')]
import strong_lensing_geometry_revision as geometry
import strong_lensing_observed_macro as old
import strong_lensing_macro_revision as local
from fit_he0435_macro_revision import setup
from strong_lensing_robustness import ProfileR1Lens


def main():
    started=time.perf_counter()
    base=json.loads((ROOT/'artifacts/strong_lensing_revision_20261020/revised_conventional_training.json').read_text())
    p0=np.array(base['parameter_vector'])
    g=json.loads((ROOT/'artifacts/strong_lensing_revision_20261020/gaia/gaia_astrometry_score.json').read_text())
    y=np.array(g['association']['relative_arcsec'])
    covariance=np.array(g['models']['conventional']['0.0']['covariance_mas2'])/1e6
    L=np.linalg.cholesky(covariance)
    a=old.load_astrometry()
    observed=dict(a)
    observed['observed_images_light_centered_arcsec']=np.vstack([np.zeros(2),y.reshape(3,2)])+a['observed_images_light_centered_arcsec'][0]
    params,dcov,_,light_likelihood=setup()
    dmean=np.array(params.delays_days[1:]); dL=np.linalg.cholesky(dcov)
    indices=[0,1,2,3,5,6,7,8]
    center=np.r_[p0[indices],11.]
    scales=np.array([.03,.03,.01,.01,.05,.05,.003,.003,.3])
    lower=np.r_[local.LOWER[indices],7.]
    upper=np.r_[local.UPPER[indices],13.]
    lo=(lower-center)/scales; hi=(upper-center)/scales
    def unpack(q):
        pp=p0.copy(); raw=center+scales*q;pp[indices]=raw[:-1]
        return pp,10**raw[-1]
    def apply(q,model='conventional',frozen_extra=None):
        pp,mass=unpack(q)
        return geometry.macro_application(pp,[-2.585,-3.637],perturber_mass_Msun=mass,
                                           astrometry=observed,model=model,frozen_extra=frozen_extra)
    def astro_terms(q,model='conventional',frozen_extra=None):
        pp,_=unpack(q); app=apply(q,model,frozen_extra)
        r=np.r_[np.linalg.solve(L,app['relative_astrometry_arcsec'][:6]-y),
                (app['relative_astrometry_arcsec'][6:]-a['observed_relative_arcsec'][6:])/a['sigma_arcsec'][6:],
                pp[7:9]/.003]
        return r,app
    def residual(q):
        try:
            r,app=astro_terms(q)
            return np.r_[r,np.linalg.solve(dL,app['delays_days']-dmean)]
        except (ArithmeticError,ValueError):
            return np.full(13,1e8)
    out={'scope':'Consumed Gaia development fit; no unused data consumed. kappa fixed to inherited value to remove mass-sheet degree; local parameters, G1 positive mass fitted.',
         'rear_galaxy_observed_arcsec':[-2.585,-3.637], 'rear_redshift':.7821,
         'compact_mass_Msun':1e8,'compact_scale_kpc':.03,'perturber_scale_kpc':.2,
         'Gaia_error_inflation':False, 'initializer_trials':[]}
    fits=[]
    for logM in [10.,11.,12.]:
        q=np.zeros(9);q[-1]=(logM-center[-1])/scales[-1]
        fit=least_squares(residual,q,bounds=(lo,hi),max_nfev=1800,ftol=1e-12,xtol=1e-12,gtol=2e-8,diff_step=1e-4,x_scale='jac')
        fits.append(fit)
        row={'logM_initial':logM,'chi2':float(fit.fun@fit.fun),'Gaia_chi2':float(fit.fun[:6]@fit.fun[:6]),
             'parameters':unpack(fit.x)[0].tolist(),'G1mass':unpack(fit.x)[1], 'success':bool(fit.success),
             'message':str(fit.message),'nfev':int(fit.nfev),'optimality':float(fit.optimality)}
        out['initializer_trials'].append(row);print(json.dumps(row),flush=True)
        Path(__file__).with_name('positive_g1_development_fit.json').write_text(json.dumps(out,indent=2)+'\n')
    q=min(fits,key=lambda f:f.fun@f.fun).x
    all_endpoints=[]
    for model in ['conventional','r1']:
        extra=None
        trace=[]
        iterations=1 if model=='conventional' else 3
        for iteration in range(iterations):
            if model=='r1':
                pp,_=unpack(q)
                screened=ProfileR1Lens(local.body_from_parameters(pp,1e8),old.frozen_coefficients())
                extra=old.FrozenExtraProjection(screened,count=400)
            def objective(v):
                try:
                    r,app=astro_terms(v,model,extra)
                    return float(r@r-2*light_likelihood(app['delays_days']))
                except (ArithmeticError,ValueError):
                    return 1e15
            fit=minimize(objective,q,method='L-BFGS-B',jac='3-point',bounds=list(zip(lo,hi)),
                         options={'maxiter':600,'ftol':1e-13,'gtol':1e-6,'maxls':40,'finite_diff_rel_step':1e-4})
            q=fit.x; pp,mass=unpack(q)
            trace.append({'iteration':iteration,'objective':float(fit.fun),'success':bool(fit.success),'message':str(fit.message),'nfev':int(fit.nfev),
                          'parameter_vector':pp.tolist(),'perturber_mass_Msun':float(mass)})
            print(json.dumps({'model':model,**trace[-1]}),flush=True)
        r,app=astro_terms(q,model,None)
        pp,mass=unpack(q)
        exact_ll=float(light_likelihood(app['delays_days']))
        def observables(v):
            return apply(v,model,extra)['relative_astrometry_arcsec'][:6]
        jac=[]
        for j in range(9):
            step=2e-5;delta=np.eye(9)[j]*step
            jac.append((observables(q+delta)-observables(q-delta))/(2*step))
        s=np.linalg.svd(np.linalg.solve(L,np.array(jac).T),compute_uv=False)
        row={'model':model,'parameter_vector':pp.tolist(),'parameter_names':list(old.PARAMETER_NAMES),
             'perturber_mass_Msun':mass,'predicted_relative_astrometry_arcsec':app['relative_astrometry_arcsec'].tolist(),
             'predicted_images_light_centered_arcsec':app['images_light_centered_arcsec'].tolist(),
             'predicted_delays_days':app['delays_days'].tolist(),'magnifications':app['magnifications'].tolist(),
             'Gaia_residual_mas':(1000*(y-app['relative_astrometry_arcsec'][:6])).tolist(),
             'Gaia_chi2':float(r[:6]@r[:6]),'rms_coordinate_residual_mas':float(1000*np.sqrt(np.mean((y-app['relative_astrometry_arcsec'][:6])**2))),
             'NIC2_lenslight_chi2':float(r[6:8]@r[6:8]),'masscenter_prior_chi2':float(r[8:]@r[8:]),
             'development_light_log_likelihood':exact_ll,'exact_training_objective':float(r@r-2*exact_ll),
             'max_lens_residual_rad':app['max_lens_residual_rad'], 'max_arrival_stationarity_Mpc':app['max_arrival_stationarity_Mpc'],
             'Gaia_jacobian_singular_values_scaled_coordinates':s.tolist(),'optimization_trace':trace,
             'active_bounds':{str(index):float(value) for index,value,l,u in zip(indices+[9],center+scales*q,lower,upper) if min(abs(value-l),abs(value-u))<1e-5},
             'rear_center_Mpc':app['system']['perturber_center_Mpc'].tolist(),'physical_admission':False}
        if model=='r1':
            impacts=np.linalg.norm(app['images_light_centered_arcsec']*ARCSEC*app['system']['distances']['ol']-np.asarray(app['system']['lens'].centre_Mpc),axis=1)
            row['original_source_checks']=app['system']['lens'].checks(impacts)
            row['original_section_source_applications']=local.source_section_checks(app['system']['lens'],[.00003,.0002,.001,.006])
        out[model]=row
        all_endpoints.append((q.copy(),row))
        Path(__file__).with_name('positive_g1_development_fit.json').write_text(json.dumps(out,indent=2)+'\n')
    cp=np.array(out['conventional']['parameter_vector']);cm=out['conventional']['perturber_mass_Msun']
    ca=geometry.macro_application(cp,[-2.585,-3.637],perturber_mass_Msun=cm)
    ra=geometry.macro_application(cp,[-2.585,-3.637],perturber_mass_Msun=cm,model='r1')
    out['fixed_body_R1_activation']={'image_offsets_microarcsec':((ra['images_light_centered_arcsec']-ca['images_light_centered_arcsec'])*1e6).tolist(),
                                   'delay_offsets_seconds':((ra['delays_days']-ca['delays_days'])*86400).tolist()}
    out['seconds']=time.perf_counter()-started
    Path(__file__).with_name('positive_g1_development_fit.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    from strong_lensing_comparison import ARCSEC
    main()
