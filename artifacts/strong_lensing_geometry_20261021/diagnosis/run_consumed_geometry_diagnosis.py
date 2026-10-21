"""Development diagnostics on already-consumed HE0435 observations only."""
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'code'))
import strong_lensing_observed_macro as old
import strong_lensing_macro_revision as m


def main():
    base = json.loads((ROOT/'artifacts/strong_lensing_revision_20261020/revised_conventional_training.json').read_text())
    p = np.array(base['parameter_vector'])
    g = json.loads((ROOT/'artifacts/strong_lensing_revision_20261020/gaia/gaia_astrometry_score.json').read_text())
    y = np.array(g['association']['relative_arcsec'])
    covariance = np.array(g['models']['conventional']['0.0']['covariance_mas2'])/1e6
    L = np.linalg.cholesky(covariance)
    a = old.load_astrometry()
    old_app = m.macro_application(p, compact_mass_Msun=1e8)
    out = {'scope':'DEVELOPMENT ONLY; Gaia and NIC2 have already been consumed. No unused observation is read.',
           'model_minus_NIC2_mas':(1000*(old_app['relative_astrometry_arcsec'][:6]-a['observed_relative_arcsec'][:6])).tolist(),
           'Gaia_minus_NIC2_mas':(1000*(y-a['observed_relative_arcsec'][:6])).tolist(),
           'trials':[]}
    # Reproduce the known nearby branches with exact original lens application.
    gaia_astrometry = dict(a)
    gaia_astrometry['observed_images_light_centered_arcsec'] = np.vstack([np.zeros(2),y.reshape(3,2)])+a['observed_images_light_centered_arcsec'][0]
    for name, use, include_lens, include_prior, delay_sigma in [
        ('Gaia_only_nine_free', list(range(9)),False,False,None),
        ('Gaia_six_source_screen', [2,3,4,5,6,7],False,False,None),
        ('Gaia_plus_NIC2_lightcenter',list(range(9)),True,True,None),
        ('Gaia_plus_NIC2_lightcenter_old_delays_1day',list(range(9)),True,True,1.),
        ('Gaia_fixed_center_and_screenkappa', [0,1,2,3,5,6],False,False,None),
    ]:
        q0=p[use].copy(); scales=np.array([.03,.03,.01,.01,.1,.05,.05,.003,.003])[use]
        qlo=(m.LOWER[use]-q0)/scales; qhi=(m.UPPER[use]-q0)/scales
        def parameters(q):
            r=p.copy(); r[use]=q0+scales*q; return r
        def residual(q):
            r=parameters(q)
            try:
                app=m.macro_application(r,compact_mass_Msun=1e8,astrometry=gaia_astrometry)
                v=[np.linalg.solve(L,app['relative_astrometry_arcsec'][:6]-y)]
                if include_lens:
                    v.append((app['relative_astrometry_arcsec'][6:]-a['observed_relative_arcsec'][6:])/a['sigma_arcsec'][6:])
                if include_prior:
                    v.append(r[7:9]/.003)
                if delay_sigma is not None:
                    v.append((app['delays_days']-old_app['delays_days'])/delay_sigma)
                return np.concatenate(v)
            except (ArithmeticError, ValueError):
                return np.full(6+2*include_lens+2*include_prior+3*(delay_sigma is not None),1e8)
        t=time.perf_counter()
        fit=least_squares(residual,np.zeros(len(use)),bounds=(qlo,qhi),max_nfev=2000,
                          ftol=2e-13,xtol=2e-13,gtol=2e-9,diff_step=1e-4,x_scale='jac')
        pp=parameters(fit.x); app=m.macro_application(pp,compact_mass_Msun=1e8,astrometry=gaia_astrometry)
        u,s,v=np.linalg.svd(fit.jac,full_matrices=False)
        row={'name':name,'free_parameter_indices':use,'parameter_vector':pp.tolist(),
             'chi2':float(fit.fun@fit.fun),'Gaia_chi2':float(fit.fun[:6]@fit.fun[:6]),
             'residual_mas':(1000*(y-app['relative_astrometry_arcsec'][:6])).tolist(),
             'rms_mas':float(1000*np.sqrt(np.mean((y-app['relative_astrometry_arcsec'][:6])**2))),
             'delays_days':app['delays_days'].tolist(),'masscenter_mas':(1000*pp[7:9]).tolist(),
             'lensresidual_rad':app['max_lens_residual_rad'], 'success':bool(fit.success),'optimality':float(fit.optimality),
             'message':str(fit.message),'nfev':int(fit.nfev),'jacobian_singular_values_scaled_coordinates':s.tolist(),
             'seconds':time.perf_counter()-t}
        out['trials'].append(row)
        print(json.dumps(row),flush=True)
        Path(__file__).with_name('consumed_geometry_diagnosis.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
