"""Field-only STPSF color/diffusion development calibration, with LOO checks.

No lens parameters or quasar pixels enter this calculation.  The native ERR
arrays and verified usable-slope mask are unchanged.  This is development,
not independent observational validation.  All writes remain beside this file.
"""
from pathlib import Path
import os, sys, json, time, hashlib, logging

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
os.environ['STPSF_PATH'] = str(ROOT/'stpsf_reference/stpsf-data')
os.environ['MPLCONFIGDIR'] = str(ROOT/'matplotlib_cache')
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import numpy as np
from astropy.io import fits
from astropy.config.paths import set_temp_cache, set_temp_config
from scipy.ndimage import gaussian_filter
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize
from scipy.stats import chi2

sys.path.insert(0, str(ROOT/'mask_patch/code'))
sys.path.insert(0, str(ROOT/'repo/code'))
from strong_lensing_jwst_rate_mask import fit_jwst_native_rate_point

TBOUNDS = (2500., 12000.)
SBOUNDS = (2., 12.)  # milli-arcseconds, an explicitly bounded exploratory family
OVERSAMPLE = 5
FOV = 25
NWAVES = 9
OPD = ROOT/'stpsf_trial/wss_opds/R2022101203-NRCA3_FP1-1.fits'
PROTO = {
    'status': 'FIELD_ONLY_DEVELOPMENT_PROTOCOL',
    'star_temperature_K_bounds': TBOUNDS,
    'charge_diffusion_sigma_mas_bounds': SBOUNDS,
    'fixed_jitter_sigma_mas': .8,
    'fixed_WSS_file': str(OPD),
    'WSS_sha256': hashlib.sha256(OPD.read_bytes()).hexdigest(),
    'detector': 'NRCB3', 'filter': 'F115W',
    'fft_oversample': OVERSAMPLE, 'detector_oversample': OVERSAMPLE,
    'fov_pixels': FOV, 'wavelength_bins': NWAVES,
    'fitted_global_parameters': ['one stellar blackbody temperature',
                                 'one charge diffusion width shared by four exposures'],
    'exposure_nuisances': ['centroid_x', 'centroid_y', 'positive flux',
                          'affine background constant', 'affine background dx',
                          'affine background dy'],
    'baseline': {'temperature_K': 4000., 'charge_diffusion_sigma_mas': 6.2},
    'coarse_grid': {'temperature_K': [2500., 4000., 6500., 12000.],
                    'diffusion_sigma_mas': [2., 4., 6.2, 9., 12.]},
    'optimization': 'Bounded Nelder-Mead in log temperature and sigma_mas; profile exposure nuisances with unchanged archived fitter',
    'mask': 'Audited native CAL2.0.1 any-group SAT semantics; no DNU or invalid ERR restoration',
    'stamp': 'Existing isolated-star association; same 17x17 stamp per exposure',
    'errors_scaled': False, 'target_pixels_used': False,
    'validation': 'Four leave-one-exposure-out field checks; held-out exposure fits only its six measurement nuisances',
    'interpretation': 'All four field exposures were already inspected; LOO is an internal development check, never new independent validation',
    'acceptance': 'Pixel chi-square p>=0.01 per held-out exposure with native ERR, successful response optimization, and valid measurement fit; also report boundaries and residual maps',
    'efficiency': 'Nine monochromatic STPSF grids per field position, exact spectral superposition, Gaussian diffusion commuting with pixel/IPC convolutions in the central fitting support; verify against direct broadband STPSF',
}
(OUT/'protocol.json').write_text(json.dumps(PROTO, indent=2)+'\n')
for d in ['astropy_cache','astropy_config']:
    (OUT/d).mkdir(exist_ok=True)

def slim(fit):
    return {k:v for k,v in fit.items() if not k.startswith('pixel_') and k != 'standardized_residual'}

def main():
    start = time.time()
    with set_temp_cache(str(OUT/'astropy_cache')), set_temp_config(str(OUT/'astropy_config')):
        import stpsf
        from synphot import SourceSpectrum
        from synphot.models import BlackBodyNorm1D
        logging.getLogger('stpsf').setLevel(logging.ERROR)
        logging.getLogger('poppy').setLevel(logging.ERROR)
        logging.getLogger('numexpr').setLevel(logging.ERROR)
        rows = json.loads((ROOT/'repo/artifacts/strong_lensing_geometry_20261021/standard_psf/field_psf_calibration.json').read_text())['field_fits']
        fields, instruments, mono_cubes, pixel_scales = [], [], [], []
        nrc0 = stpsf.NIRCam()
        nrc0.detector = 'NRCB3'; nrc0.filter = 'F115W'
        waves, weights0 = nrc0._get_weights(source=SourceSpectrum(BlackBodyNorm1D, temperature=4000), nlambda=NWAVES)
        log_temperatures = np.linspace(np.log(TBOUNDS[0]), np.log(TBOUNDS[1]), 61)
        weight_grid = []
        weight_path = OUT/'spectral_weights.npz'
        weight_configuration = json.dumps({
            'stpsf_version': stpsf.__version__, 'detector': nrc0.detector,
            'filter': nrc0.filter, 'source_model': 'BlackBodyNorm1D',
            'wavelengths_m': waves.tolist(),
            'log_temperatures': log_temperatures.tolist(),
        }, sort_keys=True)
        weights_cached = False
        if weight_path.exists():
            with np.load(weight_path) as saved:
                if 'configuration' in saved and str(saved['configuration']) == weight_configuration:
                    assert np.array_equal(saved['wavelengths_m'],waves)
                    assert np.array_equal(saved['log_temperatures'],log_temperatures)
                    weight_grid=saved['weights'].copy()
                    weights_cached = True
        if not weights_cached:
            for lt in log_temperatures:
                wv, wt = nrc0._get_weights(source=SourceSpectrum(BlackBodyNorm1D, temperature=float(np.exp(lt))), nlambda=NWAVES)
                assert np.array_equal(waves, wv)
                weight_grid.append(wt)
        weight_interp = PchipInterpolator(log_temperatures, np.array(weight_grid), axis=0)
        np.savez_compressed(OUT/'spectral_weights.npz', wavelengths_m=waves,
                            log_temperatures=log_temperatures, weights=weight_grid,
                            configuration=weight_configuration)
        for i, row in enumerate(rows):
            path = ROOT/'repo/data/strong_lensing_geometry_20261021/native'/row['exposure']
            guess = np.array(row['centroid_xy_pixel'])
            ix, iy = np.rint(guess).astype(int)
            bounds = (ix-8, iy-8, ix+9, iy+9)
            x0,y0,x1,y1 = bounds
            with fits.open(path, memmap=False) as hd:
                header = dict(hd[0].header)
                data = np.array(hd['SCI'].section[y0:y1,x0:x1], dtype=float)
                error = np.array(hd['ERR'].section[y0:y1,x0:x1], dtype=float)
                dq = np.array(hd['DQ'].section[y0:y1,x0:x1], dtype=np.uint32)
            fields.append(dict(exposure=row['exposure'], guess=guess-[x0,y0],
                               offset=np.array([x0,y0]), header=header,
                               data=data, error=error, dq=dq))
            nrc = stpsf.NIRCam(); nrc.detector='NRCB3'; nrc.filter='F115W'
            nrc.detector_position=tuple(guess)
            nrc.options['jitter_sigma'] = .0008
            nrc.options['charge_diffusion_sigma'] = 0.
            nrc.load_wss_opd(str(OPD), plot=False, verbose=False)
            pixel_scales.append(float(nrc.pixelscale)/OVERSAMPLE)
            instruments.append(nrc)
            cache = OUT/f'field_{i+1}_monochromatic.npz'
            configuration = json.dumps({
                'stpsf_version': stpsf.__version__,
                'WSS_sha256': PROTO['WSS_sha256'],
                'detector': nrc.detector, 'filter': nrc.filter,
                'detector_position': list(nrc.detector_position),
                'wavelengths_m': waves.tolist(),
                'fft_oversample': OVERSAMPLE, 'detector_oversample': OVERSAMPLE,
                'fov_pixels': FOV, 'jitter_sigma_arcsec': .0008,
                'charge_diffusion_sigma_arcsec': 0.,
                'nrc_pixelscale': float(nrc.pixelscale),
            }, sort_keys=True)
            cube = None
            if cache.exists():
                with np.load(cache) as saved:
                    if 'configuration' in saved and str(saved['configuration']) == configuration:
                        cube = saved['cube'].copy()
            if cube is None:
                basis = []
                for j, wave in enumerate(waves):
                    grid = nrc.psf_grid(num_psfs=1, all_detectors=False,
                        single_psf_centered=False, oversample=OVERSAMPLE,
                        fov_pixels=FOV, monochromatic=float(wave), verbose=False)
                    assert grid.data.shape == (1,FOV*OVERSAMPLE,FOV*OVERSAMPLE)
                    assert np.all(np.asarray(grid.oversampling) == OVERSAMPLE)
                    basis.append(np.array(grid.data[0]))
                cube = np.array(basis)
                np.savez_compressed(cache, cube=cube, wavelengths_m=waves,
                                    detector_position=nrc.detector_position,
                                    oversample=OVERSAMPLE, nrc_pixelscale=nrc.pixelscale,
                                    configuration=configuration)
            mono_cubes.append(cube)
            print('BASIS',i+1,'elapsed_seconds',round(time.time()-start,2),flush=True)

        def exact_weights(temp):
            _, w = nrc0._get_weights(source=SourceSpectrum(BlackBodyNorm1D, temperature=float(temp)), nlambda=NWAVES)
            return w

        def profile_model(field_index, logtemp, width_mas, exact=False):
            weights = exact_weights(np.exp(logtemp)) if exact else weight_interp(logtemp)
            weights = weights/weights.sum()
            base = np.einsum('i,ijk->jk', weights, mono_cubes[field_index])
            # Same Gaussian routine used by STPSF, after its field distortion.
            epsf = gaussian_filter(base, (width_mas/1000)/pixel_scales[field_index])
            return epsf

        def fit_one(i, params, exact=False):
            f=fields[i]; epsf=profile_model(i,*params,exact=exact)
            origin=((epsf.shape[1]-1)/2,(epsf.shape[0]-1)/2)
            fit=fit_jwst_native_rate_point(f['data'],f['error'],f['dq'],f['guess'],
                epsf,OVERSAMPLE,origin,primary_header=f['header'],background_degree=1)
            fit['centroid_xy_pixel']=(np.array(fit['centroid_xy_pixel'])+f['offset']).tolist()
            fit['pixel_x']=(np.asarray(fit['pixel_x'])+f['offset'][0]).tolist()
            fit['pixel_y']=(np.asarray(fit['pixel_y'])+f['offset'][1]).tolist()
            fit['stamp_bounds_xyxy']=(np.asarray(fit['stamp_bounds_xyxy'])+
                np.tile(f['offset'],2)).tolist()
            fit['coordinate_frame']='Global native detector pixels, zero-based'
            fit['exposure']=f['exposure']; fit['chi2_dof']=fit['chi2']/fit['dof']
            return fit

        # Check central-support equivalence to root's full direct baseline.
        baseline_params=(np.log(4000.),6.2)
        fast_baseline=profile_model(0,*baseline_params,exact=True)
        direct_path=ROOT/'stpsf_equal_sampling'/rows[0]['exposure'].replace('_cal.fits','_FIELD_stpsf_epsf.fits')
        direct=fits.getdata(direct_path)
        support=slice(17,108)  # Covers fitting stamp and bounded subpixel offsets.
        equivalence={
            'central_max_abs_difference':float(np.max(np.abs(fast_baseline[support,support]-direct[support,support]))),
            'central_max_difference_over_peak':float(np.max(np.abs(fast_baseline[support,support]-direct[support,support]))/direct.max()),
            'direct_baseline_file':str(direct_path),
        }
        assert equivalence['central_max_difference_over_peak'] < 1e-6, equivalence
        baseline_fits=[fit_one(i,baseline_params,exact=True) for i in range(4)]
        print('BASELINE',[round(f['chi2_dof'],5) for f in baseline_fits],flush=True)
        cache={}
        def evaluate(params):
            p=np.asarray(params,dtype=float)
            key=tuple(np.round(p,9))
            if key not in cache:
                cache[key]=[fit_one(i,p) for i in range(4)]
            return cache[key]

        bounds=[(np.log(TBOUNDS[0]),np.log(TBOUNDS[1])),SBOUNDS]
        coarse=[]
        for temp in PROTO['coarse_grid']['temperature_K']:
            for sigma in PROTO['coarse_grid']['diffusion_sigma_mas']:
                params=np.array([np.log(temp),sigma])
                fits_=evaluate(params)
                coarse.append({'params':params.tolist(), 'chi2':[f['chi2'] for f in fits_]})
        (OUT/'coarse_grid.json').write_text(json.dumps(coarse,indent=2)+'\n')
        print('COARSE COMPLETE',round(time.time()-start,2),flush=True)

        def select(indices):
            first=min(coarse,key=lambda row:sum(row['chi2'][i] for i in indices))
            def objective(p):
                return sum(evaluate(p)[i]['chi2'] for i in indices)
            opt=minimize(objective,first['params'],method='Nelder-Mead',bounds=bounds,
                         options={'xatol':2e-5,'fatol':1e-5,'maxiter':180,'maxfev':300})
            p=np.asarray(opt.x)
            fits_=[fit_one(i,p,exact=True) for i in range(4)]
            return opt,p,fits_

        opt,params,final_fits=select([0,1,2,3])
        results={
            'status':'FIELD_ONLY_PHYSICAL_RESPONSE_DEVELOPMENT_CALIBRATION',
            'protocol':PROTO, 'stpsf_version':stpsf.__version__,
            'baseline_equivalence':equivalence,
            'baseline_fits':[slim(f) for f in baseline_fits],
            'selected_temperature_K':float(np.exp(params[0])),
            'selected_diffusion_sigma_mas':float(params[1]),
            'global_optimizer_success':bool(opt.success),
            'global_optimizer_message':str(opt.message),
            'boundary_parameters':[name for k,name in enumerate(['log_temperature','diffusion_sigma'])
                if min(abs(params[k]-bounds[k][0]),abs(params[k]-bounds[k][1]))<1e-3],
            'final_field_fits':[slim(f) for f in final_fits],
            'joint_chi2':float(sum(f['chi2'] for f in final_fits)),
            'joint_dof':int(sum(f['dof'] for f in final_fits)-2),
            'errors_scaled':False,'target_pixels_used':False,
            'independent_observational_validation':False,
            'leave_one_out':[],
        }
        print('GLOBAL',results['selected_temperature_K'],results['selected_diffusion_sigma_mas'],
              [round(f['chi2_dof'],5) for f in final_fits],flush=True)
        for held in range(4):
            fit_indices=[i for i in range(4) if i!=held]
            loo_opt,lp,loo_fits=select(fit_indices)
            hf=loo_fits[held]
            results['leave_one_out'].append({
                'held_out_exposure':fields[held]['exposure'],
                'temperature_K':float(np.exp(lp[0])),
                'diffusion_sigma_mas':float(lp[1]),
                'optimizer_success':bool(loo_opt.success),
                'training_chi2':float(sum(loo_fits[i]['chi2'] for i in fit_indices)),
                'training_dof':int(sum(loo_fits[i]['dof'] for i in fit_indices)-2),
                'held_out_fit':slim(hf),
                'held_out_absolute_adequacy_pass':bool(loo_opt.success and
                    hf['valid_measurement_fit'] and hf['p_value']>=.01),
                'held_out_centroid_minus_allfit_pixel':(np.array(hf['centroid_xy_pixel'])-np.array(final_fits[held]['centroid_xy_pixel'])).tolist(),
            })
            print('LOO',held+1,float(np.exp(lp[0])),float(lp[1]),hf['chi2_dof'],flush=True)
        results['joint_p_value']=float(chi2.sf(results['joint_chi2'],results['joint_dof']))
        results['all_field_loo_adequacy_pass']=all(r['held_out_absolute_adequacy_pass'] for r in results['leave_one_out'])
        # Independent direct broadband generation checks the accelerated response.
        direct_checks=[]
        for i,nrc in enumerate(instruments):
            nrc.options['charge_diffusion_sigma']=float(params[1]/1000)
            direct_grid=nrc.psf_grid(num_psfs=1,all_detectors=False,single_psf_centered=False,
                oversample=OVERSAMPLE,fov_pixels=FOV,nlambda=NWAVES,
                source=SourceSpectrum(BlackBodyNorm1D,temperature=float(np.exp(params[0]))),verbose=False)
            da=np.array(direct_grid.data[0]);fa=profile_model(i,*params,exact=True)
            delta=float(np.max(np.abs(da[support,support]-fa[support,support]))/da.max())
            assert delta<1e-6,delta
            direct_checks.append({'exposure':fields[i]['exposure'],'max_relative_peak_difference_in_fitting_support':delta})
            hd=fits.PrimaryHDU(da);hd.header['OVERSMP']=OVERSAMPLE
            hd.header['ORIGINX']=(da.shape[1]-1)/2;hd.header['ORIGINY']=(da.shape[0]-1)/2
            hd.header['PIXINT']=True;hd.header['TEMPK']=float(np.exp(params[0]))
            hd.header['DIFFMAS']=float(params[1]);hd.header['SCOPE']='FIELD DEVELOPMENT'
            hd.writeto(OUT/f'field_{i+1}_selected_epsf.fits',overwrite=True)
        results['direct_broadband_checks']=direct_checks
        results['evaluated_parameter_points']=len(cache)
        results['elapsed_seconds']=time.time()-start
        (OUT/'results.json').write_text(json.dumps(results,indent=2,default=str)+'\n')
        (OUT/'final_field_pixel_receipts.json').write_text(json.dumps(final_fits,indent=2,default=str)+'\n')
        print('FINISHED',json.dumps({k:results[k] for k in ['selected_temperature_K','selected_diffusion_sigma_mas','joint_chi2','joint_dof','all_field_loo_adequacy_pass','elapsed_seconds']}),flush=True)

if __name__=='__main__':
    main()
