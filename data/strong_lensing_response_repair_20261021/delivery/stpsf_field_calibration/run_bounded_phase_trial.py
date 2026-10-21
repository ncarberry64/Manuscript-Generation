"""Field-only low-order pupil-phase development trial, with exact replay.

Finite-difference intensities propose bounded wavefront coefficients.  Results
are evaluated using nonlinear optical propagation; the linear intensity model
is not accepted as a physical PSF or a validation result.
"""
from pathlib import Path
import os,sys,json,time,hashlib,logging,copy
ROOT=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent/'phase_trial'
OUT.mkdir(exist_ok=True)
os.environ['STPSF_PATH']=str(ROOT/'stpsf_reference/stpsf-data')
os.environ['MPLCONFIGDIR']=str(ROOT/'matplotlib_cache')
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from astropy.io import fits
from astropy.config.paths import set_temp_cache,set_temp_config
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares
sys.path.insert(0,str(ROOT/'mask_patch/code'));sys.path.insert(0,str(ROOT/'repo/code'))
from strong_lensing_jwst_rate_mask import fit_jwst_native_rate_point

OPD=ROOT/'stpsf_trial/wss_opds/R2022101203-NRCA3_FP1-1.fits'
NOLL=list(range(4,12));STEP=2.;BOUND=10.;PRIOR=5.;OV=5;FOV=25
PROTO={
 'status':'FIELD_ONLY_PHASE_DEVELOPMENT_PROTOCOL',
 'noll_modes':NOLL,'coefficient_bounds_nm':[-BOUND,BOUND],
 'assumed_gaussian_regularizer_sigma_nm':PRIOR,
 'regularizer_scope':'Exploratory assumption, not calibrated optical uncertainty',
 'finite_difference_step_nm':STEP,'source_blackbody_K':4000.,
 'diffusion_sigma_mas':6.2,'jitter_sigma_mas':.8,
 'source_and_diffusion_fixed':'Original prescribed values; all-four color fit not reused in LOO',
 'pupil_mode_definition':'Noll Zernikes on the loaded entrance-pupil grid, piston removed and each normalized to unit RMS on illuminated pixels',
 'mode_location':'Added to measured global OTE-only pupil OPD in meters; nominal OTE and NIRCam field dependence retained',
 'WSS_file':str(OPD),'WSS_sha256':hashlib.sha256(OPD.read_bytes()).hexdigest(),
 'fft_oversample':OV,'detector_oversample':OV,'fov_pixels':FOV,'wavelength_bins':9,
 'objective':'Native-ERR pixel residuals plus assumed Gaussian phase regularizer; six measurement nuisances per training exposure',
 'identifiability_gate':'Eight projected phase-Jacobian modes must have full rank at relative1e-10 and max conditional coefficient standard error below10nm, without residual rescaling',
 'validation':'Four leave-one-exposure-out conditional field-shape checks using exact nonlinear PSFs; all data are already development',
 'errors_scaled':False,'target_pixels_used':False,'independent_observational_validation':False,
 'acceptance':'No physical-response acceptance from linear intensity; exact nonlinear per-exposure native-ERR p>=.01 and successful unconstrained measurement fit required',
}
(OUT/'protocol.json').write_text(json.dumps(PROTO,indent=2)+'\n')
for name in ['astropy_cache','astropy_config']:(OUT/name).mkdir(exist_ok=True)

def slim(f):
 return {k:v for k,v in f.items() if not k.startswith('pixel_') and k!='standardized_residual'}

def main():
 start=time.time()
 with set_temp_cache(str(OUT/'astropy_cache')),set_temp_config(str(OUT/'astropy_config')):
  import stpsf
  from poppy.zernike import zernike1
  from synphot import SourceSpectrum
  from synphot.models import BlackBodyNorm1D
  logging.getLogger('stpsf').setLevel(logging.ERROR);logging.getLogger('poppy').setLevel(logging.ERROR)
  rows=json.loads((ROOT/'repo/artifacts/strong_lensing_geometry_20261021/standard_psf/field_psf_calibration.json').read_text())['field_fits']
  fields=[];instruments=[];opds=[];bases=[];derivatives=[]
  modes=None
  for i,row in enumerate(rows):
   guess=np.array(row['centroid_xy_pixel']);ix,iy=np.rint(guess).astype(int);offset=np.array([ix-8,iy-8])
   x0,y0=offset;x1,y1=offset+17
   path=ROOT/'repo/data/strong_lensing_geometry_20261021/native'/row['exposure']
   with fits.open(path,memmap=False) as h:
    hd=dict(h[0].header)
    data=np.array(h['SCI'].section[y0:y1,x0:x1],float)
    err=np.array(h['ERR'].section[y0:y1,x0:x1],float)
    dq=np.array(h['DQ'].section[y0:y1,x0:x1],np.uint32)
   fields.append(dict(exposure=row['exposure'],guess=guess-offset,offset=offset,header=hd,data=data,error=err,dq=dq))
   nrc=stpsf.NIRCam();nrc.detector='NRCB3';nrc.filter='F115W';nrc.detector_position=tuple(guess)
   nrc.options['charge_diffusion_sigma']=.0062;nrc.options['jitter_sigma']=.0008
   nrc.load_wss_opd(str(OPD),plot=False,verbose=False)
   opd=copy.deepcopy(nrc.pupilopd);assert opd[0].header['BUNIT'] in ('m','meter')
   if modes is None:
    npix=opd[0].data.shape[0]
    mask=stpsf.utils.get_pupil_mask(npix=npix)>0
    z=[]
    for j in NOLL:
     a=zernike1(j,npix=npix,outside=0.)*mask
     a[mask]-=np.mean(a[mask]);a/=np.sqrt(np.mean(a[mask]**2))
     z.append(a)
    modes=np.array(z)
    gram=np.array([[np.mean(a[mask]*b[mask]) for b in modes] for a in modes])
    np.savez_compressed(OUT/'pupil_modes.npz',modes=modes,mask=mask,noll=NOLL,gram=gram)
   instruments.append(nrc);opds.append(opd)
  source=SourceSpectrum(BlackBodyNorm1D,temperature=4000.)
  waves,weights=instruments[0]._get_weights(source=source,nlambda=9)
  spectral={'wavelengths':waves,'weights':weights}

  def exact_psf(i,coef):
   nrc=instruments[i]
   modified=copy.deepcopy(opds[i])
   modified[0].data=np.array(opds[i][0].data)+np.einsum('i,ijk->jk',coef,modes)*1e-9
   nrc.pupilopd=modified
   g=nrc.psf_grid(num_psfs=1,all_detectors=False,single_psf_centered=False,
       oversample=OV,fov_pixels=FOV,nlambda=9,source=spectral,verbose=False)
   assert np.all(np.asarray(g.oversampling)==OV)
   return np.array(g.data[0])

  def measure(i,psf):
   f=fields[i];origin=((psf.shape[1]-1)/2,(psf.shape[0]-1)/2)
   fit=fit_jwst_native_rate_point(f['data'],f['error'],f['dq'],f['guess'],psf,OV,origin,
       primary_header=f['header'],background_degree=1)
   return fit

  def receipt(i,fit):
   record=slim(fit)
   offset=fields[i]['offset']
   record['centroid_xy_pixel']=(np.array(fit['centroid_xy_pixel'])+offset).tolist()
   record['stamp_bounds_xyxy']=(np.array(fit['stamp_bounds_xyxy'])+np.tile(offset,2)).tolist()
   record['coordinate_frame']='Global native detector pixels, zero-based'
   record['exposure']=fields[i]['exposure']
   return record

  signatures={'protocol_sha256':hashlib.sha256((OUT/'protocol.json').read_bytes()).hexdigest(),
              'stpsf_version':stpsf.__version__,'source_wavelengths':waves.tolist(),'weights':weights.tolist(),
              'positions':[list(n.detector_position) for n in instruments]}
  config=json.dumps(signatures,sort_keys=True)
  for i in range(4):
   cache=OUT/f'field_{i+1}_phase_basis.npz'
   if cache.exists():
    c=np.load(cache);assert str(c['configuration'])==config
    base=c['base'];deriv=c['derivatives']
   else:
    base=exact_psf(i,np.zeros(8));deriv=[]
    for j in range(8):
     c=np.zeros(8);c[j]=STEP
     plus=exact_psf(i,c);minus=exact_psf(i,-c)
     deriv.append((plus-minus)/(2*STEP))
    deriv=np.array(deriv)
    np.savez_compressed(cache,base=base,derivatives=deriv,configuration=config)
   bases.append(base);derivatives.append(deriv)
   print('PHASE_BASIS',i+1,'elapsed',round(time.time()-start,2),flush=True)
  baseline=[measure(i,bases[i]) for i in range(4)]
  interps=[]
  for i in range(4):
   f=fields[i];bf=baseline[i]
   f['x']=np.array(bf['pixel_x']);f['y']=np.array(bf['pixel_y'])
   f['raw']=np.array(bf['pixel_data']);f['sigma']=np.array(bf['pixel_error'])
   f['bg']=np.column_stack([np.ones(len(f['x'])),f['x']-f['guess'][0],f['y']-f['guess'][1]])
   f['initial']=np.r_[bf['centroid_xy_pixel'],np.log(bf['flux_in_SCI_pixel_units']),bf['background_plane_SCI']]
   interps.append([RectBivariateSpline(np.arange(FOV*OV),np.arange(FOV*OV),a.T,kx=3,ky=3,s=0)
                   for a in [bases[i],*derivatives[i]]])

  def field_arrays(i,coef,nuisance):
   f=fields[i];interp=interps[i]
   xx=(f['x']-nuisance[0])*OV+(FOV*OV-1)/2
   yy=(f['y']-nuisance[1])*OV+(FOV*OV-1)/2
   values=np.array([g.ev(xx,yy) for g in interp])
   dpx=np.array([g.ev(xx,yy,dx=1) for g in interp])
   dpy=np.array([g.ev(xx,yy,dy=1) for g in interp])
   inside=(xx>=0)&(yy>=0)&(xx<=FOV*OV-1)&(yy<=FOV*OV-1)
   values[:,~inside]=0.;dpx[:,~inside]=0.;dpy[:,~inside]=0.
   comb=np.r_[1.,coef]
   p=comb@values;flux=np.exp(nuisance[2])
   pred=flux*p+f['bg']@nuisance[3:]
   jp=flux*values[1:].T/f['sigma'][:,None]
   jn=np.column_stack([-flux*OV*(comb@dpx),-flux*OV*(comb@dpy),flux*p,f['bg']])/f['sigma'][:,None]
   residual=(pred-f['raw'])/f['sigma']
   return residual,jp,jn,p

  def identify(indices):
   projected=[]
   for i in indices:
    _,jp,jn,_=field_arrays(i,np.zeros(8),fields[i]['initial'])
    projected.append(jp-jn@np.linalg.lstsq(jn,jp,rcond=1e-12)[0])
   a=np.vstack(projected);_,sing,vh=np.linalg.svd(a,full_matrices=False)
   if sing[-1]>sing[0]*1e-10:
    cov=(vh.T/(sing*sing))@vh
    stderr=np.sqrt(np.diag(cov))
   else:
    stderr=np.full(8,np.inf)
   return {'singular_values_per_nm':sing.tolist(),'rank':int(np.sum(sing>sing[0]*1e-10)),
           'condition_number':float(sing[0]/sing[-1]),'conditional_coefficient_stderr_nm':stderr.tolist(),
           'gate_pass':bool(sing[-1]>sing[0]*1e-10 and stderr.max()<BOUND)}

  def solve(indices):
   ident=identify(indices)
   if not ident['gate_pass']:
    return {'identifiability':ident,'status':'STOPPED_NOT_IDENTIFIABLE'}
   initial=np.r_[np.zeros(8),np.concatenate([fields[i]['initial'] for i in indices])]
   lower=list(np.full(8,-BOUND));upper=list(np.full(8,BOUND))
   for i in indices:
    x,y=fields[i]['guess'];lower.extend([x-5,y-5,np.log(1e-8),-np.inf,-np.inf,-np.inf])
    upper.extend([x+5,y+5,np.log(1e8),np.inf,np.inf,np.inf])
   def unpack(p):
    return [field_arrays(i,p[:8],p[8+6*k:14+6*k]) for k,i in enumerate(indices)]
   def residual(p):
    return np.r_[np.concatenate([r[0] for r in unpack(p)]),p[:8]/PRIOR]
   def jacobian(p):
    arrays=unpack(p);n=sum(len(a[0]) for a in arrays)
    jac=np.zeros((n+8,len(p)));offset=0
    for k,(_,jp,jn,_) in enumerate(arrays):
     m=len(jp);jac[offset:offset+m,:8]=jp;jac[offset:offset+m,8+6*k:14+6*k]=jn;offset+=m
    jac[n:,:8]=np.eye(8)/PRIOR
    return jac
   fit=least_squares(residual,initial,jac=jacobian,bounds=(lower,upper),x_scale='jac',
       ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=150)
   arrays=unpack(fit.x)
   return {'status':'BOUNDED_LINEARIZED_PHASE_PROPOSAL','identifiability':ident,
           'coefficients_nm':fit.x[:8].tolist(),'optimizer_success':bool(fit.success),
           'optimizer_message':str(fit.message),'nfev':int(fit.nfev),
           'active_bound_indices':np.flatnonzero(fit.active_mask).tolist(),
           'linearized_training_chi2':float(sum(a[0]@a[0] for a in arrays)),
           'phase_regularizer_chi2':float(np.sum((fit.x[:8]/PRIOR)**2)),
           'linearized_negative_profile_sample_count':int(sum(np.sum(a[3]<0) for a in arrays))}

  results={'protocol':PROTO,'signatures':signatures,'baseline_field_fits':[receipt(i,f) for i,f in enumerate(baseline)],
           'all_fields':solve([0,1,2,3]),'leave_one_out':[],'target_pixels_used':False,
           'native_errors_scaled':False,'independent_observational_validation':False}
  global_fit=results['all_fields']
  print('PHASE_SOLUTION',json.dumps(global_fit),flush=True)
  if global_fit['status']=='STOPPED_NOT_IDENTIFIABLE':
   (OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n');return
  exact_global=[]
  coeff=np.array(global_fit['coefficients_nm'])
  perturbation=np.einsum('i,ijk->jk',coeff,modes)
  results['global_perturbation_pupil_rms_nm']=float(np.sqrt(np.mean(perturbation[mask]**2)))
  for i in range(4):
   psf=exact_psf(i,coeff);f=measure(i,psf)
   f['exposure']=fields[i]['exposure'];f['chi2_dof']=f['chi2']/f['dof']
   f['exact_nonlinear_phase']=True
   exact_global.append(receipt(i,f))
   fits.PrimaryHDU(psf).writeto(OUT/f'field_{i+1}_global_exact_epsf.fits',overwrite=True)
  results['global_exact_field_fits']=exact_global
  print('EXACT_GLOBAL',[f['chi2_dof'] for f in exact_global],flush=True)
  for held in range(4):
   train=[i for i in range(4) if i!=held]
   record=solve(train);record['held_out_exposure']=fields[held]['exposure']
   if record['status']!='STOPPED_NOT_IDENTIFIABLE':
    psf=exact_psf(held,np.array(record['coefficients_nm']));f=measure(held,psf)
    f['exposure']=fields[held]['exposure'];f['chi2_dof']=f['chi2']/f['dof']
    record['held_out_exact_fit']=receipt(held,f)
    record['absolute_adequacy_pass']=bool(f['p_value']>=.01 and f['valid_measurement_fit'] and record['optimizer_success'])
    print('PHASE_LOO',held+1,record['coefficients_nm'],f['chi2_dof'],flush=True)
   results['leave_one_out'].append(record)
  results['all_loo_absolute_adequacy_pass']=all(r.get('absolute_adequacy_pass',False) for r in results['leave_one_out'])
  results['elapsed_seconds']=time.time()-start
  (OUT/'results.json').write_text(json.dumps(results,indent=2)+'\n')
  print('PHASE_FINISHED',round(results['elapsed_seconds'],2),results['all_loo_absolute_adequacy_pass'],flush=True)

if __name__=='__main__':main()
