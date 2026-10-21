"""Development test of dated physical telescope response on the isolated star."""
from pathlib import Path
import os,json,sys,hashlib,traceback,time
ROOT=Path(__file__).resolve().parent
os.environ['STPSF_PATH']=str(ROOT/'stpsf_reference/stpsf-data')
os.environ['MPLCONFIGDIR']=str(ROOT/'matplotlib_cache')
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from astropy.io import fits
from astropy.config.paths import set_temp_cache,set_temp_config
sys.path.insert(0,str(ROOT/'mask_patch/code'));sys.path.insert(0,str(ROOT/'repo/code'))
from strong_lensing_jwst_rate_mask import fit_jwst_native_rate_point
OUT=ROOT/'stpsf_trial';OUT.mkdir(exist_ok=True)
for folder in ['astropy_cache','astropy_config','wss_opds']:(OUT/folder).mkdir(exist_ok=True)
with set_temp_cache(str(OUT/'astropy_cache')),set_temp_config(str(OUT/'astropy_config')):
 import stpsf
 from astroquery.mast import Mast
 Mast.cache_location=OUT/'astropy_cache'
 from synphot import SourceSpectrum
 from synphot.models import BlackBodyNorm1D
 base=ROOT/'repo/artifacts/strong_lensing_geometry_20261021'
 cal=json.loads((base/'standard_psf/field_psf_calibration.json').read_text())
 report={'status':'DEVELOPMENT_PHYSICAL_PSF_CALIBRATION_TRIAL','stpsf_version':stpsf.__version__,
         'spectral_assumption':'4000K blackbody field star, prescribed rather than measured NIR spectrum',
         'native_errors_changed':False,'target_pixels_used_for_response':False,
         'oversampling':5,'fft_oversampling':2,'fov_pixels':25,'nlambda':9,'records':[]}
 try:
  for row in cal['field_fits'][:1]:
   p=ROOT/'repo/data/strong_lensing_geometry_20261021/native'/row['exposure']
   with fits.open(p) as h:
    head=dict(h[0].header);data=np.array(h['SCI'].data);error=np.array(h['ERR'].data);dq=np.array(h['DQ'].data)
   guess=np.array(row['centroid_xy_pixel'])
   nrc=stpsf.NIRCam();nrc.detector='NRCB3';nrc.filter='F115W';nrc.detector_position=tuple(guess)
   print('STPSF',stpsf.__version__,'position',nrc.detector_position,flush=True)
   nrc.load_wss_opd_by_date(head['DATE-OBS']+'T'+head['TIME-OBS'],choice='closest',plot=False,
                            output_path=str(OUT/'wss_opds'))
   spec=SourceSpectrum(BlackBodyNorm1D,temperature=4000)
   start=time.time()
   grid=nrc.psf_grid(num_psfs=1,all_detectors=False,single_psf_centered=False,
       fft_oversample=2,detector_oversample=5,fov_pixels=25,nlambda=9,source=spec,verbose=False)
   arr=np.array(grid.data[0]);origin=((arr.shape[1]-1)/2,(arr.shape[0]-1)/2)
   path=OUT/row['exposure'].replace('_cal.fits','_FIELD_stpsf_epsf.fits')
   hd=fits.PrimaryHDU(arr);hd.header['OVERSMP']=5;hd.header['ORIGINX']=origin[0];hd.header['ORIGINY']=origin[1]
   hd.header['PIXINT']=True;hd.header['SOURCE']='STPSF dated WSS,4000K';hd.writeto(path,overwrite=True)
   fit=fit_jwst_native_rate_point(data,error,dq,guess,arr,5,origin,primary_header=head,background_degree=1)
   result={k:v for k,v in fit.items() if not k.startswith('pixel_') and k!='standardized_residual'}
   result.update({'exposure':row['exposure'],'chi2_dof':fit['chi2']/fit['dof'],'elapsed_seconds':time.time()-start,
                  'psf_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'options':dict(nrc.options),'metadata':str(grid.meta)})
   report['records'].append(result);print(json.dumps(result),flush=True)
 except Exception as exc:
  report['failure']={'type':type(exc).__name__,'message':str(exc),'traceback':traceback.format_exc()}
  print(report['failure']['traceback'],flush=True)
 (OUT/'field_trial.json').write_text(json.dumps(report,indent=2,default=str)+'\n')

