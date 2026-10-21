"""Exploratory mean-response diagnostics on the already consumed FIELD star.

No quasar array is read. Fits below are calibration development, not held-out
validation. Jitter kernels change the mean image; ERR remains unmodified.
Whole-ePSF dilation is an empirical diagnostic because pixel integration is
already included; it must not be silently identified with an optical SED model.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares
from scipy.linalg import expm
from scipy.stats import chi2
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_pixel_geometry import native_pixel_mask


def mean_response(grid,kind,parameters):
    n=grid.shape[0];o=(n-1)/2;yy,xx=np.indices(grid.shape);output=grid
    if kind in ['isotropic_dilation','affine_dilation']:
        scale=parameters[0]
        e1,e2=parameters[1:] if kind=='affine_dilation' else (0.,0.)
        A=expm(np.array([[scale+e1,e2],[e2,scale-e1]]));inverse=np.linalg.inv(A)
        coordinates=inverse@np.vstack([(xx-o).ravel(),(yy-o).ravel()])
        spline=RectBivariateSpline(np.arange(n),np.arange(n),grid.T,kx=3,ky=3)
        output=spline.ev(coordinates[0]+o,coordinates[1]+o).reshape(grid.shape)/np.linalg.det(A)
    if kind in ['isotropic_jitter','elliptic_jitter','isotropic_halo','elliptic_halo']:
        halo='halo' in kind;shape=parameters[1:] if halo else parameters
        sx=np.exp(shape[0]);sy=sx if kind in ['isotropic_jitter','isotropic_halo'] else np.exp(shape[1])
        correlation=0. if kind in ['isotropic_jitter','isotropic_halo'] else shape[2]
        padded=np.pad(output,32);freq=np.fft.fftfreq(padded.shape[0],d=.25);ky,kx=np.meshgrid(freq,freq,indexing='ij')
        kernel=np.exp(-2*np.pi**2*(sx*sx*kx*kx+sy*sy*ky*ky+2*correlation*sx*sy*kx*ky))
        blurred=np.fft.ifft2(np.fft.fft2(padded)*kernel).real[32:-32,32:-32]
        fraction=1/(1+np.exp(-parameters[0])) if halo else 1.
        output=(1-fraction)*output+fraction*blurred
    return output*(16/output.sum())


def fit_stamp(data,error,dq,guess,grid,kind):
    yy,xx=np.indices(data.shape);valid=native_pixel_mask(data,error,dq)
    x,y=xx[valid].astype(float),yy[valid].astype(float);raw=data[valid];sigma=error[valid]
    unit=np.median(sigma);d=raw/unit;err=sigma/unit;cx,cy=guess
    bg=np.median(d);flux=max(np.maximum(d-bg,0).sum(),1.)
    shape0={'standard':[],'isotropic_jitter':[np.log(.05)],'elliptic_jitter':[np.log(.05),np.log(.05),0.],
            'isotropic_dilation':[0.],'affine_dilation':[0.,0.,0.],
            'isotropic_halo':[-2.,np.log(1.5)],'elliptic_halo':[-2.,np.log(1.5),np.log(1.5),0.]}[kind]
    slo={'standard':[],'isotropic_jitter':[np.log(1e-4)],'elliptic_jitter':[np.log(1e-4),np.log(1e-4),-.95],
         'isotropic_dilation':[np.log(.85)],'affine_dilation':[np.log(.85),-.1,-.1],
         'isotropic_halo':[-10.,np.log(.1)],'elliptic_halo':[-10.,np.log(.1),np.log(.1),-.95]}[kind]
    shi={'standard':[],'isotropic_jitter':[np.log(.8)],'elliptic_jitter':[np.log(.8),np.log(.8),.95],
         'isotropic_dilation':[np.log(1.15)],'affine_dilation':[np.log(1.15),.1,.1],
         'isotropic_halo':[0.,np.log(5.)],'elliptic_halo':[0.,np.log(5.),np.log(5.),.95]}[kind]
    n=grid.shape[0];o=(n-1)/2
    def values(p):
        response=mean_response(grid,kind,p[6:])
        spline=RectBivariateSpline(np.arange(n),np.arange(n),response.T,kx=3,ky=3)
        return np.exp(p[2])*spline.ev((x-p[0])*4+o,(y-p[1])*4+o)+p[3]+p[4]*(x-cx)+p[5]*(y-cy)
    lo=[cx-2,cy-2,np.log(1e-10),-np.inf,-np.inf,-np.inf]+slo
    hi=[cx+2,cy+2,np.log(max(flux*1e4,1e8)),np.inf,np.inf,np.inf]+shi
    fit=least_squares(lambda p:(values(p)-d)/err,[cx,cy,np.log(flux),bg,0.,0.]+shape0,
        bounds=(lo,hi),jac='3-point',x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-7,max_nfev=500)
    statistic=float(fit.fun@fit.fun);rank=int(np.linalg.matrix_rank(fit.jac));dof=len(d)-rank
    model=values(fit.x)*unit
    return {'chi2':statistic,'dof':dof,'p_value':float(chi2.sf(statistic,dof)),
        'centroid_local_pixel':fit.x[:2].tolist(),'shape_parameters':fit.x[6:].tolist(),
        'all_parameters':fit.x.tolist(),'rank':rank,'optimizer_success':bool(fit.success),'nfev':fit.nfev,
        'active_bound_indices':np.flatnonzero(fit.active_mask).tolist(),'pixel_x':x.tolist(),'pixel_y':y.tolist(),
        'pixel_data':raw.tolist(),'pixel_error':sigma.tolist(),'pixel_model':model.tolist(),
        'pixel_residual':(raw-model).tolist(),'standardized_residual':(-fit.fun).tolist(),
        'SCI_DQ_counts':{str(1<<bit):int(np.count_nonzero(dq&(1<<bit))) for bit in range(32) if np.any(dq&(1<<bit))},
        'ERR_modified':False,'calibration_development_only':True}


def main():
    report=json.loads((ROOT/'artifacts/strong_lensing_geometry_20261021/standard_psf/field_psf_calibration.json').read_text())
    output={'scope':'Exploratory mean-response calibration development on the chosen Gaia FIELD star only; no quasar pixel use or held-out validation claim. Known ERR unmodified.',
            'dilation_scope':'Remaps the already pixel-integrated response empirically; not an exact optical wavelength/SED derivation.',
            'jitter_scope':'Positive Gaussian convolution of ePSF is a physically valid additional image-motion/diffusion mean response; kernel covariance changes the mean, not measurement errors.','fits':{}}
    for item in report['field_fits']:
        name=item['exposure'];cx,cy=item['centroid_xy_pixel']
        ix,iy=np.rint(item['saturated_component_association']['center_xy_pixel']).astype(int);r=8
        # This source was independently identified and its stamp checked outside
        # the protected target circle before the external-library application.
        with fits.open(ROOT/'data/strong_lensing_geometry_20261021/native'/name,memmap=True,do_not_scale_image_data=True) as hdul:
            def local_plane(label):
                hdu=hdul[label]
                raw=hdu.data[iy-r:iy+r+1,ix-r:ix+r+1].astype(float)
                return raw*hdu.header.get('BSCALE',1.)+hdu.header.get('BZERO',0.)
            data=local_plane('SCI');error=local_plane('ERR');dq=local_plane('DQ').astype(np.uint32)
        grid=fits.getdata(ROOT/'artifacts/strong_lensing_geometry_20261021/standard_psf'/name.replace('_cal.fits','_FIELD_standard_epsf.fits'))
        output['fits'][name]={}
        for kind in ['standard','isotropic_halo','elliptic_halo']:
            result=fit_stamp(data,error,dq,np.array([cx-(ix-r),cy-(iy-r)]),grid,kind)
            output['fits'][name][kind]=result
            print(name,kind,json.dumps({k:result[k] for k in ['chi2','dof','shape_parameters','centroid_local_pixel','active_bound_indices']}),flush=True)
            (Path(__file__).with_name('field_PSF_halo_response_diagnostic.json')).write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
