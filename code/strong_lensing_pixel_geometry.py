"""Frozen native-pixel point geometry, with measurement adequacy kept explicit.

The Moffat family is an instrumental approximation, not a lens equation. Its
known-ERR residual test can fail. Formal covariance is never multiplied by the
observed residual variance, and no fitted astrometric floor is introduced.
"""
from __future__ import annotations

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from scipy.interpolate import RectBivariateSpline
from scipy.stats import chi2

DO_NOT_USE = 1
SATURATED = 2


def native_pixel_mask(data, error, dq):
    data, error, dq = np.asarray(data), np.asarray(error), np.asarray(dq)
    if data.shape != error.shape or data.shape != dq.shape:
        raise ValueError('SCI, ERR and DQ shapes must agree')
    return (np.isfinite(data) & np.isfinite(error) & (error > 0)
            & ((dq.astype(np.uint64) & (DO_NOT_USE | SATURATED)) == 0))


def integrated_moffat(x, y, center, flux, alpha_x, alpha_y, angle, beta,
                      quadrature_order=7):
    """Unit-pixel area integral of a normalized elliptical Moffat intensity."""
    if flux < 0 or min(alpha_x, alpha_y) <= 0 or beta <= 1:
        raise ValueError('Moffat requires nonnegative flux, alpha>0 and beta>1')
    z, w = leggauss(quadrature_order)
    ox, oy = np.meshgrid(z * .5, z * .5)
    ww = np.outer(w, w).ravel() * .25
    dx = np.asarray(x)[..., None] - center[0] + ox.ravel()
    dy = np.asarray(y)[..., None] - center[1] + oy.ravel()
    co, si = np.cos(angle), np.sin(angle)
    r2 = ((co * dx + si * dy) / alpha_x) ** 2 + ((-si * dx + co * dy) / alpha_y) ** 2
    return flux * (beta - 1) / (np.pi * alpha_x * alpha_y) * ((1 + r2) ** -beta @ ww)


def _moffat_plane(parameters, x, y, origin, quadrature_order):
    cx, cy, logflux, logax, logay, angle, logbm1, bg, gx, gy = parameters
    return (integrated_moffat(x, y, (cx, cy), np.exp(logflux),
                            np.exp(logax), np.exp(logay), angle,
                            1 + np.exp(logbm1), quadrature_order)
            + bg + gx * (x - origin[0]) + gy * (y - origin[1]))


def fit_moffat_point(data, error, dq, guess, stamp_radius=8,
                     centroid_bound=5., quadrature_order=7):
    """Fit one fixed 17x17 native-pixel square using declared ERR, never rescaled.

    Every valid pixel in the square is used. No sigma clipping or residual-driven
    mask/PSF-family choice is applied. Three fixed angle starts address optimizer
    degeneracy within the single family. Other DQ flags are counted, not masked.
    """
    data, error, dq = np.asarray(data, float), np.asarray(error, float), np.asarray(dq)
    if data.ndim != 2 or not np.all(np.isfinite(guess)):
        raise ValueError('Finite two-dimensional image and centroid guess required')
    good = native_pixel_mask(data, error, dq)
    cx, cy = np.asarray(guess, float)
    ix, iy = int(np.rint(cx)), int(np.rint(cy))
    x0, x1 = max(0, ix-stamp_radius), min(data.shape[1], ix+stamp_radius+1)
    y0, y1 = max(0, iy-stamp_radius), min(data.shape[0], iy+stamp_radius+1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    valid = good[y0:y1, x0:x1]
    if valid.sum() < 50:
        raise ValueError('Fewer than50 usable native pixels in declared stamp')
    x, y = xx[valid].astype(float), yy[valid].astype(float)
    raw, sigma = data[y0:y1, x0:x1][valid], error[y0:y1, x0:x1][valid]
    # Parameter conditioning only; the known standardized residual is unchanged.
    unit = float(np.median(sigma))
    d, err = raw/unit, sigma/unit
    bg = float(np.median(d))
    flux0 = max(float(np.maximum(d-bg, 0).sum()), 1.)
    lo = [cx-centroid_bound, cy-centroid_bound, np.log(1e-8), np.log(.3), np.log(.3),
          -np.pi, np.log(.1), -np.inf, -np.inf, -np.inf]
    hi = [cx+centroid_bound, cy+centroid_bound, np.log(max(flux0*1e4, 1e8)),
          np.log(6.), np.log(6.), np.pi, np.log(11.), np.inf, np.inf, np.inf]
    fits = []
    for theta in [0., np.pi/4, -np.pi/4]:
        p0 = [cx, cy, np.log(flux0), np.log(1.25), np.log(1.05), theta,
              np.log(2.), bg, 0., 0.]
        result = least_squares(lambda p: (_moffat_plane(p,x,y,(cx,cy),quadrature_order)-d)/err,
                               p0, bounds=(lo,hi), jac='3-point', x_scale='jac',
                               ftol=1e-10, xtol=1e-10, gtol=1e-8, max_nfev=350)
        fits.append(result)
    fit = min(fits, key=lambda r: float(r.fun @ r.fun))
    j = fit.jac
    # Profile nuisance directions. An unidentified circular-PSF angle does not
    # justify discarding an otherwise identified two-dimensional centroid.
    jc, jn = j[:,:2], j[:,2:]
    profiled = jc - jn @ np.linalg.lstsq(jn,jc,rcond=1e-12)[0]
    information = profiled.T @ profiled
    if np.linalg.matrix_rank(information) != 2:
        raise ValueError('Centroid Fisher information is singular')
    covariance = np.linalg.inv(information)
    rank = int(np.linalg.matrix_rank(j))
    dof = int(len(d)-rank)
    statistic = float(fit.fun @ fit.fun)
    model = _moffat_plane(fit.x,x,y,(cx,cy),quadrature_order)*unit
    p = fit.x
    flags = dq[y0:y1,x0:x1].astype(np.uint64)
    counts = {str(1 << bit): int(np.count_nonzero(flags & (np.uint64(1) << np.uint64(bit))))
              for bit in range(32) if np.any(flags & (np.uint64(1) << np.uint64(bit)))}
    active = [int(i) for i in np.flatnonzero(fit.active_mask)]
    return {'centroid_xy_pixel':p[:2].tolist(),
            'centroid_covariance_pixel2':covariance.tolist(),
            'flux_in_SCI_pixel_units':float(np.exp(p[2])*unit),
            'alpha_x_pixel':float(np.exp(p[3])), 'alpha_y_pixel':float(np.exp(p[4])),
            'angle_rad':float(p[5]), 'beta':float(1+np.exp(p[6])),
            'background_plane_SCI':[float(v*unit) for v in p[7:]],
            'chi2':statistic, 'dof':dof, 'p_value':float(chi2.sf(statistic,dof)),
            'usable_pixels':int(len(d)), 'jacobian_rank':rank,
            'optimizer_success':bool(fit.success),'optimizer_message':str(fit.message),
            'active_bound_parameter_indices':active,
            'valid_measurement_fit':bool(fit.success and not active),
            'centroid_covariance_scale':'Known ERR inverse information; no chi2/dof scaling',
            'dq_flag_counts_in_stamp':counts,
            'stamp_bounds_xyxy':[x0,y0,x1,y1],
            'pixel_x':x.tolist(),'pixel_y':y.tolist(),'pixel_data':raw.tolist(),
            'pixel_error':sigma.tolist(),'pixel_model':model.tolist(),
            'pixel_residual_data_minus_model':(raw-model).tolist(),
            'standardized_residual':(-fit.fun).tolist(),
            'all_start_chi2':[float(r.fun @ r.fun) for r in fits]}


def empirical_psf_values(grid,oversampling,origin,x,y,center):
    """Sample a pixel-integrated ePSF, matching Photutils ImagePSF convention.

    The field calibration supplies the normalization and origin. This is already
    a pixel-integrated response; integrating it over another pixel would blur it
    twice. Negative estimated wings are retained rather than clipped.
    """
    grid=np.asarray(grid,float)
    if grid.ndim!=2 or min(grid.shape)<4 or not np.all(np.isfinite(grid)):
        raise ValueError('Finite two-dimensional empirical PSF grid required')
    if oversampling<=0:raise ValueError('Positive oversampling required')
    px=(np.asarray(x)-center[0])*oversampling+origin[0]
    py=(np.asarray(y)-center[1])*oversampling+origin[1]
    # ImagePSF's interpolator uses x/y order with a transposed input data array.
    spline=RectBivariateSpline(np.arange(grid.shape[1]),np.arange(grid.shape[0]),grid.T,kx=3,ky=3,s=0)
    values=spline.ev(px,py)
    outside=(px<0)|(py<0)|(px>grid.shape[1]-1)|(py>grid.shape[0]-1)
    return np.where(outside,0.,values)


def fit_empirical_point(data,error,dq,guess,psf_grid,oversampling,psf_origin,
                        stamp_radius=8,centroid_bound=5.,background_degree=1):
    """Known-ERR ePSF fit with a fixed affine or quadratic diffuse background.

    Degree1 is the original six-parameter API. A separately frozen degree2
    application adds x²,xy,y² smooth nuisance terms; no PSF or ERR is modified.
    """
    if background_degree not in [1,2]:raise ValueError('Background degree must be1 or2')
    data,error,dq=np.asarray(data,float),np.asarray(error,float),np.asarray(dq)
    good=native_pixel_mask(data,error,dq)
    cx,cy=np.asarray(guess,float)
    ix,iy=int(np.rint(cx)),int(np.rint(cy))
    x0,x1=max(0,ix-stamp_radius),min(data.shape[1],ix+stamp_radius+1)
    y0,y1=max(0,iy-stamp_radius),min(data.shape[0],iy+stamp_radius+1)
    yy,xx=np.mgrid[y0:y1,x0:x1];valid=good[y0:y1,x0:x1]
    if valid.sum()<50:raise ValueError('Fewer than50 usable native pixels in declared stamp')
    x,y=xx[valid].astype(float),yy[valid].astype(float)
    raw,sigma=data[y0:y1,x0:x1][valid],error[y0:y1,x0:x1][valid]
    unit=float(np.median(sigma));d,err=raw/unit,sigma/unit
    bg=float(np.median(d));flux0=max(float(np.maximum(d-bg,0).sum()),1.)
    # Cache exactly the same cubic interpolator for all optimizer evaluations.
    grid=np.asarray(psf_grid,float)
    if grid.ndim!=2 or min(grid.shape)<4 or not np.all(np.isfinite(grid)) or oversampling<=0:
        raise ValueError('Invalid field-calibrated empirical PSF')
    interp=RectBivariateSpline(np.arange(grid.shape[1]),np.arange(grid.shape[0]),grid.T,kx=3,ky=3,s=0)
    dx,dy=x-cx,y-cy
    background=np.column_stack([np.ones_like(x),dx,dy] +
                               ([dx*dx,dx*dy,dy*dy] if background_degree==2 else []))
    def predict(p):
        px=(x-p[0])*oversampling+psf_origin[0]
        py=(y-p[1])*oversampling+psf_origin[1]
        profile=interp.ev(px,py)
        profile=np.where((px>=0)&(py>=0)&(px<=grid.shape[1]-1)&(py<=grid.shape[0]-1),profile,0.)
        return np.exp(p[2])*profile+background@p[3:]
    nbackground=background.shape[1]
    lo=[cx-centroid_bound,cy-centroid_bound,np.log(1e-8)]+[-np.inf]*nbackground
    hi=[cx+centroid_bound,cy+centroid_bound,np.log(max(flux0*1e4,1e8))]+[np.inf]*nbackground
    result=least_squares(lambda p:(predict(p)-d)/err,[cx,cy,np.log(flux0),bg]+[0.]*(nbackground-1),
                         bounds=(lo,hi),jac='3-point',x_scale='jac',
                         ftol=1e-10,xtol=1e-10,gtol=1e-8,max_nfev=350)
    j=result.jac;jc,jn=j[:,:2],j[:,2:]
    profiled=jc-jn@np.linalg.lstsq(jn,jc,rcond=1e-12)[0]
    information=profiled.T@profiled
    if np.linalg.matrix_rank(information)!=2:raise ValueError('Centroid Fisher information is singular')
    covariance=np.linalg.inv(information);rank=int(np.linalg.matrix_rank(j));dof=len(d)-rank
    statistic=float(result.fun@result.fun);model=predict(result.x)*unit
    active=np.flatnonzero(result.active_mask).astype(int).tolist()
    flags=dq[y0:y1,x0:x1].astype(np.uint64)
    counts={str(1<<b):int(np.count_nonzero(flags&(np.uint64(1)<<np.uint64(b)))) for b in range(32)
            if np.any(flags&(np.uint64(1)<<np.uint64(b)))}
    return {'centroid_xy_pixel':result.x[:2].tolist(),'centroid_covariance_pixel2':covariance.tolist(),
            'flux_in_SCI_pixel_units':float(np.exp(result.x[2])*unit),
            'background_plane_SCI':(result.x[3:6]*unit).tolist(),
            'background_polynomial_SCI':(result.x[3:]*unit).tolist(),
            'background_polynomial_basis':['1','dx','dy']+(['dx2','dxdy','dy2'] if background_degree==2 else []),
            'background_degree':int(background_degree),
            'chi2':statistic,'dof':int(dof),'p_value':float(chi2.sf(statistic,dof)),
            'usable_pixels':int(len(d)),'jacobian_rank':rank,
            'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
            'active_bound_parameter_indices':active,'valid_measurement_fit':bool(result.success and not active),
            'centroid_covariance_scale':'Known ERR inverse information conditional on field-frozen ePSF; no chi2/dof scaling',
            'psf_oversampling':int(oversampling),'psf_origin_grid_pixel':list(psf_origin),
            'dq_flag_counts_in_stamp':counts,'stamp_bounds_xyxy':[x0,y0,x1,y1],
            'pixel_x':x.tolist(),'pixel_y':y.tolist(),'pixel_data':raw.tolist(),'pixel_error':sigma.tolist(),
            'pixel_model':model.tolist(),'pixel_residual_data_minus_model':(raw-model).tolist(),
            'standardized_residual':(-result.fun).tolist()}


def intrinsic_geometry(points):
    """Five coordinates invariant to translation/proper rotation, retaining scale."""
    p = np.asarray(points,float)
    if p.shape != (4,2) or not np.all(np.isfinite(p)):
        raise ValueError('Four finite ordered A/B/C/D coordinates required')
    b = p[1]-p[0]
    length = np.linalg.norm(b)
    if length <= 0:
        raise ValueError('A and B must be distinct')
    u = b/length
    v = np.array([-u[1],u[0]])
    c,d = p[2]-p[0],p[3]-p[0]
    return np.array([length,c@u,c@v,d@u,d@v])


def intrinsic_jacobian(points):
    p = np.asarray(points,float)
    intrinsic_geometry(p)
    b = p[1]-p[0]; length = np.linalg.norm(b); u = b/length
    turn = np.array([[0.,-1.],[1.,0.]])
    v = turn@u
    derivative_u = (np.eye(2)-np.outer(u,u))/length
    j = np.zeros((5,8))
    j[0,:2],j[0,2:4] = -u,u
    for index,vector in [(2,p[2]-p[0]),(3,p[3]-p[0])]:
        row = 1+2*(index-2)
        for delta,basis,db in [(0,u,vector@derivative_u),
                               (1,v,vector@turn@derivative_u)]:
            j[row+delta,:2] = -basis-db
            j[row+delta,2:4] = db
            j[row+delta,2*index:2*index+2] = basis
    return j


def intrinsic_with_covariance(points,covariance):
    cov = np.asarray(covariance,float)
    if cov.shape != (8,8) or not np.all(np.isfinite(cov)):
        raise ValueError('Finite full8x8 joint centroid covariance required')
    j = intrinsic_jacobian(points)
    return intrinsic_geometry(points), j@cov@j.T


def gaussian_score(residual,covariance):
    r,c = np.asarray(residual,float),np.asarray(covariance,float)
    if c.shape != (len(r),len(r)) or not np.allclose(c,c.T):
        raise ValueError('Covariance shape/symmetry mismatch')
    chol = np.linalg.cholesky(c)
    z = np.linalg.solve(chol,r)
    q = float(z@z)
    return {'residual':r.tolist(),'chi2':q,'dof':int(len(r)),
            'p_value':float(chi2.sf(q,len(r))),
            'gaussian_logdensity':float(-.5*(q+len(r)*np.log(2*np.pi)
                                            +2*np.log(np.diag(chol)).sum()))}


def combine_exposures(values,covariances):
    """Known-covariance inverse-information mean and unreweighted scatter test."""
    values,covs = np.asarray(values,float),np.asarray(covariances,float)
    if values.ndim != 2 or len(values) == 0 or covs.shape != (len(values),values.shape[1],values.shape[1]):
        raise ValueError('At least one exposure and matching full covariances required')
    precisions = np.array([np.linalg.inv(c) for c in covs])
    mean_cov = np.linalg.inv(precisions.sum(axis=0))
    mean = mean_cov @ np.einsum('nij,nj->i',precisions,values)
    scatter = float(sum((v-mean)@w@(v-mean) for v,w in zip(values,precisions)))
    dof = int(values.shape[1]*(len(values)-1))
    return {'mean':mean.tolist(),'covariance':mean_cov.tolist(),
            'exposure_consistency_chi2':scatter,'exposure_consistency_dof':dof,
            'exposure_consistency_p_value':float(chi2.sf(scatter,dof)) if dof else None,
            'covariance_rescaled':False}


def tangent_plane_arcsec(ra_dec_deg,origin_ra_dec_deg):
    """Exact gnomonic east/north coordinates about a fixed declared ICRS origin."""
    a = np.deg2rad(np.asarray(ra_dec_deg,float)); o = np.deg2rad(origin_ra_dec_deg)
    dra = (a[...,0]-o[0]+np.pi)%(2*np.pi)-np.pi
    denom = np.sin(o[1])*np.sin(a[...,1])+np.cos(o[1])*np.cos(a[...,1])*np.cos(dra)
    if np.any(denom <= 0):raise ValueError('Point outside tangent hemisphere')
    x = np.cos(a[...,1])*np.sin(dra)/denom
    y = (np.cos(o[1])*np.sin(a[...,1])-np.sin(o[1])*np.cos(a[...,1])*np.cos(dra))/denom
    return np.stack([x,y],axis=-1)*206264.80624709636


def gwcs_centroid(wcs,xy,covariance,origin_ra_dec_deg,step=.001):
    """Apply official instrument GWCS, retaining its nonlinear local Jacobian."""
    xy,cov = np.asarray(xy,float),np.asarray(covariance,float)
    def transform(p):
        sky = wcs(float(p[0]),float(p[1]))
        return tangent_plane_arcsec(np.asarray(sky[:2],float),origin_ra_dec_deg)
    jac = np.column_stack([(transform(xy+np.eye(2)[i]*step)-transform(xy-np.eye(2)[i]*step))/(2*step)
                           for i in range(2)])
    return transform(xy),jac@cov@jac.T,jac
