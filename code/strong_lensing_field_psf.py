"""Field-only empirical effective PSFs; no target selection or data access.

The ePSF already includes the native pixel response. Its samples must not be
integrated over a pixel a second time. All statistical checks retain the input
ERR values; finite calibration-star uncertainty remains explicit and conditional.
"""
from __future__ import annotations

from itertools import combinations
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.distance import cdist
from scipy.stats import chi2

from strong_lensing_pixel_geometry import native_pixel_mask, tangent_plane_arcsec


DEFAULT_ALGORITHM = {
    'version': 1,
    'target_center_icrs_deg': [69.562083333, -12.287333333],
    'exclude_target_radius_arcsec': 10.,
    'stamp_size_native_pixel': 25,
    'fit_stamp_radius_native_pixel': 8,
    'oversampling': 4,
    'dao_fwhm_native_pixel': 1.5,
    'dao_detection_sigma': 8.,
    'dao_sharpness_range': [.1, 1.],
    'dao_roundness_range': [-.5, .5],
    'minimum_peak_snr': 20.,
    'minimum_aperture_snr': 100.,
    'aperture_radius_native_pixel': 4.,
    'minimum_valid_stamp_fraction': .95,
    'central_valid_square_size': 5,
    'minimum_isolation_native_pixel': 31.,
    'minimum_detector_edge_native_pixel': 13,
    'sky_match_radius_arcsec': .12,
    'heldout_every_nth_physical_star': 4,
    'minimum_training_physical_stars': 3,
    'minimum_heldout_physical_stars': 1,
    'minimum_training_stamps_per_exposure': 3,
    'allow_common_pool_for_this_visit': False,
    'maximum_common_pool_dither_native_pixel': 64.,
    'background_annulus_native_pixel': [10., 12.],
    'epsf_smoothing_kernel': 'quartic',
    'epsf_maxiters': 12,
    'epsf_refinement_iters': 5,
    'epsf_fit_shape': [17, 17],
    'epsf_sigma_clip': 3.,
    'epsf_sigma_clip_maxiters': 10,
    'epsf_center_accuracy_native_pixel': .001,
    'epsf_converged_fraction': .95,
    'heldout_centroid_bound_native_pixel': 2.,
    'heldout_minimum_absolute_fit_p_value': .001,
    'negative_wings_policy': 'preserve learned finite values; no clipping to zero',
    'normalization': 'sum(epsf_grid)/oversampling**2 = 1',
    'covariance': 'Known ERR inverse information; no residual variance scaling, floor or inflation',
    'fallback': 'No quasar or target pixel may become a calibration source',
}


def field_mask_from_gwcs(shape, wcs, center, radius_arcsec=10., row_chunk=128):
    """Compute exclusion using WCS coordinates only, before examining SCI/ERR/DQ."""
    mask = np.zeros(shape, dtype=bool)
    for start in range(0, shape[0], row_chunk):
        yy, xx = np.mgrid[start:min(start+row_chunk, shape[0]), :shape[1]]
        sky = wcs(xx, yy)
        tangent = tangent_plane_arcsec(np.stack(sky[:2], axis=-1), center)
        separation = np.linalg.norm(tangent, axis=-1)
        mask[start:start+len(yy)] = np.isfinite(separation) & (separation > radius_arcsec)
    return mask


def guarded_field_arrays(data, error, dq, field_mask):
    """Copy only allowed field entries; excluded SCI values are never inspected."""
    field = np.asarray(field_mask, bool)
    if data.shape != field.shape or error.shape != field.shape or dq.shape != field.shape:
        raise ValueError('Data arrays and geometry-only field mask must have equal shapes')
    safe_data = np.zeros(field.shape, dtype=float)
    safe_error = np.ones(field.shape, dtype=float)
    safe_dq = np.ones(field.shape, dtype=np.uint32)
    safe_data[field] = data[field]
    safe_error[field] = error[field]
    safe_dq[field] = dq[field]
    return safe_data, safe_error, safe_dq


def _stamp_bounds(shape, center, size):
    radius = size//2
    ix, iy = np.rint(center).astype(int)
    bounds = (ix-radius, iy-radius, ix+radius+1, iy+radius+1)
    if min(bounds[:2]) < 0 or bounds[2] > shape[1] or bounds[3] > shape[0]:
        raise ValueError('A calibration cutout would extend outside the detector')
    return bounds


def select_field_stars(data, error, dq, field_mask, wcs, algorithm=None):
    """Fixed DAO morphology/SNR/isolation selection on sanitized field arrays."""
    from photutils.detection import DAOStarFinder
    cfg = DEFAULT_ALGORITHM if algorithm is None else algorithm
    field = np.asarray(field_mask, bool)
    # This assertion is a data-flow guard, not an inspection of excluded input:
    # callers must first pass through guarded_field_arrays.
    if np.any(data[~field] != 0) or np.any(error[~field] != 1) or np.any(dq[~field] != 1):
        raise ValueError('Unsanitized target pixels must not enter field selection')
    valid = field & native_pixel_mask(data, error, dq)
    if not np.any(valid):
        return [], {'status':'NO_VALID_FIELD_PIXELS', 'candidate_count':0}
    background = float(np.median(data[valid]))
    score = np.zeros(data.shape, float)
    score[valid] = (data[valid]-background)/error[valid]
    finder = DAOStarFinder(cfg['dao_detection_sigma'], cfg['dao_fwhm_native_pixel'],
                           sharpness_range=tuple(cfg['dao_sharpness_range']),
                           roundness_range=tuple(cfg['dao_roundness_range']),
                           exclude_border=True)
    table = finder(score, mask=~valid)
    if table is None:
        return [], {'status':'NO_FIELD_POINT_CANDIDATES', 'candidate_count':0}
    centers = np.c_[table['x_centroid'], table['y_centroid']]
    distances = cdist(centers, centers)
    np.fill_diagonal(distances, np.inf)
    accepted, rejected = [], []
    for index, center in enumerate(centers):
        reason = None
        if distances[index].min() < cfg['minimum_isolation_native_pixel']:
            reason = 'nonisolated'
        edge = cfg['minimum_detector_edge_native_pixel']
        if not (edge <= center[0] < data.shape[1]-edge and edge <= center[1] < data.shape[0]-edge):
            reason = 'detector_edge'
        try:
            x0,y0,x1,y1 = _stamp_bounds(data.shape, center, cfg['stamp_size_native_pixel'])
        except ValueError:
            reason = 'stamp_edge'
            x0=y0=x1=y1=0
        if reason is None and not np.all(field[y0:y1,x0:x1]):
            reason = 'target_exclusion_intersection'
        if reason is None:
            use = valid[y0:y1,x0:x1]
            r = cfg['stamp_size_native_pixel']//2
            c = cfg['central_valid_square_size']//2
            if use.mean() < cfg['minimum_valid_stamp_fraction'] or not use[r-c:r+c+1,r-c:r+c+1].all():
                reason = 'DQ_or_ERR_mask'
        if reason is None:
            yy,xx = np.mgrid[y0:y1,x0:x1]
            aperture = (xx-center[0])**2+(yy-center[1])**2 <= cfg['aperture_radius_native_pixel']**2
            aperture &= use
            peak_snr = float(np.max(score[y0:y1,x0:x1][aperture]))
            flux = float(np.sum((data[y0:y1,x0:x1]-background)[aperture]))
            snr = flux/np.sqrt(np.sum(error[y0:y1,x0:x1][aperture]**2))
            if peak_snr < cfg['minimum_peak_snr'] or snr < cfg['minimum_aperture_snr']:
                reason = 'low_SNR'
        if reason is not None:
            rejected.append({'candidate_index':index,'xy':center.tolist(),'reason':reason})
            continue
        sky = wcs(float(center[0]), float(center[1]))
        accepted.append({'candidate_index':index,'xy':center.tolist(),
                         'ra_dec_deg':[float(sky[0]),float(sky[1])],
                         'peak_snr':peak_snr,'aperture_snr':float(snr),
                         'valid_stamp_fraction':float(use.mean()),'stamp_bounds_xyxy':[x0,y0,x1,y1]})
    return accepted, {'status':'FIELD_SELECTION_EXECUTED','candidate_count':len(table),
                      'accepted_count':len(accepted),'rejected':rejected,
                      'field_median_background_SCI':background}


def assign_physical_stars(exposure_catalogues, algorithm=None):
    """Sky crossmatch and deterministic object-level training/holdout split."""
    cfg = DEFAULT_ALGORITHM if algorithm is None else algorithm
    groups = []
    for exposure, catalogue in sorted(exposure_catalogues.items()):
        for row in sorted(catalogue,key=lambda r:tuple(r['ra_dec_deg'])):
            candidates = []
            for i, group in enumerate(groups):
                if any(item['exposure']==exposure for item in group['detections']):
                    continue
                delta=tangent_plane_arcsec(np.asarray(row['ra_dec_deg']),group['ra_dec_deg'])
                if np.linalg.norm(delta) <= cfg['sky_match_radius_arcsec']:
                    candidates.append((np.linalg.norm(delta),i))
            detection={**row,'exposure':exposure}
            if candidates:
                groups[min(candidates)[1]]['detections'].append(detection)
            else:
                groups.append({'ra_dec_deg':row['ra_dec_deg'],'detections':[detection]})
    groups.sort(key=lambda group:tuple(group['ra_dec_deg']))
    for index, group in enumerate(groups):
        group['physical_star_id']=f'field{index+1:04d}'
        group['split']=('heldout' if index%cfg['heldout_every_nth_physical_star']==0 else 'training')
    return groups


def common_pool_dither_diagnostic(training_groups, exposure_names, maximum_dither):
    """Missing matched exposure pairs never establish a small dither."""
    pairs={tuple(pair):[] for pair in combinations(sorted(exposure_names),2)}
    for group in training_groups:
        detected={row['exposure']:np.asarray(row['xy'],float) for row in group['detections']}
        for pair in pairs:
            if all(name in detected for name in pair):
                pairs[pair].append(float(np.linalg.norm(detected[pair[0]]-detected[pair[1]])))
    report={'|'.join(pair):{'matched_stars':len(values),'maximum_dither_native_pixel':max(values) if values else None}
            for pair,values in pairs.items()}
    complete=bool(pairs) and all(pairs.values())
    maximum=max((v for values in pairs.values() for v in values),default=None)
    return {'all_exposure_pairs_have_matched_training_star':bool(complete),
            'maximum_dither_native_pixel':maximum,'exposure_pairs':report,
            'small_dither_condition_pass':bool(complete and maximum<=maximum_dither)}


def epsf_model(grid, oversampling=4):
    """Official Photutils cubic-spline evaluator; ePSF is already pixel integrated."""
    from photutils.psf import ImagePSF
    array=np.asarray(grid,float)
    if array.ndim!=2 or min(array.shape)<4 or not np.isfinite(array).all():
        raise ValueError('Finite two-dimensional empirical PSF grid required')
    if oversampling < 1 or int(oversampling)!=oversampling:
        raise ValueError('Positive integer oversampling required')
    return ImagePSF(array,oversampling=int(oversampling),fill_value=0.)


def build_empirical_psf(stamps, algorithm=None):
    """Build from pre-extracted field-only stamps, with fixed background subtraction."""
    from astropy.stats import SigmaClip
    from photutils.psf import EPSFBuilder, EPSFStar, EPSFStars
    cfg=DEFAULT_ALGORITHM if algorithm is None else algorithm
    stars=[]
    for item in stamps:
        data,error,dq=np.asarray(item['data'],float),np.asarray(item['error'],float),np.asarray(item['dq'])
        center=np.asarray(item['center'],float)
        yy,xx=np.indices(data.shape); valid=native_pixel_mask(data,error,dq)
        rr=np.hypot(xx-center[0],yy-center[1]);r0,r1=cfg['background_annulus_native_pixel']
        outer=valid & (rr>=r0) & (rr<=r1)
        if outer.sum()<20:
            raise ValueError('Insufficient valid field background-annulus pixels')
        basis=np.c_[np.ones(outer.sum()),(xx-center[0])[outer],(yy-center[1])[outer]]
        plane=np.linalg.lstsq(basis/error[outer,None],data[outer]/error[outer],rcond=None)[0]
        background=plane[0]+plane[1]*(xx-center[0])+plane[2]*(yy-center[1])
        weights=np.zeros(data.shape);weights[valid]=1/error[valid]
        clean=np.zeros(data.shape);clean[valid]=(data-background)[valid]
        stars.append(EPSFStar(clean,weights=weights,cutout_center=center,id_label=item['id']))
    if len(stars)<cfg['minimum_training_stamps_per_exposure']:
        raise ValueError('Insufficient training field-star stamps')
    builder=EPSFBuilder(oversampling=cfg['oversampling'],smoothing_kernel=cfg['epsf_smoothing_kernel'],
        sigma_clip=SigmaClip(sigma=cfg['epsf_sigma_clip'],maxiters=cfg['epsf_sigma_clip_maxiters']),
        maxiters=cfg['epsf_maxiters'],refinement_iters=cfg['epsf_refinement_iters'],
        center_accuracy=cfg['epsf_center_accuracy_native_pixel'],converged_fraction=cfg['epsf_converged_fraction'],
        fit_shape=tuple(cfg['epsf_fit_shape']),constrain_fluxes=True,progress_bar=False)
    result=builder(EPSFStars(stars))
    retained=len(stars)-int(result.n_excluded_stars)
    if retained<cfg['minimum_training_stamps_per_exposure']:
        raise ValueError('EPSFBuilder excluded too many training stars: fewer than the frozen minimum remain')
    grid=np.asarray(result.epsf.data,float)
    if not np.isfinite(grid).all() or grid.sum()<=0:
        raise ValueError('Empirical PSF building produced an invalid grid')
    grid=grid*(cfg['oversampling']**2/grid.sum())
    return grid, {'iterations':int(result.iterations),'converged':bool(result.converged),
        'final_center_accuracy_native_pixel':float(result.final_center_accuracy),
        'final_converged_fraction':float(result.final_converged_fraction),
        'n_excluded_stars':int(result.n_excluded_stars),
        'n_retained_training_stars':retained,
        'excluded_star_indices':[int(i) for i in result.excluded_star_indices],
        'grid_shape':list(grid.shape),'negative_grid_entries':int(np.sum(grid<0)),
        'normalization_sum_over_oversampling_squared':float(grid.sum()/cfg['oversampling']**2),
        'background_calibration_uncertainty_marginalized':False}


def fit_heldout_field_star(data,error,dq,guess,grid,algorithm=None):
    """A fixed six-parameter fit and absolute known-ERR adequacy test."""
    cfg=DEFAULT_ALGORITHM if algorithm is None else algorithm
    size=2*cfg['fit_stamp_radius_native_pixel']+1
    x0,y0,x1,y1=_stamp_bounds(data.shape,guess,size)
    yy,xx=np.mgrid[y0:y1,x0:x1];valid=native_pixel_mask(data,error,dq)[y0:y1,x0:x1]
    x,y=xx[valid].astype(float),yy[valid].astype(float)
    values,sigma=np.asarray(data[y0:y1,x0:x1])[valid],np.asarray(error[y0:y1,x0:x1])[valid]
    if len(values)<50:raise ValueError('Insufficient valid held-out field pixels')
    unit=float(np.median(sigma));values=values/unit;sigma=sigma/unit
    model=epsf_model(grid,cfg['oversampling']);cx,cy=guess
    bg=float(np.median(values));flux=max(float(np.maximum(values-bg,0).sum()),1.)
    def evaluate(p):
        return model.evaluate(x,y,np.exp(p[2]),p[0],p[1])+p[3]+p[4]*(x-cx)+p[5]*(y-cy)
    bound=cfg['heldout_centroid_bound_native_pixel']
    lo=[cx-bound,cy-bound,np.log(1e-10),-np.inf,-np.inf,-np.inf]
    hi=[cx+bound,cy+bound,np.log(max(flux*1e4,1e8)),np.inf,np.inf,np.inf]
    fit=least_squares(lambda p:(evaluate(p)-values)/sigma,[cx,cy,np.log(flux),bg,0.,0.],
        bounds=(lo,hi),jac='3-point',x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-8,max_nfev=300)
    rank=int(np.linalg.matrix_rank(fit.jac));dof=len(values)-rank;stat=float(fit.fun@fit.fun)
    jc,jn=fit.jac[:,:2],fit.jac[:,2:]
    profiled=jc-jn@np.linalg.lstsq(jn,jc,rcond=1e-12)[0];info=profiled.T@profiled
    cov=np.linalg.inv(info)
    good=bool(fit.success and not np.any(fit.active_mask) and rank==6 and chi2.sf(stat,dof)>=cfg['heldout_minimum_absolute_fit_p_value'])
    return {'centroid_xy_pixel':fit.x[:2].tolist(),'centroid_covariance_pixel2':cov.tolist(),
        'chi2':stat,'dof':int(dof),'p_value':float(chi2.sf(stat,dof)),
        'optimizer_success':bool(fit.success),'active_bound_parameter_indices':np.flatnonzero(fit.active_mask).tolist(),
        'usable_pixels':len(values),'rank':rank,'absolute_adequacy_pass':good,
        'covariance_rescaled':False,'conditional_on_fixed_empirical_PSF':True}
