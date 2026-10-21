"""Independent field-star identity and spatial standard-ePSF applications."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import label


def stellar_identity(row, minimum_probability=.9, minimum_significance=5.,maximum_ruwe=1.4):
    """Catalogue classification/astrometric evidence, without any PSF-fit result."""
    def number(key):
        try:return float(row[key])
        except (KeyError,ValueError,TypeError):return np.nan
    probability=number('classprob_dsc_combmod_star')
    parallax,perr=number('parallax'),number('parallax_error')
    pm=np.array([number('pmra'),number('pmdec')]);errors=np.array([number('pmra_error'),number('pmdec_error')])
    correlation=number('pmra_pmdec_corr');significance=None
    if np.isfinite(pm).all() and np.isfinite(errors).all() and min(errors)>0 and np.isfinite(correlation) and abs(correlation)<1:
        covariance=np.outer(errors,errors)*np.array([[1.,correlation],[correlation,1.]])
        significance=float(np.sqrt(pm@np.linalg.solve(covariance,pm)))
    parallax_significance=float(parallax/perr) if np.isfinite(parallax) and perr>0 else None
    independent_evidence=bool((np.isfinite(probability) and probability>=minimum_probability)
        or (significance is not None and significance>=minimum_significance)
        or (parallax_significance is not None and parallax_significance>=minimum_significance))
    quality=(str(row.get('duplicated_source','')).lower()=='false'
             and np.isfinite(number('ruwe')) and number('ruwe')<maximum_ruwe)
    return {'accepted':bool(independent_evidence and quality),
            'star_probability':float(probability) if np.isfinite(probability) else None,
            'proper_motion_significance':significance,'parallax_significance':parallax_significance,
            'catalogue_quality_accepted':bool(quality)}


def propagate_gaia_field_position(row,epoch_jyear):
    """First-order ICRS proper-motion transport; RA proper motion is alpha-star."""
    ra,dec=float(row['ra']),float(row['dec']);dt=float(epoch_jyear)-float(row['ref_epoch'])
    pmra,pmdec=float(row['pmra']),float(row['pmdec'])
    if not np.isfinite([ra,dec,pmra,pmdec,dt]).all() or abs(dec)>=89:
        raise ValueError('Finite nonpolar Gaia field coordinates and proper motions required')
    return np.array([ra+pmra*dt/(3600000*np.cos(np.deg2rad(dec))),dec+pmdec*dt/3600000])


def associate_saturated_field_core(dq,field_mask,nominal_xy,search_radius=60.,minimum_pixels=2):
    """Identify the nearest saturated component from DQ only, outside the guard."""
    if not np.isfinite(nominal_xy).all():raise ValueError('Finite projected catalogue position required')
    x,y=map(float,nominal_xy);ix,iy=np.rint([x,y]).astype(int);r=int(np.ceil(search_radius))
    x0,x1=max(0,ix-r),min(dq.shape[1],ix+r+1);y0,y1=max(0,iy-r),min(dq.shape[0],iy+r+1)
    yy,xx=np.mgrid[y0:y1,x0:x1]
    use=(field_mask[y0:y1,x0:x1] & (((xx-x)**2+(yy-y)**2)<=search_radius**2))
    saturated=((np.asarray(dq[y0:y1,x0:x1],np.uint64)&2)!=0)&use
    components,count=label(saturated,structure=np.ones((3,3),int));options=[]
    for index in range(1,count+1):
        py,px=np.nonzero(components==index)
        if len(px)<minimum_pixels:continue
        center=np.array([px.mean()+x0,py.mean()+y0]);distance=float(np.linalg.norm(center-[x,y]))
        options.append({'center_xy_pixel':center.tolist(),'pixels':len(px),'distance_to_catalogue_pixel':distance})
    options.sort(key=lambda row:row['distance_to_catalogue_pixel'])
    if not options:raise ValueError('No eligible saturated field component near the selected catalogue star')
    return {**options[0],'candidate_components':options,'association_uses_SCI':False}


def sample_spatial_standard_psf(model,source_xy,oversampling=4,shape=(101,101)):
    """Sample the official gridded ePSF at a detector position, without integration."""
    cx,cy=map(float,source_xy)
    xx=(np.arange(shape[1])-(shape[1]-1)/2)/oversampling+cx
    yy=(np.arange(shape[0])-(shape[0]-1)/2)/oversampling+cy
    x,y=np.meshgrid(xx,yy)
    grid=np.asarray(model.evaluate(x,y,1.,cx,cy),float)
    if grid.shape!=tuple(shape) or not np.isfinite(grid).all() or grid.sum()<=0:
        raise ValueError('Invalid spatial standard PSF response')
    return grid*(oversampling**2/grid.sum())


def registered_reference_guesses(wcs,reference_sky,nominal_field_xy,fitted_field_xy):
    """Field-only translation changes association guesses, never a lens geometry."""
    raw=np.array([wcs.invert(float(a),float(d),with_bounding_box=False)[:2] for a,d in reference_sky],float)
    delta=np.asarray(fitted_field_xy,float)-np.asarray(nominal_field_xy,float)
    if raw.shape!=(4,2) or delta.shape!=(2,) or not np.isfinite(raw).all() or not np.isfinite(delta).all():
        raise ValueError('Finite four-reference positions and field-only translation required')
    return raw+delta,delta


def apply_minimum_pixel_gate(fit,minimum_valid_pixels):
    """Consume the frozen pixel-count requirement without overriding fit failures."""
    minimum=int(minimum_valid_pixels)
    if minimum<1 or minimum!=minimum_valid_pixels:raise ValueError('Positive integer pixel minimum required')
    result=dict(fit)
    result['minimum_usable_pixels_required']=minimum
    result['minimum_pixel_count_pass']=bool(int(result['usable_pixels'])>=minimum)
    result['absolute_adequacy_pass']=bool(result['absolute_adequacy_pass'] and result['minimum_pixel_count_pass'])
    return result
