"""Synthetic calibration checks without opening a native observation."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import erf, j1

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_field_psf import (DEFAULT_ALGORITHM,field_mask_from_gwcs,
    guarded_field_arrays,select_field_stars,assign_physical_stars,epsf_model,
    common_pool_dither_diagnostic,build_empirical_psf,fit_heldout_field_star)


def simple_wcs(x,y):
    # A synthetic tangent map about the declared field center, 0.031arcsec/pixel.
    ra0,dec0=DEFAULT_ALGORITHM['target_center_icrs_deg']
    return (ra0+(np.asarray(x)-128)*.031/(3600*np.cos(np.deg2rad(dec0))),
            dec0+(np.asarray(y)-128)*.031/3600)


def gaussian_prf(x,y,cx,cy,sigma=.7):
    return (.5*(erf((x-cx+.5)/(2**.5*sigma))-erf((x-cx-.5)/(2**.5*sigma)))
            *.5*(erf((y-cy+.5)/(2**.5*sigma))-erf((y-cy-.5)/(2**.5*sigma))))


def gaussian_grid():
    v=(np.arange(101)-50)/4
    x,y=np.meshgrid(v,v)
    grid=gaussian_prf(x,y,0.,0.)
    return grid*16/grid.sum()


def test_geometry_mask_precedes_field_values_and_target_never_enters_copy():
    shape=(256,256)
    field=field_mask_from_gwcs(shape,simple_wcs,DEFAULT_ALGORITHM['target_center_icrs_deg'],1.)
    assert not field[128,128] and field[128,180]
    data=np.full(shape,3.);error=np.ones(shape);dq=np.zeros(shape,np.uint32)
    data[~field]=np.nan;error[~field]=-1;dq[~field]=np.uint32(2**31)
    safe=guarded_field_arrays(data,error,dq,field)
    assert np.all(safe[0][~field]==0) and np.all(safe[1][~field]==1) and np.all(safe[2][~field]==1)
    assert np.all(safe[0][field]==3)


def test_automatic_selection_is_invariant_to_arbitrary_protected_target_data():
    cfg=copy.deepcopy(DEFAULT_ALGORITHM);cfg['exclude_target_radius_arcsec']=1.
    yy,xx=np.indices((256,256));data=np.zeros(xx.shape)
    for x,y in [(35.2,35.7),(100.3,35.2),(180.1,40.4),(210.4,210.2)]:
        data+=10000*gaussian_prf(xx,yy,x,y)
    error=np.ones(data.shape);dq=np.zeros(data.shape,np.uint32)
    field=field_mask_from_gwcs(data.shape,simple_wcs,cfg['target_center_icrs_deg'],1.)
    selected=[]
    for hidden in [-1e30,1e30]:
        modified=data.copy();modified[~field]=hidden
        safe=guarded_field_arrays(modified,error,dq,field)
        rows,_=select_field_stars(*safe,field,simple_wcs,cfg)
        selected.append(rows)
    assert selected[0]==selected[1] and len(selected[0])==4
    for row in selected[0]:
        x0,y0,x1,y1=row['stamp_bounds_xyxy'];assert field[y0:y1,x0:x1].all()
    with pytest.raises(ValueError,match='Unsanitized'):
        select_field_stars(data,error,dq,field,simple_wcs,cfg)


def test_physical_object_split_never_leaks_repeat_exposures_into_both_sets():
    catalogues={}
    for exposure,shift in [('e2',.01),('e1',0.)]:
        catalogues[exposure]=[{'xy':[i*50+shift,40.],
            'ra_dec_deg':[69.55+i*.001+shift/3600,-12.28]} for i in range(8)]
    groups=assign_physical_stars(catalogues)
    assert len(groups)==8
    assert sum(g['split']=='heldout' for g in groups)==2
    assert all(len(g['detections'])==2 for g in groups)
    assert groups==assign_physical_stars(dict(reversed(list(catalogues.items()))))


def test_empirical_psf_is_pixel_response_without_second_integration():
    grid=gaussian_grid();model=epsf_model(grid)
    yy,xx=np.mgrid[-4:5,-4:5]
    actual=model.evaluate(xx,yy,1.,.23,-.17)
    expected=gaussian_prf(xx,yy,.23,-.17)
    # Four-times sampling and cubic interpolation have a small finite error;
    # integrating this already integrated response again differs far more.
    np.testing.assert_allclose(actual,expected,atol=3e-5,rtol=0)
    assert np.sum(grid)/16==pytest.approx(1.)


def test_missing_crossmatches_cannot_silently_authorize_common_pool():
    singleton=[{'detections':[{'exposure':name,'xy':[10.,10.]}]} for name in ['a','b','c','d']]
    result=common_pool_dither_diagnostic(singleton,['a','b','c','d'],64.)
    assert not result['small_dither_condition_pass']
    assert result['maximum_dither_native_pixel'] is None
    matched=[{'detections':[{'exposure':name,'xy':[i*270.,20.]} for i,name in enumerate(['a','b','c','d'])]}]
    result=common_pool_dither_diagnostic(matched,['a','b','c','d'],64.)
    assert result['all_exposure_pairs_have_matched_training_star'] and not result['small_dither_condition_pass']
    for row in matched[0]['detections']:row['xy'][0]/=100
    assert common_pool_dither_diagnostic(matched,['a','b','c','d'],64.)['small_dither_condition_pass']


def test_known_ERR_centroid_covariance_is_not_residual_scaled_and_bad_PSF_fails():
    yy,xx=np.indices((41,41));grid=gaussian_grid();model=epsf_model(grid)
    data=model.evaluate(xx,yy,10000.,20.23,19.83)+2.+.01*(xx-20)
    error=np.ones(data.shape);dq=np.zeros(data.shape,np.uint32)
    fit=fit_heldout_field_star(data,error,dq,[20.,20.],grid)
    doubled=fit_heldout_field_star(data,error*2,dq,[20.,20.],grid)
    np.testing.assert_allclose(fit['centroid_xy_pixel'],[20.23,19.83],atol=1e-7)
    np.testing.assert_allclose(doubled['centroid_covariance_pixel2'],np.array(fit['centroid_covariance_pixel2'])*4,rtol=2e-6)
    assert fit['absolute_adequacy_pass'] and not fit['covariance_rescaled']
    bad=data.copy();bad[19,22]+=1000
    failed=fit_heldout_field_star(bad,error,dq,[20.,20.],grid)
    assert not failed['absolute_adequacy_pass'] and failed['p_value']<.001
    assert not failed['covariance_rescaled']


def test_epsf_builder_can_represent_diffraction_ring_calibration_stamps():
    # Independent pixel quadrature of an Airy intensity, rather than a Gaussian
    # or the ePSF interpolator used by the fitted model.
    yy,xx=np.indices((25,25));z,w=np.polynomial.legendre.leggauss(5)
    def airy(cx,cy):
        result=np.zeros(xx.shape)
        for ox,wx in zip(z/2,w/2):
            for oy,wy in zip(z/2,w/2):
                rr=np.hypot(xx-cx+ox,yy-cy+oy)*2.4
                value=np.ones(rr.shape);use=rr!=0;value[use]=(2*j1(rr[use])/rr[use])**2
                result+=wx*wy*value
        return result/result.sum()
    stamps=[]
    rng=np.random.default_rng(927)
    for i in range(12):
        center=np.array([12.,12.])+rng.uniform(-.45,.45,2)
        values=10000*airy(*center)+rng.normal(0,.1,(25,25))
        stamps.append({'id':str(i),'data':values,'error':np.full((25,25),.1),
                       'dq':np.zeros((25,25),np.uint32),'center':center})
    grid,diagnostic=build_empirical_psf(stamps)
    assert np.isfinite(grid).all() and min(grid.shape)>=97
    assert grid.sum()/16==pytest.approx(1.)
    # A real diffraction ring survives in the ePSF beyond the core.
    model=epsf_model(grid)
    assert model.evaluate(3.5,0.,1.,0.,0.)>model.evaluate(3.0,0.,1.,0.,0.)
    assert diagnostic['iterations']<=DEFAULT_ALGORITHM['epsf_maxiters']


def test_execution_gate_rejects_before_any_protocol_or_native_path_access():
    spec=importlib.util.spec_from_file_location('field_calibration_runner',ROOT/'scripts/calibrate_he0435_jwst_psf.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    with pytest.raises(PermissionError,match='no native data opened'):
        module.validate_execution(Path('does-not-exist-protocol.json'),'invalid',False)


def test_epsf_internal_exclusions_cannot_bypass_minimum_training_stars(monkeypatch):
    import photutils.psf
    monkeypatch.setattr(photutils.psf,'EPSFBuilder',lambda **kwargs:
        lambda stars:SimpleNamespace(n_excluded_stars=1))
    yy,xx=np.indices((25,25))
    stamps=[{'id':str(i),'data':10000*gaussian_prf(xx,yy,12.,12.),
             'error':np.ones((25,25)),'dq':np.zeros((25,25),np.uint32),
             'center':[12.,12.]} for i in range(3)]
    with pytest.raises(ValueError,match='frozen minimum'):
        build_empirical_psf(stamps)


def test_frozen_binding_includes_transitively_consumed_pixel_geometry(tmp_path):
    spec=importlib.util.spec_from_file_location('field_calibration_runner_binding',ROOT/'scripts/calibrate_he0435_jwst_psf.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    algorithm='preregistration/he0435_field_psf_algorithm_20261021.json'
    binding={'algorithm_file':algorithm,'algorithm_sha256':module.sha(ROOT/algorithm),
             'source_sha256':{p:module.sha(ROOT/p) for p in module.SOURCES},
             'native_products':[f'jw01198007001_02101_0000{i}_nrcb3_cal.fits' for i in range(1,5)]}
    assert 'code/strong_lensing_pixel_geometry.py' in binding['source_sha256']
    path=tmp_path/'protocol.json';path.write_text(json.dumps({'field_psf_calibration':binding}),encoding='utf-8')
    cfg,products=module.validate_execution(path,module.sha(path),True)
    assert len(products)==4 and not cfg['allow_common_pool_for_this_visit']
    binding['source_sha256']['code/strong_lensing_pixel_geometry.py']='0'*64
    path.write_text(json.dumps({'field_psf_calibration':binding}),encoding='utf-8')
    with pytest.raises(ValueError,match='source changed'):
        module.validate_execution(path,module.sha(path),True)


def test_field_receipt_numpy_encoding_preserves_finite_json_guard(tmp_path):
    spec=importlib.util.spec_from_file_location('field_calibration_runner_json',ROOT/'scripts/calibrate_he0435_jwst_psf.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    path=tmp_path/'receipt.json'
    module.dump(path,{'bounds':[np.int64(12),np.int64(25)],'grid':np.array([[np.float64(.2)]])})
    assert json.loads(path.read_text())=={'bounds':[12,25],'grid':[[.2]]}
    with pytest.raises(ValueError,match='Out of range'):
        module.dump(path,{'invalid':np.array([np.nan])})
