"""Independent library calibration tests; no external PSF or target array read."""
import importlib.util
from pathlib import Path
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_standard_psf import (stellar_identity,propagate_gaia_field_position,
    associate_saturated_field_core,sample_spatial_standard_psf,registered_reference_guesses,apply_minimum_pixel_gate)
from strong_lensing_field_psf import fit_heldout_field_star,epsf_model


def test_catalogue_identity_and_alpha_star_proper_motion_are_independent_of_PSF():
    row={'classprob_dsc_combmod_star':'.999','parallax':'.2','parallax_error':'.2',
         'pmra':'6.','pmdec':'-4.','pmra_error':'.4','pmdec_error':'.4',
         'pmra_pmdec_corr':'0.','duplicated_source':'false','ruwe':'1.0',
         'ra':'69.5','dec':'-12.3','ref_epoch':'2016.0'}
    assert stellar_identity(row)['accepted']
    out=propagate_gaia_field_position(row,2022.)
    assert (out[0]-69.5)*3600000*np.cos(np.deg2rad(-12.3))==pytest.approx(36.,abs=1e-4)
    assert (out[1]+12.3)*3600000==pytest.approx(-24.,abs=1e-4)
    row['duplicated_source']='true';assert not stellar_identity(row)['accepted']


def test_saturated_component_association_uses_DQ_only_and_enforces_field_exclusion():
    dq=np.zeros((100,100),np.uint32);field=np.ones((100,100),bool)
    dq[49:52,49:52]=2;field[48:53,48:53]=False
    dq[60:63,55:58]=2;dq[80:83,80:83]=2
    result=associate_saturated_field_core(dq,field,[50.,50.],40.)
    np.testing.assert_allclose(result['center_xy_pixel'],[56.,61.])
    assert result['pixels']==9 and not result['association_uses_SCI']
    with pytest.raises(ValueError,match='No eligible'):
        associate_saturated_field_core(dq,field,[50.,50.],5.)


def test_spatial_grid_sampling_changes_response_without_translating_PSFor_reintegrating():
    class SpatialPSF:
        def evaluate(self,x,y,flux,x_0,y_0):
            sigma=.7+x_0*.0001
            return flux*np.exp(-((x-x_0)**2+(y-y_0)**2)/(2*sigma**2))/(2*np.pi*sigma**2)
    a=sample_spatial_standard_psf(SpatialPSF(),[100.,200.])
    b=sample_spatial_standard_psf(SpatialPSF(),[1500.,1000.])
    assert np.unravel_index(np.argmax(a),a.shape)==(50,50)
    assert np.unravel_index(np.argmax(b),b.shape)==(50,50)
    assert a.sum()/16==pytest.approx(1.) and b.sum()/16==pytest.approx(1.)
    assert a[50,50]>b[50,50]


def test_field_translation_changes_only_initial_association_and_preserves_quad_differences():
    class LinearWCS:
        def invert(self,a,d,**kwargs):return a*10,d*10
    sky=np.array([[1.,2.],[1.1,2.2],[1.3,2.1],[1.4,1.8]])
    corrected,delta=registered_reference_guesses(LinearWCS(),sky,[50.,60.],[32.,31.])
    np.testing.assert_allclose(delta,[-18.,-29.])
    np.testing.assert_allclose(corrected[1:]-corrected[0],(sky[1:]-sky[0])*10)


def test_official_photutils_gridded_PSF_interface_with_synthetic_library():
    from astropy.nddata import NDData
    from photutils.psf import GriddedPSFModel
    v=(np.arange(101)-50)/4;x,y=np.meshgrid(v,v)
    grids=[]
    locations=[(0.,0.),(2048.,0.),(0.,2048.),(2048.,2048.)]
    for cx,cy in locations:
        sigma=.7+cx*.0001;grid=np.exp(-(x*x+y*y)/(2*sigma**2));grids.append(grid*16/grid.sum())
    standard=GriddedPSFModel(NDData(np.array(grids),meta={'grid_xypos':locations,'oversampling':4}))
    a=sample_spatial_standard_psf(standard,[100.,900.])
    b=sample_spatial_standard_psf(standard,[1900.,900.])
    assert a.shape==(101,101) and b.shape==a.shape
    assert a.sum()/16==pytest.approx(1.) and a[50,50]>b[50,50]


def test_saturated_core_values_never_enter_known_ERR_wing_fit():
    v=(np.arange(101)-50)/4;x,y=np.meshgrid(v,v);grid=np.exp(-(x*x+y*y)/(2*.8**2));grid*=16/grid.sum()
    yy,xx=np.indices((41,41));model=epsf_model(grid)
    data=model.evaluate(xx,yy,10000.,20.1,19.8)+1.
    error=np.ones(data.shape);dq=np.zeros(data.shape,np.uint32);dq[18:22,18:22]=2
    outputs=[]
    for value in [0.,1e20]:
        changed=data.copy();changed[dq==2]=value
        outputs.append(fit_heldout_field_star(changed,error,dq,[20.,20.],grid))
    np.testing.assert_allclose(outputs[0]['centroid_xy_pixel'],[20.1,19.8],atol=1e-6)
    assert outputs[0]==outputs[1] and outputs[0]['absolute_adequacy_pass']
    assert not outputs[0]['covariance_rescaled']


def test_external_execution_gate_precedes_any_file_or_array_access():
    spec=importlib.util.spec_from_file_location('standard_psf_runner',ROOT/'scripts/calibrate_he0435_standard_psf.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    with pytest.raises(PermissionError,match='field-only'):
        module.frozen_inputs(Path('nonexistent.json'),'invalid',False)


def test_frozen_pixel_count_gate_cannot_be_ignored_or_rescue_failed_PSF():
    assert not apply_minimum_pixel_gate({'usable_pixels':49,'absolute_adequacy_pass':True},50)['absolute_adequacy_pass']
    assert apply_minimum_pixel_gate({'usable_pixels':50,'absolute_adequacy_pass':True},50)['absolute_adequacy_pass']
    assert not apply_minimum_pixel_gate({'usable_pixels':289,'absolute_adequacy_pass':False},50)['absolute_adequacy_pass']
    assert not apply_minimum_pixel_gate({'usable_pixels':100,'absolute_adequacy_pass':True},101)['absolute_adequacy_pass']
