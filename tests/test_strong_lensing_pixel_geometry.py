from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import tempfile
import types
import json
import hashlib

import numpy as np
from scipy.special import ndtr

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_pixel_geometry import (
    native_pixel_mask, integrated_moffat, fit_moffat_point, intrinsic_geometry,
    intrinsic_jacobian, intrinsic_with_covariance, combine_exposures,
    gaussian_score, tangent_plane_arcsec, gwcs_centroid,empirical_psf_values,
    fit_empirical_point)
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from run_he0435_pixel_geometry import score_models,apply_native_files
import run_he0435_pixel_geometry as pixel_runner


class PixelGeometryTests(unittest.TestCase):
    def setUp(self):
        self.points=np.array([[.2,.1],[-1.27,.65],[-2.26,-.5],[-.74,-1.52]])

    def test_intrinsic_retains_scale_rotation_translation_and_chirality(self):
        angle=.43; r=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        q=intrinsic_geometry(self.points)
        np.testing.assert_allclose(intrinsic_geometry(self.points@r.T+[3.,-2.]),q,atol=1e-14)
        np.testing.assert_allclose(intrinsic_geometry(self.points*2),q*2)
        reflected=intrinsic_geometry(self.points*[1,-1])
        np.testing.assert_allclose(reflected,q*[1,1,-1,1,-1])

    def test_intrinsic_jacobian_and_correlated_covariance(self):
        j=intrinsic_jacobian(self.points); eps=1e-6
        numerical=np.column_stack([(intrinsic_geometry((self.points.ravel()+np.eye(8)[i]*eps).reshape(4,2))
              -intrinsic_geometry((self.points.ravel()-np.eye(8)[i]*eps).reshape(4,2)))/(2*eps) for i in range(8)])
        np.testing.assert_allclose(j,numerical,atol=5e-10)
        q,c=intrinsic_with_covariance(self.points,np.eye(8)*1e-6)
        self.assertEqual(np.linalg.matrix_rank(c),5)
        self.assertGreater(abs(c[1,3]),0.)
        np.testing.assert_allclose(j@np.tile([1,0],4),0.,atol=1e-15)

    def test_only_declared_dq_flags_are_masked(self):
        d=np.ones((2,3)); e=np.ones_like(d); dq=np.array([[0,1,2],[4,8,0]])
        d[1,2]=np.nan
        np.testing.assert_array_equal(native_pixel_mask(d,e,dq),[[True,False,False],[True,True,False]])
        with self.assertRaises(ValueError):native_pixel_mask(d,e,np.zeros((2,2)))

    def test_pixel_integral_and_phase_centroid_recovery(self):
        yy,xx=np.indices((25,25)); center=np.array([12.24,11.73]); sigma=.15
        image=integrated_moffat(xx,yy,center,2000.,1.1,1.4,.27,2.7)+3+.03*(xx-12)-.02*(yy-12)
        high=integrated_moffat(xx,yy,center,2000.,1.1,1.4,.27,2.7,quadrature_order=11)
        np.testing.assert_allclose(image-3-.03*(xx-12)+.02*(yy-12),high,atol=3e-4)
        fit=fit_moffat_point(image,np.full_like(image,sigma),np.zeros_like(image,dtype='uint32'),[12,12])
        np.testing.assert_allclose(fit['centroid_xy_pixel'],center,atol=1e-6)
        self.assertLess(fit['chi2'],1e-7)
        self.assertTrue(fit['valid_measurement_fit'])
        cov=np.array(fit['centroid_covariance_pixel2'])
        self.assertTrue(np.all(np.linalg.eigvalsh(cov)>0))
        self.assertGreater(np.trace(cov),0.) # Never becomes zero at zero residual.

    def test_known_covariance_scatter_fails_without_covariance_inflation(self):
        values=np.array([[0.,0.],[1.,0.]])
        c=np.array([np.eye(2)*.01]*2)
        result=combine_exposures(values,c)
        np.testing.assert_allclose(result['covariance'],np.eye(2)*.005)
        self.assertAlmostEqual(result['exposure_consistency_chi2'],50.)
        self.assertLess(result['exposure_consistency_p_value'],1e-9)
        self.assertFalse(result['covariance_rescaled'])

    def test_gaussian_score_retains_normalization_and_failure(self):
        result=gaussian_score([1.,2.],np.eye(2))
        self.assertAlmostEqual(result['chi2'],5.)
        self.assertAlmostEqual(result['gaussian_logdensity'],-2.5-np.log(2*np.pi))
        with self.assertRaises(np.linalg.LinAlgError):gaussian_score([0.,0.],np.zeros((2,2)))

    def test_gwcs_mapping_preserves_nonlinear_local_derivative(self):
        origin=[69.562,-12.287]
        def fake_wcs(x,y):return origin[0]+(x+.001*x*x)/3600/np.cos(np.deg2rad(origin[1])),origin[1]+y/3600
        p,c,j=gwcs_centroid(fake_wcs,[2.,3.],np.diag([.01,.04]),origin)
        self.assertAlmostEqual(p[0],2.004,places=4)
        self.assertAlmostEqual(j[0,0],1.004,places=4)
        self.assertAlmostEqual(c[0,0],.01*1.004**2,places=6)
        np.testing.assert_allclose(tangent_plane_arcsec([origin],origin),[[0.,0.]],atol=1e-10)

    def test_degenerate_domain_rejected(self):
        p=self.points.copy();p[1]=p[0]
        with self.assertRaises(ValueError):intrinsic_geometry(p)
        with self.assertRaises(ValueError):intrinsic_with_covariance(self.points,np.eye(6))
        with self.assertRaises(ValueError):combine_exposures([],[])

    def test_empirical_psf_is_not_integrated_twice_and_recovers_phase(self):
        over=8;gsize=153;origin=[76.,76.];sigma=.85
        gy,gx=np.indices((gsize,gsize));px=(gx-origin[0])/over;py=(gy-origin[1])/over
        profile=lambda z:ndtr((z+.5)/sigma)-ndtr((z-.5)/sigma)
        grid=profile(px)*profile(py)
        yy,xx=np.indices((25,25));center=[12.23,11.81]
        true=1500*profile(xx-center[0])*profile(yy-center[1])+1+.02*xx-.01*yy
        fit=fit_empirical_point(true,np.ones_like(true)*.1,np.zeros_like(true,dtype='uint32'),
                                 [12,12],grid,over,origin)
        np.testing.assert_allclose(fit['centroid_xy_pixel'],center,atol=2e-5)
        self.assertLess(fit['chi2'],.001)
        self.assertTrue(fit['valid_measurement_fit'])
        values=empirical_psf_values(grid,over,origin,[12.23],[11.81],center)
        np.testing.assert_allclose(values,[profile(0.)**2],atol=1e-14)
        self.assertEqual(empirical_psf_values(grid,over,origin,[100.],[100.],center)[0],0.)

    def test_revision_improvement_is_frozen_twenty_five_percent(self):
        old=(self.points[1:]-self.points[0]).ravel()
        improved=old.copy();improved[0]+=.0003
        models={'old':{'relative_B_C_D_to_A_east_north_arcsec':old.tolist()},
                'revised':{'relative_B_C_D_to_A_east_north_arcsec':improved.tolist(),'baseline_model':'old'}}
        truth=intrinsic_geometry(np.vstack([[0,0],improved.reshape(3,2)]))
        # An observation nearly midway gives <25% reduction despite improvement.
        oldq=intrinsic_geometry(np.vstack([[0,0],old.reshape(3,2)]))
        observed=oldq+.55*(truth-oldq)
        result=score_models(observed,np.eye(5),models,.05,.25)['revised']
        self.assertTrue(result['improves_declared_baseline_raw_q_RMS'])
        self.assertFalse(result['meets_declared_baseline_raw_q_RMS_reduction'])
        result=score_models(truth,np.eye(5),models,.05,.25)['revised']
        self.assertTrue(result['meets_declared_baseline_raw_q_RMS_reduction'])

    def test_quadratic_background_recovers_centroid_bias_without_changing_PSF_or_ERR(self):
        over=8;origin=[76.,76.];yyg,xxg=np.indices((153,153));sigma=.9
        profile=lambda z:ndtr((z+.5)/sigma)-ndtr((z-.5)/sigma)
        grid=profile((xxg-76)/over)*profile((yyg-76)/over)
        yy,xx=np.indices((25,25));dx,dy=xx-12,yy-12;center=[12.34,11.73]
        # Independent analytic native-pixel PRF plus smooth diffuse curvature.
        data=1000*profile(xx-center[0])*profile(yy-center[1])+5+.04*dx-.06*dy+.5*dx*dx+.3*dx*dy+.2*dy*dy
        err=np.full_like(data,.2);dq=np.zeros(data.shape,dtype='uint32')
        affine=fit_empirical_point(data,err,dq,[12,12],grid,over,origin)
        curved=fit_empirical_point(data,err,dq,[12,12],grid,over,origin,background_degree=2)
        self.assertGreater(np.linalg.norm(np.asarray(affine['centroid_xy_pixel'])-center),.004)
        self.assertLess(affine['p_value'],.001)
        np.testing.assert_allclose(curved['centroid_xy_pixel'],center,atol=2e-6)
        self.assertGreater(curved['p_value'],.999)
        self.assertEqual(curved['jacobian_rank'],9)
        np.testing.assert_allclose(curved['background_polynomial_SCI'],[5,.04,-.06,.5,.3,.2],atol=1e-5)
        self.assertTrue(np.all(np.linalg.eigvalsh(curved['centroid_covariance_pixel2'])>0))

    def test_nested_per_label_PSF_hashes_and_all_sixteen_measurement_routes(self):
        from astropy.io import fits
        origin=np.array([69.562,-12.287]);scale=.031
        positions=np.array([[35.,35.],[75.,35.],[35.,75.],[75.,75.]])
        class FakeWCS:
            def __call__(self,x,y):return (origin[0]+np.asarray(x)*scale/3600/np.cos(np.deg2rad(origin[1])),origin[1]+np.asarray(y)*scale/3600)
            def invert(self,ra,dec,with_bounding_box=False):return ((ra-origin[0])*3600*np.cos(np.deg2rad(origin[1]))/scale,(dec-origin[1])*3600/scale)
        class FakeImage:
            def __init__(self):
                self.meta=types.SimpleNamespace(wcs=FakeWCS(),instrument=types.SimpleNamespace(filter='F115W',detector='NRCB3'))
                self.data=np.zeros((120,120));self.err=np.ones((120,120));self.dq=np.zeros((120,120),np.uint32)
            def __enter__(self):return self
            def __exit__(self,*args):pass
        refs=np.array([FakeWCS()(*p) for p in positions]);sky=tangent_plane_arcsec(refs,origin)
        truth=(sky[1:]-sky[0]).ravel();old=truth.copy();old[0]+=.005
        models={'old':{'relative_B_C_D_to_A_east_north_arcsec':old.tolist()},
                'revised':{'relative_B_C_D_to_A_east_north_arcsec':truth.tolist(),'baseline_model':'old'}}
        fake=types.ModuleType('stdatamodels.jwst');fake.datamodels=types.SimpleNamespace(open=lambda f:FakeImage())
        seen=[]
        def extract(data,error,dq,guess,grid,*args):
            seen.append(float(grid[0,0]))
            return {'centroid_xy_pixel':np.asarray(guess).tolist(),'centroid_covariance_pixel2':(np.eye(2)*1e-4).tolist(),
                    'p_value':.5,'valid_measurement_fit':True}
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);mapping={};hashes={};native=[]
            for i in range(4):
                path=root/f'cal{i}';path.write_bytes(b'SYNTHETICONLY');native.append(path);mapping[path.name]={}
                for k,label in enumerate('ABCD'):
                    relative=f'psf{i}{label}.fits';hdu=fits.PrimaryHDU(np.ones((9,9))*(i*4+k+1))
                    hdu.header['OVERSMP']=4;hdu.header['ORIGINX']=4.;hdu.header['ORIGINY']=4.;hdu.writeto(root/relative)
                    mapping[path.name][label]=relative;hashes[relative]=hashlib.sha256((root/relative).read_bytes()).hexdigest()
            calibration={'overall_adequacy_pass':True,'psf_files_by_exposure':mapping,'file_sha256':hashes}
            calpath=root/'calibration.json';calpath.write_text(json.dumps(calibration),encoding='utf-8')
            protocol={'psf_model':'field_empirical','background_degree':2,
                      'psf_calibration':{'path':'calibration.json','sha256':hashlib.sha256(calpath.read_bytes()).hexdigest()},
                      'reference_image_radec_deg':refs.tolist(),'tangent_origin_radec_deg':origin.tolist()}
            with patch.object(pixel_runner,'ROOT',root),patch.object(pixel_runner,'verify_frozen') as verified:
                checked=pixel_runner.frozen_psf_calibration(protocol,'synthetic-commit')
                self.assertEqual(verified.call_count,17) # Manifest and every distinct grid.
                with patch.dict(sys.modules,{'stdatamodels':types.ModuleType('stdatamodels'),'stdatamodels.jwst':fake,'asdf':types.ModuleType('asdf')}),patch.object(pixel_runner,'fit_empirical_point',side_effect=extract):
                    result=apply_native_files(native,protocol,models,root/'result',checked)
                self.assertEqual(seen,list(np.arange(1,17,dtype=float)))
                self.assertTrue(result['models']['revised']['conditional_geometry_validation_pass'])
                # A known failed field calibration is rejected by default. Its
                # explicit diagnostic mode cannot pass, even with perfect q.
                calibration['overall_adequacy_pass']=False
                calpath.write_text(json.dumps(calibration),encoding='utf-8')
                protocol['psf_calibration']['sha256']=hashlib.sha256(calpath.read_bytes()).hexdigest()
                with self.assertRaisesRegex(RuntimeError,'targets must remain sealed'):
                    pixel_runner.frozen_psf_calibration(protocol,'synthetic-commit')
                diagnostic_protocol={**protocol,'allow_inadequate_field_calibration_for_diagnostic':True}
                failed_calibration=pixel_runner.frozen_psf_calibration(diagnostic_protocol,'synthetic-commit')
                with patch.dict(sys.modules,{'stdatamodels':types.ModuleType('stdatamodels'),'stdatamodels.jwst':fake,'asdf':types.ModuleType('asdf')}),patch.object(pixel_runner,'fit_empirical_point',side_effect=extract):
                    diagnostic=apply_native_files(native,diagnostic_protocol,models,root/'diagnostic-result',failed_calibration)
                self.assertTrue(diagnostic['all_pixel_models_adequate'])
                self.assertTrue(diagnostic['exposures_consistent'])
                self.assertTrue(diagnostic['models']['revised']['formal_geometry_compatible_at_alpha'])
                self.assertTrue(diagnostic['models']['revised']['meets_declared_baseline_raw_q_RMS_reduction'])
                self.assertFalse(diagnostic['calibration_adequacy_pass'])
                self.assertTrue(diagnostic['diagnostic_mode'])
                self.assertEqual(diagnostic['measurement_status'],'CALIBRATION_INADEQUATE_DIAGNOSTIC_ONLY')
                self.assertTrue(all(not score['conditional_geometry_validation_pass'] for score in diagnostic['models'].values()))
                for relative in diagnostic['exposure_result_files']:
                    exposure=json.loads((root/'diagnostic-result'/relative).read_text(encoding='utf-8'))
                    self.assertTrue(exposure['diagnostic_mode'])
                    self.assertFalse(exposure['calibration_adequacy_pass'])
                    self.assertEqual(exposure['measurement_status'],'CALIBRATION_INADEQUATE_DIAGNOSTIC_ONLY')
                    self.assertTrue(all(not score['conditional_geometry_validation_pass'] for score in exposure['model_scores'].values()))
                (root/'psf3D.fits').write_bytes(b'CORRUPTED')
                with self.assertRaisesRegex(RuntimeError,'grid hash changed'):
                    pixel_runner.frozen_psf_calibration(diagnostic_protocol,'synthetic-commit')
            with self.assertRaises(ValueError):pixel_runner.psf_paths_by_label({'psf_files_by_exposure':{'cal0':{'A':'a'}}},'cal0')
            self.assertEqual(pixel_runner.psf_paths_by_label({'psf_files_by_exposure':{'cal0':'legacy.fits'}},'cal0'),{x:'legacy.fits' for x in 'ABCD'})

    def test_frozen_field_registered_initial_association_does_not_use_sky_fit(self):
        class WrongPointing:
            def invert(self,*args,**kwargs):raise AssertionError('No Gaia inverse WCS call after frozen field registration')
        reference=np.array([[30.,40.],[50.,40.],[30.,60.],[50.,60.]])
        protocol={'reference_guesses_xy_pixel_by_exposure':{'native.fits':reference.tolist()}}
        actual=pixel_runner.reference_pixel_guesses(protocol,'native.fits',WrongPointing(),np.zeros((4,2)))
        np.testing.assert_array_equal(actual,reference)
        np.testing.assert_allclose(intrinsic_geometry(actual+[17.,-23.]),intrinsic_geometry(actual))
        with self.assertRaises(KeyError):pixel_runner.reference_pixel_guesses(protocol,'not-frozen.fits',WrongPointing(),np.zeros((4,2)))

    def test_native_identity_guard_checks_exact_four_hashes_before_array_application(self):
        with tempfile.TemporaryDirectory() as temp:
            paths=[];hashes={}
            for i in range(4):
                path=Path(temp)/f'cal{i}.fits';path.write_bytes(f'NOTANIMAGE{i}'.encode())
                paths.append(path);hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
            protocol={'primary_observation':{'products':[p.name for p in paths]}}
            pixel_runner.verify_native_inputs(paths,protocol,{'input_sha256':hashes})
            with self.assertRaisesRegex(RuntimeError,'four declared'):
                pixel_runner.verify_native_inputs(paths[:3],protocol,{'input_sha256':hashes})
            paths[3].write_bytes(b'CHANGED')
            with self.assertRaisesRegex(RuntimeError,'hash missing or changed'):
                pixel_runner.verify_native_inputs(paths,protocol,{'input_sha256':hashes})

    def test_synthetic_native_runner_preserves_full_covariance_and_gates(self):
        origin=np.array([69.562,-12.287]);scale=.031
        relative=np.array([[0.,0.],[-1.474,.552],[-2.466,-.602],[-.938,-1.616]])
        class FakeWCS:
            def __init__(self,phase):self.phase=phase
            def __call__(self,x,y):
                return (origin[0]+(x-150-self.phase)*scale/3600/np.cos(np.deg2rad(origin[1])),
                        origin[1]+(y-150-self.phase)*scale/3600)
            def invert(self,ra,dec,with_bounding_box=False):
                return ((ra-origin[0])*3600*np.cos(np.deg2rad(origin[1]))/scale+150+self.phase,
                        (dec-origin[1])*3600/scale+150+self.phase)
        class FakeImage:
            def __init__(self,phase):
                self.meta=types.SimpleNamespace(wcs=FakeWCS(phase),
                      instrument=types.SimpleNamespace(filter='F115W',detector='NRCB3'))
                yy,xx=np.indices((220,220));self.data=np.zeros((220,220))
                for point in relative/scale+150+phase:
                    self.data+=integrated_moffat(xx,yy,point,2000.,1.2,1.1,.2,2.8)
                self.data+=1+.001*xx-.002*yy
                self.err=np.ones_like(self.data)*.2
                self.dq=np.zeros_like(self.data,dtype='uint32')
            def __enter__(self):return self
            def __exit__(self,*args):pass
        refs=np.array([FakeWCS(0)(*(v/scale+150)) for v in relative])
        points=tangent_plane_arcsec(refs,origin)
        truth=(points[1:]-points[0]).ravel()
        old=truth.copy();old[0]+=.005
        models={'old':{'relative_B_C_D_to_A_east_north_arcsec':old.tolist()},
                'revised':{'relative_B_C_D_to_A_east_north_arcsec':truth.tolist(),'baseline_model':'old'}}
        fake=types.ModuleType('stdatamodels.jwst')
        fake.datamodels=types.SimpleNamespace(open=lambda f:FakeImage(int(Path(f).stem[-1])*.11))
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp)/'result';paths=[]
            for i in range(4):
                path=Path(temp)/f'fakecal{i}';path.write_bytes(b'SYNTHETICONLY');paths.append(path)
            protocol={'reference_image_radec_deg':refs.tolist(),'tangent_origin_radec_deg':origin.tolist()}
            with patch.dict(sys.modules,{'stdatamodels':types.ModuleType('stdatamodels'),'stdatamodels.jwst':fake,'asdf':types.ModuleType('asdf')}):
                result=apply_native_files(paths,protocol,models,out)
            self.assertEqual(result['complete_exposures'],4)
            self.assertTrue(result['all_pixel_models_adequate'])
            self.assertTrue(result['exposures_consistent'])
            self.assertFalse(result['combined']['covariance_rescaled'])
            self.assertEqual(result['combined']['exposure_consistency_dof'],15)
            self.assertTrue(result['models']['revised']['conditional_geometry_validation_pass'])
            self.assertFalse(result['models']['old']['conditional_geometry_validation_pass'])


if __name__=='__main__':unittest.main()
