from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import tempfile
import types

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
            with patch.dict(sys.modules,{'stdatamodels':types.ModuleType('stdatamodels'),'stdatamodels.jwst':fake}):
                result=apply_native_files(paths,protocol,models,out)
            self.assertEqual(result['complete_exposures'],4)
            self.assertTrue(result['all_pixel_models_adequate'])
            self.assertTrue(result['exposures_consistent'])
            self.assertFalse(result['combined']['covariance_rescaled'])
            self.assertEqual(result['combined']['exposure_consistency_dof'],15)
            self.assertTrue(result['models']['revised']['conditional_geometry_validation_pass'])
            self.assertFalse(result['models']['old']['conditional_geometry_validation_pass'])


if __name__=='__main__':unittest.main()
