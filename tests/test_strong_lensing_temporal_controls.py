from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from strong_lensing_temporal_controls import affine_free_temporal_contrasts,condition_noise_mean


def test_seven_contrasts_remove_arbitrary_static_affine_light_at_actual_times():
    times=np.array([16.105155,69.789005,123.472855,177.156705,230.840555,284.524405,338.208255,391.892105,445.575955])
    p=affine_free_temporal_contrasts(times)
    np.testing.assert_allclose(p@p.T,np.eye(7),atol=4e-16)
    np.testing.assert_allclose(p@(418.+13.2*times),0.,atol=2e-12)


def test_conditional_mean_is_invariant_to_contrast_basis_and_reduces_covariance():
    # A common latent electronic mode gives an exact, nonzero conditional mean.
    ct=np.diag([2.,3.])+np.outer([1.,2.],[1.,2.])
    cc=np.diag([.5,.7])+np.outer([.3,.4],[.3,.4])
    cross=np.outer([1.,2.],[.3,.4]);z=np.array([.2,-.1])
    mean,cov,transfer=condition_noise_mean(ct,cross,cc,z)
    expected=np.array([1.,2.])*(np.array([.3,.4])@np.linalg.solve(cc,z))
    np.testing.assert_allclose(mean,expected,rtol=1e-14,atol=1e-16)
    angle=.7;o=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    rotated,rotcov,_=condition_noise_mean(ct,cross@o.T,o@cc@o.T,o@z)
    np.testing.assert_allclose(rotated,mean,atol=1e-16)
    np.testing.assert_allclose(rotcov,cov,atol=1e-15)
    assert np.linalg.eigvalsh(cov).min()>0
    assert np.linalg.eigvalsh(ct-cov).min()>-1e-15
