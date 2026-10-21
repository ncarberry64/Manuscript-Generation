import sys
from pathlib import Path
import numpy as np
from scipy.signal import convolve2d
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"code"))
from strong_lensing_ipc_response import (
    inverse_deconvolution_kernel, convolution_matrix, photon_to_raw_dn_jacobian)


def test_positive_inverse_reconstructs_stored_reference_without_normalization():
    g = np.array([[-.000217629422, -.00576349022, -.000217629422],
                  [-.00576349022, 1.02392447, -.00576349022],
                  [-.000217629422, -.00576349022, -.000217629422]])
    r = inverse_deconvolution_kernel(g, 8)
    assert r["kernel"].min() >= 0
    replay = convolve2d(g, r["kernel"])
    expected = np.zeros_like(replay); expected[9, 9] = 1.
    assert np.max(abs(replay-expected)) < 2*r["tail_l1_bound"]
    assert abs(r["inverse_sum"]-1/g.sum()) < r["tail_l1_bound"]*2
    assert r["renormalized"] is False


def test_source_gain_and_nonlinearity_order_in_photon_covariance():
    # Nonuniform gains and derivatives make H*D distinguishable from D*H.
    h = convolution_matrix(np.array([[0., .1, 0.], [.1, .6, .1], [0., .1, 0.]]), (3, 3))
    gain = np.arange(1., 10.)
    derivative = np.linspace(.8, 1.4, 9)
    j = photon_to_raw_dn_jacobian(h, gain, derivative)
    q = np.arange(10., 19.)
    covariance = (j@sparse.diags(q)@j.T).toarray()
    assert np.linalg.eigvalsh(covariance).min() >= -1e-14
    assert covariance[4, 5] > 0
    # This is a photon impulse; electronic post-IPC noise is not reconvolved.
    impulse = np.zeros(9); impulse[4] = 1.
    assert np.allclose(j@impulse, h.toarray()[:, 4]*derivative[4]/gain[4])
    assert not np.allclose(j.toarray(), sparse.diags(derivative/gain)@h.toarray())


def test_nonsymmetric_convolution_orientation_and_explicit_exterior():
    kernel = np.array([[0., 0., 0.], [0., .8, .2], [0., 0., 0.]])
    field = np.arange(20.).reshape(4, 5)
    matrix = convolution_matrix(kernel, field.shape)
    assert np.allclose((matrix@field.ravel()).reshape(field.shape),
                       convolve2d(field, kernel, mode="same", boundary="fill"))
    assert matrix.toarray()[0].sum() == .8  # outside charge was explicitly set zero
