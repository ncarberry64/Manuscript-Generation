"""Conditional signal-independent IPC inverse and its photon-count map.

CRDS IPC files contain a DECONVOLUTION G in raw DN, not a positive photon
redistribution kernel.  H=G^{-1} below is a conditional inverse of that exact
stored reference; it does not identify nonlinear/signal-dependent IPC.
"""
import numpy as np
from scipy.signal import convolve2d
from scipy import sparse


def inverse_deconvolution_kernel(reference, order=8):
    """Positive Neumann inverse, retaining the reference's actual DC gain.

    G=c*I-A with A>=0 and rho=sum(A)/c<1. The returned finite kernel is
    (1/c)*sum_{n=0}^{order}(A/c)^n. Its omitted l1 mass is bounded exactly
    by rho**(order+1)/(c*(1-rho)) on the infinite translation-invariant grid.
    No detector/padded-context boundary condition is silently introduced.
    """
    g = np.asarray(reference, float)
    if (g.ndim != 2 or any(n % 2 != 1 for n in g.shape)
            or not np.isfinite(g).all() or order < 0 or int(order) != order):
        raise ValueError("Finite odd-sized kernel and nonnegative integer order required")
    center = tuple(n//2 for n in g.shape)
    c = float(g[center])
    a = -g.copy()
    a[center] = 0.
    if c <= 0 or (a < 0).any():
        raise ValueError("Positive center and nonpositive off-diagonal deconvolution required")
    rho = float(a.sum()/c)
    if rho >= 1:
        raise ValueError("Strict diagonal dominance required")
    radius = tuple(int(order)*(n//2) for n in g.shape)
    shape = tuple(2*r+1 for r in radius)
    h = np.zeros(shape)
    power = np.ones((1, 1))
    for n in range(int(order)+1):
        begin = tuple((s-p)//2 for s, p in zip(shape, power.shape))
        h[begin[0]:begin[0]+power.shape[0], begin[1]:begin[1]+power.shape[1]] += power/c
        if n < order:
            power = convolve2d(power, a/c, mode="full")
    return {"kernel": h, "reference_sum": float(g.sum()), "inverse_sum": float(h.sum()),
        "tail_l1_bound": rho**(int(order)+1)/(c*(1-rho)), "rho": rho,
        "order": int(order), "renormalized": False}


def convolution_matrix(kernel, image_shape):
    """Row-output/column-input convolution with explicitly ZERO exterior.

    This finite matrix is for algebra/tests or a padded conditional context.
    It is not a claim that unobserved sky outside a detector cutout is zero.
    """
    k = np.asarray(kernel, float)
    ny, nx = map(int, image_shape)
    if k.ndim != 2 or min(ny, nx) < 1 or any(s % 2 != 1 for s in k.shape):
        raise ValueError("Odd 2D kernel and positive image shape required")
    rows, cols, values = [], [], []
    yy, xx = np.mgrid[:ny, :nx]
    for j, i in np.argwhere(k != 0):
        sy = yy-(int(j)-k.shape[0]//2)
        sx = xx-(int(i)-k.shape[1]//2)
        valid = (sy >= 0) & (sy < ny) & (sx >= 0) & (sx < nx)
        rows.append((yy[valid]*nx+xx[valid]).ravel())
        cols.append((sy[valid]*nx+sx[valid]).ravel())
        values.append(np.full(valid.sum(), k[j, i]))
    return sparse.csr_matrix((np.concatenate(values),
        (np.concatenate(rows), np.concatenate(cols))), shape=(ny*nx, ny*nx))


def photon_to_raw_dn_jacobian(forward, gain_e_per_dn, inverse_linearity_derivative):
    """d(rawDN)/dQ_e = H*diag((L^-1)'/gain), BEFORE group averaging.

    This order follows the reference provenance: native L was derived after
    IPC deconvolution. Gain belongs to the SOURCE electron coordinate. Any
    SIRS/group-average/native-L/OLS/flat/PHOTMJSR Jacobian acts on the LEFT.
    Do not apply this map again to measured post-SIRS dark/electronic noise.
    """
    h = sparse.csr_matrix(forward)
    gain = np.asarray(gain_e_per_dn, float).ravel()
    deriv = np.asarray(inverse_linearity_derivative, float).ravel()
    if (h.shape[1] != len(gain) or deriv.shape != gain.shape or (gain <= 0).any()
            or not np.isfinite(gain).all() or not np.isfinite(deriv).all()):
        raise ValueError("Finite source gains and inverse-L derivatives required")
    return h @ sparse.diags(deriv/gain)
