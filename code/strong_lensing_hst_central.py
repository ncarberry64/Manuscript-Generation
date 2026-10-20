"""Pixel tests of HE0435's frozen central-image prediction.

This deliberately separates a smooth-light fit from the much less restrictive
nonnegative-light upper bound.  A macro magnification is not automatically an
observed flux ratio: extinction, source variability, and microlensing remain
explicit physical conditions of the test.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares, lsq_linear, minimize
from scipy.signal import fftconvolve
from scipy.stats import norm


def robust_sigma(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    return float(1.4826022185 * np.median(np.abs(a - np.median(a))))


def empirical_psf(stamp):
    """Background-subtracted, positive, unit-aperture empirical stellar PSF."""
    a = np.asarray(stamp, float)
    yy, xx = np.indices(a.shape)
    cy, cx = (np.array(a.shape) - 1) / 2
    rr = np.hypot(xx - cx, yy - cy)
    bg = float(np.median(a[rr > min(a.shape) * 0.38]))
    z = np.maximum(a - bg, 0)
    core = rr < 3
    cx = float(np.sum(xx[core] * z[core]) / np.sum(z[core]))
    cy = float(np.sum(yy[core] * z[core]) / np.sum(z[core]))
    # Discard noisy far wings in the template consistently for A and the centre.
    z[np.hypot(xx - cx, yy - cy) > min(a.shape) * 0.35] = 0
    return z / z.sum(), (cx, cy), bg


def point_template(shape, xy, psf, psf_center):
    yy, xx = np.indices(shape)
    # Cubic interpolation can ring negative; retaining a positive intensity is
    # essential for the separate positivity inequality.
    p = map_coordinates(psf, [yy - xy[1] + psf_center[1],
                             xx - xy[0] + psf_center[0]], order=3,
                        mode="constant", cval=0, prefilter=True)
    return np.maximum(p, 0)


def fit_point(data, sigma, guess, psf, psf_center, radius_pixels):
    yy, xx = np.indices(data.shape)
    mask = np.hypot(xx - guess[0], yy - guess[1]) < radius_pixels
    x = xx[mask] - guess[0]
    y = yy[mask] - guess[1]
    d = data[mask]
    err = sigma[mask]

    def pred(q):
        p = point_template(data.shape, q[:2], psf, psf_center)[mask]
        return q[2] * p + q[3] + q[4] * x + q[5] * y

    q0 = [*guess, float(np.sum(np.maximum(d - np.median(d), 0))),
          float(np.median(d)), 0, 0]
    lo = [guess[0] - 1.5, guess[1] - 1.5, 0, -np.inf, -np.inf, -np.inf]
    hi = [guess[0] + 1.5, guess[1] + 1.5, np.inf, np.inf, np.inf, np.inf]
    res = least_squares(lambda q: (pred(q) - d) / err, q0,
                        bounds=(lo, hi), xtol=1e-10, ftol=1e-10, gtol=1e-9)
    return {"xy": res.x[:2], "flux": float(res.x[2]),
            "background_plane": res.x[3:],
            "chi2": float(np.sum(res.fun**2)), "npixel": int(mask.sum()),
            "rms_electrons_per_second": float(np.sqrt(np.mean((pred(res.x)-d)**2)))}


def registered_light_center(fitted_xy, sky_xy_arcsec, scale):
    """Fit translation plus small rotation with the declared east/north signs."""
    p = np.asarray(sky_xy_arcsec, float) / scale
    p[:, 0] *= -1  # drizzled north-up arrays: increasing column is west
    obs = np.asarray(fitted_xy, float)
    q = least_squares(lambda z: (p @ np.array([[np.cos(z[2]), -np.sin(z[2])],
                    [np.sin(z[2]), np.cos(z[2])]]).T + z[:2] - obs).ravel(),
                    [*np.mean(obs - p, axis=0), 0])
    pred = p @ np.array([[np.cos(q.x[2]), -np.sin(q.x[2])],
                         [np.sin(q.x[2]), np.cos(q.x[2])]]).T + q.x[:2]
    return q.x[:2], float(q.x[2]), float(np.max(np.linalg.norm(pred-obs, axis=1))*scale)


def sky_offset_to_pixel(offset_arcsec, scale, rotation_rad=0.):
    """Apply the same east/north orientation and rotation used by registration."""
    x,y=np.asarray(offset_arcsec,float)
    co,si=np.cos(rotation_rad),np.sin(rotation_rad)
    return np.array([co*(-x)-si*y,si*(-x)+co*y])/scale


def positivity_upper(data, sigma, template, alpha=0.05, support_fraction=0.03):
    """Simultaneous upper bound allowing arbitrary nonnegative nuisance light.

    For each pixel d_i = f P_i + B_i + e_i, B_i>=0. A union bound on
    Gaussian marginal errors gives f <= min_i (d_i+z sigma_i)/P_i with at
    least 1-alpha coverage, irrespective of pixel correlations. This is a
    conditional photometric bound, not a chi-square detection significance.
    """
    p = np.asarray(template, float)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one")
    use = np.isfinite(data) & np.isfinite(sigma) & (np.asarray(sigma)>0) & np.isfinite(p) & (p > support_fraction*p.max())
    n = int(use.sum())
    if n == 0:
        raise ValueError("Upper bound requires positive template support and finite positive marginal errors")
    z = float(norm.ppf(1-alpha/n))
    values = (np.asarray(data)[use] + z*np.asarray(sigma)[use])/p[use]
    return float(np.min(values)), {"npixel": n, "bonferroni_z": z,
                                  "alpha": alpha,
                                  "support_fraction_of_peak": support_fraction}


def nonnegative_profile(data, sigma, template, flux):
    """Profile any galaxy/host/nucleus intensity independently at every pixel."""
    nuisance = np.maximum(data - flux*template, 0)
    residual = (data - nuisance - flux*template)/sigma
    return float(np.sum(residual**2)), nuisance, residual


def assess_central_image_system(data, sigma, shape, light_center_pixel, scale,
                                positions_arcsec, flux_ratios_relative_A,
                                A_flux, psf, psf_center, centroid_offset=(0,0),
                                normalization_fraction=1.0, rotation_rad=0.):
    """Assess the SUM of every unresolved additional image, with shared shift.

    The returned scale upper limit multiplies the *whole predicted image
    system*. A value below one excludes its nonattenuated flux in this PSF,
    centroid and marginal-noise application. This does not apply a single-image
    flux threshold separately to overlapping images.
    """
    positions=np.asarray(positions_arcsec,float).reshape(-1,2)
    ratios=np.asarray(flux_ratios_relative_A,float)
    if len(positions)!=len(ratios) or np.any(ratios<0):
        raise ValueError("Positions require matching nonnegative absolute flux ratios")
    total=np.zeros(shape,float)
    for xy,r in zip(positions,ratios):
        pixel=np.asarray(light_center_pixel)+sky_offset_to_pixel(xy,scale,rotation_rad)+np.asarray(centroid_offset)/scale
        total+=r*point_template(shape,pixel,psf,psf_center)
    predicted=A_flux*normalization_fraction*total
    if not np.any(predicted>0):
        return {"upper_scale_95":None,"profile_excess_chi2":0.,"n_support_pixel":0},predicted
    upper,details=positivity_upper(data,sigma,predicted)
    use=predicted>.03*predicted.max()
    q,_,_=nonnegative_profile(data[use],sigma[use],predicted[use],1.)
    q0,_,_=nonnegative_profile(data[use],sigma[use],predicted[use],0.)
    return {"upper_scale_95":upper,"profile_excess_chi2":q-q0,
            "profile_chi2_is_not_a_detection_significance":True,
            "n_support_pixel":int(use.sum()),"upper_bound_details":details,
            "summed_macro_flux_relative_A":float(ratios.sum()),
            "number_of_additional_images":int(len(ratios))},predicted


def smooth_light_profile(data, sigma, scale, psf, psf_center, center,
                         central_xy_offset_arcsec, flux, minimum_width=0.08,
                         nuclear_point=False, start=None, rotation_rad=0.):
    """Free-centroid elliptical positive Gaussian-mixture galaxy and nucleus.

    Uses only r<0.65 arcsec. The outer quasar/host system enters as a free
    background plane here, rather than pretending to reconstruct an Einstein
    ring from a central cutout. The separate positivity fit allows every host
    and galaxy pixel to vary independently and is the conservative test.
    """
    yy, xx = np.indices(data.shape)
    use = np.hypot(xx-center[0], yy-center[1])*scale < 0.65
    d, err = data[use], sigma[use]
    widths = np.unique(np.r_[minimum_width, 0.08, 0.15, 0.3, 0.6, 1.2])
    cxy = center + sky_offset_to_pixel(central_xy_offset_arcsec,scale,rotation_rad)
    central = point_template(data.shape, cxy, psf, psf_center)
    # Centred convolution kernel, independent of a star's pixel phase.
    kshape = (25, 25)
    kernel = point_template(kshape, [12, 12], psf, psf_center)
    kernel /= kernel.sum()

    def solve(q, return_model=False):
        cx, cy = center + q[:2]/scale
        co, si = np.cos(q[3]), np.sin(q[3])
        x = ((xx-cx)*co + (yy-cy)*si)*scale
        y = (-(xx-cx)*si + (yy-cy)*co)*scale
        templates = []
        for width in widths:
            g = np.exp(-0.5*((x/width)**2+(y/(width*q[2]))**2))
            templates.append(fftconvolve(g, kernel, mode="same"))
        if nuclear_point:
            templates.append(point_template(data.shape, [cx, cy], psf, psf_center))
        npositive = len(templates)
        templates += [np.ones(data.shape), (xx-center[0])*scale,
                      (yy-center[1])*scale]
        a = np.array([t[use]/err for t in templates]).T
        fit = lsq_linear(a, (d-flux*central[use])/err,
                         bounds=(np.r_[np.zeros(npositive),[-np.inf]*3], np.inf),
                         tol=1e-8, max_iter=80)
        resid = a@fit.x-(d-flux*central[use])/err
        if return_model:
            model = sum(c*t for c,t in zip(fit.x,templates))+flux*central
            return model, resid, fit.x
        return float(np.sum(resid**2))

    q0 = np.asarray(start if start is not None else [0,0,0.8,0.0], float)
    opt = minimize(solve, q0, method="Powell",
                   bounds=[(-0.05,0.05),(-0.05,0.05),(0.45,1),( -np.pi,np.pi)],
                   options={"maxiter":35,"ftol":1e-6,"xtol":2e-4})
    model, resid, amplitudes = solve(opt.x, True)
    return {"chi2": float(np.sum(resid**2)), "npixel": int(use.sum()),
            "parameters": opt.x.tolist(), "success": bool(opt.success),
            "minimum_intrinsic_width_arcsec": minimum_width,
            "nuclear_point_allowed": nuclear_point,
            "amplitudes": amplitudes.tolist(),
            "rms_electrons_per_second": float(np.sqrt(np.mean((model[use]-d)**2)))}, model
