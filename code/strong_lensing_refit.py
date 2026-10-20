"""Finite inverse fits of the declared three-image encounter.

The injected R1 observables are read from the existing calculation. A fit
changes conventional nuisance parameters, not the injected state or its
physical light preparation. Numerical row scales are not an observing
covariance or an assumed statistical quantum state.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, least_squares

from strong_lensing_comparison import ARCSEC, lens_application
from strong_lensing_local_family import TolmanVIILens


@dataclass(frozen=True)
class ConventionalSheetLens:
    body: TolmanVIILens
    physical_sheet_per_Mpc: float

    @property
    def b(self):
        return self.body.b

    def application(self, x):
        potential, alpha, hessian = self.body.application(x)
        offset = np.asarray(x, dtype=float)-np.asarray(self.body.centre_Mpc)
        scale = self.physical_sheet_per_Mpc
        return (potential+.5*scale*(offset@offset), alpha+scale*offset,
                hessian+scale*np.eye(2))


def model_from_coordinates(coordinates, geometry, mass_Msun, radius_Mpc, source_Mpc):
    """Coordinates: log M, log R, source shift and center in arcsec, kappa.

    The sheet is part of the fitted encounter once; it is not also inserted
    into a propagation map. kappa=q D_ol sheet is the lensing convergence
    for this particular source plane.
    """
    if abs(source_Mpc[1]) > 1e-20 or any(
        not np.allclose(matrix, matrix[0, 0]*np.eye(2), atol=1e-12, rtol=0)
        for matrix in (geometry.D_ol, geometry.D_os, geometry.D_ls)
    ):
        raise ValueError('This inverse application owns the axial scalar-map domain only')
    log_mass, log_radius, source_shift, center, kappa = coordinates
    dol, dos, q = geometry.D_ol[0, 0], geometry.D_os[0, 0], geometry.G[0, 0]
    body = TolmanVIILens(mass_Msun*np.exp(log_mass), radius_Mpc*np.exp(log_radius),
                        (dol*ARCSEC*center, 0.))
    lens = ConventionalSheetLens(body, kappa/(dol*q))
    source = np.asarray(source_Mpc, dtype=float)+np.array([dos*ARCSEC*source_shift, 0.])
    return lens, source


def original_rows(coordinates, geometry, mass_Msun, radius_Mpc, source_Mpc,
                  image_angles, relative_delays_days):
    """Evaluate original lens equations and relative Fermat times at the data."""
    lens, source = model_from_coordinates(coordinates, geometry, mass_Msun, radius_Mpc, source_Mpc)
    applications = [lens_application([theta, 0.], geometry, lens, source) for theta in image_angles]
    image_rows = np.array([row['residual'][0] for row in applications])/ARCSEC
    arrivals = np.array([row['arrival_days'] for row in applications])
    delay_rows = arrivals[1:]-arrivals[0]-relative_delays_days
    return np.r_[image_rows, delay_rows]


def independent_forward_solution(coordinates, geometry, mass_Msun, radius_Mpc,
                                 source_Mpc, image_seeds):
    """Re-solve all image roots instead of reusing their fitted coordinates."""
    lens, source = model_from_coordinates(coordinates, geometry, mass_Msun, radius_Mpc, source_Mpc)
    images = []
    for seed in image_seeds:
        def residual(theta):
            return lens_application([theta, 0.], geometry, lens, source)['residual'][0]
        images.append(brentq(residual, seed-.01*ARCSEC, seed+.01*ARCSEC,
                             xtol=2e-20, rtol=1e-14))
    applications = [lens_application([theta, 0.], geometry, lens, source) for theta in images]
    arrival = np.array([row['arrival_days'] for row in applications])
    return dict(images_rad=np.array(images), relative_delays_days=arrival[1:]-arrival[0],
                max_original_lens_residual_rad=max(abs(row['residual'][0]) for row in applications),
                max_stationarity_Mpc=max(np.max(abs(row['stationarity'])) for row in applications))


def fit_conventional(geometry, mass_Msun, radius_Mpc, source_Mpc, image_angles,
                     relative_delays_days, *, free_center=True, initial=None):
    """Fit the five original nonlinear rows in a bounded local neighborhood."""
    data_angles = np.asarray(image_angles, dtype=float)
    data_delays = np.asarray(relative_delays_days, dtype=float)
    owned = np.arange(5) if free_center else np.array([0, 1, 2, 4])
    scales = np.r_[np.full(3, .001), np.full(2, .1)]

    def unpack(x):
        value = np.zeros(5)
        value[owned] = x
        return value

    def residual(x):
        return original_rows(unpack(x), geometry, mass_Msun, radius_Mpc, source_Mpc,
                             data_angles, data_delays)/scales

    # Absolute steps make derivatives meaningful at the zero initial point.
    def jacobian(x):
        h = 1e-4
        return np.column_stack([(residual(x+np.eye(len(x))[j]*h)
                                 -residual(x-np.eye(len(x))[j]*h))/(2*h)
                                for j in range(len(x))])

    lower = np.array([-.2, -.2, -.02, -.02, -.2])[owned]
    upper = -lower
    fit = least_squares(residual, np.zeros(len(owned)) if initial is None else initial,
                        jac=jacobian, bounds=(lower, upper), xtol=1e-14,
                        ftol=1e-14, gtol=1e-13, max_nfev=100)
    coordinates = unpack(fit.x)
    forward = independent_forward_solution(coordinates, geometry, mass_Msun, radius_Mpc,
                                          source_Mpc, data_angles)
    singular = np.linalg.svd(jacobian(fit.x), compute_uv=False)
    return dict(coordinates=coordinates, optimizer_success=bool(fit.success),
                optimizer_message=fit.message, evaluations=int(fit.nfev),
                original_rows=original_rows(coordinates, geometry, mass_Msun, radius_Mpc,
                                            source_Mpc, data_angles, data_delays),
                scaled_residual_norm=float(np.linalg.norm(fit.fun)),
                numerical_Jacobian_singular_values=singular,
                forward=forward,
                image_errors_microarcsec=(forward['images_rad']-data_angles)/ARCSEC*1e6,
                delay_errors_seconds=(forward['relative_delays_days']-data_delays)*86400)
