"""Bounded independent review of the shared GR 2PM thin-plane composition.

This checks the supplied finite-body coefficient, its potential primitive,
and the change of the already frozen R1 comparison when the same GR term is
added to both models. It does not perform metric integration or claim scalar
post-Newtonian or finite-distance cosmological completeness.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "code"))

from strong_lensing_comparison import (
    ARCSEC, MPC_OVER_C_DAYS, LensGeometry, complete_map_response, lens_application,
)
from strong_lensing_local_family import TolmanVIILens
from strong_lensing_r1_screened import R1Coefficients, R1ScreenedLens
from strong_lensing_gr2 import (
    coefficient_I, coefficient_I_prime, coefficient_J,
    TolmanVII2PMCorrection, TolmanVII2PMLens,
)


def solve_inherited_branches(geometry, lens, source, seeds):
    """Use one identical axial root and arrival calculation for all four cases."""
    seeds = np.asarray(seeds, dtype=float)
    images = []
    for seed in seeds:
        residual = lambda t: lens_application(
            np.array([t, 0.]), geometry, lens, source)["residual"][0]
        root = brentq(residual, seed[0]-.001*ARCSEC,
                      seed[0]+.001*ARCSEC, xtol=2e-20, rtol=1e-14)
        images.append([root, 0.])
    images = np.asarray(images)
    applications = [lens_application(t, geometry, lens, source) for t in images]
    arrival = np.array([a["arrival_days"] for a in applications])
    return dict(images_rad=images, images_arcsec=images/ARCSEC,
                arrival_days=arrival, relative_delays_days=arrival[1:]-arrival[0],
                max_lens_residual_rad=max(np.max(abs(a["residual"])) for a in applications))


def change(new, old):
    return dict(images_microarcsec=(new["images_arcsec"]-old["images_arcsec"])*1e6,
                relative_delays_seconds=(new["relative_delays_days"]-
                                        old["relative_delays_days"])*86400)


def encode(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def main(output=None):
    saved = json.loads((REPO / "artifacts/strong_lensing_local_family_20261020/results.json").read_text())
    inputs = saved["physical_inputs"]
    geometry = LensGeometry(*(np.array(inputs["physical_maps_Mpc"][name])
                              for name in ("D_ol", "D_os", "D_ls")), inputs["z_l"])
    source = np.array(inputs["source_transverse_Mpc"])
    seeds = np.array(saved["optical_encounter"]["solution"]["images_rad"])
    previous = json.loads((REPO / "artifacts/strong_lensing_r1_screened_20261020/results.json").read_text())
    coefficients = R1Coefficients(**previous["coefficients"])
    base = TolmanVIILens(inputs["mass_Msun"], inputs["areal_radius_Mpc"])
    correction = TolmanVII2PMCorrection(base)
    corrected = TolmanVII2PMLens(base)
    screened_old = R1ScreenedLens(base, coefficients)
    screened_new = R1ScreenedLens(corrected, coefficients)
    impact = np.linalg.norm(seeds @ geometry.D_ol.T, axis=1)

    # Independent quadrature of I checks the primitive used by arrival times.
    primitive_checks = []
    for y in (0., 1e-5, .12549630368161785, .5, .8, .97, 1., 1.004275925146641, 1.2761280266175101, 3.):
        integrated = quad(coefficient_I, 1., y, epsabs=5e-13,
                          epsrel=5e-13, points=[1.] if y>1 else None)[0]
        primitive_checks.append(dict(y=y, analytic=coefficient_J(y),
                                     integrated=integrated,
                                     absolute_difference=abs(integrated-coefficient_J(y))))
    derivative_checks = []
    for y in (.0001, .12549630368161785, .5, .97, 1.004275925146641, 1.2761280266175101):
        h = min(1e-5, y*1e-3)
        independent = (-coefficient_I(y+2*h)+8*coefficient_I(y+h)
                       -8*coefficient_I(y-h)+coefficient_I(y-2*h))/(12*h)
        derivative_checks.append(dict(y=y, analytic=coefficient_I_prime(y),
                                      finite_difference=independent,
                                      relative_difference=abs(independent-coefficient_I_prime(y))/max(abs(independent),1e-20)))

    models = dict(one_pm=base, one_pm_plus_r1=screened_old,
                  two_pm=corrected, two_pm_plus_r1=screened_new)
    solutions = {name: solve_inherited_branches(geometry, lens, source, seeds)
                 for name, lens in models.items()}
    original_signal = change(solutions["one_pm_plus_r1"], solutions["one_pm"])
    corrected_signal = change(solutions["two_pm_plus_r1"], solutions["two_pm"])
    drift = {key: corrected_signal[key]-original_signal[key] for key in original_signal}

    zero = np.zeros((3, 2, 2))
    linear_gr = complete_map_response(solutions["one_pm"]["images_rad"], geometry,
        base, source, zero, intrinsic_encounter=lambda x: correction.extra_application(x)[:2])
    linear_r1 = complete_map_response(solutions["two_pm"]["images_rad"], geometry,
        corrected, source, zero, intrinsic_encounter=screened_new.intrinsic_extra)
    gr_change = change(solutions["two_pm"], solutions["one_pm"])
    # Independently predict the four-case difference from the mixed derivative
    # of the stationary arrival action, without re-solving perturbed roots.
    mixed_arrival = []
    for theta, gr_motion in zip(solutions["one_pm"]["images_rad"], linear_gr["images_arcsec"]):
        extra_alpha = screened_old.extra_application(geometry.D_ol @ theta)[1]
        mixed_arrival.append(-(1+geometry.z_l)*extra_alpha @ geometry.D_ol @
                             (gr_motion*ARCSEC)*MPC_OVER_C_DAYS*86400)
    mixed_arrival = np.asarray(mixed_arrival)
    mixed_delays = mixed_arrival[1:]-mixed_arrival[0]

    result = dict(
        scope="Shared ordinary-GR 2PM robustness test inside the same isolated-scattering thin-plane convention; no scalar PN completion.",
        compactness=correction.epsilon,
        checks=dict(
            primitive=primitive_checks, derivative=derivative_checks,
            small_y_I_over_y=coefficient_I(1e-7)/1e-7,
            expected_small_y_limit=-6/7,
            surface_I=coefficient_I(1.), expected_surface_I=15*np.pi/8,
            surface_I_prime=coefficient_I_prime(1.), expected_surface_I_prime=-15*np.pi/4,
            central_correction=correction.extra_application([0.,0.]),
            exterior_coefficients=[dict(impact_Mpc=float(b), computed=correction.extra_alpha(b),
                expected=15*np.pi/4*(base.mass_length_Mpc/b)**2) for b in impact if b>base.radius_Mpc]),
        old_image_impact_Mpc=impact,
        gr_extra_alpha_rad=np.array([correction.extra_alpha(b) for b in impact]),
        gr_extra_potential_Mpc=np.array([correction.extra_potential(b) for b in impact]),
        solutions=solutions,
        common_gr_changes=gr_change,
        original_screened_difference=original_signal,
        corrected_screened_difference=corrected_signal,
        screened_difference_drift=drift,
        screened_delay_difference_fractional_drift=drift["relative_delays_seconds"]/original_signal["relative_delays_seconds"],
        mixed_stationary_prediction_of_delay_drift_seconds=mixed_delays,
        four_case_minus_mixed_stationary_prediction_seconds=drift["relative_delays_seconds"]-mixed_delays,
        first_order_common_gr_delays_seconds=linear_gr["relative_delays_days"]*86400,
        first_order_screened_delays_on_corrected_base_seconds=linear_r1["relative_delays_days"]*86400,
        common_gr_finite_minus_first_order_seconds=gr_change["relative_delays_seconds"]-linear_gr["relative_delays_days"]*86400,
        screened_finite_minus_first_order_seconds=corrected_signal["relative_delays_seconds"]-linear_r1["relative_delays_days"]*86400,
        qualification="Tiny numerical drift is conditional on this shared GR correction; it is not a bound on omitted scalar-GR cross terms, local time evolution, finite-distance matching, or profile degeneracies.")
    output = output or REPO / "artifacts/strong_lensing_r1_screened_20261020/review/common_gr_composition_review_results.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, default=encode, indent=2)+"\n", encoding='utf-8', newline='\n')
    print(json.dumps({key:result[key] for key in (
        "gr_extra_alpha_rad", "gr_extra_potential_Mpc", "common_gr_changes",
        "original_screened_difference", "corrected_screened_difference",
        "screened_difference_drift", "screened_delay_difference_fractional_drift",
        "common_gr_finite_minus_first_order_seconds", "screened_finite_minus_first_order_seconds")}, default=encode, indent=2))
    print("max primitive discrepancy", max(v["absolute_difference"] for v in primitive_checks))
    print("max derivative relative discrepancy", max(v["relative_difference"] for v in derivative_checks))
    print("output", output)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, help='Write a separate replay receipt; default preserves the original command behavior.')
    args = parser.parse_args()
    main(args.output)
