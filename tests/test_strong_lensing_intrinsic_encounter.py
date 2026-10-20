"""Independent image solves test the physical fixed-impact encounter terms."""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from strong_lensing_comparison import (
    LensGeometry, GalaxyPotential, solve_images, complete_map_response,
    observable_vector, response_vector, MPC_OVER_C_DAYS,
)


def test_simultaneous_intrinsic_source_and_map_variation():
    g = LensGeometry(900*np.eye(2), 1700*np.eye(2), 1000*np.eye(2), .3)
    p = GalaxyPotential(6.7e-6, core_Mpc=5e-5)
    source = np.array([.0013, .00012])
    solution = solve_images(g, p, source)
    dmaps = np.array([.07*g.D_ol, -.03*g.D_os, .02*g.D_ls])
    ds = np.array([.00008, -.00003])
    strength = .13

    def intrinsic(x):
        potential, alpha, _ = p.application(x)
        return strength*potential, strength*alpha

    response = complete_map_response(solution['images_rad'], g, p, source,
        dmaps, d_source_Mpc=ds, intrinsic_encounter=intrinsic)
    h = 2e-5
    plus = solve_images(g.displaced(dmaps, h), replace(p, b=p.b*(1+h*strength)),
        source+h*ds, seeds=solution['images_rad'])
    minus = solve_images(g.displaced(dmaps, -h), replace(p, b=p.b*(1-h*strength)),
        source-h*ds, seeds=solution['images_rad'])
    np.testing.assert_allclose(response_vector(response),
        (observable_vector(plus)-observable_vector(minus))/(2*h), rtol=2e-8, atol=2e-8)
    expected = [-(1+g.z_l)*intrinsic(g.D_ol@theta)[0]*MPC_OVER_C_DAYS
        for theta in solution['images_rad']]
    np.testing.assert_allclose(response['arrival_response_terms_days'][:, 3], expected)
    assert np.linalg.norm(response['arrival_response_terms_days'][:, 3]) > 1


def test_intrinsic_potential_gauge_cancels_only_from_relative_delays():
    g = LensGeometry(900*np.eye(2), 1700*np.eye(2), 1000*np.eye(2), .3)
    p = GalaxyPotential(6.7e-6)
    source = np.array([.0013, 0])
    solution = solve_images(g, p, source)
    response = complete_map_response(solution['images_rad'], g, p, source,
        np.zeros((3, 2, 2)), intrinsic_encounter=lambda x: (1e-10, np.zeros(2)))
    np.testing.assert_array_equal(response['images_arcsec'], 0)
    np.testing.assert_array_equal(response['relative_delays_days'], 0)
    assert np.all(response['arrival_days'] < 0)


@pytest.mark.parametrize('bad', [(np.nan, [0., 0.]), (0., [0.]), (0., [0., np.inf])])
def test_invalid_intrinsic_application_is_rejected(bad):
    g = LensGeometry(900*np.eye(2), 1700*np.eye(2), 1000*np.eye(2), .3)
    p = GalaxyPotential(6.7e-6)
    source = np.array([.0013, 0])
    solution = solve_images(g, p, source)
    with pytest.raises(ValueError, match='intrinsic encounter'):
        complete_map_response(solution['images_rad'], g, p, source,
            np.zeros((3, 2, 2)), intrinsic_encounter=lambda x: bad)
