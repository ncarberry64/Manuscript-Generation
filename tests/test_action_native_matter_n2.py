from __future__ import annotations

import math
import numpy as np
import pytest

import background as bg
import closed_horndeski_n2 as owner
import action_native_matter_n2 as an


@pytest.fixture(scope="module")
def sol():
    return an.frozen_background()


ZS = [2.1, 1.5, 1.0, 0.5, 0.1, 0.0]


def test_full_hessian_is_symmetric(sol):
    for z in ZS:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        H = an.full_hessian(N, y)

        assert np.max(np.abs(H - H.T)) < 1e-13


def test_unreduced_polynomial_matches_hessian(sol):
    vec = np.array([
        .13, -.07, .11, -.05,
        .17, -.09, .06,
        -.03, .08, -.04,
    ])

    for z in ZS:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        H = an.full_hessian(N, y)

        u = vec[:4]
        xd = vec[4:7]
        x = vec[7:]

        direct = an.unreduced_lagrangian(
            N, y, u, x, xd
        )

        matrix = 0.5 * vec @ H @ vec

        assert abs(direct - matrix) < 2e-13


def test_zero_fluid_limit_recovers_exact_n2_owner(sol):
    for z in ZS:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)
        a = math.exp(N)

        M, _, _ = an.reduced_blocks(
            N,
            y,
            include_fluids=False,
        )

        G = owner.gs_n2(
            y[1],
            y[2],
        )

        assert np.isclose(
            M[0, 0],
            2.0 * G,
            rtol=1e-11,
            atol=1e-13,
        )

        zeta = 0.17
        zetadot = -0.031

        native = an.reconstruct(
            N,
            y,
            np.array([zeta]),
            np.array([zetadot]),
            include_fluids=False,
        )

        exact = np.asarray(
            owner.constraint_solution_n2(
                a,
                y[1],
                y[2],
                zeta,
                zetadot,
            )
        )

        np.testing.assert_allclose(
            native,
            exact,
            rtol=0.0,
            atol=5e-13,
        )


def test_action_constraints_close_on_full_phase_space(sol):
    for z in ZS:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        vel, aux = an.canonical_maps(
            N, y
        )

        for state in np.eye(6):

            r = an.constraint_residuals(
                N,
                y,
                state[:3],
                vel @ state,
                aux @ state,
            )

            assert np.max(np.abs(r)) < 5e-10


def test_schutz_fluid_euler_lagrange_equations_close(sol):
    for z in ZS:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        for state in np.eye(6):

            r = an.fluid_residuals(
                N,
                y,
                state,
            )

            assert np.max(np.abs(r)) < 2e-9


def test_reduced_kinetic_matrix_positive_on_R1(sol):
    minimum = float("inf")

    for N in np.linspace(
        an.N_START,
        0.0,
        401,
    ):
        y = np.asarray(
            sol.sol(N),
            dtype=float,
        )

        M, _, _ = an.reduced_blocks(
            N,
            y,
        )

        eig = np.linalg.eigvalsh(M)

        minimum = min(
            minimum,
            float(np.min(eig)),
        )

    assert minimum > 0.002


def test_no_extra_force_parameters_present():
    import inspect

    src = inspect.getsource(
        an.action_blocks
    )

    forbidden = [
        "fitted_coupling",
        "q2J",
        "fifth_force",
        "extra_force",
    ]

    for term in forbidden:
        assert term not in src
