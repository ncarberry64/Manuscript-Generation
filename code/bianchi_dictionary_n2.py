"""Round-5 ADM/EFT dictionary and differential-algebraic consistency audit.

Key correction:
The 00 and 0i equations are constraints. They must not be treated as the
complete evolution system and then expected to make both the ij-trace and
scalar equations vanish for an arbitrary pointwise state.

Correct procedure:
1. Use the ij-trace and scalar-pi equations to solve for
   (Phi_ddot, pi_ddot).
2. Evolve Phi, pi and the ideal-fluid perturbations.
3. Impose the 00 and 0i constraints on the initial state.
4. Monitor preservation of those constraints.

Conventions:
- Gambino-Pace/Gleyzes braiding normalization is used inside all GP equations.
- alpha_B^GP = -alpha_B^BS / 2 for this constant-Planck-mass branch.
- topographic reduced response is OFF throughout this audit.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.integrate import solve_ivp

import background
import matter_sourced_n2 as ms


def adm_to_newtonian(
    delta_n: float,
    chi: float,
    chi_dot: float,
    zeta: float,
    H: float,
) -> tuple[float, float, float]:
    pi = -chi
    Phi = delta_n + chi_dot
    Psi = -zeta - H * chi
    return Phi, Psi, pi


def newtonian_to_adm(
    Phi: float,
    Psi: float,
    pi: float,
    pi_dot: float,
    H: float,
) -> tuple[float, float, float, float]:
    chi = -pi
    chi_dot = -pi_dot
    delta_n = Phi + pi_dot
    zeta = H * pi - Psi
    return delta_n, chi, chi_dot, zeta


def _bg(sol, N: float) -> np.ndarray:
    return np.asarray(sol.sol(N), dtype=float)


def _central_N(func, N: float, h: float = 2.0e-5) -> float:
    return float((func(N + h) - func(N - h)) / (2.0 * h))


def _H(sol, N: float) -> float:
    return float(_bg(sol, N)[2])


def _Hdot(sol, N: float) -> float:
    y = _bg(sol, N)
    H = float(y[2])
    return H * float(background.rhs(N, y)[2])


def _aB(sol, N: float) -> float:
    _, v, H = [float(x) for x in _bg(sol, N)]
    return ms.gp_alpha_quantities(v, H)["alpha_B"]


def _aK(sol, N: float) -> float:
    _, v, H = [float(x) for x in _bg(sol, N)]
    return ms.gp_alpha_quantities(v, H)["alpha_K"]


def background_derivatives(sol, N: float) -> dict[str, float]:
    H = _H(sol, N)
    Hdot = _Hdot(sol, N)
    aB = _aB(sol, N)
    aK = _aK(sol, N)

    aBdot = H * _central_N(lambda n: _aB(sol, n), N)
    aKdot = H * _central_N(lambda n: _aK(sol, n), N)

    HaBdot = H * _central_N(
        lambda n: _H(sol, n) * _aB(sol, n), N
    )
    HaKdot = H * _central_N(
        lambda n: _H(sol, n) * _aK(sol, n), N
    )
    Hddot = H * _central_N(lambda n: _Hdot(sol, n), N)

    # Exact GP Eq. (12): H * d/dt(H * dot(alpha_B)).
    H_aBdot_dot = H * _central_N(
        lambda n: (
            _H(sol, n)
            * (
                _H(sol, n)
                * _central_N(lambda nn: _aB(sol, nn), n)
            )
        ),
        N,
    )

    return {
        "H": H,
        "Hdot": Hdot,
        "Hddot": Hddot,
        "aB": aB,
        "aK": aK,
        "aBdot": aBdot,
        "aKdot": aKdot,
        "HaBdot": HaBdot,
        "HaKdot": HaKdot,
        "H_aBdot_dot": H_aBdot_dot,
    }


def acceleration_matrix(sol, N: float) -> np.ndarray:
    bd = background_derivatives(sol, N)
    H = bd["H"]
    aB = bd["aB"]
    aK = bd["aK"]

    # Rows: GP ij-trace Eq. (10), scalar Eq. (12).
    return np.array(
        [
            [2.0, -2.0 * H * aB],
            [6.0 * H * aB, H * H * aK],
        ],
        dtype=float,
    )


def acceleration_determinant_identity(
    sol, N: float
) -> tuple[float, float]:
    y = _bg(sol, N)
    _, v, H = [float(x) for x in y]
    q = ms.gp_alpha_quantities(v, H)
    direct = float(np.linalg.det(acceleration_matrix(sol, N)))
    analytic = 2.0 * H * H * q["D_kin"]
    return direct, analytic


def nonacceleration_terms(
    sol,
    N: float,
    y8: np.ndarray,
) -> tuple[float, float]:
    # y8 = [Phi, pi, Phi_dot, pi_dot, delta_m, v_m, delta_r, v_r]
    Phi, pi, u, w, delta_m, v_m, delta_r, v_r = [
        float(x) for x in y8
    ]

    bg = _bg(sol, N)
    a = math.exp(N)
    bd = background_derivatives(sol, N)

    H = bd["H"]
    Hdot = bd["Hdot"]
    Hddot = bd["Hddot"]
    aB = bd["aB"]
    aK = bd["aK"]
    aBdot = bd["aBdot"]
    aKdot = bd["aKdot"]
    HaBdot = bd["HaBdot"]
    HaKdot = bd["HaKdot"]
    H_aBdot_dot = bd["H_aBdot_dot"]

    rho_m, rho_r, rho_plus_p, _ = ms.background_sources(a)
    delta_p = rho_r * delta_r / 3.0
    p_dot = -4.0 * H * rho_r / 3.0
    kappa = background.KCURV / (a * a)
    shifted = -5.0 * kappa

    # Exact GP Eq. (10), alpha_M=0, Phi=Psi, sigma=0,
    # excluding the two acceleration terms.
    trace_rest = (
        6.0 * H * u
        + 2.0 * (1.0 + aB) * H * u
        + 2.0
        * (
            Hdot
            + 0.5 * rho_plus_p
            - HaBdot
            - 3.0 * aB * H * H
        )
        * w
        + 2.0
        * (
            Hdot
            - 0.5 * rho_plus_p
            + HaBdot
            + 3.0 * (1.0 + aB) * H * H
        )
        * Phi
        + 2.0
        * (
            3.0 * H * Hdot
            + 0.5 * p_dot
            + Hddot
        )
        * pi
        - 6.0 * kappa * (w + H * pi)
        - delta_p
    )

    # Exact GP Eq. (12), alpha_M=alpha_T=alpha_H=0, Phi=Psi,
    # excluding H^2 alpha_K pi_ddot + 6 H alpha_B Phi_ddot.
    scalar_rest = (
        H * H * (6.0 * aB - aK) * u
        + (
            ((3.0 * H * H + Hdot) * aK + HaKdot)
            * H
            * w
        )
        + 2.0
        * shifted
        * (
            Hdot
            + 0.5 * rho_plus_p
            + H * H * aB
            + HaBdot
        )
        * pi
        + 6.0
        * (
            Hdot * (Hdot + 0.5 * rho_plus_p)
            + H * aBdot * (3.0 * H * H + Hdot)
            + H * H_aBdot_dot
        )
        * pi
        + 6.0
        * (
            Hdot
            + 0.5 * rho_plus_p
            + 3.0 * H * H * aB
            + HaBdot
        )
        * u
        - 6.0
        * kappa
        * (H * H * aB + HaBdot)
        * pi
        + (
            6.0 * (Hdot + 0.5 * rho_plus_p)
            + 3.0 * H * H * (6.0 * aB - aK)
            + 2.0 * (9.0 * aB - aK) * Hdot
            + H * (6.0 * aBdot - aKdot)
        )
        * H
        * Phi
        - 2.0 * shifted * H * aB * Phi
    )

    return float(trace_rest), float(scalar_rest)


def solve_accelerations(
    sol,
    N: float,
    y8: np.ndarray,
) -> tuple[float, float]:
    tr, sc = nonacceleration_terms(sol, N, y8)
    rhs = -np.array([tr, sc], dtype=float)
    acc = np.linalg.solve(acceleration_matrix(sol, N), rhs)
    return float(acc[0]), float(acc[1])


def dae_rhs_N(sol, N: float, y8: np.ndarray) -> np.ndarray:
    Phi, pi, u, w, delta_m, v_m, delta_r, v_r = [
        float(x) for x in y8
    ]
    bg = _bg(sol, N)
    H = float(bg[2])

    Phi_ddot, pi_ddot = solve_accelerations(sol, N, y8)

    state = ms.SourceState(
        Phi=Phi,
        pi=pi,
        delta_m=delta_m,
        v_m=v_m,
        delta_r=delta_r,
        v_r=v_r,
    )
    fluids = ms.ideal_fluid_rates(N, bg, state, u)

    return np.array(
        [
            u / H,
            w / H,
            Phi_ddot / H,
            pi_ddot / H,
            fluids["delta_m_dot"] / H,
            fluids["v_m_dot"] / H,
            fluids["delta_r_dot"] / H,
            fluids["v_r_dot"] / H,
        ],
        dtype=float,
    )


def constraint_vector(
    sol,
    N: float,
    y8: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    Phi, pi, u, w, delta_m, v_m, delta_r, v_r = [
        float(x) for x in y8
    ]
    bg = _bg(sol, N)
    state = ms.SourceState(
        Phi=Phi,
        pi=pi,
        delta_m=delta_m,
        v_m=v_m,
        delta_r=delta_r,
        v_r=v_r,
    )
    _, vbg, H = [float(x) for x in bg]
    A = ms.constraint_matrix(vbg, H)
    rhs = ms.constraint_rhs(N, bg, state)
    rates = np.array([u, w], dtype=float)
    return A @ rates - rhs, np.array(
        [
            abs(A[0, 0] * u)
            + abs(A[0, 1] * w)
            + abs(rhs[0]),
            abs(A[1, 0] * u)
            + abs(A[1, 1] * w)
            + abs(rhs[1]),
        ],
        dtype=float,
    )


def normalized_constraint_residuals(
    sol,
    N: float,
    y8: np.ndarray,
) -> np.ndarray:
    residual, scale = constraint_vector(sol, N, y8)
    return np.abs(residual) / np.maximum(scale, 1.0e-30)


def constraint_operator(sol, N: float) -> np.ndarray:
    """Return B(N) such that C_00,0i = B(N) y for the 8-state DAE."""
    cols = []
    for j in range(8):
        e = np.zeros(8, dtype=float)
        e[j] = 1.0
        residual, _ = constraint_vector(sol, N, e)
        cols.append(residual)
    return np.column_stack(cols)


def dae_matrix_N(sol, N: float) -> np.ndarray:
    """Return F(N) for y_N = F(N)y."""
    cols = []
    for j in range(8):
        e = np.zeros(8, dtype=float)
        e[j] = 1.0
        cols.append(dae_rhs_N(sol, N, e))
    return np.column_stack(cols)


def constraint_tangency_operator(
    sol,
    N: float,
    h: float = 2.0e-5,
) -> np.ndarray:
    """Return D(N) such that dC/dN = D(N)y on the DAE flow."""
    B0 = constraint_operator(sol, N)
    Bp = constraint_operator(sol, N + h)
    Bm = constraint_operator(sol, N - h)
    BN = (Bp - Bm) / (2.0 * h)
    F = dae_matrix_N(sol, N)
    return BN + B0 @ F


def constraint_consistent_initial_state(
    sol,
    z_start: float = 2.1,
) -> tuple[float, np.ndarray]:
    """Construct index-consistent initial data.

    The seed must satisfy both the constraints and their first tangency
    conditions. The extra relations select one small adiabatic diagnostic
    mode; they are not observational calibration.
    """
    N = math.log(1.0 / (1.0 + z_start))
    B = constraint_operator(sol, N)
    D = constraint_tangency_operator(sol, N)

    rows = [B[0], B[1], D[0], D[1]]

    adiabatic_density = np.zeros(8, dtype=float)
    adiabatic_density[6] = 1.0
    adiabatic_density[4] = -4.0 / 3.0
    rows.append(adiabatic_density)

    adiabatic_velocity = np.zeros(8, dtype=float)
    adiabatic_velocity[7] = 1.0
    adiabatic_velocity[5] = -1.0
    rows.append(adiabatic_velocity)

    pi_row = np.zeros(8, dtype=float)
    pi_row[1] = 1.0
    rows.append(pi_row)

    phi_row = np.zeros(8, dtype=float)
    phi_row[0] = 1.0
    rows.append(phi_row)

    M = np.vstack(rows)
    rhs = np.zeros(8, dtype=float)
    rhs[-1] = -1.0e-5

    y0, *_ = np.linalg.lstsq(M, rhs, rcond=None)

    mismatch = M @ y0 - rhs
    if np.max(np.abs(mismatch)) > 5.0e-10:
        raise RuntimeError(
            "Could not construct index-consistent initial data: "
            f"max mismatch={np.max(np.abs(mismatch))}"
        )

    c0 = normalized_constraint_residuals(sol, N, y0)
    tangent0 = D @ y0
    tangent_scale = max(1.0e-30, np.linalg.norm(y0))
    tangent_norm = float(np.linalg.norm(tangent0) / tangent_scale)

    if np.max(c0) > 2.0e-10:
        raise RuntimeError(f"Initial constraints not closed: {c0}")
    if tangent_norm > 2.0e-7:
        raise RuntimeError(
            f"Initial constraint tangency not closed: {tangent_norm}"
        )

    return N, np.asarray(y0, dtype=float)

def integrate_dae(
    sol,
    z_start: float = 2.1,
):
    N0, y0 = constraint_consistent_initial_state(sol, z_start)
    result = solve_ivp(
        lambda N, y: dae_rhs_N(sol, N, y),
        (N0, 0.0),
        y0,
        method="DOP853",
        rtol=2.0e-10,
        atol=2.0e-13,
        max_step=0.02,
        dense_output=True,
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result


def constraint_preservation_report(
    sol,
    pert_sol,
) -> dict:
    zs = [2.1, 1.8, 1.5, 1.2, 1.0, 0.8, 0.5, 0.3, 0.0]
    rows = []
    maxima = np.zeros(2, dtype=float)

    for z in zs:
        N = math.log(1.0 / (1.0 + z))
        y = np.asarray(pert_sol.sol(N), dtype=float)
        rn = normalized_constraint_residuals(sol, N, y)
        maxima = np.maximum(maxima, rn)
        rows.append(
            {
                "z": z,
                "constraint_00_normalized": float(rn[0]),
                "constraint_0i_normalized": float(rn[1]),
            }
        )

    return {
        "rows": rows,
        "max_constraint_00_normalized": float(maxima[0]),
        "max_constraint_0i_normalized": float(maxima[1]),
    }


def main() -> None:
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    det_rows = []
    for z in [2.1, 1.5, 1.0, 0.5, 0.0]:
        N = math.log(1.0 / (1.0 + z))
        direct, analytic = acceleration_determinant_identity(sol, N)
        det_rows.append(
            {
                "z": z,
                "acceleration_det_direct": direct,
                "acceleration_det_identity": analytic,
            }
        )

    pert = integrate_dae(sol)
    report = constraint_preservation_report(sol, pert)

    print(
        json.dumps(
            {
                "scope": "DAE evolution + 00/0i constraint preservation",
                "topographic_response": "off",
                "acceleration_determinant": det_rows,
                **report,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
