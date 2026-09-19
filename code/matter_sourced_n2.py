"""Matter/radiation sourced n=2 closed-Horndeski constraint block.

This module specializes the Newtonian-gauge 00, 0i and traceless-ij equations
of Gambino & Pace (arXiv:2412.01781) to the manuscript cubic-Galileon/KGB
subclass:

    alpha_M = alpha_T = alpha_H = 0
    M^2 = M_Pl^2 = 1

and to the first physical closed scalar harmonic

    manuscript n = 2  <->  literature N = 3,
    D^2 Q = -8 K Q,
    (D^2 + 3K) Q = -5 K Q.

For ideal dust+radiation with vanishing anisotropic stress, the traceless
equation gives Phi = Psi.

Scope:
- exact local sourced-constraint/continuity transfer block;
- no claim yet of a globally integrated evolution system;
- ij-trace/scalar evolution equations still require a convention dictionary
  to the Round-3 Akama-Kobayashi ADM normalization and a Bianchi/
  constraint-preservation audit;
- extra topographic reduced-response operator is OFF here by design.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

import numpy as np

import background
import closed_horndeski_n2 as ch

OMEGA_M0 = 0.31
OMEGA_R0 = 9.0e-5
MPL2 = 1.0


def gp_alpha_quantities(v: float, H: float) -> dict[str, float]:
    """Return the alpha functions in the Gambino-Pace/Gleyzes convention.

    The repository's closed_horndeski_n2.alpha_quantities() stores the
    Bellini-Sawicki braiding alpha_B^BS, for which
        D_kin = alpha_K + (3/2) (alpha_B^BS)^2.

    The curvature-aware Gambino-Pace equations use the convention
        alpha = alpha_K + 6 (alpha_B^GP)^2,

    so for constant Planck mass:
        alpha_B^GP = -alpha_B^BS / 2.
    """
    q_bs = ch.alpha_quantities(v, H)
    return {
        "alpha_B": -0.5 * q_bs["alpha_B"],
        "alpha_K": q_bs["alpha_K"],
        "D_kin": q_bs["D_kin"],
        "alpha_B_BS": q_bs["alpha_B"],
    }


@dataclass(frozen=True)
class SourceState:
    Phi: float
    pi: float
    delta_m: float
    v_m: float
    delta_r: float
    v_r: float


def background_sources(a: float) -> tuple[float, float, float, float]:
    rho_m = 3.0 * OMEGA_M0 * a**-3
    rho_r = 3.0 * OMEGA_R0 * a**-4
    rho_plus_p = rho_m + (4.0 / 3.0) * rho_r
    pressure = rho_r / 3.0
    return rho_m, rho_r, rho_plus_p, pressure


def harmonic_scales(a: float) -> tuple[float, float]:
    kappa = background.KCURV / (a * a)
    k2_phys = 8.0 * kappa
    shifted_laplacian = -5.0 * kappa
    return k2_phys, shifted_laplacian


def hdot_from_background(N: float, y: np.ndarray) -> float:
    dy_dN = background.rhs(N, np.asarray(y, dtype=float))
    H = float(y[2])
    return H * float(dy_dN[2])


def weighted_momentum_source(
    rho_m: float,
    rho_r: float,
    v_m: float,
    v_r: float,
) -> float:
    return rho_m * v_m + (4.0 / 3.0) * rho_r * v_r


def density_source(
    rho_m: float,
    rho_r: float,
    delta_m: float,
    delta_r: float,
) -> float:
    return rho_m * delta_m + rho_r * delta_r


def constraint_matrix(v: float, H: float) -> np.ndarray:
    q = gp_alpha_quantities(v, H)
    alpha_b = q["alpha_B"]
    alpha_k = q["alpha_K"]
    return np.array(
        [
            [
                6.0 * (1.0 + alpha_b) * H,
                (alpha_k - 6.0 * alpha_b) * H * H,
            ],
            [2.0, -2.0 * H * alpha_b],
        ],
        dtype=float,
    )


def determinant_identity(v: float, H: float) -> tuple[float, float]:
    q = gp_alpha_quantities(v, H)
    alpha_b = q["alpha_B"]
    alpha_k = q["alpha_K"]
    dkin = q["D_kin"]

    direct = float(np.linalg.det(constraint_matrix(v, H)))
    analytic = -2.0 * H * H * (
        alpha_k + 6.0 * alpha_b * alpha_b
    )
    via_dkin = -2.0 * H * H * dkin

    if not math.isclose(
        analytic, via_dkin, rel_tol=2e-14, abs_tol=2e-14
    ):
        raise AssertionError("GP/BS alpha_B convention identity failed")

    return direct, analytic


def constraint_rhs(
    N: float,
    y_bg: np.ndarray,
    state: SourceState,
) -> np.ndarray:
    a = math.exp(N)
    _, v, H = [float(x) for x in y_bg]
    Hdot = hdot_from_background(N, y_bg)

    rho_m, rho_r, rho_plus_p, _ = background_sources(a)
    _, shifted = harmonic_scales(a)
    kappa = background.KCURV / (a * a)

    q = gp_alpha_quantities(v, H)
    alpha_b = q["alpha_B"]
    alpha_k = q["alpha_K"]

    delta_rho = density_source(
        rho_m, rho_r, state.delta_m, state.delta_r
    )
    momentum = weighted_momentum_source(
        rho_m, rho_r, state.v_m, state.v_r
    )

    # Gambino-Pace Eq. (8), exact source form:
    # (6 - alpha_K + 12 alpha_B) H^2 Phi
    # - 2(D^2+3K)/a^2 Phi, after Phi=Psi.
    phi_coeff = (
        (6.0 - alpha_k + 12.0 * alpha_b) * H * H
        - 2.0 * shifted
    )
    pi_coeff = 6.0 * H * (
        rho_plus_p / (2.0 * MPL2)
        + Hdot * (1.0 + alpha_b)
        + kappa * (1.0 - alpha_b)
        + (alpha_b / 3.0) * shifted
    )

    rhs_00 = (
        -delta_rho / MPL2
        - phi_coeff * state.Phi
        - pi_coeff * state.pi
    )

    rhs_0i = (
        -momentum / MPL2
        - 2.0 * (1.0 + alpha_b) * H * state.Phi
        - (2.0 * Hdot + rho_plus_p / MPL2) * state.pi
    )

    return np.array([rhs_00, rhs_0i], dtype=float)


def solve_metric_scalar_rates(
    N: float,
    y_bg: np.ndarray,
    state: SourceState,
) -> tuple[float, float]:
    _, v, H = [float(x) for x in y_bg]
    mat = constraint_matrix(v, H)
    rhs = constraint_rhs(N, y_bg, state)
    Phi_dot, pi_dot = np.linalg.solve(mat, rhs)
    return float(Phi_dot), float(pi_dot)


def constraint_residuals(
    N: float,
    y_bg: np.ndarray,
    state: SourceState,
    Phi_dot: float,
    pi_dot: float,
) -> np.ndarray:
    _, v, H = [float(x) for x in y_bg]
    return (
        constraint_matrix(v, H)
        @ np.array([Phi_dot, pi_dot], dtype=float)
        - constraint_rhs(N, y_bg, state)
    )


def no_slip_residual(
    Phi: float,
    Psi: float,
    sigma_total: float = 0.0,
) -> float:
    return Phi - Psi + sigma_total / MPL2


def ideal_fluid_rates(
    N: float,
    y_bg: np.ndarray,
    state: SourceState,
    Phi_dot: float,
) -> dict[str, float]:
    """Covariant conservation for ideal dust and ideal radiation.

    Velocity convention:
        delta T^0_i = (rho+p) D_i v.

    With D^2 Q = -k^2 Q:
      delta_m_dot = (k^2/a^2) v_m + 3 Phi_dot
      v_m_dot     = -Phi

      delta_r_dot = (4/3)(k^2/a^2) v_r + 4 Phi_dot
      v_r_dot     = H v_r - Phi - delta_r/4
    """
    a = math.exp(N)
    H = float(y_bg[2])
    k2_phys, _ = harmonic_scales(a)

    return {
        "delta_m_dot": k2_phys * state.v_m + 3.0 * Phi_dot,
        "v_m_dot": -state.Phi,
        "delta_r_dot": (
            (4.0 / 3.0) * k2_phys * state.v_r
            + 4.0 * Phi_dot
        ),
        "v_r_dot": (
            H * state.v_r
            - state.Phi
            - 0.25 * state.delta_r
        ),
    }


def local_transfer_matrix(
    N: float,
    y_bg: np.ndarray,
) -> np.ndarray:
    basis = np.eye(6, dtype=float)
    cols = []

    for j in range(6):
        s = SourceState(*basis[:, j])
        Phi_dot, pi_dot = solve_metric_scalar_rates(N, y_bg, s)
        f = ideal_fluid_rates(N, y_bg, s, Phi_dot)
        cols.append(
            np.array(
                [
                    Phi_dot,
                    pi_dot,
                    f["delta_m_dot"],
                    f["v_m_dot"],
                    f["delta_r_dot"],
                    f["v_r_dot"],
                ],
                dtype=float,
            )
        )

    return np.column_stack(cols)


def diagnostics(sol, z: float) -> dict:
    N = math.log(1.0 / (1.0 + z))
    y_bg = np.asarray(sol.sol(N), dtype=float)
    _, v, H = [float(x) for x in y_bg]
    q = gp_alpha_quantities(v, H)

    direct, analytic = determinant_identity(v, H)
    T = local_transfer_matrix(N, y_bg)

    return {
        "z": z,
        "alpha_B": q["alpha_B"],
        "alpha_K": q["alpha_K"],
        "D_kin": q["D_kin"],
        "constraint_det_direct": direct,
        "constraint_det_identity": analytic,
        "constraint_condition_number": float(
            np.linalg.cond(constraint_matrix(v, H))
        ),
        "transfer_spectral_radius": float(
            max(abs(np.linalg.eigvals(T)))
        ),
    }


def main() -> None:
    v0 = background.shoot_v0()
    sol = background.integrate(v0)

    payload = {
        "model": "R1",
        "harmonic": "manuscript n=2 / literature N=3",
        "topographic_response": "off",
        "anisotropic_stress": "zero",
        "scope": (
            "local matter/radiation sourced constraint+continuity transfer"
        ),
        "global_evolution_certified": False,
        "rows": [
            diagnostics(sol, z)
            for z in [10.0, 3.0, 2.1, 1.5, 1.0, 0.5, 0.0]
        ],
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
