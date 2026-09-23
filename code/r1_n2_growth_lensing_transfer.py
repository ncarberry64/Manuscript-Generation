from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import background
import closed_horndeski_n2 as ch
import growth
import r1_n2_reduced_propagator as rp


Z_SAMPLES = rp.Z_SAMPLES
N_START = rp.N_START


def bg(N: float) -> np.ndarray:
    return np.asarray(growth.R1.sol(N), dtype=float)


def topographic_state(
    prop,
    N: float,
    column: int,
) -> tuple[float, float, float, float, float]:
    """
    Returns (zeta, zeta_N, delta_n, chi, H).
    """
    y_bg = bg(N)
    _, v, H = [float(x) for x in y_bg]
    a = math.exp(N)

    U = np.asarray(
        prop.sol(N),
        dtype=float,
    ).reshape(2, 2)

    zeta = float(U[0, column])
    zeta_N = float(U[1, column])
    zeta_dot = H * zeta_N

    delta_n, chi = ch.constraint_solution_n2(
        a,
        v,
        H,
        zeta,
        zeta_dot,
    )
    return zeta, zeta_N, delta_n, chi, H


def chi_value(prop, N: float, column: int) -> float:
    return topographic_state(prop, N, column)[3]


def dchi_dN(prop, N: float, column: int) -> float:
    """
    Fifth-order-safe local differentiation in the interior, with
    second-order one-sided stencils near the integration endpoints.
    """
    h = 2.0e-5

    if N - 2.0*h >= N_START and N + 2.0*h <= 0.0:
        fm2 = chi_value(prop, N - 2.0*h, column)
        fm1 = chi_value(prop, N - h, column)
        fp1 = chi_value(prop, N + h, column)
        fp2 = chi_value(prop, N + 2.0*h, column)
        return (
            fm2 - 8.0*fm1 + 8.0*fp1 - fp2
        ) / (12.0*h)

    if N + 2.0*h <= 0.0:
        f0 = chi_value(prop, N, column)
        f1 = chi_value(prop, N + h, column)
        f2 = chi_value(prop, N + 2.0*h, column)
        return (-3.0*f0 + 4.0*f1 - f2) / (2.0*h)

    f0 = chi_value(prop, N, column)
    f1 = chi_value(prop, N - h, column)
    f2 = chi_value(prop, N - 2.0*h, column)
    return (3.0*f0 - 4.0*f1 + f2) / (2.0*h)


def bardeen_potentials(
    prop,
    N: float,
    column: int,
) -> dict[str, float]:
    """
    Unitary-gauge scalar metric:
      ds^2=-(1+2 delta_n)dt^2
           +2 D_i chi dt dx^i
           +a^2(1+2 zeta)gamma_ij dx^i dx^j.

    In the De Felice-Koyama-Tsujikawa sign convention:
      Psi = delta_n + dot(chi)
      Phi = zeta + H chi
      Phi_eff = (Psi-Phi)/2

    Thus Psi-Phi is the standard lensing/Weyl combination in their
    convention.
    """
    zeta, _, delta_n, chi, H = topographic_state(
        prop,
        N,
        column,
    )
    chi_dot = H * dchi_dN(prop, N, column)

    psi_gi = delta_n + chi_dot
    phi_gi = zeta + H*chi
    phi_eff = 0.5*(psi_gi - phi_gi)

    return {
        "Psi_DFKT": psi_gi,
        "Phi_DFKT": phi_gi,
        "Phi_eff": phi_eff,
        "Weyl_combination": psi_gi - phi_gi,
        "delta_n": delta_n,
        "chi": chi,
        "chi_dot": chi_dot,
    }


def dust_rhs(
    N: float,
    y: np.ndarray,
    prop,
    column: int,
) -> np.ndarray:
    """
    Gauge-invariant dust response driven directly by the exact n=2
    reduced metric constraints.

    For w=c_m^2=0 in unitary gauge:
      dot(v_m) = -delta_n

      dot(delta_m)+3 dot(H v_m)
        = -(3 dot(zeta) + k^2/a^2 chi)
          + k^2/a^2 v_m

    where delta_m is the gauge-invariant comoving matter density
    contrast of De Felice-Koyama-Tsujikawa.
    """
    delta_m, v_m = [float(x) for x in y]

    y_bg = bg(N)
    _, _, H = [float(x) for x in y_bg]
    a = math.exp(N)

    zeta, zeta_N, delta_n, chi, _ = topographic_state(
        prop,
        N,
        column,
    )
    del zeta

    dH_dN = float(background.rhs(N, y_bg)[2])
    Hdot = H * dH_dN
    zeta_dot = H * zeta_N

    k2_over_a2 = 8.0 * background.KCURV / (a*a)

    v_dot = -delta_n
    delta_dot = (
        -3.0*Hdot*v_m
        + 3.0*H*delta_n
        - 3.0*zeta_dot
        - k2_over_a2*chi
        + k2_over_a2*v_m
    )

    return np.array(
        [
            delta_dot/H,
            v_dot/H,
        ],
        dtype=float,
    )


def integrate_dust_response(prop, column: int):
    sol = solve_ivp(
        lambda N, y: dust_rhs(N, y, prop, column),
        (N_START, 0.0),
        np.zeros(2, dtype=float),
        method="DOP853",
        rtol=2e-11,
        atol=1e-13,
        max_step=0.002,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def dust_transfer_matrix(prop, dust_solutions, z: float) -> np.ndarray:
    N = math.log(1.0/(1.0+z))
    cols = [
        np.asarray(sol.sol(N), dtype=float)
        for sol in dust_solutions
    ]
    return np.column_stack(cols)


def dust_deltaN_kernel(
    prop,
    dust_solutions,
    z: float,
) -> np.ndarray:
    N = math.log(1.0/(1.0+z))
    out = []

    for column, sol in enumerate(dust_solutions):
        state = np.asarray(sol.sol(N), dtype=float)
        rhs = dust_rhs(N, state, prop, column)
        out.append(float(rhs[0]))

    return np.asarray(out, dtype=float)


def reference_growth_prime(z: float) -> float:
    if not hasattr(reference_growth_prime, "_sol"):
        reference_growth_prime._sol = growth.integrate_growth(
            growth.growth_rhs_r1
        )
    N = math.log(1.0/(1.0+z))
    return float(reference_growth_prime._sol.sol(N)[1])


def fractional_fsigma8_kernel(
    prop,
    dust_solutions,
    z: float,
) -> np.ndarray:
    """
    Since f sigma_8 is proportional to dD/dN for fixed primordial
    normalization, this is the directional response per unit initial
    n=2 state, before the single SN calibration selects a state vector.
    """
    Dp = reference_growth_prime(z)
    return dust_deltaN_kernel(
        prop,
        dust_solutions,
        z,
    ) / Dp


def weyl_kernel(prop, z: float) -> np.ndarray:
    N = math.log(1.0/(1.0+z))
    return np.asarray(
        [
            bardeen_potentials(prop, N, j)[
                "Weyl_combination"
            ]
            for j in range(2)
        ],
        dtype=float,
    )


def numerical_derivative_check(
    prop,
    dust_solutions,
    ngrid: int = 24,
) -> float:
    """
    Check integrated dust delta_N against a numerical derivative of
    dense-output delta_m.
    """
    worst = 0.0
    h = 1.0e-5

    for N in np.linspace(N_START + 3*h, -3*h, ngrid):
        for j, sol in enumerate(dust_solutions):
            num = (
                float(sol.sol(N+h)[0])
                - float(sol.sol(N-h)[0])
            )/(2*h)
            rhs = float(
                dust_rhs(
                    float(N),
                    np.asarray(sol.sol(N), dtype=float),
                    prop,
                    j,
                )[0]
            )
            worst = max(worst, abs(num-rhs))

    return worst


def diagnostics() -> dict:
    prop = rp.integrate_propagator()
    dust_solutions = [
        integrate_dust_response(prop, 0),
        integrate_dust_response(prop, 1),
    ]

    rows = []
    for z in Z_SAMPLES:
        N = math.log(1.0/(1.0+z))

        dust = dust_transfer_matrix(
            prop,
            dust_solutions,
            z,
        )
        ddelta_dN = dust_deltaN_kernel(
            prop,
            dust_solutions,
            z,
        )
        f8 = fractional_fsigma8_kernel(
            prop,
            dust_solutions,
            z,
        )
        wk = weyl_kernel(prop, z)

        potentials = [
            bardeen_potentials(prop, N, j)
            for j in range(2)
        ]

        rows.append(
            {
                "z": z,
                "dust_transfer_delta_m_v_m": dust.tolist(),
                "delta_m_dN_kernel": ddelta_dN.tolist(),
                "fractional_fsigma8_kernel_per_unit_initial_n2_state": (
                    f8.tolist()
                ),
                "weyl_kernel_Psi_minus_Phi": wk.tolist(),
                "basis_potentials": potentials,
            }
        )

    return {
        "schema": "r1-physical-n2-growth-lensing-transfer-v1",
        "status": (
            "DIRECT_COVARIANT_MATTER_AND_WEYL_TRANSFER_EVALUATED"
        ),
        "model": "R1",
        "anchor": {
            "z": rp.Z_START,
            "n2_state_basis": [
                "zeta",
                "d_zeta_dN",
            ],
            "dust_initial_response": [
                0.0,
                0.0,
            ],
        },
        "matter_equations": {
            "velocity": "dot(v_m)=-delta_n",
            "continuity": (
                "dot(delta_m)+3 dot(H v_m)="
                "-(3 dot(zeta)+k2/a2 chi)+k2/a2 v_m"
            ),
            "closed_S3_scalar_eigenvalue": "k2/a2=8 K/a2",
        },
        "gauge_invariant_potentials": {
            "Psi_DFKT": "delta_n+dot(chi)",
            "Phi_DFKT": "zeta+H chi",
            "Weyl_combination": "Psi_DFKT-Phi_DFKT",
        },
        "dust_dense_output_derivative_residual_max": (
            numerical_derivative_check(
                prop,
                dust_solutions,
            )
        ),
        "rows": rows,
        "claim_boundary": {
            "derived_now": [
                "dust gauge-invariant response to physical n=2 mode",
                "state-space fractional f_sigma8 response kernel",
                "gauge-invariant local Weyl/lensing-potential kernel",
            ],
            "still_open": [
                "single-calibration state vector from full D_L light cone",
                "line-of-sight convergence/shear kernel",
                "full luminosity-distance kernel",
                "full-redshift H observable kernel",
            ],
            "not_used": (
                "the curvature-defective Newtonian-gauge DAE is not "
                "used in this calculation"
            ),
        },
    }


def write_tex(report: dict, path: Path) -> None:
    row0 = report["rows"][-1]
    f8 = row0[
        "fractional_fsigma8_kernel_per_unit_initial_n2_state"
    ]
    wk = row0["weyl_kernel_Psi_minus_Phi"]

    lines = [
        r"\section{Direct $n=2$ matter-growth and Weyl transfer}",
        r"\label{sec:r1-n2-growth-lensing-transfer}",
        "",
        (
            "The reduced physical $n=2$ propagator can be coupled "
            "directly to minimally coupled matter through covariant "
            "energy--momentum conservation, without invoking the "
            "curvature-defective Newtonian-gauge DAE."
        ),
        "",
        "For pressureless matter in unitary gauge,",
        r"\begin{align}",
        r"\dot v_m&=-\delta n,\\",
        r"\dot\delta_m+3\frac{d}{dt}(Hv_m)"
        r"&=-\left(3\dot\zeta+\frac{k_2^2}{a^2}\chi\right)"
        r"+\frac{k_2^2}{a^2}v_m,",
        r"\end{align}",
        (
            "with $k_2^2/a^2=8K/a^2$ for the physical "
            "manuscript $n=2$ harmonic."
        ),
        "",
        (
            "In the same unitary-gauge convention the gauge-invariant "
            "potentials are"
        ),
        r"\begin{equation}",
        r"\Psi=\delta n+\dot\chi,\qquad"
        r"\Phi=\zeta+H\chi,",
        r"\end{equation}",
        (
            "so the local Weyl/lensing combination in this sign "
            "convention is $\Psi-\Phi$."
        ),
        "",
        (
            "The response is integrated for both basis columns of the "
            "physical state $(\\zeta,d\\zeta/dN)^T$, with zero "
            "matter response at the $z=2.1$ anchor.  Thus the result "
            "is a row-vector kernel on the same two-dimensional "
            "physical state space."
        ),
        "",
        r"\subsection{Present-epoch state-space kernels}",
        r"\begin{align}",
        (
            r"\mathcal F^{\rm state}_{f\sigma_8}(0)"
            f"&=({f8[0]:.9e},\,{f8[1]:.9e}),\\\\"
        ),
        (
            r"\mathcal F^{\rm state}_{\rm Weyl}(0)"
            f"&=({wk[0]:.9e},\,{wk[1]:.9e})."
        ),
        r"\end{align}",
        "",
        (
            "These are not two new amplitudes.  They are the two "
            "components of linear functionals acting on the same "
            "initial $n=2$ state.  The single-calibration theorem "
            "collapses them to one observable prediction once the "
            "luminosity-distance light-cone functional fixes that "
            "state."
        ),
        "",
        r"\subsection{Claim boundary}",
        (
            "The $f\\sigma_8$ object above is a state-space "
            "fractional response kernel, using the R1 baseline "
            "$dD/dN$ normalization.  It is not yet the final "
            "single-number directional prediction because the full "
            "$D_L$ calibration functional has not yet selected the "
            "initial-state vector."
        ),
        (
            "Likewise $\Psi-\Phi$ is the local Weyl potential "
            "combination; the integrated convergence/shear "
            "line-of-sight kernel remains to be evaluated."
        ),
    ]

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--tex", type=Path)
    args = parser.parse_args()

    report = diagnostics()

    if args.artifact is not None:
        args.artifact.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        args.artifact.write_text(
            json.dumps(
                report,
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )

    if args.tex is not None:
        args.tex.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        write_tex(report, args.tex)

    print(
        json.dumps(
            report,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
