from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp

import background
import closed_horndeski_n2 as ch
import growth


Z_START = 2.1
Z_SAMPLES = (2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0)
N_START = math.log(1.0 / (1.0 + Z_START))


def background_state(N: float) -> np.ndarray:
    return np.asarray(growth.R1.sol(N), dtype=float)


def r_and_dr_dN(N: float, y: np.ndarray) -> tuple[float, float]:
    _, v, E = [float(x) for x in y]
    _, dv_dN, dE_dN = background.rhs(N, y)

    theta, sigma = ch.theta_sigma(v, E)
    dtheta_dN = (
        dE_dN
        - 3.0 * v * v * dv_dN / background.LAMBDA
    )
    dsigma_dN = (
        v * dv_dN
        + 12.0 / background.LAMBDA
        * (
            dE_dN * v**3
            + 3.0 * E * v * v * dv_dN
        )
        - 6.0 * E * dE_dN
    )

    r = sigma / theta**2
    dr_dN = (
        dsigma_dN / theta**2
        - 2.0 * sigma * dtheta_dN / theta**3
    )
    return float(r), float(dr_dN)


def dln_gs_dN(N: float, y: np.ndarray) -> float:
    r, dr_dN = r_and_dr_dN(N, y)
    G = 5.0 * (r + 3.0) / (r + 8.0)
    dG_dN = 25.0 * dr_dN / (r + 8.0) ** 2
    return float(dG_dN / G)


def coefficients(N: float) -> dict[str, float]:
    y = background_state(N)
    _, v, E = [float(x) for x in y]
    a = math.exp(N)

    G = float(ch.gs_n2(v, E))
    F = float(ch.formal_fs_n2_gravity_sector(N, y))
    dlnH_dN = float(background.rhs(N, y)[2] / E)
    dlnG_dN = dln_gs_dN(N, y)

    # For manuscript n=2 / literature N=3:
    # -D^2 = 8 K, hence -(D^2+3K) = 5 K.
    shifted_eigenvalue_over_H2 = (
        5.0
        * background.KCURV
        * a**-2
        / (E * E)
    )

    friction = 3.0 + dlnH_dN + dlnG_dN
    omega2 = shifted_eigenvalue_over_H2 * (F / G)

    return {
        "a": a,
        "E": E,
        "G_S_n2": G,
        "F_S_n2_gravity_scalar": F,
        "c_s2_reduced": F / G,
        "dlnH_dN": dlnH_dN,
        "dlnG_dN": dlnG_dN,
        "friction": friction,
        "omega2": omega2,
    }


def system_matrix(N: float) -> np.ndarray:
    c = coefficients(N)
    return np.array(
        [
            [0.0, 1.0],
            [-c["omega2"], -c["friction"]],
        ],
        dtype=float,
    )


def propagator_rhs(N: float, flat_u: np.ndarray) -> np.ndarray:
    U = np.asarray(flat_u, dtype=float).reshape(2, 2)
    return (system_matrix(N) @ U).reshape(-1)


def integrate_propagator():
    sol = solve_ivp(
        propagator_rhs,
        (N_START, 0.0),
        np.eye(2, dtype=float).reshape(-1),
        method="DOP853",
        rtol=2e-11,
        atol=1e-13,
        max_step=0.002,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def transfer_matrix(sol, z: float) -> np.ndarray:
    N = math.log(1.0 / (1.0 + z))
    return np.asarray(sol.sol(N), dtype=float).reshape(2, 2)


def constraint_transfer(sol, z: float) -> np.ndarray:
    N = math.log(1.0 / (1.0 + z))
    y_bg = background_state(N)
    _, v, E = [float(x) for x in y_bg]
    a = math.exp(N)
    U = transfer_matrix(sol, z)

    cols = []
    for j in range(2):
        zeta = float(U[0, j])
        zeta_N = float(U[1, j])
        zeta_dot = E * zeta_N
        delta_n, chi = ch.constraint_solution_n2(
            a,
            v,
            E,
            zeta,
            zeta_dot,
        )
        cols.append(
            np.array([delta_n, chi], dtype=float)
        )
    return np.column_stack(cols)


def constraint_residual_max(sol, ngrid: int = 80) -> float:
    residual = 0.0
    for N in np.linspace(N_START, 0.0, ngrid):
        y_bg = background_state(float(N))
        _, v, E = [float(x) for x in y_bg]
        a = math.exp(float(N))
        U = np.asarray(
            sol.sol(float(N)),
            dtype=float,
        ).reshape(2, 2)

        for j in range(2):
            zeta = float(U[0, j])
            zeta_dot = E * float(U[1, j])
            c1, c2 = ch.constraint_residuals_n2(
                a,
                v,
                E,
                zeta,
                zeta_dot,
            )
            residual = max(
                residual,
                abs(float(c1)),
                abs(float(c2)),
            )
    return residual


def liouville_residual(sol) -> float:
    U0 = transfer_matrix(sol, 0.0)
    det_numeric = float(np.linalg.det(U0))

    integral, _ = quad(
        lambda N: coefficients(float(N))["friction"],
        N_START,
        0.0,
        epsabs=1e-11,
        epsrel=1e-11,
        limit=300,
    )
    det_expected = math.exp(-integral)
    return abs(det_numeric - det_expected) / abs(det_expected)


def stability_scan(ngrid: int = 400) -> dict[str, float]:
    rows = [
        coefficients(float(N))
        for N in np.linspace(N_START, 0.0, ngrid)
    ]
    return {
        "min_G_S_n2": min(r["G_S_n2"] for r in rows),
        "min_F_S_n2_gravity_scalar": min(
            r["F_S_n2_gravity_scalar"] for r in rows
        ),
        "min_c_s2_reduced": min(
            r["c_s2_reduced"] for r in rows
        ),
        "max_c_s2_reduced": max(
            r["c_s2_reduced"] for r in rows
        ),
        "min_omega2": min(r["omega2"] for r in rows),
        "max_omega2": max(r["omega2"] for r in rows),
    }


def diagnostics() -> dict:
    sol = integrate_propagator()
    stab = stability_scan()

    rows = []
    for z in Z_SAMPLES:
        c = coefficients(
            math.log(1.0 / (1.0 + z))
        )
        U = transfer_matrix(sol, z)
        C = constraint_transfer(sol, z)
        rows.append(
            {
                "z": z,
                "G_S_n2": c["G_S_n2"],
                "F_S_n2_gravity_scalar": (
                    c["F_S_n2_gravity_scalar"]
                ),
                "c_s2_reduced": c["c_s2_reduced"],
                "friction_N": c["friction"],
                "omega2_N": c["omega2"],
                "U_z_from_z2p1": U.tolist(),
                "constraint_transfer_delta_n_chi": C.tolist(),
            }
        )

    return {
        "schema": "r1-physical-n2-reduced-propagator-v1",
        "status": (
            "REDUCED_GRAVITY_SCALAR_N2_PROPAGATOR_EVALUATED"
        ),
        "model": "R1",
        "harmonic": {
            "manuscript_n": 2,
            "literature_N": 3,
            "minus_D2_over_K": 8.0,
            "minus_D2_plus_3K_over_K": 5.0,
        },
        "state": [
            "zeta",
            "d_zeta_dN",
        ],
        "anchor_redshift": Z_START,
        "ode": (
            "zeta_NN + "
            "(3+dlnH/dN+dlnG_S/dN) zeta_N + "
            "[5 K/(a^2 H^2)](F_S/G_S) zeta = 0"
        ),
        "stability_scan_z_0_to_2p1": stab,
        "constraint_residual_max_abs": (
            constraint_residual_max(sol)
        ),
        "liouville_relative_residual": (
            liouville_residual(sol)
        ),
        "rows": rows,
        "claim_boundary": {
            "exact": [
                "R1 background",
                "closed-Horndeski physical n=2 G_S",
                "closed-Horndeski physical n=2 lapse/shift constraints",
                "reduced gravity+scalar F_S expression",
                "two-by-two phase-space fundamental matrix",
            ],
            "not_yet_claimed": [
                "complete matter+radiation scalar stability matrix",
                "Newtonian-gauge observable transfer",
                "directional f_sigma8 kernel",
                "lensing kernel",
                "full luminosity-distance line-of-sight kernel",
            ],
        },
    }


def write_tex(report: dict, path: Path) -> None:
    final_U = np.asarray(
        report["rows"][-1]["U_z_from_z2p1"],
        dtype=float,
    )
    s = report["stability_scan_z_0_to_2p1"]

    lines = [
        r"\section{First numerical R1 physical-$n=2$ reduced propagator}",
        r"\label{sec:r1-n2-reduced-propagator}",
        "",
        (
            "The representation-invariant prediction track can now be "
            "evaluated at the reduced gravity--scalar level. "
            "For the physical $n=2$ harmonic, "
            "$-(\\mathcal D^2+3K)=5K$."
        ),
        "",
        "The exact reduced equation integrated here is",
        r"\begin{equation}",
        r"\boxed{",
        r"\zeta_{,NN}",
        r"+\left(3+\frac{d\ln H}{dN}"
        r"+\frac{d\ln\mathcal G_{S,2}}{dN}\right)\zeta_{,N}",
        r"+\frac{5K}{a^2H^2}"
        r"\frac{\mathcal F_{S,2}}{\mathcal G_{S,2}}\zeta=0.",
        r"}",
        r"\end{equation}",
        "",
        (
            "The fundamental matrix is normalized to the identity at "
            "$z=2.1$ in the state "
            "$(\\zeta,d\\zeta/dN)^T$."
        ),
        "",
        r"\begin{equation}",
        r"U(0,2.1)=",
        r"\begin{pmatrix}",
        (
            f"{final_U[0,0]:.9f} & "
            f"{final_U[0,1]:.9f}\\\\"
        ),
        (
            f"{final_U[1,0]:.9f} & "
            f"{final_U[1,1]:.9f}"
        ),
        r"\end{pmatrix}.",
        r"\end{equation}",
        "",
        (
            "Across $0\\le z\\le2.1$ the numerical scan gives "
            f"$\\min\\mathcal G_{{S,2}}={s['min_G_S_n2']:.6e}$ "
            "and "
            f"$\\min\\mathcal F_{{S,2}}={s['min_F_S_n2_gravity_scalar']:.6e}$."
        ),
        (
            "The maximum absolute residual of the exact lapse/shift "
            "constraints over both fundamental solutions is "
            f"${report['constraint_residual_max_abs']:.3e}$, "
            "and the Liouville determinant identity is reproduced with "
            f"relative residual ${report['liouville_relative_residual']:.3e}$."
        ),
        "",
        r"\subsection{Claim boundary}",
        (
            "This is the first numerical phase-space transfer kernel for "
            "the physical $n=2$ R1 mode.  It uses the exact reduced "
            "closed-Horndeski $n=2$ kinetic coefficient and constraints."
        ),
        (
            "The current $\\mathcal F_{S,2}$ implementation is explicitly "
            "the gravity--scalar expression; the complete "
            "matter+radiation scalar stability/response system is not "
            "yet included.  Consequently this result is not yet the "
            "directional $f\\sigma_8$, lensing, or luminosity-distance "
            "observable kernel."
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
