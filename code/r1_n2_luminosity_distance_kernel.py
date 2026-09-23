from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

import background
import growth
import r1_n2_growth_lensing_transfer as gl
import r1_n2_reduced_propagator as rp


K = background.KCURV
SQRT_K = math.sqrt(K)

CAL_Z = 0.02
A_MU = -0.0412
SIGMA_A_MU = 0.009

Z_ROWS = (
    0.01,
    0.02,
    0.03,
    0.10,
    0.30,
    0.50,
    0.80,
    1.00,
    1.50,
    2.10,
)


def N_of_z(z: float) -> float:
    return math.log(1.0/(1.0+z))


def E_of_z(z: float) -> float:
    return float(growth.R1.sol(N_of_z(z))[2])


def integrate_background_lightcone():
    sol = solve_ivp(
        lambda z, y: np.array([1.0/E_of_z(z)]),
        (0.0, rp.Z_START),
        np.array([0.0]),
        method="DOP853",
        rtol=2e-12,
        atol=1e-14,
        max_step=0.002,
        dense_output=True,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def S_K(chi: float) -> float:
    return math.sin(SQRT_K*chi)/SQRT_K


def C_K(chi: float) -> float:
    return math.cos(SQRT_K*chi)


def dipole_radial(chi: float) -> float:
    """
    Pure A_4i manuscript-n=2 observer-shell dipole per unit mu:
        Q_2,dip = sin(2 psi) mu,
        psi = sqrt(K) chi.
    """
    return math.sin(2.0*SQRT_K*chi)


def dipole_radial_prime(chi: float) -> float:
    return 2.0*SQRT_K*math.cos(2.0*SQRT_K*chi)


def potential_rows(prop, N: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    psi = []
    phi = []
    chi_shift = []

    for j in range(2):
        p = gl.bardeen_potentials(prop, N, j)
        psi.append(float(p["Psi_DFKT"]))
        phi.append(float(p["Phi_DFKT"]))
        chi_shift.append(float(p["chi"]))

    return (
        np.asarray(psi, dtype=float),
        np.asarray(phi, dtype=float),
        np.asarray(chi_shift, dtype=float),
    )


def W_row(prop, N: float) -> np.ndarray:
    psi, phi, _ = potential_rows(prop, N)
    return psi - phi


def dW_dN_row(prop, N: float) -> np.ndarray:
    h = 2.0e-5
    lo = rp.N_START
    hi = 0.0

    if N - 2*h >= lo and N + 2*h <= hi:
        return (
            W_row(prop, N - 2*h)
            - 8.0*W_row(prop, N - h)
            + 8.0*W_row(prop, N + h)
            - W_row(prop, N + 2*h)
        )/(12.0*h)

    if N + 2*h <= hi:
        return (
            -3.0*W_row(prop, N)
            + 4.0*W_row(prop, N + h)
            - W_row(prop, N + 2*h)
        )/(2.0*h)

    return (
        3.0*W_row(prop, N)
        - 4.0*W_row(prop, N - h)
        + W_row(prop, N - 2*h)
    )/(2.0*h)


def newtonian_velocity_row(
    prop,
    dust_solutions,
    N: float,
) -> np.ndarray:
    """
    The unitary-gauge matter potential transforms to Newtonian gauge as
        v_N = v_unitary + chi
    for the same time shift that yields
        Psi = delta_n + dot(chi),
        Phi = zeta + H chi.
    """
    _, _, chi_shift = potential_rows(prop, N)
    v_unitary = np.asarray(
        [
            float(sol.sol(N)[1])
            for sol in dust_solutions
        ],
        dtype=float,
    )
    return v_unitary + chi_shift


def lightcone_rows(
    prop,
    dust_solutions,
    chi_sol,
    z_s: float,
) -> dict[str, np.ndarray | float]:
    if not (0.0 < z_s <= rp.Z_START):
        raise ValueError(
            f"z_s must lie in (0,{rp.Z_START}]"
        )

    N_s = N_of_z(z_s)
    a_s = 1.0/(1.0+z_s)
    H_s = E_of_z(z_s)
    Hconf_s = a_s*H_s
    chi_s = float(chi_sol.sol(z_s)[0])

    S_s = S_K(chi_s)
    if abs(S_s) < 1e-14:
        raise ValueError("source distance too close to zero")

    q_s = dipole_radial(chi_s)
    psi_s, phi_s, _ = potential_rows(prop, N_s)

    V_o = (
        newtonian_velocity_row(
            prop,
            dust_solutions,
            0.0,
        )
        * dipole_radial_prime(0.0)
    )
    V_s = (
        newtonian_velocity_row(
            prop,
            dust_solutions,
            N_s,
        )
        * dipole_radial_prime(chi_s)
    )

    def isw_integrand(z: float) -> np.ndarray:
        N = N_of_z(z)
        chi = float(chi_sol.sol(z)[0])
        # prime is conformal-time derivative at fixed position:
        # W' dchi = a H W_N * dz/H = a W_N dz.
        return (
            dW_dN_row(prop, N)
            * dipole_radial(chi)
            / (1.0+z)
        )

    I_isw, _ = quad_vec(
        isw_integrand,
        0.0,
        z_s,
        epsabs=2e-9,
        epsrel=2e-9,
    )

    delta_z = (
        -V_o
        - psi_s*q_s
        + V_s
        - I_isw
    )

    def radial_integrand(z: float) -> np.ndarray:
        N = N_of_z(z)
        chi = float(chi_sol.sol(z)[0])
        return (
            W_row(prop, N)
            * dipole_radial(chi)
            / E_of_z(z)
        )

    I_rad, _ = quad_vec(
        radial_integrand,
        0.0,
        z_s,
        epsabs=2e-9,
        epsrel=2e-9,
    )

    delta_chi = -delta_z/Hconf_s + I_rad

    def lensing_integrand(z: float) -> np.ndarray:
        N = N_of_z(z)
        chi = float(chi_sol.sol(z)[0])

        if chi < 1e-10:
            # Q_dip ~ 2 sqrt(K) chi and S_K(chi) ~ chi.
            geometric_times_q = (
                2.0*SQRT_K
                * S_K(chi_s-chi)
                / S_s
            )
        else:
            geometric_times_q = (
                S_K(chi_s-chi)
                / (S_s*S_K(chi))
                * dipole_radial(chi)
            )

        # hat-Laplacian(mu) = -2 mu.
        # Hence hat-Laplacian[(Psi-Phi)/2]
        # = -(Psi-Phi) for the sky dipole.
        return (
            -geometric_times_q
            * W_row(prop, N)
            / E_of_z(z)
        )

    I_lens, _ = quad_vec(
        lensing_integrand,
        0.0,
        z_s,
        epsabs=2e-9,
        epsrel=2e-9,
    )

    kappa = -V_o + I_lens

    cot_K = C_K(chi_s)/S_s

    # Curved-FLRW extension of the gauge-invariant geometric result
    # delta D_L = phi_s + delta z
    #             + (S_K'/S_K) delta chi - kappa.
    F_DL = (
        phi_s*q_s
        + delta_z
        + cot_K*delta_chi
        - kappa
    )

    return {
        "z": z_s,
        "chi": chi_s,
        "S_K": S_s,
        "dipole_radial": q_s,
        "observer_velocity": V_o,
        "source_velocity": V_s,
        "isw": np.asarray(I_isw, dtype=float),
        "radial_integral": np.asarray(I_rad, dtype=float),
        "delta_z": np.asarray(delta_z, dtype=float),
        "delta_chi": np.asarray(delta_chi, dtype=float),
        "kappa": np.asarray(kappa, dtype=float),
        "phi_source": np.asarray(phi_s*q_s, dtype=float),
        "F_DL": np.asarray(F_DL, dtype=float),
        "F_DA": np.asarray(F_DL, dtype=float),
    }


def calibration_fraction() -> tuple[float, float]:
    f = math.log(10.0)/5.0
    return f*A_MU, f*SIGMA_A_MU


def unit_null(row: np.ndarray) -> np.ndarray:
    row = np.asarray(row, dtype=float)
    n = np.array([-row[1], row[0]], dtype=float)
    norm = float(np.linalg.norm(n))
    if norm == 0.0:
        raise RuntimeError("calibration row is identically zero")
    return n/norm


def row_det(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0]*b[1] - a[1]*b[0])


def diagnostics() -> dict:
    prop = rp.integrate_propagator()
    dust = [
        gl.integrate_dust_response(prop, 0),
        gl.integrate_dust_response(prop, 1),
    ]
    chi_sol = integrate_background_lightcone()

    rows = []
    for z in Z_ROWS:
        lc = lightcone_rows(
            prop,
            dust,
            chi_sol,
            z,
        )
        fg = np.asarray(
            gl.fractional_fsigma8_kernel(
                prop,
                dust,
                z,
            ),
            dtype=float,
        )
        fw = np.asarray(
            gl.weyl_kernel(prop, z),
            dtype=float,
        )

        rows.append(
            {
                "z": z,
                "chi": float(lc["chi"]),
                "dipole_radial": float(
                    lc["dipole_radial"]
                ),
                "F_DL": lc["F_DL"].tolist(),
                "F_DA": lc["F_DA"].tolist(),
                "F_kappa": lc["kappa"].tolist(),
                "F_fractional_fsigma8": fg.tolist(),
                "F_local_Weyl": fw.tolist(),
            }
        )

    cal = lightcone_rows(
        prop,
        dust,
        chi_sol,
        CAL_Z,
    )
    Fcal = np.asarray(cal["F_DL"], dtype=float)
    null = unit_null(Fcal)
    obs, obs_sigma = calibration_fraction()

    rank_rows = []
    for z in (0.30, 0.50, 0.80, 1.00, 1.50, 2.10):
        fg = np.asarray(
            gl.fractional_fsigma8_kernel(
                prop,
                dust,
                z,
            ),
            dtype=float,
        )
        fw = np.asarray(
            gl.weyl_kernel(prop, z),
            dtype=float,
        )
        rank_rows.append(
            {
                "z": z,
                "det_DL_fsigma8": row_det(
                    Fcal,
                    fg,
                ),
                "det_DL_Weyl": row_det(
                    Fcal,
                    fw,
                ),
                "fsigma8_response_along_calibration_null": (
                    float(fg @ null)
                ),
                "Weyl_response_along_calibration_null": (
                    float(fw @ null)
                ),
            }
        )

    return {
        "schema": "r1-n2-gauge-invariant-DL-state-kernel-v1",
        "status": "GAUGE_INVARIANT_DL_STATE_SPACE_KERNEL_EVALUATED",
        "model": "R1",
        "geometry": {
            "K": K,
            "S_K": "sin(sqrt(K) chi)/sqrt(K)",
            "Q_2_dip_per_mu": "sin(2 sqrt(K) chi)",
            "angular_laplacian_Q_dip": "-2 Q_dip",
        },
        "lightcone_formula": {
            "redshift": (
                "delta_z=(psi-V)_o-(psi-V)_s"
                "-integral (psi-phi)' dchi"
            ),
            "radial": (
                "delta_chi=-delta_z/Hconf_s"
                "+integral(psi-phi)dchi"
            ),
            "convergence": (
                "kappa=-V_o+curved f_K lensing integral"
            ),
            "distance": (
                "deltaD_L/D_L=phi_s+delta_z"
                "+(C_K/S_K)delta_chi-kappa"
            ),
            "distance_duality": "F_DA=F_DL",
        },
        "calibration": {
            "z_star": CAL_Z,
            "A_mu": A_MU,
            "A_mu_sigma": SIGMA_A_MU,
            "observed_deltaDL_over_DL": obs,
            "observed_deltaDL_over_DL_sigma": obs_sigma,
            "F_DL_row_at_z_star": Fcal.tolist(),
            "calibration_equation": (
                "F_DL(z_star) dot X_anchor = observed_deltaDL_over_DL"
            ),
            "unit_null_direction": null.tolist(),
            "null_residual": float(Fcal @ null),
        },
        "calibration_rank_audit": {
            "phase_space_dimension": 2,
            "independent_scalar_calibrations": 1,
            "rank": 1,
            "remaining_state_direction_dimension": 1,
            "conclusion": (
                "single amplitude calibration fixes the physical amplitude "
                "only after the normalized mode direction Xhat_2 is fixed "
                "by theory or a preregistered branch condition"
            ),
            "cross_observable_rows": rank_rows,
        },
        "rows": rows,
        "claim_boundary": {
            "derived_now": [
                "gauge-invariant luminosity-distance state-space row",
                "closed-S3 distance geometry",
                "curved weak-lensing kernel for pure observer-sky dipole",
                "distance-duality equality F_DA=F_DL",
                "calibration-rank audit",
            ],
            "not_yet_claimed": [
                "unique calibrated physical state vector",
                "final scalar f_sigma8 anisotropy amplitude",
                "final scalar lensing amplitude",
                "full-redshift radial H observable",
            ],
            "important_logic": (
                "the earlier single-calibration theorem remains valid "
                "conditional on a fixed normalized mode trajectory; the "
                "SN amplitude alone does not determine that trajectory"
            ),
        },
    }


def write_tex(report: dict, path: Path) -> None:
    c = report["calibration"]
    f = c["F_DL_row_at_z_star"]
    n = c["unit_null_direction"]

    lines = [
        r"\section{Gauge-invariant luminosity-distance kernel and calibration-rank audit}",
        r"\label{sec:r1-n2-DL-kernel}",
        "",
        (
            "The remaining light-cone observable is evaluated with the "
            "gauge-invariant geometric luminosity-distance construction."
            r"\footnote{J. Yoo and F. Scaccabarozzi, JCAP 09 (2016) 046, "
            r"arXiv:1606.08453.}"
        ),
        (
            "For positive spatial curvature the background transverse "
            "distance and lensing efficiency use "
            "$S_K(\\chi)=\\sin(\\sqrt K\\chi)/\\sqrt K$, with the "
            "corresponding curved-FLRW line-of-sight kernel."
            r"\footnote{See also T. Pyne and M. Birkinshaw, MNRAS 348 "
            r"(2004) 581 for the linear luminosity distance in curved "
            r"perturbed FLRW spacetimes.}"
        ),
        "",
        r"For the pure $A_{4i}$ physical $n=2$ mode, per unit "
        r"$\mu=\hat{\mathbf n}\cdot\hat{\mathbf p}$,",
        r"\begin{equation}",
        r"Q_{2,\rm dip}(\chi)=\sin(2\sqrt K\,\chi),",
        r"\qquad \widehat\nabla^2 Q_{2,\rm dip}=-2Q_{2,\rm dip}.",
        r"\end{equation}",
        "",
        r"The fixed-observed-redshift distance fluctuation is",
        r"\begin{equation}",
        r"\boxed{",
        r"\frac{\delta D_L}{D_L}",
        r"=\phi_s+\delta z+\frac{C_K(\chi_s)}{S_K(\chi_s)}"
        r"\delta\chi-\kappa,",
        r"}",
        r"\end{equation}",
        r"with $C_K(\chi)=\cos(\sqrt K\chi)$ and the redshift, radial "
        r"and convergence distortions evaluated along the same R1 "
        r"physical-$n=2$ metric response.",
        "",
        r"\subsection{State-space calibration row}",
        (
            f"At the reference bin center $z_\\star={c['z_star']:.2f}$, "
            "the resulting luminosity-distance row on the "
            "$(\\zeta,d\\zeta/dN)$ state anchored at $z=2.1$ is"
        ),
        r"\begin{equation}",
        r"\boxed{",
        (
            r"\mathcal F_{D_L}(z_\star)="
            f"({f[0]:.9e},\,{f[1]:.9e})."
        ),
        r"}",
        r"\end{equation}",
        "",
        (
            "The carried-forward SN calibration imposes one scalar "
            "constraint on this two-dimensional state:"
        ),
        r"\begin{equation}",
        (
            r"\mathcal F_{D_L}(z_\star)\cdot X_{2,\rm anchor}"
            f"={c['observed_deltaDL_over_DL']:.9e}"
            r"\pm"
            f"{c['observed_deltaDL_over_DL_sigma']:.9e}."
        ),
        r"\end{equation}",
        "",
        r"The corresponding unit null direction is",
        r"\begin{equation}",
        (
            r"n_{\rm cal}="
            f"({n[0]:.9e},\,{n[1]:.9e}),"
            r"\qquad "
            r"\mathcal F_{D_L}\cdot n_{\rm cal}=0."
        ),
        r"\end{equation}",
        "",
        r"\subsection{Calibration-rank consequence}",
        (
            "This numerical light-cone calculation exposes a logical "
            "condition implicit in the single-calibration theorem: "
            "one amplitude calibration determines the model only after "
            "the normalized physical mode trajectory "
            "$\widehat X_2$ has been fixed independently."
        ),
        (
            "The luminosity-distance row has rank one on a "
            "two-dimensional phase space.  The growth and Weyl rows "
            "are not generally parallel to it, so motion along "
            "$n_{\\rm cal}$ leaves the SN calibration unchanged while "
            "altering the other observables."
        ),
        "",
        r"\boxed{\text{The next prediction-closure object is the theory-selected "
        r"normalized state direction }\widehat X_2.}",
        "",
        r"\subsection{Claim boundary}",
        (
            "No second observational calibration is introduced.  "
            "Until the branch direction is fixed by the retained action "
            "or by an explicit preregistered theoretical branch "
            "condition, the manuscript reports the observable kernels "
            "as state-space rows rather than unique scalar amplitudes."
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
