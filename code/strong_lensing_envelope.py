"""Action-owned envelope projections, with the coordinate domains kept explicit.

The available R1 projection is in ADM auxiliaries (alpha, chi), per a**3 and
one normalized n=2 harmonic. It is NOT the local (Phi,Psi)/phase block.
The static spatial target-phase null is a separate conditional application.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import action_native_matter_n2 as an
import bhsm_optical_current_gate as current_gate


SOURCE_COORDINATES = (
    "zeta_dot", "D_m_dot", "D_r_dot", "zeta", "D_m", "D_r"
)


def explicit_two_coordinate_response(metric: np.ndarray, mixed: np.ndarray) -> np.ndarray:
    """Exact -A^-1 B; the caller owns the interpretation of its coordinates."""
    A, B = np.asarray(metric, float), np.asarray(mixed, float)
    if A.shape != (2, 2) or B.shape[0] != 2:
        raise ValueError("requires a 2x2 metric block and two mixed rows")
    if not np.allclose(A, A.T, rtol=0, atol=1e-14):
        raise ValueError("the action block must be symmetric")
    a, b, c = A[0, 0], A[0, 1], A[1, 1]
    det = a*c-b*b
    if det == 0:
        raise ValueError("singular action block")
    return np.stack((b*B[1]-c*B[0], b*B[0]-a*B[1]))/det


def adm_projected_blocks(N: float, y: np.ndarray) -> dict[str, np.ndarray]:
    """Eliminate the two fluid velocity auxiliaries from the original action.

    q=(xdot,x). L/a^3=.5 m A_m m + m J_m q + .5 q C_m q.
    All entries come from an.action_blocks, not from optical fitting.
    """
    A, B, D, K0, E = an.action_blocks(N, y)
    J = np.column_stack((B, D))
    C = np.block([[K0, np.zeros_like(K0)], [np.zeros_like(K0), E]])
    Av, Amv = A[2:, 2:], A[:2, 2:]
    Jv = J[2:]
    metric = A[:2, :2]-Amv@np.linalg.solve(Av, Amv.T)
    mixed = J[:2]-Amv@np.linalg.solve(Av, Jv)
    state = C-Jv.T@np.linalg.solve(Av, Jv)
    return {"metric": metric, "mixed": mixed, "state": state}


def adm_explicit_blocks(N: float, y: np.ndarray) -> dict[str, np.ndarray]:
    """Independent analytic polynomial coefficients after v_I=chi+Ddot_I/k."""
    theta, sigma = an.owner.theta_sigma(y[1], y[2])
    _, kap, k, s, h = an.scales(N)
    hm, hr = h
    metric = np.array([[2*sigma, 2*k*theta],
                       [2*k*theta, k*(hm+hr-2*kap)]])
    mixed = np.array([[6*theta, 0, 0, 2*s+3*(hm+hr), -hm, -hr],
                      [-2*k, hm, hr, 0, 0, 0]])
    kinetic = np.diag([-6., hm/k, hr/k])
    state_potential = np.diag([2*s, 0., 0.])
    radiation_density = np.array([-3., 0., 1.])
    state_potential -= (hr/3)*np.outer(radiation_density, radiation_density)
    state = np.block([[kinetic, np.zeros((3, 3))],
                      [np.zeros((3, 3)), state_potential]])
    return {"metric": metric, "mixed": mixed, "state": state}


def adm_projection_snapshot(z: float = .3, background=None) -> dict:
    """Numerical native Hessian, source columns, Schur and original residuals.

    Unit source coordinates are derivative probes, not realized amplitudes.
    The one harmonic obeys integral sqrt(gamma)*Q**2=1, H0=Mpl=1.
    """
    if not 0 <= z <= 2.1:
        raise ValueError("outside the inherited audited R1 domain")
    sol = an.frozen_background() if background is None else background
    N = -np.log1p(z)
    y = sol.sol(N)
    blocks = adm_projected_blocks(N, y)
    explicit = adm_explicit_blocks(N, y)
    A, J, C = (blocks[key] for key in ("metric", "mixed", "state"))
    response = explicit_two_coordinate_response(A, J)
    schur = C-J.T@np.linalg.solve(A, J)
    M, L, V = an.reduced_blocks(N, y)
    native_schur = np.block([[M, L], [L.T, V]])
    reconstruction_errors, original_residuals = [], []
    for q, m in zip(np.eye(6), response.T):
        xd, x = q[:3], q[3:]
        _, _, k, _, _ = an.scales(N)
        velocities = m[1]+xd[1:]/k
        u = np.r_[m, velocities]
        reconstruction_errors.append(np.linalg.norm(u-an.reconstruct(N, y, x, xd)))
        original_residuals.append(an.constraint_residuals(N, y, x, xd, u).tolist())
    zeta_column = J[:, 3]
    return {
        "status": "EXECUTED_ACTION_NATIVE_ADM_PROJECTION_NOT_LOCAL_PHASE_EXTRACTION",
        "z": float(z), "N": float(N), "background_phi_v_H": y.tolist(),
        "units": "H0=Mpl=1; action density divided by a^3",
        "normalization": "one real closed-S3 n=2 harmonic: integral sqrt(gamma) Q^2=1",
        "metric_basis_normalization": (
            "unit lapse harmonic amplitude alpha; unit covariant shift potential "
            "harmonic amplitude chi with N_i=D_i chi; chi measured in H0^-1"
        ),
        "metric_coordinates": ["alpha_ADM_lapse", "chi_ADM_shift"],
        "source_coordinates": list(SOURCE_COORDINATES),
        "a_ADM": float(A[0, 0]), "b_ADM": float(A[0, 1]), "c_ADM": float(A[1, 1]),
        "metric_block": A.tolist(), "mixed_columns": J.tolist(),
        "source_state_self_block": C.tolist(),
        "determinant": float(np.linalg.det(A)),
        "metric_condition_number": float(np.linalg.cond(A)),
        "response_columns": response.tolist(),
        "response_residual_norm": float(np.linalg.norm(A@response+J)),
        "explicit_inverse_vs_solve_norm": float(np.linalg.norm(response+np.linalg.solve(A, J))),
        "analytic_vs_projected_max": float(max(np.max(np.abs(blocks[k]-explicit[k])) for k in blocks)),
        "original_constraints_by_unit_column": original_residuals,
        "original_constraints_max": float(np.max(np.abs(original_residuals))),
        "simultaneous_aux_reconstruction_max_norm": float(max(reconstruction_errors)),
        "schur_state_block": schur.tolist(),
        "schur_vs_native_reduction_norm": float(np.linalg.norm(schur-native_schur)),
        "c_zeta_before_metric_elimination": float(C[3, 3]),
        "zeta_schur_curvature": float(C[3, 3]-zeta_column@np.linalg.solve(A, zeta_column)),
        "kinetic_eigenvalues": np.linalg.eigvalsh(M).tolist(),
        "legendre_gate": bool(np.min(np.linalg.eigvalsh(M)) > 0),
        "common_domain": "same N, R1 background, normalized n=2 harmonic, endpoint-vanishing temporal variations",
        "local_phase_extraction_dependency": (
            "ADM lapse/shift entries cannot be renamed local Phi/Psi entries. "
            "Phi=alpha+chi_dot; Psi=-zeta-H chi is a derivative/state map, "
            "not an instantaneous invertible 2x2 change of metric coordinates."
        ),
    }


def static_spatial_phase_null() -> dict:
    """Execute the two exact null channels in the retained original bilinear.

    This algebraic application holds pointwise on every admissible normalized
    domain and for every finite F',F'',A_ij,g_ij. The displayed finite tensor
    is a representative identity check, not a physical profile selection.
    """
    fp, fpp = np.array([1.25]), np.array([-.3])
    tensor = np.array([np.diag([.7, .2, 1.1])])
    g = np.array([np.eye(3)])
    nonzero = np.array([[.4, -.2, .1]])
    zero = np.zeros((1, 3))
    applications = {}
    for name, j, grad in (("zero_spatial_current", zero, nonzero),
                          ("constant_local_texture", nonzero, zero)):
        source = current_gate.spatial_mixed_source_tensor(fp, fpp, j, grad, tensor, g)
        uv = current_gate.project_mixed_coefficients(np.ones(1), tensor, g, source)
        applications[name] = {"pointwise_source_tensor": source[0].tolist(),
                              "u": float(uv[0]), "v": float(uv[1]),
                              "source_norm": float(np.linalg.norm(source)),
                              "scope": "exact static spatial phase bilinear identity, conditional on stated channel"}
    return {
        "status": "DERIVED_EXACT_ZERO_IN_SPECIFIED_STATIC_SPATIAL_PHASE_CHANNELS",
        "applications": applications,
        "metric_response": "zero if the same-domain physical metric block is invertible",
        "schur_curvature": "c_phi (because u=v=0); not determined by the current-null condition",
        "phase_invariant_FX_sector_self_entry": (
            "integral dmu*w*[4 F''(X)(j.dphi)^2+2 F'(X)|T eta|^2 (dphi)^2]"
        ),
        "constant_texture_FX_sector_self_entry": 0.0,
        "constant_texture_self_entry_scope": (
            "only the phase-invariant F(X) sector on a domain admitting the constant; "
            "additional action and boundary terms must be included before claiming a full c_phi"
        ),
        "not_inferred": "D0 eta=0 and j0=0 alone do not imply ji=0 or zero full lensing",
    }


def checkpoint_snapshot(z: float = .3, background=None) -> dict:
    return {"schema": "strong-lensing-envelope-action-projection-v1",
            "adm_native_projection": adm_projection_snapshot(z, background),
            "local_static_phase": static_spatial_phase_null(),
            "seam_used": "no seam amplitude inserted; static export remains symbolic",
            "moving_seam_status": "kinematic proposal only; not selected as an action-derived source",
            "muon_input_consumed": False,
            "realized_physical_envelope_amplitude_claimed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z", type=float, default=.3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = checkpoint_snapshot(args.z)
    content = json.dumps(data, indent=2, allow_nan=False)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content, encoding="utf8")
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
