#!/usr/bin/env python3
"""
Action-native closed-S^3 gradient / dispersion audit for the R1 cosmology.

Purpose
-------
This is a paper-level stability audit of the existing action-native
matter+radiation quadratic system.  It does NOT change the frozen background,
fit any observational amplitude, invoke Gate 7, or retune the cosmology.

The script:
  1. imports code/action_native_matter_n2.py as the authoritative n=2 action;
  2. verifies that the generalized harmonic continuation reproduces n=2
     exactly at q=n(n+2)=8;
  3. constructs the fully constraint-reduced quadratic action for physical
     scalar harmonics n>=2;
  4. checks positivity of the reduced kinetic matrix;
  5. constructs the instantaneous second-order dispersion pencil

         P(s) = s^2 M + s Gamma + Omega,

     including background-time derivatives:
         Gamma = 3 H M + Mdot + L - L^T,
         Omega = 3 H L + Ldot - V;

  6. solves the generalized quadratic eigenvalue problem without explicitly
     inverting M;
  7. separates finite-harmonic growth from genuine high-k gradient behavior;
  8. extracts high-k characteristic-speed-squared limits and checks that the
     radiation branch approaches c_s^2=1/3;
  9. records the pure gravity/scalar Akama-Kobayashi F_S/G_S diagnostic as an
     independent cross-check.

Outputs
-------
artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_V1.json
docs/R1_ACTION_NATIVE_GRADIENT_DISPERSION_V1.md

Interpretation boundary
-----------------------
A pressureless-dust growing mode is not, by itself, a gradient instability.
The relevant high-k diagnostic is whether Re(s)/sqrt(k_phys^2) approaches a
nonzero positive constant as k increases.  The script therefore reports both
finite-harmonic roots and a high-k principal-symbol extrapolation.

The continuation from n=2 to n>=3 changes only the S^3 scalar-harmonic
eigenvalues in the SAME quadratic action:
    q_n = n(n+2),
    k_phys^2 = q_n K/a^2,
    shifted curvature factor = (q_n-3) K/a^2.
It is used here as a stability/principal-symbol audit, not as a new
observational prediction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.linalg import eig


# ---------------------------------------------------------------------------
# Repository discovery
# ---------------------------------------------------------------------------

def discover_repo(explicit: Path | None) -> Path:
    candidates: list[Path] = []
    if explicit is not None:
        candidates.append(explicit.expanduser().resolve())

    cwd = Path.cwd().resolve()
    candidates.extend([cwd, *cwd.parents])

    home = Path.home()
    candidates.extend([
        home / "Manuscript-Generation",
        home / "Documents" / "Manuscript-Generation",
        home / "OneDrive" / "Documents" / "Manuscript-Generation",
    ])

    seen: set[Path] = set()
    for p in candidates:
        try:
            p = p.resolve()
        except Exception:
            continue
        if p in seen:
            continue
        seen.add(p)
        if (
            (p / "code" / "action_native_matter_n2.py").exists()
            and (p / "code" / "closed_horndeski_n2.py").exists()
            and (p / "preregistration" / "prediction_manifest.json").exists()
        ):
            return p

    raise FileNotFoundError(
        "Could not locate Manuscript-Generation. Run from the repository or "
        "pass --repo C:\\Users\\carbe\\Manuscript-Generation"
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Generic harmonic continuation of the already-derived action
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HarmonicPoint:
    n: int
    q: float
    k2: float
    shifted: float


def harmonic_point(N: complex, n: int, bg) -> HarmonicPoint:
    if n < 2:
        raise ValueError("Physical scalar audit begins at n=2.")
    a = np.exp(N)
    kappa = bg.KCURV / a**2
    q = float(n * (n + 2))
    return HarmonicPoint(
        n=n,
        q=q,
        k2=q * kappa,
        shifted=(q - 3.0) * kappa,
    )


def action_blocks_harmonic(N, y, n: int, an, bg, include_fluids: bool = True):
    """
    Same ADM + Schutz-Sorkin quadratic action as action_native_matter_n2,
    with only the S^3 harmonic eigenvalues generalized.

    x=(zeta,D_m,D_r), u=(alpha,chi,v_m,v_r).
    """
    _, v, H = y
    theta, sigma = an.owner.theta_sigma(v, H)

    hp = harmonic_point(N, n, bg)
    a = np.exp(N)
    h_all = np.array(
        [3 * bg.OMEGA_M0 / a**3, 4 * bg.OMEGA_R0 / a**4]
    )
    h = h_all if include_fluids else h_all[:0]
    w = an.W if include_fluids else an.W[:0]

    nx = 1 + len(h)
    nu = 2 + len(h)
    dtype = np.result_type(N, y, float)

    A = np.zeros((nu, nu), dtype=dtype)
    B = np.zeros((nu, nx), dtype=dtype)
    D = np.zeros((nu, nx), dtype=dtype)
    K0 = np.zeros((nx, nx), dtype=dtype)
    E = np.zeros((nx, nx), dtype=dtype)

    k = hp.k2
    kap = bg.KCURV / a**2
    s = hp.shifted

    # Gravity/scalar block: identical formula to the n=2 owner.
    A[0, 0] = 2 * sigma
    A[0, 1] = A[1, 0] = 2 * k * theta
    A[1, 1] = -2 * k * kap

    B[0, 0] = 6 * theta
    B[1, 0] = -2 * k

    D[0, 0] = 2 * s + 3 * np.sum(h)
    K0[0, 0] = -6
    E[0, 0] = 2 * s

    # Exact ideal-fluid auxiliaries in the same IBP convention.
    for i, (hi, wi) in enumerate(zip(h, w)):
        j = i + 1
        u = i + 2

        A[u, u] = -hi * k
        A[1, u] = A[u, 1] = hi * k

        B[u, j] = hi
        D[0, j] = -hi

        delta = np.zeros(nx, dtype=dtype)
        delta[0] = -3
        delta[j] = 1
        E -= hi * wi * np.outer(delta, delta)

    return A, B, D, K0, E


def reduced_blocks_harmonic(N, y, n: int, an, bg, include_fluids: bool = True):
    A, B, D, K0, E = action_blocks_harmonic(
        N, y, n, an, bg, include_fluids
    )
    AB = np.linalg.solve(A, B)
    AD = np.linalg.solve(A, D)
    M = K0 - B.T @ AB
    L = -B.T @ AD
    V = E - D.T @ AD
    return M, L, V


def along_background_derivative_matrix(function, N, y, bg, step=1e-25):
    """
    Complex-step derivative d/dN along the frozen background trajectory.
    This mirrors action_native_matter_n2.along_background_derivative.
    """
    y_real = np.asarray(y, dtype=float)
    yN = bg.rhs(float(N), y_real)
    value = function(
        N + 1j * step,
        y_real.astype(complex) + 1j * step * yN,
    )
    return np.imag(value) / step


def exact_eom_matrices(N, y, n: int, an, bg):
    """
    Return M, Gamma, Omega for

        M xddot + Gamma xdot + Omega x = 0

    obtained from d/dt[a^3(M xdot + L x)] -
    a^3(L^T xdot + V x)=0.
    """
    M, L, V = reduced_blocks_harmonic(N, y, n, an, bg, True)

    MN = along_background_derivative_matrix(
        lambda nn, yy: reduced_blocks_harmonic(nn, yy, n, an, bg, True)[0],
        N,
        y,
        bg,
    )
    LN = along_background_derivative_matrix(
        lambda nn, yy: reduced_blocks_harmonic(nn, yy, n, an, bg, True)[1],
        N,
        y,
        bg,
    )

    H = float(np.real(y[2]))
    Mdot = H * MN
    Ldot = H * LN

    Gamma = 3.0 * H * M + Mdot + L - L.T
    Omega = 3.0 * H * L + Ldot - V
    return (
        np.asarray(np.real_if_close(M), dtype=float),
        np.asarray(np.real_if_close(Gamma), dtype=float),
        np.asarray(np.real_if_close(Omega), dtype=float),
        np.asarray(np.real_if_close(L), dtype=float),
        np.asarray(np.real_if_close(V), dtype=float),
    )


def dispersion_roots(M, Gamma, Omega):
    """
    Generalized linearization of
        det(s^2 M + s Gamma + Omega)=0.

    Avoids an explicit M^{-1}, which is important because the fluid-coordinate
    kinetic scales become strongly harmonic-dependent.
    """
    d = M.shape[0]
    Z = np.zeros((d, d))
    I = np.eye(d)

    A = np.block([
        [Z, I],
        [-Omega, -Gamma],
    ])
    B = np.block([
        [I, Z],
        [Z, M],
    ])

    roots = eig(A, B, right=False, check_finite=True)
    roots = np.asarray(roots, dtype=complex)

    # Reject infinite generalized roots explicitly.
    if not np.all(np.isfinite(roots.real)) or not np.all(np.isfinite(roots.imag)):
        raise FloatingPointError("Non-finite generalized dispersion root.")

    return roots


def encode_complex(values: Iterable[complex]):
    return [[float(complex(v).real), float(complex(v).imag)] for v in values]


def sort_roots(roots):
    return np.array(
        sorted(
            roots,
            key=lambda z: (
                round(float(z.real), 13),
                round(float(z.imag), 13),
            ),
        ),
        dtype=complex,
    )


def positive_imag_speeds(roots, k2, imag_floor=1e-9):
    scale = math.sqrt(float(k2))
    vals = [
        abs(float(r.imag)) / scale
        for r in roots
        if float(r.imag) > imag_floor * max(1.0, scale)
    ]
    vals.sort()
    return vals


# ---------------------------------------------------------------------------
# Audits
# ---------------------------------------------------------------------------

def n2_owner_match(sol, an, bg):
    rows = []
    worst_action = 0.0
    worst_reduced = 0.0

    for z in [2.1, 1.5, 1.0, 0.5, 0.1, 0.0]:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        original = an.action_blocks(N, y, True)
        generic = action_blocks_harmonic(N, y, 2, an, bg, True)

        action_err = max(
            float(np.max(np.abs(a - b)))
            for a, b in zip(original, generic)
        )

        r0 = an.reduced_blocks(N, y, True)
        r1 = reduced_blocks_harmonic(N, y, 2, an, bg, True)
        reduced_err = max(
            float(np.max(np.abs(a - b)))
            for a, b in zip(r0, r1)
        )

        worst_action = max(worst_action, action_err)
        worst_reduced = max(worst_reduced, reduced_err)

        rows.append({
            "z": z,
            "action_block_max_abs_difference": action_err,
            "reduced_block_max_abs_difference": reduced_err,
        })

    return {
        "rows": rows,
        "worst_action_block_difference": worst_action,
        "worst_reduced_block_difference": worst_reduced,
        "exact_to_float_tolerance": bool(
            worst_action < 5e-12 and worst_reduced < 5e-12
        ),
    }


def gravity_sector_owner_scan(sol, owner):
    rows = []
    min_g = math.inf
    min_f = math.inf
    min_c2 = math.inf

    for z in np.linspace(0.0, 2.1, 85):
        N = -math.log1p(float(z))
        y = np.asarray(sol.sol(N), dtype=float)
        v, H = float(y[1]), float(y[2])

        G = float(owner.gs_n2(v, H))
        F = float(owner.formal_fs_n2_gravity_sector(N, y))
        c2 = F / G

        min_g = min(min_g, G)
        min_f = min(min_f, F)
        min_c2 = min(min_c2, c2)

        rows.append({
            "z": float(z),
            "G_S_n2": G,
            "F_S_n2_gravity_sector_only": F,
            "c_s2_gravity_sector_only": c2,
        })

    return {
        "min_G_S_n2": min_g,
        "min_F_S_n2_gravity_sector_only": min_f,
        "min_c_s2_gravity_sector_only": min_c2,
        "rows": rows,
        "scope": (
            "Independent pure gravity+scalar Akama-Kobayashi diagnostic; "
            "not substituted for the full fluid-coupled dispersion pencil."
        ),
    }


def kinetic_scan(sol, an, bg, ns):
    rows = []
    global_min = math.inf
    global_owner = None

    for z in np.linspace(0.0, 2.1, 85):
        N = -math.log1p(float(z))
        y = np.asarray(sol.sol(N), dtype=float)

        for n in ns:
            M, _, _ = reduced_blocks_harmonic(N, y, n, an, bg, True)
            ev = np.linalg.eigvalsh(np.asarray(M, dtype=float))
            mn = float(ev[0])

            if mn < global_min:
                global_min = mn
                global_owner = {
                    "z": float(z),
                    "n": int(n),
                    "eigenvalues": ev.tolist(),
                }

            rows.append({
                "z": float(z),
                "n": int(n),
                "kinetic_eigenvalues": ev.tolist(),
            })

    return {
        "minimum_eigenvalue": global_min,
        "minimum_owner": global_owner,
        "all_positive": bool(global_min > 0.0),
        "rows": rows,
    }


def finite_dispersion_scan(sol, an, bg, ns):
    rows = []

    for z in [2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.1, 0.0]:
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)
        H = float(y[2])

        for n in ns:
            hp = harmonic_point(N, n, bg)
            M, Gamma, Omega, _, _ = exact_eom_matrices(
                N, y, n, an, bg
            )
            roots = sort_roots(dispersion_roots(M, Gamma, Omega))

            rows.append({
                "z": z,
                "n": int(n),
                "q": hp.q,
                "k_phys2": float(np.real(hp.k2)),
                "kinetic_eigenvalues": np.linalg.eigvalsh(M).tolist(),
                "roots_s": encode_complex(roots),
                "max_Re_s_over_H": float(np.max(roots.real) / H),
                "max_abs_Im_s_over_H": float(np.max(np.abs(roots.imag)) / H),
                "positive_imag_speeds_over_sqrt_k2":
                    positive_imag_speeds(roots, hp.k2),
            })

    return rows


def principal_symbol_scan(sol, an, bg, high_ns):
    rows = []
    z_summaries = []

    for z in np.linspace(0.0, 2.1, 22):
        z = float(z)
        N = -math.log1p(z)
        y = np.asarray(sol.sol(N), dtype=float)

        seq = []
        for n in high_ns:
            hp = harmonic_point(N, n, bg)
            M, Gamma, Omega, _, _ = exact_eom_matrices(
                N, y, n, an, bg
            )
            roots = dispersion_roots(M, Gamma, Omega)
            scale = math.sqrt(float(np.real(hp.k2)))
            normalized = roots / scale
            speeds = positive_imag_speeds(roots, hp.k2)

            item = {
                "z": z,
                "n": int(n),
                "k_phys2": float(np.real(hp.k2)),
                "normalized_roots_s_over_sqrt_k2":
                    encode_complex(sort_roots(normalized)),
                "max_positive_normalized_growth":
                    float(max(0.0, np.max(normalized.real))),
                "positive_imag_speeds":
                    [float(v) for v in speeds],
                "positive_imag_speed2":
                    [float(v * v) for v in speeds],
            }
            seq.append(item)
            rows.append(item)

        # Fit the normalized growth to an intercept in x=1/sqrt(k^2).
        xs = np.array([
            1.0 / math.sqrt(r["k_phys2"]) for r in seq
        ])
        growth = np.array([
            r["max_positive_normalized_growth"] for r in seq
        ])
        growth_fit = np.polyfit(xs, growth, 1)
        growth_intercept = float(growth_fit[1])

        # Two nonzero propagating scalar/radiation branches are expected
        # at high k.  Fit speed^2 branch-wise, sorted by speed.
        branch_c2 = []
        if all(len(r["positive_imag_speed2"]) >= 2 for r in seq):
            for j in range(2):
                ys = np.array([
                    sorted(r["positive_imag_speed2"])[j] for r in seq
                ])
                fit = np.polyfit(xs, ys, 1)
                branch_c2.append(float(fit[1]))

        radiation_c2 = None
        scalar_c2 = None
        radiation_error = None

        if len(branch_c2) == 2:
            distances = [abs(v - 1.0 / 3.0) for v in branch_c2]
            ir = int(np.argmin(distances))
            radiation_c2 = branch_c2[ir]
            scalar_c2 = branch_c2[1 - ir]
            radiation_error = abs(radiation_c2 - 1.0 / 3.0)

        # A pressureless matter growth exponent is O(H), so after division
        # by sqrt(k^2) it must decay.  Genuine gradient instability instead
        # approaches a nonzero positive normalized growth.
        monotone_tail = bool(
            seq[-1]["max_positive_normalized_growth"]
            <= seq[-2]["max_positive_normalized_growth"] * 1.05 + 1e-12
        )

        z_summaries.append({
            "z": z,
            "high_n_sequence": seq,
            "growth_intercept_fit": growth_intercept,
            "growth_tail_is_nonincreasing": monotone_tail,
            "branch_c2_intercepts": branch_c2,
            "radiation_c2_intercept": radiation_c2,
            "radiation_c2_abs_error_from_1_over_3": radiation_error,
            "scalar_c2_intercept": scalar_c2,
        })

    growth_intercepts = [
        s["growth_intercept_fit"] for s in z_summaries
    ]
    scalar_c2_values = [
        s["scalar_c2_intercept"]
        for s in z_summaries
        if s["scalar_c2_intercept"] is not None
    ]
    radiation_errors = [
        s["radiation_c2_abs_error_from_1_over_3"]
        for s in z_summaries
        if s["radiation_c2_abs_error_from_1_over_3"] is not None
    ]

    return {
        "high_ns": [int(n) for n in high_ns],
        "rows": rows,
        "z_summaries": z_summaries,
        "max_growth_intercept_fit":
            float(max(growth_intercepts)) if growth_intercepts else None,
        "min_scalar_c2_intercept":
            float(min(scalar_c2_values)) if scalar_c2_values else None,
        "max_radiation_c2_error":
            float(max(radiation_errors)) if radiation_errors else None,
        "all_tail_growth_nonincreasing":
            bool(all(s["growth_tail_is_nonincreasing"] for s in z_summaries)),
    }


def numerical_classification(kinetic, principal):
    """
    Conservative numerical classification, not a theorem label.

    Thresholds are intentionally loose and are emitted into the artifact.
    They detect an O(k) instability, not ordinary pressureless-matter growth.
    """
    thresholds = {
        "kinetic_min_gt": 0.0,
        "max_high_k_growth_intercept_lt": 5e-3,
        "min_scalar_c2_gt": 1e-4,
        "radiation_c2_error_lt": 5e-2,
    }

    if not kinetic["all_positive"]:
        status = "FAIL_REDUCED_KINETIC_MATRIX_NOT_POSITIVE"
    elif principal["min_scalar_c2_intercept"] is None:
        status = "INCONCLUSIVE_COULD_NOT_ISOLATE_TWO_PROPAGATING_HIGH_K_BRANCHES"
    elif (
        principal["max_growth_intercept_fit"] is not None
        and principal["max_growth_intercept_fit"]
        >= thresholds["max_high_k_growth_intercept_lt"]
    ):
        status = "FAIL_HIGH_K_SCALED_GROWTH_PERSISTS"
    elif (
        principal["min_scalar_c2_intercept"]
        <= thresholds["min_scalar_c2_gt"]
    ):
        status = "FAIL_SCALAR_PRINCIPAL_C2_NONPOSITIVE_OR_NEAR_ZERO"
    elif (
        principal["max_radiation_c2_error"] is not None
        and principal["max_radiation_c2_error"]
        >= thresholds["radiation_c2_error_lt"]
    ):
        status = "INCONCLUSIVE_RADIATION_BRANCH_NOT_CLEANLY_RECOVERED"
    elif not principal["all_tail_growth_nonincreasing"]:
        status = "INCONCLUSIVE_HIGH_K_GROWTH_TAIL_NOT_CONVERGED"
    else:
        status = "PASS_NUMERICAL_PRINCIPAL_GRADIENT_DISPERSION_AUDIT"

    return status, thresholds


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def markdown_report(data):
    p = data["principal_symbol"]
    k = data["kinetic"]
    g = data["gravity_sector_owner"]

    lines = [
        "# R1 action-native gradient / dispersion audit v1",
        "",
        f"**Status:** `{data['status']}`",
        "",
        "This audit uses the existing action-native ADM + Schutz--Sorkin "
        "quadratic system. It performs no observational fit, no parameter "
        "retuning, and no Gate-7 operation.",
        "",
        "## Owner-reproduction gate",
        "",
        f"- Worst n=2 action-block difference: "
        f"`{data['n2_owner_match']['worst_action_block_difference']:.6e}`",
        f"- Worst n=2 reduced-block difference: "
        f"`{data['n2_owner_match']['worst_reduced_block_difference']:.6e}`",
        "",
        "## Kinetic gate",
        "",
        f"- Minimum reduced kinetic eigenvalue over the scan: "
        f"`{k['minimum_eigenvalue']:.12e}`",
        f"- All scanned kinetic matrices positive: `{k['all_positive']}`",
        "",
        "## Full fluid-coupled principal dispersion",
        "",
        f"- High-n sequence: `{p['high_ns']}`",
        f"- Maximum extrapolated Re(s)/sqrt(k^2) intercept: "
        f"`{p['max_growth_intercept_fit']}`",
        f"- Minimum extrapolated non-radiation propagating c^2: "
        f"`{p['min_scalar_c2_intercept']}`",
        f"- Maximum radiation-branch |c^2-1/3|: "
        f"`{p['max_radiation_c2_error']}`",
        f"- All high-k growth tails nonincreasing: "
        f"`{p['all_tail_growth_nonincreasing']}`",
        "",
        "A pressureless-dust growing mode is not classified as a gradient "
        "instability unless its growth scales as sqrt(k^2) at high k.",
        "",
        "## Independent pure gravity/scalar owner diagnostic",
        "",
        f"- min G_S,2: `{g['min_G_S_n2']:.12e}`",
        f"- min F_S,2 (gravity/scalar only): "
        f"`{g['min_F_S_n2_gravity_sector_only']:.12e}`",
        f"- min F_S,2/G_S,2: "
        f"`{g['min_c_s2_gravity_sector_only']:.12e}`",
        "",
        "The F_S/G_S result is reported as an independent owner diagnostic; "
        "the manuscript-level stability classification comes from the full "
        "matter+radiation dispersion pencil.",
        "",
        "## Claim boundary",
        "",
        "- This is a linear/quadratic stability audit on frozen R1.",
        "- It is not nonlinear fluid closure.",
        "- The n>=3 continuation is used to expose the principal spatial "
        "symbol; it is not a new high-n observational fit.",
        "- A numerical PASS supports the paper-level statement that no "
        "ghost or high-k gradient instability was detected over the audited "
        "redshift/harmonic domain. It is not a microscopic BHSM theorem.",
        "",
    ]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--repo",
        type=Path,
        default=None,
        help="Path to Manuscript-Generation repository.",
    )
    ap.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional JSON output path.",
    )
    ap.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional Markdown report path.",
    )
    ap.add_argument(
        "--high-ns",
        type=int,
        nargs="+",
        default=[40, 80, 160, 320],
        help="Harmonics used for principal-symbol extrapolation.",
    )
    args = ap.parse_args()

    repo = discover_repo(args.repo)
    code = repo / "code"
    if str(code) not in sys.path:
        sys.path.insert(0, str(code))

    import action_native_matter_n2 as an
    import background as bg
    import closed_horndeski_n2 as owner

    output = (
        args.output
        if args.output is not None
        else repo / "artifacts" / "R1_ACTION_NATIVE_GRADIENT_DISPERSION_V1.json"
    )
    report = (
        args.report
        if args.report is not None
        else repo / "docs" / "R1_ACTION_NATIVE_GRADIENT_DISPERSION_V1.md"
    )

    print("=" * 96)
    print("R1 ACTION-NATIVE GRADIENT / DISPERSION AUDIT")
    print("=" * 96)
    print("repo =", repo)
    print("NO OBSERVATIONAL FIT")
    print("NO PARAMETER RETUNING")
    print("NO GATE-7 DEPENDENCY")
    print()

    sol = an.frozen_background()

    print("[1/5] Verifying generalized action reproduces physical n=2 owner...")
    match = n2_owner_match(sol, an, bg)
    print(
        "  worst action difference =",
        f"{match['worst_action_block_difference']:.3e}",
    )
    print(
        "  worst reduced difference =",
        f"{match['worst_reduced_block_difference']:.3e}",
    )
    if not match["exact_to_float_tolerance"]:
        raise RuntimeError(
            "Generalized harmonic action does not reproduce n=2 owner."
        )

    print("[2/5] Scanning reduced kinetic matrices...")
    kinetic_ns = [2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320]
    kinetic = kinetic_scan(sol, an, bg, kinetic_ns)
    print(
        "  minimum kinetic eigenvalue =",
        f"{kinetic['minimum_eigenvalue']:.12e}",
    )
    print("  owner =", kinetic["minimum_owner"])

    print("[3/5] Solving finite-harmonic full dispersion pencils...")
    finite_ns = [2, 3, 4, 5, 8, 12, 20, 40]
    finite = finite_dispersion_scan(sol, an, bg, finite_ns)
    print("  finite dispersion rows =", len(finite))

    print("[4/5] Extracting high-k principal symbol...")
    principal = principal_symbol_scan(sol, an, bg, args.high_ns)
    print(
        "  max growth intercept Re(s)/sqrt(k2) =",
        principal["max_growth_intercept_fit"],
    )
    print(
        "  min scalar c^2 intercept =",
        principal["min_scalar_c2_intercept"],
    )
    print(
        "  max radiation |c^2-1/3| =",
        principal["max_radiation_c2_error"],
    )

    print("[5/5] Evaluating independent gravity/scalar F_S/G_S owner...")
    gravity = gravity_sector_owner_scan(sol, owner)
    print(
        "  min G_S2 =",
        f"{gravity['min_G_S_n2']:.12e}",
    )
    print(
        "  min F_S2 =",
        f"{gravity['min_F_S_n2_gravity_sector_only']:.12e}",
    )
    print(
        "  min F_S2/G_S2 =",
        f"{gravity['min_c_s2_gravity_sector_only']:.12e}",
    )

    status, thresholds = numerical_classification(kinetic, principal)

    data = {
        "schema": "R1-action-native-gradient-dispersion-v1",
        "status": status,
        "repo": str(repo),
        "source_hashes": {
            "action_native_matter_n2.py":
                sha256(repo / "code" / "action_native_matter_n2.py"),
            "closed_horndeski_n2.py":
                sha256(repo / "code" / "closed_horndeski_n2.py"),
            "background.py":
                sha256(repo / "code" / "background.py"),
            "prediction_manifest.json":
                sha256(repo / "preregistration" / "prediction_manifest.json"),
        },
        "observational_fits": 0,
        "parameter_retuning": False,
        "gate7_dependency": False,
        "harmonic_rule": {
            "scalar_laplacian_eigenvalue_q_n": "n(n+2)",
            "k_phys2": "q_n*K/a^2",
            "shifted_scalar_curvature_factor": "(q_n-3)*K/a^2",
            "physical_minimum_n": 2,
        },
        "dispersion_equation": (
            "det[s^2 M + s Gamma + Omega]=0; "
            "Gamma=3HM+Mdot+L-L^T; Omega=3HL+Ldot-V"
        ),
        "classification_thresholds": thresholds,
        "n2_owner_match": match,
        "kinetic": kinetic,
        "finite_dispersion": finite,
        "principal_symbol": principal,
        "gravity_sector_owner": gravity,
        "claim_boundary": [
            "Linear/quadratic R1 stability audit only.",
            "Pressureless-dust growth is not called a gradient instability "
            "unless it scales as sqrt(k_phys2).",
            "n>=3 is used to expose the principal spatial symbol, not to "
            "create or fit a new cosmological prediction.",
            "No local-source normalization or Gate-7 result enters.",
            "A numerical PASS is a paper-level stability result, not a "
            "complete microscopic BHSM proof.",
        ],
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    report.parent.mkdir(parents=True, exist_ok=True)

    output.write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    report.write_text(markdown_report(data), encoding="utf-8")

    print()
    print("=" * 96)
    print("RESULT")
    print("=" * 96)
    print("STATUS =", status)
    print("WROTE:", output)
    print("WROTE:", report)
    print("DONE")


if __name__ == "__main__":
    main()
