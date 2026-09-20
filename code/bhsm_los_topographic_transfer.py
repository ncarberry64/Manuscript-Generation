from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad

import background
import growth


# Existing carried-forward low-z calibration.  This normalizes the LOCAL
# topographic response in this exploratory module; it no longer calibrates
# the global n=2 mode at the same time.
Z_REF = 0.02
A_MU_REF = -0.0412
SIGMA_A_MU_REF = 0.009

# Illustrative local-coherence cutoff only.  It is NOT a fitted/frozen BHSM
# parameter.  A foreground catalog / action-derived source field must replace it.
Z_COHERENCE_ILLUSTRATIVE = 0.03

K = background.KCURV
SQRT_K = math.sqrt(K)


def E_of_z(z: float) -> float:
    if z < 0.0:
        raise ValueError("z must be nonnegative")
    N = math.log(1.0 / (1.0 + z))
    return float(growth.R1.sol(N)[2])


def chi_of_z(z: float) -> float:
    if z < 0.0:
        raise ValueError("z must be nonnegative")
    if z == 0.0:
        return 0.0
    value, _ = quad(
        lambda zp: 1.0 / E_of_z(zp),
        0.0,
        z,
        epsabs=2e-12,
        epsrel=2e-12,
        limit=200,
    )
    return float(value)


def C_over_S(chi: float, Kcurv: float = K) -> float:
    """Closed-S3 C_K/S_K = sqrt(K) cot(sqrt(K) chi)."""
    if chi <= 0.0:
        raise ValueError("chi must be positive")
    if Kcurv <= 0.0:
        # Flat-limit continuation used only as a numerical guard.
        return 1.0 / chi
    root = math.sqrt(Kcurv)
    x = root * chi
    return root * math.cos(x) / math.sin(x)


def localized_accumulation(
    chi_source: float,
    chi_coherence: float,
) -> float:
    """Constant effective LOS source inside one coherent local region."""
    if chi_source < 0.0 or chi_coherence <= 0.0:
        raise ValueError("invalid comoving distances")
    return min(chi_source, chi_coherence)


def normalized_local_transfer(
    z: float,
    z_ref: float = Z_REF,
    z_coherence: float = Z_COHERENCE_ILLUSTRATIVE,
) -> float:
    """Normalized closed-S3 radial response for a compact local source.

    The effective radial perturbation accumulates only while the ray remains
    in the coherent source region.  Once it exits, delta-chi saturates while
    C_K/S_K continues to decrease with source distance.
    """
    if z <= 0.0 or z_ref <= 0.0 or z_coherence <= 0.0:
        raise ValueError("redshifts must be positive")

    chi = chi_of_z(z)
    chi_ref = chi_of_z(z_ref)
    chi_c = chi_of_z(z_coherence)

    raw = C_over_S(chi) * localized_accumulation(chi, chi_c)
    raw_ref = (
        C_over_S(chi_ref)
        * localized_accumulation(chi_ref, chi_c)
    )
    return float(raw / raw_ref)


def magnitude_dipole_local(
    z: float,
    A_mu_ref: float = A_MU_REF,
    z_ref: float = Z_REF,
    z_coherence: float = Z_COHERENCE_ILLUSTRATIVE,
) -> float:
    return float(
        A_mu_ref
        * normalized_local_transfer(
            z,
            z_ref=z_ref,
            z_coherence=z_coherence,
        )
    )


def path_average_rms_ratio(
    chi: float,
    correlation_length: float,
) -> float:
    """RMS of an exponentially correlated path average / point RMS.

    C(Delta chi)=sigma^2 exp(-|Delta chi|/L_c).
    """
    if chi <= 0.0 or correlation_length <= 0.0:
        raise ValueError("positive chi and correlation_length required")
    L = correlation_length
    variance_ratio = (
        2.0 / (chi * chi)
        * (
            L * chi
            - L * L * (1.0 - math.exp(-chi / L))
        )
    )
    return float(math.sqrt(max(variance_ratio, 0.0)))


def seam_radius(
    GM: float,
    R_H: float,
    xi_c: float = 1.0,
    c: float = 1.0,
) -> float:
    """Leading spherical BHSM seam radius candidate."""
    if GM <= 0.0 or R_H <= 0.0 or xi_c <= 0.0 or c <= 0.0:
        raise ValueError("all inputs must be positive")
    return math.sqrt(GM * R_H / (xi_c * c * c))


def seam_charge_speed2(GM: float, r_star: float) -> float:
    """V_Sigma^2 = GM/r_star in the spherical weak-field candidate."""
    if GM <= 0.0 or r_star <= 0.0:
        raise ValueError("GM and r_star must be positive")
    return GM / r_star


def exterior_topographic_acceleration_magnitude(
    r: float,
    r_star: float,
    V_sigma2: float,
) -> float:
    """Stored leading exterior neighborhood-flow magnitude.

    Sign/orientation is deliberately not fixed here because the retained
    source-export convention and g_T=-Q F_T convention still require a
    common action-level orientation choice.
    """
    if r < r_star or r_star <= 0.0 or V_sigma2 < 0.0:
        raise ValueError("requires r >= r_star > 0 and V_sigma2 >= 0")
    return V_sigma2 * (1.0 / r - r_star / (r * r))


def diagnostics() -> dict:
    z_rows = (0.02, 0.10, 0.30, 0.50, 0.80, 1.00, 1.50, 2.10)
    rows = []
    for z in z_rows:
        rows.append(
            {
                "z": z,
                "chi_H0_units": chi_of_z(z),
                "normalized_local_transfer": normalized_local_transfer(z),
                "A_mu_local_mag": magnitude_dipole_local(z),
            }
        )

    # Dimensionless algebraic seam witness; not a physical mass choice.
    GM = 2.0
    RH = 10.0
    xi = 1.0
    c = 1.0
    rs = seam_radius(GM, RH, xi, c)
    V2 = seam_charge_speed2(GM, rs)

    return {
        "schema": "bhsm-los-topographic-transfer-candidate-v1",
        "status": "EXPLORATORY_SOURCE_EXPORT_AND_LOS_TRANSFER_CANDIDATE",
        "geometry": {
            "parent": "S^3(R_H)",
            "local_metric_ansatz": "gamma_ij=a^2 exp(2 psi) gamma_bar_ij",
            "local_volume_rate": "H_local=H+dot(psi)",
            "linear_curvature_response": (
                "delta R3=-(4/a^2)(bar_Laplacian+3/R_H^2) psi"
            ),
            "n2_shifted_factor": 5,
        },
        "bhsm_source_provenance": {
            "static_flux": "F_T=grad(T)-B grad(Laplacian(T))",
            "static_operator": "div(F_T)=Laplacian(T)-B Laplacian^2(T)=S",
            "candidate_acceleration": "g_T=-Q F_T",
            "seam_charge": (
                "V_Sigma^2=(c^2/8pi) integral_Sigma "
                "H_Sigma DeltaK_uu dA"
            ),
            "spherical_weak_field": "V_Sigma^2=GM/r_star",
            "neighborhood_profile": (
                "|g_T|=V_Sigma^2(1/r-r_star/r^2), r>=r_star"
            ),
            "source_export_orientation": "OPEN_ACTION_LEVEL_SIGN_CONVENTION",
            "psi_source_normalization": "OPEN_ACTION_ATTACHMENT",
        },
        "calibration_policy": {
            "low_z_SN_normalizes": "LOCAL_TOPOGRAPHIC_RESPONSE_ONLY",
            "z_ref": Z_REF,
            "A_mu_ref": A_MU_REF,
            "A_mu_ref_sigma": SIGMA_A_MU_REF,
            "global_n2_amplitude": (
                "NOT_SIMULTANEOUSLY_INFERRED_FROM_THE_SAME_SCALAR_CALIBRATION"
            ),
        },
        "illustrative_coherence_proxy": {
            "z_coherence": Z_COHERENCE_ILLUSTRATIVE,
            "status": (
                "ILLUSTRATIVE_NOT_FITTED_NOT_FROZEN;"
                " replace with foreground/source map"
            ),
        },
        "rows": rows,
        "seam_identity_witness": {
            "GM": GM,
            "R_H": RH,
            "xi_c": xi,
            "r_star": rs,
            "V_sigma2": V2,
            "V_sigma4": V2 * V2,
            "expected_V_sigma4": xi * GM * c * c / RH,
        },
        "claim_boundary": {
            "derived_or_reused": [
                "closed-S3 radial C_K/S_K geometry",
                "localized-source saturation mechanism",
                "BHSM seam charge candidate",
                "BHSM spherical neighborhood-flow shape",
                "finite-correlation path decorrelation formula",
            ],
            "not_yet_claimed": [
                "action-normalized de-encapsulation-to-psi source map",
                "foreground-catalog prediction",
                "SMBH-specific source law",
                "unique global-n2/local-topography amplitude split",
                "formal BHSM cosmology closure",
            ],
        },
    }


def write_tex(report: dict, path: Path) -> None:
    rows = report["rows"]
    lines = [
        r"\section{Local BHSM topography on the hyperspherical parent}",
        r"\label{sec:bhsm-los-topography}",
        "",
        (
            "The smooth $S^3(R_H)$ background and the physical $n=2$ "
            "response do not exhaust the geometry sampled by a photon. "
            "We therefore separate the global parent from a local "
            "deformation field $\\psi(\\mathbf x,t)$:"
        ),
        r"\begin{equation}",
        r"\gamma_{ij}=a^2 e^{2\psi}\bar\gamma^{S^3}_{ij}.",
        r"\end{equation}",
        (
            "To first order its intrinsic-curvature perturbation is"
        ),
        r"\begin{equation}",
        r"\delta{}^{(3)}R=-\frac{4}{a^2}"
        r"\left(\bar\nabla^2+\frac{3}{R_H^2}\right)\psi,",
        r"\end{equation}",
        (
            "so the physical manuscript $n=2$ harmonic carries the "
            "shifted factor $8-3=5$.  The local volume rate is"
        ),
        r"\begin{equation}",
        r"H_{\rm local}=H+\dot\psi.",
        r"\end{equation}",
        "",
        (
            "The retained BHSM topographic flux candidate is "
            r"$F_T=\nabla T-B\nabla(\nabla^2T)$ with "
            r"$\nabla\cdot F_T=S$.  The existing seam-charge candidate is"
        ),
        r"\begin{equation}",
        r"V_\Sigma^2=\frac{c^2}{8\pi}"
        r"\int_{\Sigma_\star}H_\Sigma\Delta K_{uu}\,dA,",
        r"\end{equation}",
        (
            r"giving $V_\Sigma^2=GM/r_\star$ at leading spherical "
            "weak-field order.  The exact action-level export from this "
            "boundary charge into the dynamical $\\psi$ source remains "
            "open; no normalization is invented here."
        ),
        "",
        r"\subsection{Localized line-of-sight response}",
        (
            "For a minimal coherent local-source proxy the accumulated "
            "radial perturbation grows only to a coherence distance "
            r"$\chi_c$ and then saturates.  The closed-$S^3$ radial "
            "distance term therefore has the normalized shape"
        ),
        r"\begin{equation}",
        r"\mathcal F_{\rm loc}(z)="
        r"\frac{[C_K/S_K](\chi)\min(\chi,\chi_c)}"
        r"{[C_K/S_K](\chi_\star)\min(\chi_\star,\chi_c)}.",
        r"\end{equation}",
        (
            "The numerical example below uses $z_c=0.03$ only as an "
            "illustrative local-coherence proxy, not as a fitted or "
            "frozen model parameter."
        ),
        "",
        r"\begin{center}",
        r"\begin{tabular}{r r}",
        r"$z$ & $A_\mu^{\rm loc}$ [mag]\\",
        r"\hline",
    ]
    for row in rows:
        lines.append(f"{row['z']:.2f} & {row['A_mu_local_mag']:.6f}\\\\")
    lines += [
        r"\end{tabular}",
        r"\end{center}",
        "",
        r"\paragraph{Calibration boundary.}",
        (
            "The carried-forward low-redshift SN scalar amplitude "
            "normalizes the local topographic response in this exploratory "
            "decomposition.  It cannot simultaneously determine an "
            "independent global $n=2$ amplitude.  The global/local split "
            "must ultimately be fixed by the common BHSM action or an "
            "independent preregistered observable."
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path)
    parser.add_argument("--tex", type=Path)
    args = parser.parse_args()

    report = diagnostics()
    if args.artifact is not None:
        args.artifact.parent.mkdir(parents=True, exist_ok=True)
        args.artifact.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.tex is not None:
        args.tex.parent.mkdir(parents=True, exist_ok=True)
        write_tex(report, args.tex)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
