#!/usr/bin/env python3
"""Generate the final-scrutiny scientific figures for Geometry Before Fields.

Inputs:
- exact analytic n=2 geometry;
- artifacts/provenance/BHSM_R1_COUPLED_ENVIRONMENTAL_STATE_AUDIT_V1.json;
- artifacts/R1_ACTION_NATIVE_GRADIENT_DISPERSION_CONVERGENCE_V1.json.

No observational fit, parameter retuning, branch selection, or source normalization
is performed here.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures"
COUPLED = ROOT / "artifacts" / "provenance" / "BHSM_R1_COUPLED_ENVIRONMENTAL_STATE_AUDIT_V1.json"
STABILITY = ROOT / "artifacts" / "R1_ACTION_NATIVE_GRADIENT_DISPERSION_CONVERGENCE_V1.json"

def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

def n2_shell_geometry():
    psi = np.linspace(0.0, np.pi, 800)
    c, s = np.cos(psi), np.sin(psi)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(psi/np.pi, c*c-s*s/3.0, label=r"$A_{44}$ monopole radial factor")
    ax.plot(psi/np.pi, np.sin(2.0*psi), label=r"$A_{4i}$ dipole radial factor $\sin(2\psi)$")
    ax.plot(psi/np.pi, s*s, label=r"$S_{ij}$ quadrupole radial factor")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel(r"Geodesic angle $\psi/\pi$")
    ax.set_ylabel("Normalized radial factor")
    ax.set_title(r"Exact observer-shell factors of the physical $n=2$ scalar harmonic")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "n2_observer_shell_radial_basis")

def n2_sky_dipole():
    lon = np.linspace(-np.pi, np.pi, 361)
    lat = np.linspace(-np.pi/2, np.pi/2, 181)
    LON, LAT = np.meshgrid(lon, lat)
    mu = np.cos(LAT) * np.cos(LON)
    fig = plt.figure(figsize=(7.2, 4.2))
    ax = fig.add_subplot(111, projection="mollweide")
    mesh = ax.pcolormesh(LON, LAT, mu, shading="auto", rasterized=True)
    ax.grid(True, linewidth=0.5)
    ax.set_title(r"Pure $A_{4i}$ observer-sky factor $(\hat n\cdot\hat p)$")
    cb = fig.colorbar(mesh, ax=ax, orientation="horizontal", pad=0.08, fraction=0.06)
    cb.set_label(r"$\mu=\hat n\cdot\hat p$")
    save(fig, "n2_observer_sky_dipole")

def coupled_environment():
    if not COUPLED.exists():
        raise FileNotFoundError(
            f"Missing {COUPLED}. Bring provenance commit e32b13ad86d1b66889a4a0483fed9ef20ae84561 "
            "onto the scrutiny branch before generating figures."
        )
    obj = json.loads(COUPLED.read_text(encoding="utf-8"))
    rows = obj["temporal_result"]["rows"]
    z = np.array([float(r["z"]) for r in rows])
    s1 = np.array([float(r["s1"]) for r in rows])
    s2 = np.array([float(r["s2"]) for r in rows])
    det = np.array([float(r["matter_det"]) for r in rows])

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(z, s1, marker="o", label=r"$s_1(U_{XE})$")
    ax.plot(z, s2, marker="o", label=r"$s_2(U_{XE})$")
    ax.set_yscale("log")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Singular value")
    ax.set_title(r"Environmental $\rightarrow$ topographic temporal transfer is rank two")
    ax.invert_xaxis()
    ax.legend(frameon=False)
    save(fig, "coupled_environment_singular_values")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(z, det, marker="o")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\det U_{XE}^{(\delta_m,v_m)}$")
    ax.set_title("Matter density and velocity alone span both temporal components")
    ax.invert_xaxis()
    save(fig, "coupled_environment_matter_determinant")

def principal_stability():
    obj = json.loads(STABILITY.read_text(encoding="utf-8"))
    summaries = obj["principal_symbol"]["z_summaries"]
    z = np.array([float(r["z"]) for r in summaries])
    radiation = np.array([float(r["radiation_c2_intercept"]) for r in summaries])
    scalar = np.array([float(r["scalar_c2_intercept"]) for r in summaries])

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(z, scalar, marker="o", label=r"non-radiation propagating $c^2$")
    ax.plot(z, radiation, marker="o", label=r"radiation branch $c^2$")
    ax.axhline(1.0/3.0, linewidth=0.8, linestyle="--", label=r"$1/3$")
    ax.axhline(0.0, linewidth=0.8)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"High-$k$ characteristic speed squared")
    ax.set_title("Frozen R1 principal-symbol characteristic speeds")
    ax.legend(frameon=False, fontsize=9)
    save(fig, "r1_principal_characteristic_speeds")

    # Use the two endpoint redshifts to show the positive high-k growth tail.
    endpoints = [summaries[0], summaries[-1]]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for row in endpoints:
        seq = row["high_n_sequence"]
        n = np.array([int(r["n"]) for r in seq])
        g = np.array([float(r["max_positive_normalized_growth"]) for r in seq])
        ax.plot(n, g, marker="o", label=rf"$z={float(row['z']):g}$")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel("Harmonic index n")
    ax.set_ylabel(r"max positive $\mathrm{Re}(s)/\sqrt{k_{\rm phys}^2}$")
    ax.set_title("High-k positive-growth tail decreases with harmonic scale")
    ax.legend(frameon=False)
    save(fig, "r1_highk_growth_convergence")

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    n2_shell_geometry()
    n2_sky_dipole()
    coupled_environment()
    principal_stability()
    print(f"Final scrutiny figures written to {OUT}")

if __name__ == "__main__":
    main()
