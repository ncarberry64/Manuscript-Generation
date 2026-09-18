"""Generate manuscript figures from the canonical R1 implementation.

No observational comparison data are read by this script.
"""
from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import background
import growth
import perturbations

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "manuscript" / "figures"
OUT.mkdir(parents=True, exist_ok=True)


def save_alpha_b():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)
    zs = np.linspace(0.0, 10.0, 300)
    vals = [background.diagnostics(sol, v0, float(z))["alpha_B"] for z in zs]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\alpha_B$")
    ax.set_title("Reference Branch R1: late-time kinetic braiding")
    ax.set_yscale("log")
    fig.tight_layout()
    fig.savefig(OUT / "r1_alphaB.png", dpi=220)
    plt.close(fig)


def save_omega_t():
    v0 = background.shoot_v0()
    sol = background.integrate(v0)
    zs = np.linspace(0.0, 10.0, 300)
    vals = [background.diagnostics(sol, v0, float(z))["Omega_T"] for z in zs]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\Omega_T$")
    ax.set_title("Reference Branch R1: scalar energy fraction")
    fig.tight_layout()
    fig.savefig(OUT / "r1_omegaT.png", dpi=220)
    plt.close(fig)


def save_growth():
    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)
    zs = np.linspace(0.0, 2.1, 250)
    vals = []
    for z in zs:
        N = math.log(1.0 / (1.0 + float(z)))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        vals.append(float(dpr / dpl))

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(zs, vals)
    ax.axhline(1.0, linewidth=1.0)
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$f\sigma_{8,\mathrm{R1}}/f\sigma_{8,\Lambda\mathrm{CDM}}$")
    ax.set_title("R1 prospective growth ratio")
    fig.tight_layout()
    fig.savefig(OUT / "r1_growth_ratio.png", dpi=220)
    plt.close(fig)


def save_filter():
    ell = 0.32
    ns = np.arange(1, 11)
    f = np.array([1.0 - ell * ell * perturbations.nu(int(n)) for n in ns])

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(ns, f, marker="o")
    ax.axhline(0.0, linewidth=1.0)
    ax.set_xlabel("S3 harmonic index n")
    ax.set_ylabel(r"$1-\ell^2 n(n+2)$")
    ax.set_title(r"Physical harmonic filter sign for $\ell=0.32$")
    ax.set_xticks(ns)
    fig.tight_layout()
    fig.savefig(OUT / "harmonic_filter.png", dpi=220)
    plt.close(fig)


def main():
    save_alpha_b()
    save_omega_t()
    save_growth()
    save_filter()
    print(f"Figures written to {OUT}")


if __name__ == "__main__":
    main()
