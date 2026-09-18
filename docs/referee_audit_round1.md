
# Scientific referee audit â€” Round 1

Date: 2026-09-18

Scope: claim strength, closed-FLRW harmonic consistency, inherited Horndeski/KGB stability language, observational provenance, and current-data framing.

This audit intentionally leaves the frozen R1 numerical manifest unchanged.

## Major issue: exceptional first nonconstant scalar harmonic

The manuscript uses scalar eigenvalues `n(n+2)`. Closed-universe perturbation literature often uses `N^2-1`, so `N=n+1`; manuscript `n=1` maps to literature `N=2`.

Kiefer & Vardanyan's gauge-invariant closed-FLRW Hamiltonian contains factors `1/(N^2-4 K)`, singular for `N=2` at positive unit curvature, and their harmonic appendix identifies the lowest modes as exceptional in the usual scalar-derived tensor-harmonic construction.

Action taken:
- preserve the exact S3 scalar monopole+dipole geometry;
- stop treating that geometry alone as proof of a physical metric/light-cone dipole;
- make all n=1 metric observables conditional on a dedicated constrained gauge-invariant derivation;
- add absence of such a physical metric response as an explicit falsification route.

Source:
- Claus Kiefer & Tatevik Vardanyan, *General Relativity and Gravitation* 54, 30 (2022), arXiv:2111.07835, DOI 10.1007/s10714-022-02918-3.

## Stability language

Bellini & Sawicki's `D = alpha_K + 3 alpha_B^2/2` is an appropriate inherited no-ghost/kinetic diagnostic for the Horndeski sector, but their standard perturbation formulation assumes spatial flatness.

Action taken:
- retain D_kin;
- remove any implication that D_kin alone proves full stability of the closed model;
- distinguish positivity of the reduced spatial kernel from complete kinetic/gradient stability.

Source:
- Emilio Bellini & Ignacy Sawicki, JCAP 07 (2014) 050, arXiv:1404.3713.

## Kinetic-gravity-braiding background formulas

The energy density, pressure, shift current, and scalar-metric braiding interpretation are supported by the KGB literature.

Action taken:
- add the primary KGB citation next to those equations.

Source:
- Cedric Deffayet, Oriol Pujolas, Ignacy Sawicki, Alexander Vikman, JCAP 10 (2010) 026, arXiv:1008.0048.

## Harmonic-level versus sky-quadrupole claim

The compact S3 n=2 sector is not synonymous with the observed sky l=2 sector.

Action taken:
- rename the former "quadrupole null" as an n=2 harmonic-level response bound;
- require observer-shell projection/source weights before turning it into a sky-quadrupole prediction.

## Radial coordinate consistency

The old draft used chi both as a dimensionless S3 angle and as a dimensional comoving distance.

Action taken:
- use dimensionless `psi` in the S3 metric;
- define `D_C(z)` as dimensional comoving radial distance and `psi(z)=D_C/R_H`.

## Low-z supernova calibration provenance

Action taken:
- mark the numerical low-z amplitude and frozen direction as carried-forward inputs from the earlier analysis;
- cite Pantheon+ dataset/cosmology papers;
- state that the current repository does not rerun that likelihood.

## Curvature forecast versus current data

DESI DR2 reports its BAO results as well described by flat LambdaCDM while also discussing evidence for extensions when additional datasets are combined.

Action taken:
- retain `Omega_k=-0.018 Â± 0.004` unchanged as a frozen prospective model forecast;
- explicitly state it is not presented as a current empirical best fit.

Source:
- DESI Collaboration, *Physical Review D* 112, 083515 (2025), arXiv:2503.14738.

## Next referee calculation

The highest-priority derivation after this claims audit is the exceptional closed-FLRW scalar-metric problem:
1. construct the constrained gauge-invariant physical variables at the manuscript n=1 / literature N=2 level;
2. determine whether an independent physical scalar-metric degree of freedom survives after constraints/gauge quotient when the topographic field is included;
3. derive its kinetic and gradient matrices;
4. derive the luminosity-distance and velocity source terms;
5. only then promote the scalar sky pattern into a metric/light-cone prediction.

MDPI template conversion is intentionally deferred to a separate formatting pass so it cannot obscure the scientific diff.
