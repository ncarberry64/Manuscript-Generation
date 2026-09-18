# Scientific referee audit â€” Round 3
## Exact closed-Horndeski constraint kernel for the physical n=2 mode

Date: 2026-09-18

### Literature basis

Akama & Kobayashi, Phys. Rev. D 99, 043522 (2019), arXiv:1810.01863,
derive the general quadratic action for scalar perturbations in open and closed
FLRW backgrounds in Horndeski theory.

For positive curvature their propagating scalar harmonics begin at their
index N=3. This is manuscript n=2 because N=n+1.

### Specialization

For G2=X-V, G3=2X/Lambda^3, G4=M_Pl^2/2, G5=0:

- F_T = G_T = M_Pl^2
- Theta_K = 0
- Sigma_K = 0
- Theta = H M_Pl^2 - dot(T)^3/Lambda^3
- Sigma = dot(T)^2/2 + 12 H dot(T)^3/Lambda^3 - 3 M_Pl^2 H^2

For manuscript n=2 / literature N=3, k^2=8K.

The lapse and shift constraints are solved exactly in
`code/closed_horndeski_n2.py`.

### Exact kinetic identity

Let r=M_Pl^2 Sigma/Theta^2. Then

`G_S,2=M_Pl^2*5(r+3)/(r+8)`.

The alpha-function dictionary gives

`r+3=D_kin/[2(1-alpha_B/2)^2]`.

Therefore

`G_S,2=M_Pl^2*5D_kin/[D_kin+10(1-alpha_B/2)^2]`.

Thus D_kin>0 implies a positive gravity-scalar kinetic coefficient for the
first physical closed scalar harmonic, provided Theta does not cross zero.

### Scope limit

This is not yet the full late-time transfer matrix:
- R1 contains matter and radiation.
- Matter/radiation perturbations add dynamical/source variables.
- The reduced topographic p^2 response is an additional bridge operator.
- Luminosity distance requires a gauge-invariant light-cone projection.

The next step is to add matter/radiation perturbations to this exact closed
gravity-scalar kernel before converting to Newtonian potentials and observables.
