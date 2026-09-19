# Scientific referee audit ÃƒÆ’Ã†â€™Ãƒâ€šÃ‚Â¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡Ãƒâ€šÃ‚Â¬ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€šÃ‚Â Round 5
## ADM/EFT dictionary and DAE constraint preservation

Date: 2026-09-18

### Exact dictionary

`pi=-chi`

`Phi=delta_n-pi_dot=delta_n+chi_dot`

`Psi=-zeta+H*pi=-zeta-H*chi`

Inverse:

`chi=-pi`

`delta_n=Phi+pi_dot`

`zeta=H*pi-Psi`

### Audit correction

The original Round-5 test was structurally invalid.

The 00 and 0i field equations are constraints. Solving them pointwise for
`(Phi_dot, pi_dot)` and then demanding that both the ij-trace and scalar
equations vanish for an arbitrary state incorrectly promotes the constraints
to the complete evolution system.

Correct DAE formulation:

1. use the ij-trace and scalar equations to solve for
   `(Phi_ddot, pi_ddot)`;
2. evolve `(Phi, pi, Phi_dot, pi_dot)` plus ideal dust/radiation;
3. choose initial rates that satisfy 00 and 0i;
4. monitor preservation of 00 and 0i during evolution.

### Exact acceleration determinant

For alpha_M=alpha_T=alpha_H=0 and Phi=Psi:

`A = [[2, -2H alpha_B], [6H alpha_B, H^2 alpha_K]]`.

Using the GP/Gleyzes braiding normalization,

`det A = 2 H^2 (alpha_K + 6 alpha_B^2)`
`      = 2 H^2 D_kin`.

Therefore the acceleration block is nonsingular on the established
positive-D_kin branch.

### Source transcription status

The audit uses the equations as printed in Gambino-Pace:
- Eq. (8): `(6-alpha_K+12 alpha_B) H^2 Phi`;
- Eq. (10): trace equation;
- Eq. (12): scalar equation with
  `H dot(alpha_B)[H^2(3+alpha_M)+Hdot]`
  and `H d/dt[H dot(alpha_B)]`.

The factor-of-two GP/BS braiding convention bridge remains explicit.

### Certification criterion

Starting at z=2.1 with small generic linear amplitudes, the initial
`Phi_dot` and `pi_dot` are chosen from the sourced 00/0i constraints.

The trace/scalar equations are then integrated to z=0.

Certification requires both normalized 00 and 0i residuals to remain below
`2e-5` over the sampled interval. This threshold is not relaxed if the test
fails.

Topographic response remains OFF.

### Braiding sign correction

A direct convention audit against the EFT literature showed that the
Gleyzes/Gambino--Pace braiding parameter differs from the Bellini--Sawicki
parameter by a factor **-2**, not merely +2:

`alpha_B^BS = -2 alpha_B^GP`

and therefore

`alpha_B^GP = - alpha_B^BS / 2`.

The kinetic invariant is unchanged because the braiding contribution is
quadratic:

`alpha_K + 6(alpha_B^GP)^2`
`= alpha_K + 3/2(alpha_B^BS)^2`
`= D_kin`.

All curvature-aware GP equations, including every time derivative of
`alpha_B`, are rerun with the corrected sign. No numerical threshold is
relaxed.

### Radiation velocity-potential correction

The Gambino--Pace matter convention is

`delta T^0_i = (rho+p) partial_i v`.

For this covariant velocity potential, the ideal-fluid Euler equation is

`v_dot + 3H(c_s^2-w)v = -Phi - c_s^2 delta/(1+w)`.

Therefore:
- dust: `v_m_dot = -Phi`;
- ideal radiation (`w=c_s^2=1/3`):
  `v_r_dot = -Phi - delta_r/4`.

The earlier Round-4 code had incorrectly included `+H v_r`, which belongs to
a different velocity normalization. That error violates total
stress-energy conservation and therefore necessarily spoils the Bianchi
constraint-preservation test.

No Bianchi threshold is changed.

### Index-consistent initial-data correction

The previous DAE runs imposed the sourced 00 and 0i constraints at the
starting point but did not impose the first DAE consistency conditions.

For a linear constraint system

`C = B(N)y = 0`

and evolution

`y_N = F(N)y`,

a valid initial state must also satisfy

`dC/dN = [B_N + B F] y = 0`.

The former arbitrary seed lay on the constraint surface but was not tangent to
it. Rapid departure from the surface was therefore expected and was not a
physical-instability result.

The revised diagnostic seed satisfies:
- 00 and 0i;
- their first tangency conditions;
- adiabatic density relation `delta_r=4 delta_m/3`;
- common initial velocity potential `v_r=v_m`;
- `pi=0`;
- a harmless normalization `Phi=-1e-5`.

The normalization is not an observational fit.

The radiation Euler equation is restored to the covariant GP velocity
potential convention:
`v_r_dot = H v_r - Phi - delta_r/4`.

### Exact GP Eq. (12) source correction

A literal audit against Gambino & Pace, arXiv:2412.01781v1, Eq. (12)
identified two coupled transcription errors in the coefficient multiplying
the StÃƒÂ¼ckelberg perturbation `pi`.

The source contains

`Hdot * alpha_B * [H^2(3+alpha_M)+Hdot]`

and

`H * d/dt(Hdot * alpha_B)`.

For `alpha_M=0` these are

`Hdot * alpha_B * (3 H^2 + Hdot)`

and

`H * d/dt(Hdot * alpha_B)`.

The previous Round-5 implementation instead used

`H * alpha_B_dot * (3 H^2 + Hdot)`

and

`H * d/dt(H * alpha_B_dot)`.

Those are not algebraically equivalent. The implementation is corrected
literally to the published equation before any further dynamical
interpretation. No physical parameter, preregistered prediction, or Bianchi
threshold is changed. Topographic response remains off.

The stale module prose for the convention bridge is also corrected to
`alpha_B^GP = -alpha_B^BS / 2`; the numerical bridge had already been
changed earlier.

### Published Eq. (12) literal-form restoration

A later direct read of the arXiv PDF showed that the prior attempted
"Eq. (12) correction" had itself misread the source.

The published v1 equation contains

`H * alpha_B_dot * [H^2(3+alpha_M)+Hdot]`

and

`H * d/dt(H * alpha_B_dot)`.

The temporary replacement by

`Hdot * alpha_B * [...]`

and

`H * d/dt(Hdot * alpha_B)`

is therefore reverted. The repository is restored to the literal published
Gambino--Pace Eq. (12) before further curvature diagnostics.

No scientific parameter or certification threshold is changed.
