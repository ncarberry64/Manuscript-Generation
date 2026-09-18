# Scientific referee audit â€” Round 4
## Matter/radiation sourced constraint transfer for physical n=2

Date: 2026-09-18

### Source equations

Gambino & Pace (arXiv:2412.01781) give curvature-aware Newtonian-gauge
linear EFT equations with explicit matter density, momentum, pressure and
anisotropic-stress sources.

For the manuscript subclass:
- alpha_M = 0
- alpha_T = 0
- alpha_H = 0
- M^2 = M_Pl^2

and ideal sources with zero scalar anisotropic stress, their traceless spatial
equation gives Phi=Psi.

For manuscript n=2:
- D^2 -> -8 K
- D^2+3K -> -5 K

The sourced 00 and 0i equations become a 2x2 system for
(Phi_dot, pi_dot).

### Determinant result

`det C = -2 H^2 (alpha_K + 6 alpha_B^2)`

and therefore

`det C = -2 H^2 (D_kin + 9 alpha_B^2/2)`.

On the positive D_kin branch this local sourced constraint block is
nonsingular.

### Matter/radiation closure

For `delta T^0_i=(rho+p)D_i v`, ideal covariant conservation gives

dust:
- delta_m_dot = (k^2/a^2) v_m + 3 Phi_dot
- v_m_dot = -Phi

ideal radiation:
- delta_r_dot = (4/3)(k^2/a^2) v_r + 4 Phi_dot
- v_r_dot = H v_r - Phi - delta_r/4

with k^2/a^2=8K/a^2.

### Compatibility gate

This round intentionally does NOT integrate the global matter-era system.

The Round-3 Akama-Kobayashi ADM kernel and the Gambino-Pace EFT
Newtonian-gauge equations use different perturbation variables and
normalizations. Their local sourced constraints can be specialized and checked
directly, but the ij-trace/scalar evolution equations must be analytically
matched and then shown to preserve the 00/0i constraints under time evolution.

No topographic reduced-response term is included in this round.

### Next calculation

1. derive the exact variable/convention dictionary between the Round-3 ADM
   variables and the curvature-aware EFT Newtonian-gauge variables;
2. derive the matched ij-trace/scalar evolution block;
3. verify the Bianchi identity and numerical constraint preservation;
4. only then integrate the global n=2 matter/radiation system;
5. add the topographic response as a separate final operator.
