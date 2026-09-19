# R1 physical-n=2 growth and Weyl transfer

This checkpoint couples the exact reduced n=2 mode to pressureless matter
using covariant conservation directly in unitary gauge.

It does not use the previously identified curvature-defective Newtonian-gauge
DAE.

Dust equations:

    dot(v_m) = -delta_n

    dot(delta_m) + 3 d(H v_m)/dt
      = -(3 dot(zeta) + k2/a2 chi) + k2/a2 v_m

with k2/a2 = 8 K/a2 for manuscript n=2.

Gauge-invariant metric potentials in the adopted convention:

    Psi = delta_n + dot(chi)
    Phi = zeta + H chi

so the local Weyl/lensing combination is

    Psi - Phi.

The code propagates both basis columns of the physical state and emits:

- dust response matrix;
- delta_m,N kernel;
- fractional f sigma_8 state-space response;
- local Weyl state-space response.

These remain two-component state-space functionals until the full
luminosity-distance light-cone calibration selects the one physical initial
state allowed by the single-calibration theorem.
