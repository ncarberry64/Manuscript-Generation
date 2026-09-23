# R1 physical-n=2 reduced propagator

This checkpoint evaluates the first actual numerical phase-space transfer
kernel on the representation-invariant prediction track.

State:

    X = (zeta, d zeta/dN)

Anchor:

    z = 2.1, U = I

Equation:

    zeta_NN
    + (3 + d ln H/dN + d ln G_S,2/dN) zeta_N
    + [5 K/(a^2 H^2)] (F_S,2/G_S,2) zeta
    = 0

The factor 5 comes from the physical manuscript n=2 / literature N=3
closed-S3 harmonic:

    -D^2 = 8 K,
    -(D^2+3K) = 5 K.

The code also propagates the exact closed-Horndeski lapse/shift constraint
response for both columns of the fundamental matrix and checks those
constraint residuals directly.

Claim boundary:
the current F_S routine is explicitly gravity+scalar only. This is therefore
a reduced physical-n=2 phase-space propagator, not yet the complete
matter+radiation directional observable kernel.
