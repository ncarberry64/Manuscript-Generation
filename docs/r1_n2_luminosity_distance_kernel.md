# R1 physical-n=2 gauge-invariant luminosity-distance kernel

This checkpoint evaluates the luminosity-distance response as a row vector on
the same two-dimensional state space used by the reduced n=2 propagator.

Flat-FLRW gauge-invariant geometric structure:

    delta D_L / D_L
      = phi_s + delta z + delta r/r - kappa

is extended to the closed background by

    r -> S_K(chi)=sin(sqrt(K)chi)/sqrt(K)

so that

    delta D_L / D_L
      = phi_s + delta z
        + [C_K(chi_s)/S_K(chi_s)] delta chi
        - kappa.

The pure observer-sky dipole has

    Q_2,dip = sin(2 sqrt(K) chi) mu
    hat Laplacian Q_2,dip = -2 Q_2,dip.

The redshift, radial and convergence distortions are all evaluated from the
same gauge-invariant Psi/Phi response and transformed matter velocity.

## Calibration-rank result

At z*=0.02 the SN distance calibration gives one scalar equation

    F_DL(z*) . X_anchor = observed distance dipole.

But X_anchor is two-dimensional. The code explicitly constructs the null
direction of F_DL and checks that the growth and Weyl response rows are not
all parallel to F_DL.

Therefore the SN amplitude alone does not select the temporal phase/state
direction. The single-calibration theorem is valid only after the normalized
mode trajectory Xhat_2 is fixed independently.

No second observational fit is authorized. The next closure object is a
theory-selected or explicitly preregistered branch direction Xhat_2.
