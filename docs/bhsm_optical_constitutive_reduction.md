# BHSM optical constitutive reduction

This checkpoint reduces the two open maps from the optical-transfer theorem.

## Seam export

The retained stationary spherical neighborhood law already fixes the static
source profile:

    E0[Q_Sigma]
      = epsilon_Sigma (Q_Sigma/Q)
        Theta(r-r_star)/r^2.

Only the common-action orientation/normalization remains open.

If de-encapsulation is represented by the time variation of that exported
source, a moving seam gives

    E1
      = (epsilon_Sigma/Q)
        [ Qdot_Sigma Theta/r^2
          - Q_Sigma rdot_star
            delta(r-r_star)/r_star^2 ].

The second term is the seam-localized pulse produced purely by boundary motion.

## Local metric response

For physical scalar metric coordinates `(Phi,Psi)` with

    A_g = [[a,b],[b,c]]
    B_gpsi = [u,v]^T,

the response is exactly

    R_gpsi
      = -A_g^-1 B_gpsi
      = [(b v-c u)/Delta,
         (b u-a v)/Delta]^T,

    Delta = a c-b^2.

The Weyl/lensing coefficient is

    R_W
      = [u(b-c)+v(b-a)]
        / [2(ac-b^2)].

So the local metric response is no longer an arbitrary function once the five
projected action matrix elements are evaluated.

## Optical covariance is upstream-owned

The unresolved optical covariance is

    Q_opt
      = O_opt R_gpsi G_psi
        C_J
        G_psi^dagger R_gpsi^dagger O_opt^dagger.

This means the optical covariance should not be fit as 28 unrelated numbers.
Its rank, correlations, redshift evolution and cross-observable structure are
inherited from the much smaller BHSM source covariance.

For rare seam/de-encapsulation events with rate `lambda_Sigma`,

    C_Sigma^shot
      = lambda_Sigma <q_Sigma q_Sigma^dagger>.

That is the natural next statistical closure target.
