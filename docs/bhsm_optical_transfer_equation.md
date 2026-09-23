# BHSM optical transfer equation

## Purpose

Replace object-by-object foreground accounting with propagation of the
measurable photon-bundle state through the total geometry.

## Optical state

    Y_gamma =
      (zeta_gamma,
       kappa, p_kappa,
       gamma_1, p_gamma_1,
       gamma_2, p_gamma_2)^T

where `zeta_gamma = delta ln(1+z)` and primes denote the closed-S3 comoving
radial parameter used by the optical system.

The exact first-order closed-background equation is

    Y_gamma' = A_K(chi) Y_gamma + B_opt S_opt

with

    S_opt = (Z, F, G1, G2)^T.

`Z` is the redshift/energy forcing; `F` is the trace part of the optical tidal
perturbation; `G1,G2` are its trace-free shear components.

## BHSM closure

The optical forcing is sourced by the same metric response already used in
the cosmology bridge:

    h_2 = R_g2 X_2

and by a sourced local response

    h_topo = R_gpsi psi.

The local topographic field is written

    L_psi psi = J_grav + E_Sigma[Q_Sigma] + xi

with the existing seam charge

    Q_Sigma = V_Sigma^2
            = (c^2/8pi) integral H_Sigma DeltaK_uu dA.

The action normalization/sign of `E_Sigma` and the complete `R_gpsi` remain
open. They are now the constitutive closure problem; the optical propagation
itself is no longer open.

## Unresolved terrain

For unresolved optical forcing covariance `Q_opt`,

    C_Y' = A_K C_Y + C_Y A_K^T + B_opt Q_opt B_opt^T.

This is the key computational compression: foreground complexity is represented
by a 7x7 propagated photon-state covariance rather than an exhaustive list of
objects.

Rare nonlinear/de-encapsulation encounters can be modeled separately as a
compound jump process without changing the deterministic closed-S3 drift.

## Observable map

At fixed physical source endpoint,

    delta ln D_L = 2 zeta_gamma - kappa.

For a Hubble diagram at fixed observed redshift, the accumulated optical state
must then be passed through the existing closed-S3 gauge-invariant source
remapping,

    delta D_L/D_L
      = phi_s + delta z + (C_K/S_K) delta chi - kappa.

Thus this module supplies the physical path state; it does not replace the
existing luminosity-distance kernel.
