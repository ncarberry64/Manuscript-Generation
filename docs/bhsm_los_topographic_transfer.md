# BHSM line-of-sight topographic transfer candidate

## Result

The failed homogeneous soft-eigenline surrogate showed that state selection
alone cannot make the current R1 n=2 luminosity-distance response decay at
high redshift.  The next model therefore restores the physical distinction
between:

1. the global `S^3(R_H)` parent geometry;
2. the coherent physical n=2 response;
3. local BHSM topography encountered along each photon path.

Write

    gamma_ij = a^2 exp(2 psi) gamma_bar^S3_ij.

Then

    H_local = H + dot(psi)

and, to first order,

    delta R3 = -(4/a^2) (bar_Laplacian + 3/R_H^2) psi.

For manuscript n=2, `-bar_Laplacian -> 8/R_H^2`, so the shifted curvature
factor is exactly `8-3=5`.

## BHSM source provenance

Retained BHSM work already provides the candidate static topographic flux

    F_T = grad(T) - B grad(Laplacian(T)),
    div(F_T) = S,

the candidate acceleration `g_T=-Q F_T`, and the seam charge

    V_Sigma^2
      = (c^2/8pi) integral_Sigma H_Sigma Delta K_uu dA.

At leading spherical weak-field order,

    V_Sigma^2 = GM/r_star.

The stored exterior neighborhood-flow magnitude is

    |g_T(r)|
      = V_Sigma^2 (1/r - r_star/r^2),
      r >= r_star,

which activates continuously at the seam.

The exact common-action orientation/sign and the normalization that exports
this boundary charge into a dynamical cosmological `psi` source remain open.
This patch therefore preserves those items as explicit claim boundaries
instead of inventing a coupling.

For a time-dependent seam, the minimal de-encapsulation *shape* source is the
time variation of the same boundary charge distributed over its exported
topographic profile.  Schematically,

    J_de(x,t) propto d(V_Sigma^2)/dt * W_Sigma(x,t),

with the proportionality/field normalization to be fixed by the retained
action.  The overall normalization is not needed for the normalized
line-of-sight shape test below.

## Localized LOS consequence

A local coherent source contributes to the radial light-cone distortion only
while the ray is inside that source region.  If its accumulated `delta chi`
saturates after `chi_c`, the closed-S3 distance factor gives

    F_loc(z)
      = [(C_K/S_K)(chi) min(chi,chi_c)]
        / [(C_K/S_K)(chi_ref) min(chi_ref,chi_c)].

In the low-curvature limit, after the ray exits the local region this falls
approximately as `chi_c/chi`.

The code carries `z_c=0.03` only as an illustrative coherence proxy.  It is
not fitted or frozen.  The intended replacement is an actual foreground
source map whose weights come from matter/structure plus BHSM seam charges.

## Calibration / identifiability boundary

The low-z SN amplitude cannot calibrate an independent global n=2 amplitude
and a local-topographic amplitude simultaneously.

For this exploratory LOS test, the carried-forward low-z SN amplitude
normalizes the local component only.  The global n=2 amplitude is left
theory-owned/open until the BHSM action or an independent preregistered
observable fixes the global/local split.

This prevents the earlier mistake of forcing the entire observed SN dipole
onto one homogeneous global mode.

## Falsification target

With an actual foreground catalog, the model predicts that residuals should
correlate with shared foreground topography, and that sightlines with little
shared foreground structure should decorrelate with path length.  Failure of
that correlation after conventional peculiar-velocity, density and lensing
effects are controlled would count against the added BHSM topographic layer.
