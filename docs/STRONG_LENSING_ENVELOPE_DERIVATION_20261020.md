# Envelope action projection and the available static-channel result

The original ADM/control derivation below is retained. An evaluated
physical metric/phase Hessian, target source problem and finite encounter
now follow in [the local-family derivation](STRONG_LENSING_LOCAL_FAMILY_20261020.md).
Its restricted effective sector and remaining exported-scalar dependency
are explicit there.

This calculation extends the current action and optical construction. It does
not infer an unspecified state from a missing numerical array. The executed
projection is the lapse/shift block of the current R1 action, followed by its
exact constraint elimination. The local target-phase optical block is a
different projection and is kept distinct below.

## Source trail and the domain used

The original polynomial is `unreduced_lagrangian` in
`code/action_native_matter_n2.py`; `action_blocks` is its Hessian assembly.
`docs/R1_ACTION_NATIVE_MATTER_BACKREACTION_V1_DERIVATION.md` records its ADM
and Schutz--Sorkin derivation. The newer gradient/dispersion section preserves
the already audited effective R1 domain. This checkpoint does not repeat that
broad stability audit. `frozen_background` reads the existing manifest's
shooting value and performs no new shooting or observational calibration.

Use cosmic-time dots, H0=Mpl=1, N=ln(a), and one real n=2 harmonic satisfying
integral sqrt(gamma) Q^2=1 on closed S3. Define kappa=K/a^2, k=8 kappa,
s=5 kappa, h_I=rho_I+p_I, h=h_m+h_r, w_m=0, w_r=1/3.
The quadratic action density is divided by a^3. Temporal variations vanish at
the variation endpoints; S3 has no spatial boundary. The pointwise projection
and Schur calculation use the same N, harmonic and background throughout.
The native metric directions are unit harmonic amplitudes of the lapse alpha
and shift potential chi (N_i=D_i chi); chi is measured in H0^-1. This specifies
their coordinate normalization without declaring the ADM auxiliary quadratic
form to be a positive metric or a physical Newtonian-potential norm.

The retained original Lagrangian is

```text
L/a^3 = -3 zdot^2+s z^2+Sigma alpha^2+2k Theta alpha chi
        -2k zdot chi+6Theta alpha zdot+2s alpha z-k kappa chi^2
        +sum_I h_I [v_I Ddot_I-k v_I^2/2+k chi v_I
                    -w_I(D_I-3z)^2/2-alpha(D_I-3z)].
```

Here z denotes zeta, not redshift. The exact auxiliary equations are

```text
0 = Sigma alpha+k Theta chi+3Theta zetadot+s zeta
    -(1/2)sum_I h_I(D_I-3zeta),
0 = Theta alpha-zetadot-kappa chi+(1/2)sum_I h_I v_I,
0 = Ddot_I-k(v_I-chi).
```

The implementation evaluates these original equations independently of the
matrix projection. Its unit probes concern instantaneous source coordinates
(zetadot,Dmdot,Drdot,zeta,Dm,Dr); they are derivative columns, not a claim of six
independently prescribed physical solutions. The full canonical evolution
and physical six-coordinate anchor select actual combinations elsewhere.

## Extracted native metric and source blocks

First solve v_I=chi+Ddot_I/k. This gives

```text
m=(alpha,chi), q=(zetadot,Dmdot,Drdot,zeta,Dm,Dr),
L/a^3 = (1/2)m^T A_ADM m + m^T J_ADM q + (1/2)q^T C_ADM q,

A_ADM = [[2Sigma,2kTheta],[2kTheta,k(h-2kappa)]],
J_ADM = [[6Theta,0,0,2s+3h,-h_m,-h_r],
         [-2k,h_m,h_r,0,0,0]],
C_ADM = diag_block(diag(-6,h_m/k,h_r/k), E),
E = diag(2s,0,0)-(h_r/3)(-3,0,1)^T(-3,0,1).
```

These are actual action projections after eliminating fluid velocities.
For any symmetric 2x2 block A=[[a,b],[b,c]], Delta=ac-b^2, the explicit
inverse application is

```text
-A^-1 J = [(b J_row2-c J_row1)/Delta,
           (b J_row1-a J_row2)/Delta].
```

For a single local mixed column J=(u,v), this is precisely the manuscript's
(bv-cu,bu-av)/Delta formula. The computation checks the explicit expression
against a numerical solve and against simultaneous elimination of all four
original auxiliaries. The full same-domain Schur block is
C_ADM-J_ADM^T A_ADM^-1 J_ADM. It agrees with the original reduced block
[[M,L],[L^T,V]].

At lens redshift 0.3 the executed snapshot gives

| Native quantity | Value |
| --- | ---: |
| a_ADM | -7.96725027661529 |
| b_ADM | 0.5695448002713244 |
| c_ADM | 0.48267978497855984 |
| determinant | -4.1700119299031435 |
| condition number | 15.368657973059397 |
| norm(A_ADM R+J_ADM) | 1.41621372581761e-15 |
| explicit inverse versus solve norm | 4.849413209257403e-16 |
| maximum original constraint residual | 5.273559366969494e-16 |
| maximum simultaneous auxiliary reconstruction error | 1.4657992515888388e-15 |
| full Schur versus native reduction norm | 4.154655943942955e-16 |
| c_zeta before metric elimination | 0.301115412 |
| zeta Schur curvature | 5.0970947797763175 |

The corresponding source columns are

```text
J_ADM = [[ 7.021015782437431, 0, 0, 6.436914588,
          -2.04321, -0.001028196],
         [-0.48672, 2.04321, 0.001028196, 0, 0, 0]].
```

The native kinetic eigenvalues are
(0.004182113778467213,0.17251064482078457,0.4338635425190717),
so the Legendre gate for this numerical application is positive. A negative
auxiliary determinant is not itself a ghost criterion: these are constrained
lapse/shift variables, not the propagating kinetic form.

The machine-readable snapshot retains every response column, original
residual, self block and Schur entry. No entry above is a fit parameter.

## Why these entries are not the local phase coefficients

The current provider's conventional Newtonian-potential dictionary is

```text
Phi=alpha+chidot, Psi=-zeta-H chi.
```

It depends on the time derivative and the evolving state. It is not an
invertible instantaneous 2x2 change from (alpha,chi) to (Phi,Psi). Consequently
the extracted native a_ADM,b_ADM,c_ADM and zeta source columns cannot be
renamed a,b,c,u,v of a local target phase. The complete canonical-to-potential
rows supply the cosmological optical response on their declared state; they
do not by themselves supply the local eta/phase response.

For the latter, the retained manuscript equations explicitly prescribe

```text
a=<e_A,H_gg e_A>, b=<e_A,H_gg e_psi>, c=<e_psi,H_gg e_psi>,
u=B_gphi[e_A,phi_loc], v=B_gphi[e_psi,phi_loc],
c_phi=<phi_loc,H_phiphi phi_loc>,
B_gphi[gamma,phi]=integral dmu*w*gamma^{mu nu} [
  2F''(X)(j.dphi) A_mu_nu +2F'(X)j_(mu partial_nu)phi
  -F'(X)g_mu_nu(j.dphi)].
```

The required state contractions are j_mu=<D_mu eta,T eta> and the tensor
A_mu_nu in that same action, with F'(X),F''(X), the integration weight, and
the constrained metric embeddings e_A,e_psi. The equation requiring them is
the displayed mixed action variation, not an inclusive photon or muon
correlation function. To determine them for a local envelope one must supply
the eta/common-action background and its boundary/state prescription and
identify the projection of the exported field psi onto the normalized local
phase profile. The source equation

```text
L_psi psi=J_grav+E_Sigma[Q_Sigma]+xi
```

also needs its common-domain operator and boundary data, the normalization
and orientation of E_Sigma, the state/source contractions defining J_grav,
Q_Sigma,xi, and the matching between this psi and phi_loc. Specifying a source
amplitude alone does not fix a normalized profile or this matching. Neither
the R1 background's isotropy nor its radiation-fluid delta_r,v_r chooses the
local eta state. The incident beam and those fluid variables are distinct.

This is the precise operator/state dependency. The available cosmological
Hessian has been applied as far as its declared variables allow; identifying
zeta with the phase phi would be an extra, unestablished constitutive relation.

## Actual static spatial phase null application

The source trail is manuscript sections 04aq and 04ar, implemented in
`bhsm_optical_current_gate.spatial_mixed_source_tensor`. On the identified
static block-diagonal phase channel its original source is

```text
S_ij=2F''(j.grad phi) A_ij+2F' j_(i grad_j)phi
     -F' g_ij(j.grad phi).
```

If j_i=0, each of its three terms vanishes pointwise, for every normalized
admissible phi_loc. If grad_i phi=0 the same is true for every finite j_i.
The new executable applies this original tensor and both metric projections;
it returns S_ij=0 and **u=v=0 exactly** for each specified channel. Numerical
representative tensors in that identity check are not selected physical
profiles and set no amplitude. Since this is a pointwise identity it is
independent of quadrature and the still-undetermined metric embeddings.

For an invertible same-domain metric block, R_gphi=0 and the Schur curvature
is exactly c_phi. Neither j_i=0 alone nor the static condition j_0=0 fixes
c_phi. Also D_0 eta=0 and j_0=0 do not imply j_i=0. This checkpoint therefore
does not choose the null branch as the universal physical envelope state.
This channel-specific zero does not cancel conventional gravity, induced
fluid response, other common-action sectors or the cosmological optical maps.

The phase-invariant F(X) sector supplies a further symbolic self-projection.
If X=g^{mu nu}<D_mu eta,D_nu eta>, and the target symmetry used in the phase
specialization acts as eta_epsilon=exp(epsilon phi T)eta, write
q_T=|T eta|^2. Target invariance gives

```text
X_epsilon=X+2epsilon(j.dphi)+epsilon^2 q_T(dphi)^2,
c_phi[F(X)] = integral dmu*w*[
                 4F''(X)(j.dphi)^2+2F'(X)q_T(dphi)^2].
```

This follows by differentiating F(X_epsilon) twice and is the sector whose
metric derivative yields the retained mixed tensor. At j_i=0 a nonconstant
static texture still has self curvature 2 integral dmu*w*F' q_T|grad phi|^2.
It is not determined without the same target state and profile. For a
constant phase this sector contributes zero, provided that constant is in
the physical domain. Boundary terms or other phase-dependent action terms
must be added before reporting the complete c_phi. No new vacuum, Gaussian
factorization, muon covariance or photon interaction is chosen.

## Static seam versus moving continuation

The inherited static export is
epsilon_Sigma(Q_Sigma/Q)Theta(r-r_star)/r^2. Its normalization, common-action
orientation and matching onto the phase/source operator remain explicit.
Its time derivative contains the bulk Qdot term and the moving-seam delta
term. That derivative is a proposed kinematic continuation; this checkpoint
does not use it as an action-selected de-encapsulation law or insert an
unfixed seam amplitude. Any benchmark using smooth cosmological metric maps
counts that metric contribution once.

## Reproduction and focused checks

```powershell
python code/strong_lensing_envelope.py --z 0.3 --output artifacts/STRONG_LENSING_ENVELOPE_PROJECTION_20261020.json
python -m pytest tests/test_strong_lensing_envelope.py -q
```

Seven tests pass. They check the original polynomial action and constraint
equations, sequential versus simultaneous Schur elimination at z=0,0.3,1,2.1,
the explicit inverse, the native/local domain distinction and the restricted
static phase null. These are algebraic/constraint applications, not physical
admission of a finite envelope state. The next local computation is to supply
the same-domain eta/profile/source matching and evaluate its displayed
mixed and self projections; the conventional benchmark and conditional
cosmological transport proceed independently of that normalization.
