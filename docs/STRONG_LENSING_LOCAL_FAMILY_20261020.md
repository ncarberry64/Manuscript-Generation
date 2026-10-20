# Constructed local reference family and exact physical action projection

This calculation constructs a conditional, deterministic envelope reference
family and carries its complete gravitational encounter in the companion
finite-body calculation.  The family has a derived zero additional phase
response.  Its gravitational metric and optical encounter are nonzero.  It is
an equivalence control, rather than a selection of the state of an observed
galaxy.

## Recovered action and scope

The original physical Path-B sector is

\[
 S_{\rm color,\eta}=S_{\rm YM}[A]
 -\int d\mu_g\,w(\sigma)F(X)+S_{\rm unit},\qquad
 F(X)=\frac{\kappa_1 X}{2}+\frac{X^4}{8},\quad
 X=g^{\mu\nu}\langle D_\mu\eta,D_\nu\eta\rangle .
\]

Here \(\eta\) is the physical four-dimensional section of
\(P_{SU(3)}\times_{SU(3)}S^6\), not a second independently varied
seven-dimensional ultraviolet texture.  The ambient target is
\((s,z)\in\mathbb R\oplus\mathbb C^3\), with
\(s^2+2z^\dagger z=1\) and target metric
\(ds^2+2\operatorname{Re}(dz^\dagger dz)\).  \(A\) is the physical
\(SU(3)\) connection; \(w>0\) is the inherited scalar weight.
The constants and length conversion of this sector are conditional inputs,
not parameters obtained by fitting a deflection.

The local gravity calculation explicitly uses the retained four-dimensional
Einstein term, \(C_\partial=c^4/(16\pi G)\), with the ordinary matter
stress and the Dirichlet boundary term.  This is a conditional effective
gravity problem.  The original source separately has five-dimensional cap
equations and matcher reactions in its intrinsic metric equation.  A cap
solution and a derivation of its reduction to this local Einstein problem
are **not** supplied here.  The additional eta forcing of any inherited
scalar equation is nevertheless exactly zero, \(w'(\sigma)F(0)=0\).
This establishes no eta-dependent scalar or cap source.  The independent
local sigma and cap equations lie outside this effective truncation; a
constant inherited sigma is not thereby asserted to solve those equations.

The eight original source files are pinned in
`artifacts/strong_lensing_local_family_source_receipt_20261020.json`.
They were recovered read-only from BHSM commit
`524ed90689bd5923c249bba2e699abf627e703cd`; their repository-normalized blobs
all agree with that commit.  Exact local-byte SHA256 hashes are also retained.
The original definitions are independently accessible in:

- [physical action and current, v14.31](https://github.com/ncarberry64/Berger-Hopf-Standard-Model/blob/524ed90689bd5923c249bba2e699abf627e703cd/src/bhsm/interface/completion/path_b_foundational_action_v14_31.py),
- [physical S6 target, v14.32](https://github.com/ncarberry64/Berger-Hopf-Standard-Model/blob/524ed90689bd5923c249bba2e699abf627e703cd/src/bhsm/interface/completion/path_b_physical_topology_v14_32.py),
- [mixed inverse-metric/eta variation, v14.39](https://github.com/ncarberry64/Berger-Hopf-Standard-Model/blob/524ed90689bd5923c249bba2e699abf627e703cd/src/bhsm/interface/completion/static_eta_metric_spin4_source_v14_39.py),
- [phase self variation, v14.36](https://github.com/ncarberry64/Berger-Hopf-Standard-Model/blob/524ed90689bd5923c249bba2e699abf627e703cd/src/bhsm/interface/completion/degree_one_phase_hessian_v14_36.py),
- [Lorentzian kinetic condition, v14.87](https://github.com/ncarberry64/Berger-Hopf-Standard-Model/blob/524ed90689bd5923c249bba2e699abf627e703cd/src/bhsm/interface/completion/eta_relative_periodic_legendre_current_gate_v14_87.py).

The master-action `terms.py` and `variations.py`, pinned in the same receipt,
fix the intrinsic Einstein coefficient, oriented GHY term, and the distinction
between this conditional four-dimensional problem and cap matching.

## Actual state, equations, and zero-flux derivation

Take \(e_1=(1,0,0)\), and the meridian family

\[
 \eta_q=(\cos q,\,\sin q\,e_1/\sqrt2),\qquad
 A_\mu=0,\qquad q\text{ spatially and temporally constant}.
\]

For example \(q=\pi/4\) specifies a nonzero triplet component.  The
conditional parameter \(q\) is a field value, not an incident-light
amplitude.  Its unit constraint is exact, and its target tangent
\(\partial_q\eta_q\) has norm one.  The incoming optical pulse remains
an independent probe of the spacetime metric.

The constancy is also a solution of an explicit radial boundary problem.
Write the exact static spherical metric as
\(ds^2=-N(r)^2(cdt)^2+H(r)^2dr^2+r^2d\Omega^2\), and initially allow a
meridian angle \(\theta(r)\).  The original eta equation restricted to
this geodesic is

\[
 \frac{d}{dr}\left[
   \frac{N r^2}{H}\,w(\sigma)
   (\kappa_1+X^3)\theta'\right]=0,
 \qquad X=\frac{\theta'^2}{H^2}.
\]

Regularity at the center, finite positive \(N,H,w\), and bounded
\(\theta'\) make the constant radial flux zero.  On the retained
\(\kappa_1>0,X\ge0\) branch the bracket therefore implies
\(\theta'=0\).  Prescribe the same constant \(q\) at the outer
boundary or at infinity; the radial meridian solution is \(\theta=q\).
This argument is restricted to this source-free meridian family.  It does
not assert that every allowed nonlinear or sourced eta configuration is
constant.

Evaluation on this solution gives, pointwise on **any** matched metric,

\[
 D_\mu\eta_q=0,\quad X=F(0)=0,\quad
 A_{\mu\nu}=\langle D_\mu\eta,D_\nu\eta\rangle=0,\quad
 j_\mu=\langle D_\mu\eta,T\eta\rangle=0.
\]

The original tangent eta Euler equation, its unit multiplier, eta stress,
and momentum density all vanish.  The original color equation is

\[
 g_3^{-2}(D_\nu F^{\nu\mu})_a
 =J^\mu_{\eta,a}+J^\mu_{{\rm other},a},\qquad
 J^\mu_{\eta,a}=w(\kappa_1+X^3)K_{aI}(\eta)D^\mu\eta^I.
\]

The declared neutral perfect-fluid body has no color current, so both
sides are zero for \(A=0\).  The scalar conormal
\(w(\kappa_1+X^3)n^\mu D_\mu\eta\) is zero on both sides of every
ordinary stellar matching surface.  Thus no eta source or surface layer is
hidden at that surface.  Ordinary gravity remains sourced by the prescribed
body; it is not canceled by the phase null.

The companion body calculation uses a finite Tolman-VII density profile
with \(M=10^{11}M_\odot\), areal radius \(R=3\,\mathrm{kpc}\),
and pressure/lapse obtained from the original TOV equations.  Its pressure
vanishes at the surface.  The interior metric matches exterior Schwarzschild
with continuous induced metric and extrinsic curvature.  Hence the ordinary
surface has no Israel layer.  This is the same body in the conventional and
BHSM calculations, and its metric is counted once.

## Exact common-annulus metric and phase projection

The primary local action projection lies in the vacuum exterior, where there
is no unspecified matter second variation.  The matched body fixes the
boundary metric.  Use isotropic radius \(r\), mass length \(m=GM/c^2\),
and physical log-potential coordinates

\[
 ds^2=-e^{2A}(cdt)^2+e^{-2P}d\mathbf x^2,\qquad
 A_0=\ln\frac{1-m/(2r)}{1+m/(2r)},\qquad
 P_0=-2\ln(1+m/(2r)).
\]

The conventional weak-field potentials equal \(A,P\) at first order.
These coordinates retain their exact metric definition in this projection.
They are not a relabeling of the cosmological ADM lapse and shift.

The domain is the same for both physical metric directions and the target
phase:

\[
 r_0=\frac{R-m+\sqrt{R(R-2m)}}2,\quad
 r_0<r<L=2R,\quad \ell=L-r_0,
 \qquad
 f(r)=\frac{\ell}{\sqrt{2\pi}\,r}
       \sin\frac{\pi(r-r_0)}\ell .
\]

Both endpoint traces vanish, and
\(\int f^2d^3x/\ell^3=1\),
\(I_{\rm flat}=\int|\nabla f|^2d^3x=\pi^2\ell\).
The outward annular normals point inward at \(r_0\) and outward at
\(L\).  The zero eta conormal makes both fluxes zero.  These annular
boundaries are fixed computational variation boundaries, not physical new
interfaces.

Use \(\delta A=x_Af\), \(\delta P=x_Pf\), and the canonical
target phase \(\delta\theta=p_\theta f\).  Equivalently its inverse
metric basis has
\((e_A)^{00}=2e^{-2A_0}f\),
\((e_P)^{ij}=2e^{2P_0}f\delta^{ij}\), with other entries zero.
All three amplitudes are dimensionless.  \(T\eta_q=\partial_q\eta_q\)
has unit norm.  The phase tangent is not the optical pulse, a radiation-fluid
variable, or a reference-state covariance.

In the **static energy** convention \(E=-L_{\rm static}\), the exact
Einstein-Hilbert plus Dirichlet boundary reduction is

\[
 E_g=K\int d^3x\,e^{A-P}
       [2\nabla A\cdot\nabla P-|\nabla P|^2],
 \qquad K=\frac{c^4}{8\pi G}.
\]

This follows directly from \(N\sqrt h\,{}^{(3)}R\) for a static
conformally flat spatial metric, integrating its Laplacian term with fixed
boundary values.  Reversing to the Lorentzian static action reverses the
entire static Hessian.  The physical Lorentzian kinetic gate is evaluated
separately below, rather than inferred from the signature of this static
constraint block.

For \(p=A_0'\), \(q_r=P_0'\), \(g=f'\),
\(B_0=q_r^2-2pq_r\), and \(D=A_0-P_0\), differentiating this actual
functional twice gives

\[
\begin{aligned}
 a&=K\int d^3x\,e^D[-f^2B_0+4fq_rg],\\
 b&=K\int d^3x\,e^D[f^2B_0-4fq_rg+2fpg+2g^2],\\
 c&=K\int d^3x\,e^D[-f^2B_0+4fq_rg-4fpg-2g^2].
\end{aligned}
\]

The independently evaluated original vacuum Euler equations are

\[
 B_0+2[P_0''+2P_0'/r+(A_0'-P_0')P_0']=0,
\]
\[
 B_0+2[P_0''-A_0''+2(P_0'-A_0')/r
               +(A_0'-P_0')(P_0'-A_0')]=0.
\]

The first constraint and the homogeneous endpoint traces imply \(a=0\)
as an integrated identity.  The code evaluates the original integral rather
than inserting zero.  In the flat limit,
\(a=0,b=2KI_{\rm flat},c=-2KI_{\rm flat}\), so
\(\Delta_g=-4K^2I_{\rm flat}^2\ne0\).

There is also an exact curved simplification.  Schwarzschild obeys
\(r^2e^{A_0-P_0}A_0'=m\), so integration by parts gives
\(\int d^3x\,e^D fA_0'f'=0\).  Combining this with the first
constraint yields
\[
 A_g^{\rm loc}=2K I_{\rm curv}
       \begin{pmatrix}0&1\\1&-1\end{pmatrix},\qquad
 I_{\rm curv}=\int d^3x\,e^D|\nabla f|^2>0.
\]
The executable independently integrates the unsimplified original second
variations and checks this identity.  Thus its nonzero determinant and finite
condition number are not numerical accidents of a weak-field example.

The retained mixed-source tensor defines \(F,X,j,A_{\mu\nu},w\)
explicitly as above and is evaluated on this same annulus:

\[
 \mathcal S_{\mu\nu}
 =2F''(X)(j\cdot d\varphi)A_{\mu\nu}
 +2F'(X)j_{(\mu}\partial_{\nu)}\varphi
 -F'(X)g_{\mu\nu}(j\cdot d\varphi).
\]

Its every term is zero on the constructed eta state, so \(u=v=0\)
exactly for both physical metric directions.  The actual phase energy
self-entry is nevertheless nonzero:

\[
 c_\varphi=k_\eta\int d^3x\,e^{A_0-P_0}|\nabla f|^2,
 \qquad k_\eta=w\kappa_1>0,
 \quad e^{A_0-P_0}=1-\frac{m^2}{4r^2}.
\]

Here \(k_\eta\) denotes the dimensionful two-derivative stiffness after
the sector's chosen length/energy conversion; it has units energy/length.
\(K\) has the same units.  Thus \(a,b,c\) have units energy and
\(c_\varphi\) has units energy for dimensionless amplitudes.  If the
original density is expressed in dimensionless reference units
\(E_\eta,\ell_\eta\), then
\(k_\eta=E_\eta w\kappa_1/\ell_\eta\); this conversion remains a
conditional physical normalization.  It does not affect the derived zero
metric response.  No physical length is extracted by equating the eta
normalization to the galaxy radius.

The Lorentzian tangent kinetic matrix is \(k_\eta\) times the identity
in a local orthonormal frame.  For this coordinate-time phase mode its
positive integrated coefficient is
\[
 K_\varphi^{(t)}=\frac{k_\eta}{c^2}
 \int d^3x\,e^{-A_0-3P_0}f^2>0.
\]

The evaluated physical inputs give
\(m/R=1.595138613862614\times10^{-6}\).  In units \(KR\) for the
metric block and \(k_\eta R\) for its phase self-entry, the application
returns

| Quantity | Evaluated value |
| --- | ---: |
| \(a/(KR)\) | \(-1.27\times10^{-21}\) (analytic zero) |
| \(b/(KR)\) | \(19.73924028895931\) |
| \(c/(KR)\) | \(-19.73924028895931\) |
| \(u,v\) | \(0,0\), exactly |
| \(c_\varphi/(k_\eta R)\) | \(9.869620144479654\) |
| Maximum original vacuum Euler residual, in \(R\) units | \(3.27\times10^{-21}\) |

The metric condition number is approximately \(2.618\), the determinant
is negative and nonzero, and

\[
 \mathcal R_{g\varphi}^{\rm loc}
 =\frac1{ac-b^2}(bv-cu,\,bu-av)^T=(0,0)^T,
 \qquad k_{\rm loc}=c_\varphi>0.
\]

This inverse is justified for the retained **two-dimensional Dirichlet
projection**.  It is not a proof of invertibility or stability of the entire
metric functional.  Its indefinite static gravitational signature is not
identified as a propagating kinetic signature.  The full curved projection,
normalization, quadrature errors, kinetic norm, and original-equation residuals
are recorded in `artifacts/strong_lensing_local_hessian_20261020.json`.

## Normalized phase source and the precise constitutive boundary

The canonical target phase has a directly derived, dimensionful continuum
operator on this background:

\[
 \mathcal L_\theta\varphi
 =-\nabla\cdot[k_\eta e^{A_0-P_0}\nabla\varphi].
\]

With regularity across the body, homogeneous perturbation traces, and zero
conormal source, its solution is \(\varphi=0\).  Its projection onto
\(f\) is \(c_\varphi p_\theta=0\).  A gravitational source in this
phase equation vanishes because its mixed action block is exactly zero; this
is derived, rather than filled in from a missing array.  The ordinary stellar
matching has \(\Delta K_{uu}=0\), so the specified seam charge
\((c^2/8\pi)\int H_\Sigma\Delta K_{uu}dA\) is zero for that surface.
The specified static material-surface seam contribution is therefore zero,
independently of its otherwise unresolved orientation. This does not set the
separate exported scalar's gravitational source to zero. The declared encounter is deterministic and has no
stochastic source sector; no covariance assumption is needed to solve it.

This solves the target-phase source restriction and its zero metric response.
It does not solve the separately exported scalar in the inherited chain
\(\mathcal L_\psi\psi=J_{\rm grav}+\mathcal E_\Sigma[Q_\Sigma]+\xi\),
\(h_{\rm topo}=\mathcal R_{g\psi}\psi\), or determine that scalar's
gravitational source \(J_{\rm grav}\).  The smooth material surface has
zero geometrically defined seam charge on that surface only.  The full
metric source \(J_g\) from the finite body still produces its ordinary
gravitational metric, equivalently \(A_gh+B\varphi+J_g=0\); it is not
omitted.

For a nonzero exported topographic field, the constitutive identification
between that field and the canonical target angle is **still unresolved**.
In particular, \(\mathcal L_\theta\) is not renamed as the earlier
\(-\nabla^2+B\nabla^4\) topographic operator.  A nonzero map would
require its action coupling, normalization, and common boundary data.  In the
metric convention used here, a spatial conformal metric amplitude is
\(-\delta P\), so the evaluated canonical phase channel has zero such
amplitude.  A zero response does not determine the
off-control constitutive map or a universal absence of interface response.

Thus the construction specifies the constant eta boundary value, finite
ordinary body, and effective Einstein gravity normalization; derives the eta
equations, zero currents/flux/stress, exact curved local Hessian, and phase
null; and leaves the selection of \(q,k_\eta\), cap reduction, and a
nonzero export-to-phase coupling conditional or unresolved.  The companion
optical calculation uses the nonzero matched metric in both comparisons,
so its exact equivalence in this family is a substantive optical application.

## Focused reproduction

```powershell
python -m pytest tests/test_strong_lensing_local_hessian.py -q
python code/strong_lensing_local_hessian.py --m-over-R 1.595138613862614e-6 --output artifacts/strong_lensing_local_hessian_20261020.json
```

The six focused checks include the analytic flat limit, original curved
vacuum constraints, finite-difference second variation of the original exact
energy at appreciable compactness, endpoint-domain validity, and the
positive phase Schur entry.  The companion body runner records the same
physical constants and assembles this result with its TOV and optical data.
