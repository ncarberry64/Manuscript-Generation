# Conditional axial optical transport for the first strong-lensing comparison

This checkpoint evaluates current action-native metric responses as physical,
measured-redshift angle-to-separation maps. It supplies the cosmological
propagation for the declared thin galaxy encounter. It is conditional on a
chosen normalized spatial profile and sky axis; it is neither an arbitrary-sky
transport calculation nor a determination of the envelope's physical amplitude.
The current action-native provider and its inherited closure qualifications are
retained. This calculation does not certify a new action or repeat its broad
stability audits.

## State, normalization and the equations actually applied

The anchor is at background redshift 2.1. For the chosen profile it has six
independent physical coordinates

\[
Y=(\zeta,\Pi_2,\delta_m^{\rm ADM},\delta_r^{\rm ADM},
v_m^{\rm ADM},v_r^{\rm ADM}),\qquad
\Pi_2=2a_i^3\mathcal G_{S,2}\dot\zeta_i.
\]

The existing action has configuration variables
\(x=(\zeta,D_m,D_r)\), where
\(D_I=\delta_I^{\rm ADM}/(1+w_I)+3\zeta\), and canonical momenta
\(p=a^3(M\dot x+Lx)\). The new provider inverts six current action rows:
\(\zeta\), \(2a_i^3\mathcal G_{S,2}\dot\zeta\),
\(D_m-3\zeta\), \(4(D_r-3\zeta)/3\), and the two reconstructed ADM
velocities. This constructs \((x,p)_i=C_iY\). The applied evolution is

\[
\dot x=M^{-1}(p/a^3-Lx),\qquad
\dot p=a^3(L^T\dot x+Vx),\qquad
u=-A^{-1}(B\dot x+Dx).
\]

These are calls to `canonical_maps`, `hamiltonian_matrix`,
`integrate_propagator` and `newtonian_map` in
`code/action_native_matter_n2.py`; no legacy growth/luminosity solution is
substituted. In this background the resulting canonical momentum is

\[
p_\zeta=0.0184841134\,\zeta+\Pi_2
-0.296053883\,\delta_m^{\rm ADM}
-0.000266448495\,\delta_r^{\rm ADM}
-0.00218954758\,v_m^{\rm ADM}
-0.00000262745709\,v_r^{\rm ADM}.
\]

Thus the coefficient of \(\Pi_2\) happens to be one here, but \(p_\zeta\)
is not generally \(\Pi_2\). The historical two columns are recovered by
setting the other four physical anchor amplitudes to zero; they are not the
full state space.

The profile is

\[
Q(\chi,\vartheta)=\mathcal N\sin(2\sqrt K\chi)\cos\vartheta,
\quad \mathcal N=\sqrt{\frac{3}{\pi^2R^3}},\quad R=K^{-1/2}.
\]

For the dimensional spatial volume
\(dV_\gamma=R^3\sin^2(\chi/R)d(\chi/R)d\Omega\),
\(\int Q^2dV_\gamma=1\), and \(D^2Q=-8KQ\). At frozen
\(K=0.018\), \(\mathcal N=0.027093501352896836\) in \(H_0=1\) units.
This is one profile in the nine-dimensional spatial scalar eigenspace. Its
positive axis is an explicit conditional input, not an axis selected by an
isotropic background or by earlier supernova fits.

Write the current conventional-sign metric potentials as
\(\Phi=\phi(N)Q\), \(\Psi=\psi(N)Q\). Both rows are used, without forcing
equality or adopting the legacy positive spatial-potential sign. The metric
used for finite checks is the explicitly declared completion

\[
ds^2=-A^2d\eta^2+B^2[d\chi^2+S_K^2(\chi)d\Omega^2],\quad
A=a\sqrt{1+2\Phi},\quad B=a\sqrt{1-2\Psi}.
\]

Only its derivatives at \(Y=0\) are interpreted as the linear action-native
response. The action supplies no second-order metric, so finite trial rays are
not nonlinear physical cosmologies.

## Radial geodesic, observed redshift and matter tetrads

The observer is at \(\chi=0,N=0\); the central ray follows the positive axis
\(\vartheta=0\). Let \(q(\chi)=Q(\chi,0)\) and let \(\lambda\) increase
into the past. Set

\[
k^\eta=-E/A,\qquad k^\chi=E/B,
\]

where \(E\) is the static-tetrad frequency, not the matter-frame frequency.
Applying the original radial null geodesic equation gives

\[
\frac{dN}{d\chi}=-aH\sqrt{\frac{1-2\psi q}{1+2\phi q}},
\]
\[
\frac{d\ln E}{d\chi}=
aH\sqrt{\frac{1-2\psi q}{1+2\phi q}}
-\frac{aH\psi_Nq}{\sqrt{(1-2\psi q)(1+2\phi q)}}
-\frac{\phi q_{,\chi}}{1+2\phi q},\qquad
\frac{d\lambda}{d\chi}=B/E.
\]

The second equation also follows from
\(d(AE)/d\lambda=\tfrac12\partial_\eta g_{ab}k^ak^b\); it includes time
variation of the spatial potential. No fixed-affine Green-function kernel is
treated as a complete measured-redshift observable.

The velocity sign and scale follow directly from the action's
Schutz--Sorkin convention:
\(\partial_i\ell=\mu u_i\), \(\delta\ell=\mu v\), and current Newtonian
\(v_N=v_{\rm ADM}-\chi_{\rm ADM}\). Consequently
\(u_\chi=v_Nq_{,\chi}\). Defining

\[
s=\frac{v_Nq_{,\chi}}B,\qquad
\gamma=\sqrt{1+s^2},\qquad
v_{\rm physical}=s/\gamma,
\]

provides a normalized matter tetrad in the finite metric completion. In
particular, the first-order physical velocity is
\(v_Nq_{,\chi}/a\). The legacy expression without this factor \(1/a\)
cannot be reused as a physical velocity here. Along this past-directed ray,
the positive matter-frame frequency is

\[
\omega=u\cdot k=E(\gamma+s).
\]

The observer normalization is \(\omega_o=1\). Although \(Q_o=0\),
\(q_{,\chi}(0)=2\sqrt K\mathcal N\ne0\), so the observer's matter velocity
and aberration generally remain nonzero. Setting the observer velocity to zero
on the basis of \(Q_o=0\) would change the experiment.

Lens and source events are roots \(\omega_l=1.3\), \(\omega_s=2\). Their
\(N,\chi,\lambda\) coordinates move with the state. The code does not hold
background emission times or affine coordinates fixed while claiming fixed
measured redshifts.

### Fixed-affine frequency, fixed-radial frequency and the moving event

The original numerical event-root application above is retained. Its
coordinate endpoints require two distinct first variations:

\[
Z_\lambda=\delta\ln\omega|_\lambda,\qquad
Z_\chi=\delta\ln\omega|_\chi,\qquad
\delta\chi_{\rm ray}|_\lambda=\delta\chi|_\lambda.
\]

Since \(d\ln\bar\omega/d\chi=\mathcal H\), the chain rule is

\[
Z_\chi=Z_\lambda-\mathcal H\delta\chi_{\rm ray}|_\lambda.
\]

At a fixed measured-redshift event, observer normalization fixed,

\[
\delta\lambda_{\rm event}
=-\frac{Z_\lambda}{d\ln\bar\omega/d\lambda},
\]
\[
\delta\chi_{\rm event}
=\delta\chi_{\rm ray}|_\lambda
 +\bar\chi'_\lambda\delta\lambda_{\rm event}
=\delta\chi_{\rm ray}|_\lambda-\frac{Z_\lambda}{\mathcal H}
=-\frac{Z_\chi}{\mathcal H}.
\]

The converted background-affine shift
\(\delta\chi_{\rm affine}\equiv\bar\chi'_\lambda\delta\lambda_{\rm event}\)
has its own symbol and is not the actual radial-coordinate event shift.
For the saved unit-\(\zeta_i\) source application,
\(\partial_{\zeta_i}\chi_s=-0.005611965495588844\) whereas
\(\partial_{\zeta_i}\lambda_s=0.00002815169979075716\).
At \(z_s=1\), \(\bar\chi'_\lambda=a_s^{-2}=4\), giving
\(\partial_{\zeta_i}\chi_{{\rm affine},s}=+0.00011260679916302863\).
The chain rule consequently infers the fixed-affine ray displacement
\(-0.005724572294751873\) per unit \(\zeta_i\). The last number is
inferred from the saved moving-event responses, not a newly integrated
fixed-affine sample.

For Jacobi initialization at a moving vertex, the original conditions are
\(J_a(\lambda_a)=0\), \(J'_a(\lambda_a)=\omega_a I\). Their variation is

\[
\delta J_a|_a=-\omega_a I\delta\lambda_a,\qquad
(\delta J_a)'|_a=\Delta\omega_a I,
\quad
\Delta\omega_a=\delta\omega_a|_\lambda
+\bar\omega'_a\delta\lambda_a.
\]

The frequency on this right-hand side is the **total** variation at the
moving vertex, so \(\Delta\omega_a=0\) at fixed measured redshift.
The term \(\bar J_a''\delta\lambda_a\) vanishes because the initial
Jacobi separation vanishes. The propagation variation remains
\((\delta J_a)''=\mathcal T_0\delta J_a+\delta\mathcal T J_{a0}\),
and the output includes \(J'_{a0}(\lambda_b)\delta\lambda_b\) and
screen/frame conversion. The event-root implementation already includes
these moving redshift events and is unchanged by this notation correction.

## Angular Jacobi equation and compatible three-segment initialization

For a small screen deviation from the axial ray, let \(\vartheta\) denote
the polar angular perturbation and
\(L=B^2S_K^2\,d\vartheta/d\lambda\). The angular Euler--Lagrange equation
linearized about the actual axial ray is

\[
\frac{d\vartheta}{d\chi}=\frac{L}{BS_K^2E},\qquad
\frac{dL}{d\chi}=EBa^2q
\left(\frac\phi{A^2}+\frac\psi{B^2}\right)\vartheta.
\]

Indeed, \(\tfrac12\partial_\vartheta g_{ab}k^ak^b\) equals
\(a^2q[\phi(k^\eta)^2+\psi(k^\chi)^2]\vartheta\) at first order in the
angular variation. This is the original metric angular forcing; it is not a
substitution of only the Born transverse Weyl Hessian. With local angular
screen factor \(W=BS_K\), the physical Jacobi separation is
\(D=W\vartheta\). The two normalized screen vectors are parallel along the
axis, and the equivalent Sachs equation is

\[
\frac{d^2D}{d\lambda^2}=
\left[\frac{W''}{W}+
\frac{a^2q\{\phi(k^\eta)^2+\psi(k^\chi)^2\}}{W^2}\right]D.
\]

Primes in this equation are affine derivatives. Both the varying physical
screen scale and the angular metric forcing enter. Rotations about the axis
make both screen components identical, so the executed maps are scalar times
\(I_2\). This conditional isotropy is not a claim of an arbitrary-direction
shear null or of zero lensing response.

Observer matter-frame angles initialize
\(\vartheta_o=\gamma_o+s_o,L_o=0\); aberration then gives
\(D_o=0,D'_o=\omega_o=1\). The polar origin is initialized through its
regular series at \(\chi=10^{-7}\). For the lens-to-source map, initialization
is at the independently found, moving lens event:

\[
\vartheta_l=0,\qquad L_l=W_l\omega_l,
\quad\text{so}\quad D_l=0,\quad D'_l=\omega_l.
\]

The same globally normalized affine parameter is used for all segments. A
unit affine derivative at the lens without the \(\omega_l\) factor would
not represent a unit lens-frame angle. In particular,
\(D_{ls}\) is independently integrated and is never obtained by subtracting
observer distances.

The inherited seven-state implementation is also explicitly evaluated as an
independent first-order conformal Jacobi control. On the background ray,
\(\delta\vartheta''+2H_K\delta\vartheta'=q(\phi+\psi)/S_K^2\), where primes
now mean background \(\chi\) derivatives. Consequently
\(\kappa=\psi q-\delta\vartheta\) obeys the existing equation implemented by
`bhsm_optical_transfer.optical_rhs`, with

\[
F=(\psi q)''+2H_K(\psi q)'-\frac{q(\phi+\psi)}{S_K^2},
\qquad G_1=G_2=0.
\]

The evaluated, numerically stable harmonic form is

\[
F=q(\psi''-8K\psi)+2\psi'(q'+H_Kq)
+(\psi-\phi)q/S_K^2,
\]

with
\(N'=-aH\), \(N''=a^2H(H+H_N)\),
\(\psi'=\psi_NN'\), and
\(\psi''=\psi_{NN}(N')^2+\psi_NN''\). Observer aberration gives
\(\kappa_o=-v_{N,o}q'_o\) and
\(p_{\kappa,o}=q'_o(\psi_o-\phi_o)/2\). This control embeds the conformal
Jacobi pair in the original seven-state system; its frequency slot and source
are unused zeros. The full physical frequency is supplied by the metric ray,
so these zero control slots do not imply zero physical redshift response.
At each measured-redshift endpoint the resulting independent map control is

\[
\frac{\delta D_{ob}}{\bar D_{ob}}=
\delta N_b+H_K(\bar\chi_b)\delta\chi_b-\kappa_b.
\]

The full metric-derived endpoint shifts enter this last conversion, and all
six observer-to-lens/source unit derivatives agree with the independently
finite-differenced physical maps to maximum absolute
\(9.72\times10^{-11}\). This explicitly reuses the inherited optical forcing
and drift while retaining the complete metric ray as the primary calculation.

Reverse beams initialize \(L_b=-W_b\omega_b\) and are independently
integrated. The conserved two-solution Wronskian gives reciprocity
\(D_{ba}/D_{ab}=\omega_b/\omega_a\). The code checks this for observer--lens,
observer--source and lens--source. At the reverse polar origin the regular
momentum gives \(D_{bo}=-L(0)/E_o\); no divergent coordinate angle is
interpreted as a physical divergence.

## Interface to the separately counted galaxy encounter

The maps contain smooth cosmological propagation only. The root comparison
uses the same physical maps, in a common parallel screen, to compute

\[
\beta=D_{os}^{-1}X_s,\qquad
\beta=\theta-D_{os}^{-1}D_{ls}\hat\alpha(D_{ol}\theta),
\]
\[
T=(1+z_l)D_{ol}^{T}D_{ls}^{-1}D_{os},\qquad
ct=\tfrac12(\theta-\beta)^TT(\theta-\beta)
-(1+z_l)\hat\psi(D_{ol}\theta)+f(\beta).
\]

The galaxy potential and its deflection are defined in the local matter/lens
rest frame. The factor \(1+z_l=\omega_l/\omega_o\) converts lens-frame
proper delay to observer-frame proper delay. Because the physical lens-angle
normalization already appears in \(D_{ls}\), an additional independent
aberration or Doppler multiplier would count it twice. The arrival-time matrix
is tested for symmetry without arbitrary symmetrization.

The matched experiment fixes the observer event, sky axis, measured lens and
source redshifts, local galaxy parameters, and transverse physical source
displacement \(X_s\). It allows radial event locations, inferred angular
\(\beta\), and images to respond. “Fixed physical source” here specifically
means this transverse screen displacement; it does not additionally freeze a
three-dimensional coordinate position or a particular cosmological worldline.
The source term \(f(\beta)\) cancels between images of the same source.

### Exact nuisance factorization for the existing scalar-map controls

For the centered softened-isothermal family, write

\[
\beta=s/D_{os},\qquad\theta_E=bD_{ls}/D_{os},\qquad
\theta_c=r_c/D_{ol},\qquad
T=(1+z_l)D_{ol}D_{os}/D_{ls}.
\]

Its reduced potential is
\(\psi_{\rm red}(\theta)=\theta_E\sqrt{|\theta|^2+\theta_c^2}\).
Image coordinates depend on \((\beta,\theta_E,\theta_c)\), and relative
delays are \(T/c\) times the corresponding reduced Fermat differences.
The SIS sets \(r_c=0\). The source lies on the positive screen axis.
For a state coordinate \(a\), let
\(f_{ol,a}=\partial_a\ln D_{ol}\),
\(f_{os,a}=\partial_a\ln D_{os}\), and
\(f_{ls,a}=\partial_a\ln D_{ls}\).

**Proposition.** With fixed measured \(z_l\), dimensional transverse source
\(s\), physical \(b,r_c\), scalar propagation maps, and fixed noncaustic
image branches, the existing response satisfies

\[
B=N_{\min}R,\qquad
N_{\min}=[N_{{\rm source}\,x}\ N_{\ln b}\ N_{\ln r_c}\ N_{\rm MST}],
\]
\[
R=\begin{pmatrix}
-\beta_{\rm as}f_{os}\\
f_{ls}-f_{os}\\
-f_{ol}\\
f_{ol}+f_{os}-f_{ls}
\end{pmatrix}.
\]

Rows run across state coordinates. Delete the core-radius column and third
row for the SIS. The arcsecond normalization matches the saved source
nuisance; the other selected nuisance coordinates are dimensionless.
Differentiating the four reduced parameters gives these four rows.
At fixed maps, varying the source, physical strength and core reproduces
the first three variations. The implemented compensated mass-sheet column
has zero image response and multiplies every relative delay, reproducing
\(\delta\ln T\). Their chain rule gives the factorization.

If \(W^TW=C^{-1}\) for any positive-definite observational covariance,
and \(P\) is the orthogonal complement projector of
\(\operatorname{range}(WN_{\min})\), then \(PWB=0\) exactly for this
model family when those nuisances are unrestricted. Further line-of-sight
and \(H_0\) columns are unnecessary for this identity. These compensating
changes alter physical source/galaxy parameters; the proposition asserts
observable degeneracy, not physical equivalence of preparations.

`code/strong_lensing_degeneracy.py` applies this identity to the saved
`response_arrays.npz` and `results.json` without rerunning propagation.
The relative whitened Frobenius errors are \(2.452944603865752\times10^{-12}\)
(SIS) and \(6.377724303772925\times10^{-12}\) (smooth core), consistent
with the finite-differenced nuisance columns. The receipt is
`artifacts/strong_lensing_first_comparison_20261020/nuisance_factorization.json`.
The smooth-core covariance reconstructs the original runner's equal
\(0.001\) arcsec coordinate errors and \(0.1\) day delay errors. Its central
image has only \(2.7956547901766042\times10^{-5}\) of the positive outer
image's magnification, so this remains an idealized precision control.

The fixed-physical-source SIS identity is

\[
\Delta t=\frac{2(1+z_l)bs}{c}\frac{D_{ol}}{D_{os}}.
\]

It cancels a common \(H_0^{-1}\) distance rescaling exactly when the
dimensionless expansion/curvature history, measured redshifts, \(b,s\)
and two-image branch are fixed. It does not cancel arbitrary distance
responses or the response of a core with a fixed dimensional radius.
Holding the angular source instead gives the usual \(H_0^{-1}\) delay
scaling. A new intrinsic encounter or nonscalar map needs a separate
response and nuisance analysis; the proposition does not classify it.

Reproduce this focused saved-array check with

```powershell
python code/strong_lensing_degeneracy.py
python -m pytest tests/test_strong_lensing_degeneracy.py -q
```

Paraxial propagation about the axial central ray and an instantaneous weak
thin-plane galaxy encounter are the approximation regime. Across day-scale
delays, cosmological and galaxy parameters are frozen; evolution during those
delays, a time-dependent intrinsic galaxy potential, finite envelope
thickness, off-axis anisotropic transport, and additional interface physics
are not computed. A zero static seam export does not make these absent
constructions into demonstrated zero responses. The incoming optical beam,
cosmological radiation-fluid amplitudes, and muon/reference-state covariance
remain different objects.

## Executed numerical checks

The reusable implementation is `code/strong_lensing_transport.py`. A focused
verification command is

```powershell
$env:PYTHONPATH='code'
python -m pytest tests/test_strong_lensing_transport.py -q
```

Eight tests pass. They check dimensional harmonic normalization/eigenvalue,
all six physical anchor rows and original action constraints, historical-slice
recovery, analytic closed-FLRW maps, moving measured-redshift endpoints, all
three reverse reciprocity relations, unsymmetrized time-delay symmetry,
independent half-step responses, weak trial metrics, interpolation refinement,
the inherited seven-state conformal Jacobi application and the explicit
provider domain.

At \(z_l=0.3,z_s=1\), the zero-state maps in \(c/H_0\) are

| Map | Physical angle-to-separation length |
|---|---:|
| \(D_{ol}\) | 0.2124935890606791 |
| \(D_{os}\) | 0.3797948944897192 |
| \(D_{ls}\) | 0.24213223889447272 |

Their six physical-anchor unit derivatives are:

| Map | \(\zeta\) | \(\Pi_2\) | \(\delta_m^{ADM}\) | \(\delta_r^{ADM}\) | \(v_m^{ADM}\) | \(v_r^{ADM}\) |
|---|---:|---:|---:|---:|---:|---:|
| \(D_{ol}\) | 6.40540365e-6 | -2.35724358e-3 | -2.61199601e-5 | -7.98634076e-9 | 2.78253192e-4 | 1.18776829e-7 |
| \(D_{os}\) | 1.58658000e-5 | -5.29022437e-3 | -1.52761019e-4 | -4.40376208e-8 | 1.52101379e-3 | 1.05289388e-8 |
| \(D_{ls}\) | 2.32240515e-5 | -2.23022927e-3 | -1.84661299e-4 | -9.66363378e-8 | 1.82947417e-3 | 6.97220497e-7 |

Central finite trials use a nominal 0.002 step rescaled by each current metric
and velocity row. Independently solved half steps differ by at most
\(4.26\times10^{-11}\) in an individual map derivative; the largest relative
difference is \(7.91\times10^{-6}\), in a very small response. Across these
finite trials the largest relative metric perturbation is
\(2.15\times10^{-5}\), the largest reciprocity relative residual is
\(2.00\times10^{-15}\), and the largest fixed-redshift residual is
\(8.89\times10^{-16}\). The potential rows are used separately; their maximum
row-relative no-slip discrepancy is \(2.39\times10^{-10}\), an inherited
numerical diagnostic, not a new physical-admission test.

A mixed six-coordinate trial with maximum relative metric perturbation
\(3.77\times10^{-5}\) changes the maps by at most
\(2.78\times10^{-16}\) when the temporal interpolation grid is refined from
801 to 1601 points. Physical anchor reconstruction is verified to absolute
\(2\times10^{-12}\), and original action constraints to
\(3\times10^{-11}\). These tolerances verify the declared numerical
interfaces; they do not select a physical state, fix its normalization, or
replace any unchanged physical-admission guard in another preparation task.
