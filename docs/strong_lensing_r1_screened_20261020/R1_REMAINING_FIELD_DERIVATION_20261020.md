# Continuing the R1 field problem beyond the retained encounter

This note derives the next field and domain problem from the **same declared R1 action**. It does not add a microscopic BHSM coupling or a new lensing amplitude. The retained screened solution and common ordinary-GR 2PM control remain unchanged. The equations below show what must actually be applied to calculate scalar-dependent higher orders, a nonlinear environment, and a finite-domain cosmological continuation. They also identify why an exactly static isolated continuation of the rolling reference is not the problem already solved.

The source equations and conventions are those in `R1_LOCAL_COUPLING_DERIVATION_20261020.md`, the action in `manuscript/sections/04_scalar_metric_bridge.tex`, and `code/background.py`. The common-GR control is derived in `FINITE_BODY_GR_2PM_DERIVATION_20261020.md`. This is a bulk classical field derivation within that effective action, not a derivation of the open BHSM surface-to-R1 map.

## 1. Exact covariant applications

Use signature \((-+++)\), \(c=\hbar=1\), and \(\kappa=M_{\rm Pl}^{-2}=8\pi G\). The symbol \(\lambda\) remains the existing cubic coupling. Write

\[
 U_\mu=\nabla_\mu T,\quad Y=U_\mu U^\mu,\quad
 \Pi_{\mu\nu}=\nabla_\mu\nabla_\nu T,\quad S=\Box T,
 \quad g_3=\frac{2}{\lambda},\quad V=V_0-\mu^3T.
\]

Vary the retained action, without dropping the local gradient stress or metric time derivatives:

\[
 G_{\mu\nu}=\kappa\left(T^{m}_{\mu\nu}
 +U_\mu U_\nu-g_{\mu\nu}[Y/2+V]
 +T^{(3)}_{\mu\nu}\right),                                      \tag{1}
\]

\[
 T^{(3)}_{\mu\nu}=g_3\left[-S U_\mu U_\nu
 +U_\mu\Pi_{\nu\alpha}U^\alpha
 +U_\nu\Pi_{\mu\alpha}U^\alpha
 -g_{\mu\nu}\Pi_{\alpha\beta}U^\alpha U^\beta\right],           \tag{2}
\]

\[
 S+g_3\left[\Pi_{\mu\nu}\Pi^{\mu\nu}-S^2
 +R_{\mu\nu}U^\mu U^\nu\right]+\mu^3=0,
 \qquad \nabla_\mu T_m^{\mu\nu}=0.                              \tag{3}
\]

These equations retain the actual Ricci contraction and cubic stress; adding a GR correction to the photon potential does not execute them at the next order. In a weak-field expansion the first correction is obtained by substituting the leading metric and scalar into (1)--(3), retaining their next-order terms, and solving the resulting coupled metric/scalar equations with the same boundary/state data. The component system below supplies those applications without an unspecified preparation residual.

The scalar equation also has an exact current form:

\[
 J^\mu=U^\mu+g_3(\Pi^{\mu\nu}U_\nu-SU^\mu),\qquad
 \nabla_\mu J^\mu=V_T=-\mu^3.                                  \tag{4}
\]

Indeed,
\(\nabla_\mu(SU^\mu-\Pi^{\mu\nu}U_\nu)
=S^2-\Pi_{\mu\nu}\Pi^{\mu\nu}-R_{\mu\nu}U^\mu U^\nu\).
On the homogeneous reference, \(J^0=-[v+6Hv^2/\lambda]\); (4) gives the stored equation
\(\dot J_{\rm bg}+3HJ_{\rm bg}=\mu^3\).
Thus the source sign is fixed by the same background equation, rather than chosen in a local completion.

## 2. A concrete spherical continuation with time dependence

On a local spherical patch where the areal radius is spacelike, use polar-areal coordinates

\[
 ds^2=-N(t,r)^2dt^2+A(t,r)^2dr^2+r^2d\Omega^2,
 \quad \nu=\ln N,\quad \gamma=\ln A,\quad
 f=A^{-2}=1-\frac{2m(t,r)}r.                                   \tag{5}
\]

Here \(m\) is a mass **length**, not the physical mass. This is a local coordinate choice; a cosmological boundary must be pulled back into these coordinates. A homogeneous FLRW scalar in cosmic time need not be a function of the polar time alone.

The field unknowns are \(m,\nu,T\), together with the material number density \(n\) and radial fluid velocity \(\beta\) if the matter is evolved. Define orthonormal covariant scalar derivatives

\[
 E=\frac{T_t}{N},\qquad P=\frac{T_r}{A},\qquad Y=-E^2+P^2,
\]

and the actual Hessian applications

\[
 H_{00}=N^{-2}\left[T_{tt}-\nu_tT_t
                    -\frac{N^2}{A^2}\nu_rT_r\right],
\]
\[
 H_{01}=(NA)^{-1}\left[T_{tr}-\nu_rT_t-\gamma_tT_r\right],
\]
\[
 H_{11}=A^{-2}\left[T_{rr}-\frac{A^2}{N^2}\gamma_tT_t
                         -\gamma_rT_r\right],\qquad
 H_\Omega=A^{-2}\frac{T_r}{r}.                                  \tag{6}
\]

Then

\[
 S=-H_{00}+H_{11}+2H_\Omega,\quad
 \Pi_{\mu\nu}\Pi^{\mu\nu}
   =H_{00}^2-2H_{01}^2+H_{11}^2+2H_\Omega^2,
\]
\[
 K\equiv\Pi_{\mu\nu}U^\mu U^\nu
       =H_{00}E^2-2H_{01}EP+H_{11}P^2.                          \tag{7}
\]

Let \(\rho,p_r,p_\perp,\tau\) denote the total orthonormal components
\(T_{00},T_{11},T_{22},T_{01}\). Applying (2) gives the explicit scalar contributions

\[
 \rho_T=\tfrac12(E^2+P^2)+V
       +g_3[H_{11}Y-2H_\Omega E^2],
\]
\[
 p_{r,T}=\tfrac12(E^2+P^2)-V
       +g_3[H_{00}Y-2H_\Omega P^2],
\]
\[
 p_{\perp,T}=\tfrac12(E^2-P^2)-V-g_3K,\qquad
 \tau_T=EP+g_3[H_{01}Y-2H_\Omega EP].                           \tag{8}
\]

For a perfect material fluid with rest-frame density \(\rho_m\), pressure \(p_m\), and \(\Gamma=(1-\beta^2)^{-1/2}\), add

\[
 \rho_m^{\rm frame}=(\rho_m+p_m)\Gamma^2-p_m,\quad
 p_{r,m}^{\rm frame}=(\rho_m+p_m)\Gamma^2\beta^2+p_m,
\]
\[
 p_{\perp,m}^{\rm frame}=p_m,\qquad
 \tau_m=-(\rho_m+p_m)\Gamma^2\beta.                             \tag{9}
\]

The sign convention is explicit: a positive outward material energy flux has \(\tau_m<0\).

The spherical Einstein applications are

\[
 m_r=\frac{\kappa r^2}{2}\rho,\qquad
 \nu_r=\frac{m+\kappa r^3p_r/2}{r^2f},\qquad
 m_t=\frac{\kappa r^2}{2}\frac{N}{A}\tau,                        \tag{10}
\]

\[
 A^{-2}\left[\nu_{rr}+\nu_r^2-\nu_r\gamma_r
                  +\frac{\nu_r-\gamma_r}{r}\right]
 -N^{-2}\left[\gamma_{tt}+\gamma_t^2-\nu_t\gamma_t\right]
 =\kappa p_\perp.                                             \tag{11}
\]

The angular row is a consistency/evolution equation; which rows are independently evolved depends on the constrained evolution formulation. They must not be replaced by the exact GR-fluid rows while the scalar stress is changed.

The exact scalar application can use the Einstein trace-reversal directly:

\[
 R_{\mu\nu}U^\mu U^\nu=
 \frac{\kappa}{2}\left[
  (\rho+p_r+2p_\perp)E^2-4\tau EP
  +(\rho+p_r-2p_\perp)P^2\right].                              \tag{12}
\]

Insert (6)--(9) and (12) into

\[
 -H_{00}+H_{11}+2H_\Omega
 +g_3\left[H_{00}^2-2H_{01}^2+H_{11}^2+2H_\Omega^2
            -(-H_{00}+H_{11}+2H_\Omega)^2
            +R_{\mu\nu}U^\mu U^\nu\right]+\mu^3=0.           \tag{13}
\]

Equations (6)--(13), not a GR2 optical insertion, determine the next scalar/metric field response. The lapse row and scalar acceleration are coupled because \(p_{r,T}\) contains \(H_{00}\). The momentum row is coupled because \(\tau_T\) contains \(\gamma_t\) through \(H_{01}\). The admissibility/invertibility of that implicit system must be checked on its actual solution; positivity of the retained diagonal scalar coefficients alone is not that check.

For material evolution one explicit completion is a declared barotropic relation \(\rho_m=\rho_m(n)\), \(p_m=n\,d\rho_m/dn-\rho_m\). With this relation the actual continuity and Euler applications are

\[
 \partial_t(Ar^2n\Gamma)+\partial_r(Nr^2n\Gamma\beta)=0,
\]
\[
 (\rho_m+p_m)\left[
 \Gamma^2\left(\frac{\beta_t}{N}+\frac{\beta\beta_r}{A}\right)
 +\frac{\nu_r}{A}+\frac{\beta\gamma_t}{N}\right]
 +\frac{\beta\,\partial_t p_m}{N}+\frac{\partial_r p_m}{A}=0.    \tag{14}
\]

This relation is **not supplied by the prescribed density alone**. The original comparison intentionally fixes the same density and adjusts pressure at leading order. For a genuinely common fixed barotropic EOS and \(\delta\rho=0\), one would require \(\delta p=(dp/d\rho)\delta\rho=0\), whereas the retained scalar force gives
\(\delta p(r)=\int_r^R\rho\,\delta\Phi' ds>0\).
The next physical material sequence must therefore allow density readjustment, specify extra support, or explicitly retain a prescribed-stress comparison. This is the particular material defining relation consumed by (14); it is not evidence that the retained prescribed-density calculation is numerically incomplete.

## 3. Why freezing the rolling state does not produce an exact static extension

The current components from (4) in the same coordinates are

\[
 J^t=N^{-1}\{-E+g_3[(H_{11}+2H_\Omega)E-H_{01}P]\},
\]
\[
 J^r=A^{-1}\{P+g_3[-H_{01}E+(H_{00}-2H_\Omega)P]\},
\]

and the exact scalar balance is

\[
 \partial_t(NAr^2J^t)+\partial_r(NAr^2J^r)
       =-\mu^3NAr^2.                                          \tag{15}
\]

For an exactly static metric and \(T=vt+\pi(r)\) with constant \(v\), the left temporal derivative vanishes. Regularity then gives

\[
 J^r(r)=-\frac{\mu^3}{NAr^2}\int_0^r N(s)A(s)s^2ds.
                                                                    \tag{16}
\]

In an asymptotically flat static patch this grows as \(-\mu^3r/3\), rather than approaching the decaying isolated perturbation. Independently, \(V(T)=V_0-\mu^3vt-\mu^3\pi\) makes the canonical stress time dependent when \(\mu^3v\ne0\). The nonzero \(\tau_T\) in (8) must also satisfy the momentum row if the metric is static. The retained rolling R1 reference has nonzero \(v\), \(d\), \(H\) and \(\mu^3\). An exactly static full-action solution with static matter and the same isolated boundary data does not follow by appending higher powers of its static potential.

The homogeneous time-current balance must instead be retained and subtracted consistently, or the full coupled time-dependent problem must be evolved. The existing local solution does the former at its leading quasi-static order. Equation (16) identifies the term lost if that approximation is mistaken for an exact static solve.

The saved values quantify coefficient variation during a screening-radius light-crossing time without producing an error bar:

| Parameter from the saved reference | Value |
|---|---:|
| \(H_0r_V/c\) | \(4.766887685833153\times10^{-5}\) |
| \(H_lr_V/c\) | \(5.642138685209712\times10^{-5}\) |
| \(|d/v|r_V/c\) | \(4.981553281754921\times10^{-5}\) |
| \(|\dot b_T/b_T|r_V/c=2|d/v|r_V/c\) | \(9.963106563509841\times10^{-5}\) |

Here the last two expressions use the physical derivatives represented by the saved hatted coefficients. These are small local expansion parameters, not measured uncertainties or bounds on time-delay errors. First-order temporal effects can also cancel under a particular encounter symmetry; that cancellation must be derived from the full application, not assumed. No existing lens result is recomputed by this table.

## 4. The nonlinear environmental equation already follows from the retained action

A projected optical convergence/shear does not determine the three-dimensional scalar Hessian. To see the exact dependency at the retained quasi-static order, write
\(\pi=\pi_L+\pi_E\), \(\delta\rho=\rho_L+\rho_E\), and require the environmental field itself to solve the same equation. Subtract that environmental equation from the total equation. The local field then obeys

\[
 C\Delta\pi_L-\frac{2}{\lambda}
 [ (\Delta\pi_L)^2-\Pi_L:\Pi_L]
 -\frac{4}{\lambda}
 [ (\Delta\pi_E)(\Delta\pi_L)-\Pi_E:\Pi_L]
 =-\frac{b_T}{2M_{\rm Pl}^2}\rho_L,                            \tag{17}
\]

where \(\Pi_{L,ij}=\partial_i\partial_j\pi_L\), and similarly for \(E\). Thus the actual environmental contraction consumed by the local problem is

\[
 [(\operatorname{tr}\Pi_E)\delta^{ij}-\Pi_E^{ij}]
     \partial_i\partial_j\pi_L,\qquad
 Z_E^{ij}=C\delta^{ij}-\frac4\lambda
 [(\operatorname{tr}\Pi_E)\delta^{ij}-\Pi_E^{ij}].              \tag{18}
\]

A constant or uniform-gradient scalar environment has no Hessian and no such leading cross term. A nonzero tidal Hessian generally does. Its trace alone is insufficient, and a two-dimensional optical tide does not specify it. Once physical environmental density, Cauchy/boundary data and the same R1 coefficients are supplied, the defining relation is the joint total scalar equation together with
\(2M_{\rm Pl}^2\Delta\Phi=\rho_L+\rho_E-b_T\Delta(\pi_L+\pi_E)\).
Equation (17) is the deterministic conditional problem for an environmental continuation; no independent Gaussian or stochastic state law is needed to derive it. The nuisance tide calculation remains a distinct optical inference calculation.

## 5. Matching data and the compact-domain obstruction to a naive isolated export

Choose an overlap worldtube with \(r_V\ll r_m\ll c/H_l\). The regular branch already fixes its leading far-zone coefficients:

\[
 \pi_T'(r_m)=-\frac{b_TGM}{Cr_m^2}
       [1+O((r_V/r_m)^3)],
\]
\[
 \Phi'(r_m)=\Psi'(r_m)=\frac{GM}{r_m^2}
       [1+\epsilon_\infty+O((r_V/r_m)^3)].                     \tag{19}
\]

There is no new free local scalar charge. The corresponding \(1/r\) terms are \(\pi_T=\pi_0+b_TGM/(Cr)\) and \(\Phi=\Phi_0-(1+\epsilon_\infty)GM/r\). Their additive values and clock normalization must be fixed by the outer state. They cannot be chosen independently in a finite-distance arrival-time comparison.

For a smooth artificial matching surface, impose continuity of the induced metric and extrinsic curvature, and smooth scalar value and first normal derivative. Evaluate the normal current from (4), rather than transplanting the quasi-static radial flux into a moving boundary. Smooth matching with no surface source must satisfy the integrated current balance (15). If a physical surface stress or scalar charge is admitted instead, its action-derived metric/scalar junction conditions are required; the present note does not replace those with an assumed moving-seam law.

The outer problem is the linear or nonlinear perturbation of the **same** rolling, closed background with these worldtube data and its own initial state. A localized source generally requires more than the stored six-amplitude axial harmonic. Its temporal and homogeneous components must satisfy the full Einstein and scalar constraints.

A simple compatibility calculation explains why a flat isolated profile is not already a global solution. If the retained local Poisson equation were extended unchanged to a compact slice with a constant \(b_T\), integration would give

\[
 0=\int_{S^3}\Delta\Phi\,dV
 =\frac{1}{2M_{\rm Pl}^2}\int_{S^3}\delta\rho\,dV
   -\frac{b_T}{2M_{\rm Pl}^2}\int_{S^3}\Delta\pi_T\,dV,
\]

so that \(\int\delta\rho\,dV=0\). The prescribed positive finite excess with zero excess outside does not satisfy that extended equation. One must solve the actual curvature/time-dependent outer constraint, allow the homogeneous response, or specify compensating excess elsewhere; merely putting the isolated flat equation on a closed grid cannot supply the match. This is a compatibility test of that naive extension, **not** a theorem that the full closed Einstein system forbids positive local mass or universally requires a negative physical density. Negative *excess* would mean underdensity relative to the background, if compensation is the chosen outer state.

Similarly, integration of (4) on a closed Cauchy slice fixes the total scalar-current evolution. It is not permissible to impose zero scalar flux everywhere while discarding the background time-current source. The open global dependency is therefore a specific constrained Cauchy/matching problem, not an unavailable array.

## 6. The next optical application once that field problem is solved

For the spherical metric (5), the null Hamiltonian and its actual ray applications are

\[
 \mathscr H=\tfrac12[-p_t^2/N^2+p_r^2/A^2+L^2/r^2]=0,
\]
\[
 \dot t=-p_t/N^2,\quad \dot r=p_r/A^2,\quad \dot\phi=L/r^2,
\]
\[
 \dot p_t=-\nu_t p_t^2/N^2+\gamma_t p_r^2/A^2,
\qquad
 \dot p_r=-\nu_r p_t^2/N^2+\gamma_r p_r^2/A^2+L^2/r^3.          \tag{20}
\]

Dots here denote affine ray differentiation, not cosmic-time differentiation. Photon frequency is \(-u^\mu p_\mu\) on the actual source/observer worldlines. The same emitted marker and observer proper-clock protocol determine relative arrival times after matching the rays into the outer metric. Time dependence removes the static conserved \(p_t\), so a complete finite-distance 2PM prediction is not recovered solely by integrating an asymptotic deflection coefficient.

The existing common-GR2 exercise consistently adds one known ordinary-GR term to both leading encounters. It includes optical image-response cross effects; it does not supply the missing scalar-dependent field terms in (8), (10)--(13), their material response in (14), their environment in (17), or their outer match. Conversely, none of those unfinished terms invalidates the executed leading solution or establishes that the positive leading effective-mass equivalence is broken. Their calculation requires the same source, material and global-state prescriptions, with the equivalent scalar/effective-mass contribution accounted for once.

## Derived progress and precise remaining choices

The exact bulk component equations, conserved-current balance, environmental cross-Hessian and overlap coefficients are now explicit. A static frozen-background extrapolation fails the exact current balance as shown in (16); the correct continuation retains time dependence. An evolved finite body consumes a material relation in (14), and an actual environmental/global application consumes outer Cauchy and matching data in (17)--(19). The original leading comparison fixes neither a common material EOS nor a complete global compensation/homogeneous state. Those are specific additional physical prescriptions; no vacuum, factorization, new optical amplitude or state-population covariance has been silently selected here. The independently open microscopic BHSM source reduction remains the subject of the separate surface-to-R1 derivation.

## Algebra checks for this continuation derivation

New symbolic checks, separate from the existing numerical replay, verify the four cubic stress projections, two scalar-current components and trace-reversed Ricci contraction by direct matrix contraction. An independent Christoffel/Ricci construction verifies the spherical Einstein tt, rr, tr and angular rows. A separate fluid-acceleration projection verifies the radial Euler row. All twelve identities reduce to zero in SymPy. These checks validate the component algebra; they do not supply a material EOS, global state or evolved solution. The receipt is `artifacts/strong_lensing_bhsm_integration_20261020/remaining_field_algebra_check.json`.
