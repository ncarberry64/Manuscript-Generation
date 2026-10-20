# Common-source reduction: executed application and a nonuniqueness theorem

This continuation works out the source relation consumed by the lens calculation, its material matching, and what the retained action data can determine. It does not promote an assumed microscopic identification to a derived one. The outcome is a nonzero evaluated effective source, zero additional material layer in the specified smooth problem, and a constructive proof that the leading R1 normalization does not determine a microscopic fourth-order coefficient or cutoff.

The original finite-body inputs, incident light, clocks, maps, action parameters, delay outputs and constant-target null are unchanged. No new amplitude is fitted. The two complement models below are mathematical witnesses with explicit boundary data; neither is selected as a galaxy state or admitted as a completed BHSM preparation. No quantum vacuum or statistical kernel is assigned.

## 1. The exact contraction consumed after metric elimination

For the action already used in the manuscript,

\[
 S=\int\sqrt{-g}\left[\frac{M_{\rm Pl}^2}{2}R
 -\frac12(\nabla T)^2+\frac{(\nabla T)^2\Box T}{\lambda}
 -V(T)+\mathcal L_m\right],
\]

write \(U_\mu=\nabla_\mu T\) and \(\Pi_{\mu\nu}=\nabla_\mu\nabla_\nu T\). Integrating the variation of the cubic term twice gives

\[
 E_T=\Box T+\frac2\lambda
 [\Pi_{\mu\nu}\Pi^{\mu\nu}-(\Box T)^2+R_{\mu\nu}U^\mu U^\nu]-V_T.
\]

The cancellation uses \(\Box(U^2)=2\Pi^2+2U\cdot\nabla\Box T+2R_{\mu\nu}U^\mu U^\nu\). On the homogeneous field it gives
\(\dot J+3HJ=\mu^3\), \(J=v+6Hv^2/\lambda\), fixing the same sign as the inherited background.

The cubic stress is

\[
 T^{(3)}_{\mu\nu}=-\frac2\lambda(\Box T)U_\mu U_\nu
 +\frac4\lambda U_{(\mu}\Pi_{\nu)\alpha}U^\alpha
 -\frac2\lambda g_{\mu\nu}\Pi_{\alpha\beta}U^\alpha U^\beta.
\]

It yields \(\rho_3=6Hv^3/\lambda\), \(p_3=-2v^2d/\lambda\). The leading static perturbation has \(\delta T^{(3)}_{00}=-b_T\Delta\pi_T\), \(b_T=2v^2/\lambda\). The background cross term in the scalar Hessian invariant gives \(A=1+4(d+2Hv)/\lambda\).

To expose a possible additional surface/cap source without choosing one, **define** \(\rho_\Sigma\) as its contribution to the leading metric energy-density row and \(J_T\) as its contribution to the independent scalar Euler row, with the signs below. The original applications are then

\[
 2M_{\rm Pl}^2\Delta\Phi
 =\rho+\rho_\Sigma-b_T\Delta\pi_T,
\]
\[
 A\Delta\pi_T-\frac2\lambda
 [ (\Delta\pi_T)^2-|\operatorname{Hess}\pi_T|^2]
 +b_T\Delta\Phi=J_T.
\]

Eliminating this same metric row gives

\[
 C\Delta\pi_T-\frac2\lambda
 [ (\Delta\pi_T)^2-|\operatorname{Hess}\pi_T|^2]
 =\underbrace{J_T-\frac{b_T}{2M_{\rm Pl}^2}(\rho+\rho_\Sigma)}_{J_{\rm eff}},
 \qquad C=A-\frac{b_T^2}{2M_{\rm Pl}^2}.
 \tag{1}
\]

This is the precise consumed contraction. An exported surface charge is not simply another density or an extra lensing amplitude. Its scalar and metric variations must be projected into the two slots in (1). A leading traceless spatial stress would also require the corresponding slip equation; the retained \(\Phi=\Psi\) cannot be assumed for an arbitrary cap/interface source.

In the declared pure-cubic finite-body problem, no additional layer action occurs: \(\rho_\Sigma=J_T=0\), while the nonzero volume source is \(-b_T\rho/(2M_{\rm Pl}^2)\). This follows from that stated effective action and boundary problem, rather than from assigning an unknown microscopic source its vacuum value.

## 2. Nonlinear canonical flux and a source application

Equation (1) has divergence form with

\[
 \mathbf F_{\rm R1}=C\nabla\pi_T-\frac2\lambda
 [(\Delta\pi_T)\nabla\pi_T
       -(\operatorname{Hess}\pi_T)\nabla\pi_T],
 \qquad \nabla\cdot\mathbf F_{\rm R1}=J_{\rm eff}.
 \tag{2}
\]

For spherical symmetry this reduces to the evaluated original flux,

\[
 r^2 F_{\rm R1}^r=Cr^2\pi_T'-\frac4\lambda r(\pi_T')^2
 =\int_0^r s^2J_{\rm eff}(s)\,ds.
 \tag{3}
\]

For a proper-radial delta layer, \(\rho_\Sigma=\sigma_\Sigma\delta(r-R)\), \(J_T=j_{T,\Sigma}\delta(r-R)\), the retained weak-metric jump is

\[
 [r^2 F_{\rm R1}^r]_{R^-}^{R^+}
 =R^2\left[j_{T,\Sigma}
       -\frac{b_T}{2M_{\rm Pl}^2}\sigma_\Sigma\right].
 \tag{4}
\]

Multiplying by \(4\pi\) restores the integrated spherical charge. This is a jump of the nonlinear flux, not just a jump of \(\pi_T'\).

The new application reads the **existing frozen coefficient receipt**, evaluates its already-defined scalar sections, and integrates the density source independently. In \(M_{\rm Pl}=H_0=c=1\), with \(r\) measured in \(c/H_0\), the declared density is

\[
 \widehat\rho(r)=\frac{15\widehat m}{\widehat R^3}
 [1-r^2/\widehat R^2],\qquad
 \int4\pi r^2\widehat\rho\,dr=8\pi\widehat m,
 \quad \widehat m=\frac{GM/c^2}{c/H_0}.
\]

Actual results in `source_application.json` are:

| Application | Result |
|---|---:|
| \(\widehat R\) | \(6.744666004906636\times10^{-7}\) |
| \(\widehat m\) | \(1.0758677182033066\times10^{-12}\) |
| \(\pi_T'(R)/(M_{\rm Pl}H_0)\) | \(-2.119825778845613\times10^{-4}\) |
| \(\int_0^R r^2J_{\rm eff}\,dr=-b_T\widehat m\) | \(-1.216417408185084\times10^{-13}\) |
| Maximum relative residual of (3), at \(r/R=0.1,0.5,1,2,10\) | \(4.15\times10^{-16}\) |
| Original metric Gauss residual at those sections | zero at printed precision |
| Relative material-surface flux jump in (4) | \(6.23\times10^{-16}\) |
| Extra effective enclosed mass at \(R\) | \(506705.8798\,M_\odot\) |
| Asymptotic extra effective mass | \(1.5078862438\times10^8\,M_\odot\) |

The interior direct quadratic root and exterior computed rational root are compared separately at the surface. Thus the zero layer accompanies an evaluated nonzero volume source and scalar section.

For the unchanged exact-GR material control, same-oriented \(+dr\) limits give

\[
 K_{uu}^{\rm in}=K_{uu}^{\rm out}
 =-\frac{m_g}{R^2\sqrt{1-2m_g/R}}
 =-0.0005317137194453\ {\rm Mpc}^{-1}.
\]

Consequently its material-worldtube \(\Delta K_{uu}=0\) and the stated material \(Q_\Sigma\) functional vanishes. This is the smooth material junction, **not** the cosmological cap equator, an evaluated active encapsulation event, or an exact modified-gravity TOV solution. The R1 matching checked here is its leading scalar flux and metric Gauss matching. The null target sector and the rolling scalar remain distinct.

The R1 metric itself is evaluated separately at its retained order:
\(K_{uu}^{\rm R1}=-\Phi'(R)=-0.0005317155655079203\ {\rm Mpc}^{-1}\)
on both material limits, so its leading jump is also zero. This different
value is not replaced by the ordinary-GR extrinsic curvature above. The
statement concerns the leading weak metric and does not combine an exact
GR fluid with modified scalar stress as a purported exact solution.

## 3. What the original cap action actually fixes

The source trail does not end with the old symbolic scaffold. The retained v7.1 `master_action/reduction.py` supplies:

\[
 \mathcal R_{85}=(\pi_!,P_{\rm ret},Q_H),\quad
 d\mu_F=a_F^3\eta_1\wedge\eta_2\wedge\eta_3,\quad
 V_F=16\pi^2a_F^3,
\]

and a constrained cap critical-value/trace reduction

\[
 S_{4,\rm cap}[h,\varphi,k]
 =\operatorname{Crit}_{q_5:\operatorname{Tr}q_5=\varphi,\ C(q_5)=0,\ k\ retained}
 S_{5|4}[q_5,h].
 \tag{5}
\]

Only the closed-range complement is eliminated. Gauge and constraint directions are quotiented as owned by the domain; Lyapunov--Schmidt coordinates \(k\) remain explicit. The construction does not make the intrinsic M4 Wilson action equal to a literal copy of the fiber-pushed M8 action. A varying fiber scale also generates radion/connection terms, so it cannot be replaced silently by constant \(V_F\).

The coefficient-locked Einstein/GHY variation supplies cap canonical momentum. To fix its normalization without guessing a density convention, define

\[
 \delta S_{5,\epsilon}^{\rm on-shell}
   =\int_{B_1}\Pi_\epsilon^{ab}\,\delta\gamma_{\epsilon ab},
 \quad
 S_{\rm match}=\sum_\epsilon\int_{B_1}
 \Lambda_\epsilon^{ab}(h_{ab}-\gamma_{\epsilon ab}).
\]

The two variations give \(\Pi_\epsilon^{ab}=\Lambda_\epsilon^{ab}\) and

\[
 2C_\partial G^{ab}[h]
 =T_m^{ab}+\tau_{\rm cap}^{ab},\qquad
 \tau_{\rm cap}^{ab}=\frac{2}{\sqrt{-h}}
                 \sum_\epsilon\Pi_\epsilon^{ab}.
 \tag{6}
\]

The sign follows from varying the **covariant** \(h_{ab}\) in the displayed matcher convention. The familiar density-free Brown--York form is proportional to the oriented \(\kappa_1(K^{ab}-Kh^{ab})\); its coefficient and density must be converted using this variational definition, not guessed from a symbol. The reflection-even mean recorded in `mass_curvature_response.py` is half the cap sum and cannot be inserted for that sum without restoring its factor.

Equation (6) identifies the actual metric-source contraction \(u_a u_b\delta\tau_{\rm cap}^{ab}\) once the physical observer and cap solution are supplied. The cap scalar normal momentum is separately \(-\sqrt{-h}Z_5 n\cdot D\sigma\). It is not a canonical R1 \(T\) source merely because both fields are called scalars. The intrinsic coefficient \(C_\partial\), cap response and any absorbed local Einstein term must match \(M_{\rm Pl}^2/2\) together. Counting the same cap response both as \(\tau_{\rm cap}\) and again as an independent R1 stress would double count it.

At the round equatorial inclusion the source supplies \(K_{ab}=0\) and therefore zero background cap Brown--York response. That is a defined geometric application, not a galaxy-source solution or a proof that all perturbed cap momenta vanish. We do not insert the galaxy R1 metric into that round cap and claim the cap Einstein/scalar equations have thereby been solved.

The later owner-authorized interface action source is also explicit about the relevant defining relation:

\[
 \Delta_{\rm enc}=\Pi_c+C(F_B)^*\Pi_e
 =-[E_{q_c}(S_{\rm enc})+C(F_B)^*E_{q_e}(S_{\rm enc})].
 \tag{7}
\]

It supplies the Euler operators for a selected invariant density \(W_s\), while leaving those universal constitutive functions unselected. Equation (7) cannot be replaced with \(Q_\Sigma\), a scalar vacuum, or a numerical lensing amplitude. It must first be pulled back into the independent scalar/metric slots of (1).

## 4. Reduction changes the source as well as the operator

On an owned physical complement the original two block equations are

\[
 A\psi+Bq=J_P,\qquad B^\dagger\psi+Dq=J_Q.
\]

Solving the second application and substituting it into the first gives

\[
 (A-BD_\perp^{-1}B^\dagger)\psi
 =J_P-BD_\perp^{-1}J_Q-Bq_0.
 \tag{8}
\]

Here \(q_0\in\ker D\) stays explicit and obeys its projected compatibility equation. The inverse is not an unrestricted pseudoinverse. This source correction is required alongside the Schur Hessian in (5); a reduced kernel alone does not fix the forcing.

To identify the retained field with R1 \(T\), the same reduction must supply the physical profile/embedding \(e_T\), its kinetic pairing \(Z_T\), its cubic derivative vertex and the source pullback. If a valid scalar coordinate \(\tau\) is supplied, \(T=\sqrt{Z_T}\tau\) makes the cubic coefficient \(c_3/Z_T^{3/2}\). This is a necessary normalization identity, not a license to choose \(Z_T\), \(c_3\), or the profile to reproduce the frozen \(\lambda\). A time-dependent profile or moving domain requires the full symplectic pullback as well.

## 5. Two explicit positive lifts prove that B and the cutoff are not determined

Use the already-retained quadratic family in `03_reduced_response.tex`, on one common physical two-scalar domain. With \(q=p\ell_0^2\), take

\[
 H(q)=\begin{pmatrix}Z_0q&g\sqrt q\\g\sqrt q&M^2+cq\end{pmatrix},
 \qquad K_{\rm time}=\operatorname{diag}(D,1),
\]

where \(C=4.238858559433196\) and \(D=4.396221028708471\) are the already-frozen R1 coefficients. Choose **mathematical witnesses**, not galaxy parameters,

\[
 \begin{array}{c|cccc}
 &M^2&c&g&Z_0\\\hline
 A&1&1&1&C+1\\
 B&4&1&2&C+1
 \end{array}
\]

Both have the same leading \(Z_0-g^2/M^2=C\), the same temporal coefficient \(D\), and the same retained source \(J_P\) with explicitly unforced complement \(J_Q=0\) for this counterexample. An unforced classical boundary problem is not a quantum vacuum prescription.

The common domain is a ball of radius \(L=10\pi\ell_0\), with both real fields regular at the origin and Dirichlet at \(L\). Its normalized radial basis is

\[
 f_n(r)=\frac{\sin(n\pi r/L)}{r\sqrt{2\pi L}},\quad
 \int_0^L4\pi r^2 f_n^2\,dr=1,\quad q_n=(n/10)^2.
\]

For every \(q>0\),

\[
 D_{QQ}=M^2+cq>0,\qquad
 K_{\rm red}(q)=Cq+\frac{g^2c q^2}{M^2(M^2+cq)}>0.
\]

Thus the full spatial block is positive and the Legendre matrix is positive. The value \(q=0\) is a formal limit excluded by this Dirichlet ball. A constant mode, if admitted on a different physical domain, would require separate kernel handling rather than an inverse. At \(n=1,\ldots,20\), the minimum spatial eigenvalues are 0.0420574 and 0.0423072; the kinetic eigenvalues are \(1,D\).

Nevertheless,

\[
 \frac{B_A}{\ell_0^2}=\frac1C=0.2359125661,\qquad
 \frac{B_B}{\ell_0^2}=\frac1{4C}=0.05897814152,
\]

after dividing by the same leading spatial coefficient. The expansion conditions are respectively \(p\ell_0^2\ll1\) and \(p\ell_0^2\ll4\). These are different exact spectra/operators with identical leading R1 normalization. The common reference unit \(\ell_0\) has not been identified with a physical kpc scale.

This is a constructive underdetermination proof, not an inference from an unavailable numerical array. Leading \(C,D\) and a normalized low-energy source do not determine the complement gap, mixing, higher derivative coefficient or physical cutoff. The finite density profile also has a surface derivative discontinuity, so a physical derivative expansion would need control of its weighted high-momentum tail rather than a universal sharp maximum momentum inferred from one image impact.

## 6. Why the old static seam law cannot be identified directly with R1

The inherited neighborhood law is a distributed exterior source,

\[
 \nabla\cdot\mathbf F_{\rm top}
 =a\frac{\Theta(r-r_\star)}{r^2},\qquad
 a=\epsilon_\Sigma Q_\Sigma/Q.
\]

With its stated activation condition \(F_{\rm top}^r(r_\star)=0\), direct integration gives

\[
 r^2F_{\rm top}^r=a(r-r_\star),\qquad
 QF_{\rm top}^r=\epsilon_\Sigma Q_\Sigma
       (1/r-r_\star/r^2).
 \tag{9}
\]

It is not a delta layer: its integrated source is \(4\pi a(r-r_\star)\). The dimensions also distinguish \(Q_\Sigma\) (speed squared), a volume source, and a canonical R1 scalar charge. In contrast, outside the finite R1 body equation (3) has constant charge \(-b_TGM\), and both its scalar flux and metric correction decay as \(r^{-2}\) in the linear exterior.

For nonzero \(a\), no constant field normalization identifies the \(r^{-1}\) tail in (9) with the isolated R1 exterior on an unbounded domain. The result rejects that direct identification only. A genuinely owned finite-support export, different boundary problem or nonlocal source/metric map could still relate the sectors. None is silently chosen here, and the material radius \(R\) is not declared to be the distinct seam radius \(r_\star\).

## 7. Review conclusion and exact next derivation

The retained effective source is now explicitly reduced, applied and checked. Its material matching is evaluated. The operator/source Schur reduction and two positive same-domain lifts resolve the mathematical question of whether the leading R1 inputs uniquely fix the microscopic source normalization and cutoff: **they do not**.

The next physical input is correspondingly precise: an action-owned embedding and boundary/state prescription for the retained scalar/cap/interface sector, including the source pullback in (8), whose canonical variation reproduces the frozen R1 \(T\) kinetic and cubic terms and whose complement spectrum controls its range of validity. Equations (5)-(8) supply the actual defining relations; the physical profile and constitutive interface density required to evaluate them are not selected by the presently adopted low-energy action. This is not repaired by identifying \(T\), the target phase, and the metric coordinate, nor by appending a fitted \(B\).

The nonlinear environmental equation, scalar higher orders and global matching can be continued only within a declared common action/domain. These conclusions preserve the successful leading optical equivalence and the fixed-body offsets without claiming a complete physical BHSM source construction.

Reproduce the new application and focused checks:

```text
python scripts/run_strong_lensing_surface_source.py
python -m pytest -q tests/test_strong_lensing_surface_source.py
```

Seven focused tests pass. They test the consumed coupled source, nonzero complement forcing, regular nonlinear flux, delta-layer normalization, distributed-source support, positive same-domain nonuniqueness and the actual computed finite-body sections. The nine earlier reproduction commands are reused and are not rerun by this continuation.

## Source trail

- `manuscript/sections/04_scalar_metric_bridge.tex`: retained covariant R1 action and background normalization.
- `manuscript/sections/03_reduced_response.tex`: exact reduced block and conditional two-field complement family.
- `manuscript/sections/04ao_bhsm_optical_transfer.tex` through `04ar_bhsm_optical_current_gate.tex`: inherited source, metric Hessian and target-current distinctions.
- Original BHSM `src/bhsm/interface/master_action/reduction.py`: oriented fiber pushforward, cap trace critical value, KKT domain and source-compatible reduction.
- Original BHSM `src/bhsm/interface/master_action/mass_curvature_response.py`: oriented metric/scalar canonical momenta and cap sum/mean convention.
- Original BHSM `src/bhsm/interface/master_action/terms.py`, `variations.py`, `coefficients.py`: displayed matcher convention, intrinsic coefficient ownership and levelwise equations.
- Original BHSM `src/bhsm/interface/owner_authorized_encapsulation_interface_action.py`: actual interface Euler-source relation and unselected invariant density class.
- Original BHSM `docs/BHSM_GLOBAL_ACTION_X2_SELECTION_V1.md`: reconciliation with later retained work and the distinction between generic reduction constructions and an evaluated microscopic R1 intertwiner.

This is an AI-assisted mathematical/research review and reproducible conditional calculation, not external human peer review.
