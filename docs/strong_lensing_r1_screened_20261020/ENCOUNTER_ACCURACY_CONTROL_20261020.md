# Finite-width control of the prepared R1 light transfer

This study tests the **optical application already used in the paper**. It distributes the same finite-body potential along the three rays instead of collapsing it into a thin plane. The source, observer clock convention, measured-redshift background maps, body and frozen scalar coefficients are identical. Each exterior wavefront is derived separately from the closed-FLRW boundary data. Nothing is fitted to obtain a different signal.

The result is a conditional numerical control of the retained weak-field, paraxial action. At a common 20 Mpc half-window, distributing the encounter changes the two scalar-minus-fixed-body delay offsets by approximately **−0.619 ms and +0.0010 ms**. This is larger than the previously computed common-GR2 drift in the first pair. It leaves the seconds-level leading result intact but gives a concrete reason not to interpret all its printed digits as physical accuracy. It is not a complete higher-order scalar metric, an exact finite-distance null-ray solution, or a global BHSM match.

## 1. The physical exterior split

Let \(S_K(\chi)=\sin(\sqrt K\chi)/\sqrt K\), \(C_K=S_K'\), and use the fixed \(a_l=(1+z_l)^{-1}\), \(a_s=(1+z_s)^{-1}\). The saved physical maps give

\[
 D_{ol}=a_lS_K(\chi_l),\quad D_{os}=a_sS_K(\chi_s),\quad
 D_{ls}=a_sS_K(\chi_s-\chi_l).
\]

The declared redshifts are in the unambiguous short branch of the sine. Their inverse, with the frozen curvature, determines \(\chi_l,\chi_s\). The third map is checked independently rather than inferred silently.

Use lens-scaled conformal screen and longitudinal coordinates

\[
 x=a_l q,\quad z=a_l(\chi-\chi_l),\quad
 \kappa=K/a_l^2,\quad
 \ell_o=a_l\chi_l,\quad \ell_s=a_l(\chi_s-\chi_l),\quad
 X_s=(a_l/a_s)s.
\]

The background paraxial action between endpoints a distance \(d\) apart is

\[
 S_\kappa(x_1,x_2;d)=
 \frac{C_\kappa(d)(x_1^2+x_2^2)-2x_1x_2}{2S_\kappa(d)},
 \qquad S_\kappa(d)=\frac{\sin(\sqrt\kappa d)}{\sqrt\kappa}.
\]

Consequently, at the central plane the **individual** wavefront curvatures and source term are

\[
 k_o=\frac{C_K(\chi_l)}{D_{ol}},\quad
 k_s=\frac{C_K(\chi_s-\chi_l)}{a_lS_K(\chi_s-\chi_l)},\quad
 j_s=\frac{X_s}{S_\kappa(\ell_s)}=\frac{s}{D_{ls}}.
\]

The closed-space addition identity then gives
\(k_o+k_s=D_{os}/(D_{ol}D_{ls})\). This fixes each side from its own physical boundary, not merely a convenient split of the required sum.

The actual values are

| Boundary quantity | Value |
|---|---:|
| Observer wavefront curvature \(k_o\) | 0.0010572918625775208 Mpc\(^{-1}\) |
| Source wavefront curvature \(k_s\) | 0.0006022563206977872 Mpc\(^{-1}\) |
| Source term \(j_s\) | \(1.3534549833466666\times10^{-6}\) |
| Observer endpoint length \(\ell_o\) | 945.3793478167898 Mpc |
| Source endpoint length \(\ell_s\) | 1658.0823008477794 Mpc |
| Lens-scaled fixed source \(X_s\) | 0.0022425590267960926 Mpc |
| Local background curvature \(\kappa\) | \(1.5375795596997213\times10^{-9}\) Mpc\(^{-2}\) |

The sum agrees with the saved three-map expression within floating-point rounding. In particular \(k_o\ne1/D_{ol}\); replacing the physical split with flat effective distances would discard known background data.

## 2. The distributed action actually solved

For a common half-window \(L\), evaluate

\[
 ct=\frac1{a_l}\left\{
 S_\kappa(0,x_-;\ell_o-L)
 +\int_{-L}^{L}\left[\frac12(x'^2-\kappa x^2)
                       -2\Phi(\sqrt{x^2+z^2})\right]dz
 +S_\kappa(x_+,X_s;\ell_s-L)\right\}+\text{common constant}.
\]

The factor \(1/a_l=1+z_l\) is the same observer-clock normalization as the manuscript. \(X_s\), rather than the angular source position, is fixed. Potential constants cancel in relative arrivals because the potential is integrated over the same interval for every image.

Variation gives the original field evaluated along the computed section:

\[
 x''+\kappa x=-2\Phi'(r)x/r,\qquad r=\sqrt{x^2+z^2},
\]

\[
 x'(-L)=k_-x_-,\qquad x'(L)=j_+-k_+x_+,
\]

where
\(k_-=C_\kappa(\ell_o-L)/S_\kappa(\ell_o-L)\),
\(k_+=C_\kappa(\ell_s-L)/S_\kappa(\ell_s-L)\), and
\(j_+=X_s/S_\kappa(\ell_s-L)\).
The observer direction is
\(\theta=x_-/S_\kappa(\ell_o-L)\), not an assumed extrapolation of \(x(0)\).

The fixed-body force is \(\Phi_N'=m_g f_M(r)/r^2\). The R1 application uses exactly the saved force

\[
 \Phi_{\rm R1}'=\Phi_N'[1+F(r)],\qquad
 F(r)=\frac{2\epsilon_\infty}
 {1+\sqrt{1+r_V^3 f_M(r)/r^3}}.
\]

No GR2 potential is included in this study. The body and scalar equations are not re-solved at another field order. The background curvature term is counted in the free propagation and the finite interval once; the local force is the same excess source as before.

The comparison thin plane uses **the same window**:

\[
 \hat\alpha_L(b)=4b\int_0^L\frac{\Phi'(\sqrt{b^2+z^2})}
                                      {\sqrt{b^2+z^2}}dz,
\]
\[
 \hat\psi_L(b)-\hat\psi_L(b_0)
 =4\int_0^L[\Phi(\sqrt{b^2+z^2})-\Phi(\sqrt{b_0^2+z^2})]dz.
\]

Its stationary equation is \((k_o+k_s)x-j_s-\hat\alpha_L(x)=0\). Collapsing this same-window projection to the central plane gives the existing plane action with its window-projected potential. The same three noncritical branches are initialized from the saved images and solved again in each application.

The ray solve is nonlinear in its trajectory. It remains a **leading weak-field paraxial** action: terms such as \(\Phi x'^2\), quartic angular terms, a scalar-dependent higher-order metric, and time evolution have not been supplied by solving the full field equations. Its finite-width correction is therefore a specific optical control, not a complete post-Minkowskian prediction.

## 3. Actual three-ray results

All values below are scalar-minus-fixed-body differences in the same finite window. They do not replace the archived infinite-plane prediction.

| Half-window \(L\) (Mpc) | Thin-window delay offsets (s) | Distributed delay offsets (s) | Width/path change (microseconds) |
|---:|---:|---:|---:|
| 0.212029225 | 25.668832976; 8.537235483 | 25.668372786; 8.537237038 | −460.190; +1.555 |
| 0.848116900 | 25.939611607; 8.641457491 | 25.939033033; 8.641458638 | −578.575; +1.147 |
| 4.240584500 | 25.958471067; 8.648716431 | 25.957858682; 8.648717463 | −612.385; +1.032 |
| 20.000000000 | 25.959222770; 8.649005757 | 25.958603720; 8.649006762 | −619.050; +1.006 |

At 20 Mpc the relative changes are approximately \(-23.85\) ppm and \(+0.116\) ppm of the same-window state response. The ordinary fixed-body delays themselves change by \(-4.583872\) s and \(+1.361295\) s. These sizeable common transfer corrections and the much smaller change in the state-dependent offset are different quantities.

The final reference-relative image offsets are approximately
\((6.1232664,0.7880596,-4.3485861)\) microarcseconds. Their width/path changes are approximately \((+0.0000004,-0.0000467,-0.0000009)\) microarcseconds; those very small changes have less numerical relative precision than the delay response and are retained in the receipt rather than advertised as absolute astrometric accuracy.

A direct stationary first variation provides a stronger calculation of the state response without subtracting ordinary arrival times:

\[
 \delta\Delta t_{i0}^{\rm ray}
 =-\frac{2}{a_lc}\int_{-L}^L
   [\delta\Phi(r_i(z))-\delta\Phi(r_0(z))]dz.
\]

This evaluates the already derived extra potential on the computed fixed-body sections. It does not introduce a fitted scalar amplitude. At 20 Mpc its distributed-minus-plane response is
\((-619.0342,+1.0089)\) microseconds. To calculate the finite change without subtracting large ordinary arrivals, a numerical homotopy integrates this stationary derivative along \(\Phi_\eta=\Phi_N+\eta\delta\Phi\), from \(\eta=0\) to 1. This is an integration coordinate between the two fixed optical fields, not a fitted state parameter or a claim that intermediate fields solve another scalar theory. Gauss orders 3 and 5 give the finite thickness change \((-619.0435,+1.0064)\) microseconds and agree at the displayed precision. The direct finite-arrival subtraction agrees within 7 ns. The direct stationary and homotopy applications are the preferred numerical diagnostics; they avoid attributing precision to cancellation of the much larger ordinary arrival terms.

## 4. A rigorous, narrowly scoped exterior bound

For the declared isolated extra field, \(0\le f_M\le1\) and \(0\le F\le\epsilon_\infty\). Therefore its fixed-impact Born deflection tail satisfies

\[
 |\delta\hat\alpha_{z>L}(b)|\le
 \frac{4m_g\epsilon_\infty}{b}
 \left(1-\frac{L}{\sqrt{L^2+b^2}}\right)
 \le\frac{2m_g\epsilon_\infty b}{L^2}.
\]

Integrating this bound between two impacts gives a relative-potential and stationary-delay bound:

\[
 |\delta\hat\psi_{z>L}(b_i)-\delta\hat\psi_{z>L}(b_0)|
 \le4m_g\epsilon_\infty
 \left|\ln\frac{L+\sqrt{L^2+b_i^2}}{L+\sqrt{L^2+b_0^2}}\right|
 \le\frac{m_g\epsilon_\infty|b_i^2-b_0^2|}{L^2}.
\]

Multiply by \((1+z_l)/c\) for observer seconds. At the 20 Mpc fixed-body impacts these bounds are **35.04 and 13.47 microseconds**. Independent quadrature verifies the inequality. They bound the *isolated scalar fixed-impact Born tail*, not an arbitrary environmental state, the full nonlinear re-solved difference, or global matching. Common-body tail changes and their shifted image positions are evaluated numerically, not silently included in that analytic bound. The modest change from \(20r_V\) to 20 Mpc also remains visible in the table; no infinite-domain thickness value is asserted.

## 5. Numerical checks and approximation hierarchy

The ray equations were solved after scaling length by the material radius and slope by the fixed compactness. No physical parameter changes under that numerical scaling. The displayed solve uses collocation tolerance \(2\times10^{-11}\), with independent integrated momentum and boundary applications. Across the recorded window solves, the integrated momentum residual is below \(3.5\times10^{-16}\) rad and the scaled boundary residual is below \(10^{-15}\). A refinement to \(8\times10^{-12}\) and radial potential quadrature order 64 instead of 32 changes the direct stationary thickness response by less than **0.65 ns**, and the finite re-solved thickness change by less than **0.24 ns**. The direct-response quadrature estimates are below 0.06 ns per displayed delay response. These convergence and quadrature diagnostics concern this optical action, not omitted field physics.

Seven focused checks verify the individual physical wavefront split, zero-field propagation with the same endpoints, a separately parameterized radial potential integral, independent line-of-sight deflection, the rigorous extra-field tail, the computed section's boundary/force balance, and reproduction of the saved force ratios. They all pass. The original source artifacts remain unchanged.

The hierarchy relevant to the paper is now concrete:

- The state transfer is approximately 25.96 s and 8.65 s relative to the fixed-body reference.
- The same-window finite-width/path application changes it by approximately −0.62 ms and +0.0010 ms at 20 Mpc.
- The earlier common-GR2 **separate control** changed it by +0.0690 ms and −0.0101 ms. These controls are not added as a complete higher-order prediction.
- The isolated fixed-impact scalar tail beyond that window has the stated 35 and 13 microsecond bounds.
- Scalar-dependent higher-order stress/slip, material/EOS readjustment, background time dependence, environmental coupling and outer matching remain governed by the previously derived continuation equations. This numerical optical control supplies no bound for all of them.

Thus the new calculation tests the transfer central to the paper and resolves one explicit approximation. It supports a seconds-level conditional response while replacing an unquantified thin-plane concern with an actual same-preparation calculation. It does not promote numerical precision or a selected window to a physical admission certificate.

## Reproduction

From the repository root:

```powershell
python scripts/run_strong_lensing_encounter_accuracy.py
python -m pytest tests/test_strong_lensing_encounter_accuracy.py -q
```

The script writes only `artifacts/strong_lensing_acceptance_20261020/encounter_accuracy.json` by default; `--output` writes another receipt. Inputs are the immutable saved R1 results and frozen prediction manifest, with SHA256 identities recorded. The code is `code/strong_lensing_encounter_accuracy.py`. No broad test suite or background integration is repeated by these commands.
