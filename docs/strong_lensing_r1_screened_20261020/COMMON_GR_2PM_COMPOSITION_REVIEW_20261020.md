# Independent review: common GR 2PM correction and the screened comparison

> Scope reconciliation: the [current BHSM scope](../STRONG_LENSING_SCOPE_20261020.md)
> supersedes earlier project framing in this retained technical record. The
> governing question is surface contributions to envelope light transfer;
> conventional effective-mass equivalence is a constructive result at the
> retained order. Reference-body offsets are not discrepancies between fully
> matched representations. The full surface-to-R1 common-source derivation
> remains unfinished. Equations, numerical evidence and review chronology
> below keep their original scope.

> This is a separate AI-assisted computational verification, not an external human peer-review report.

## Finding

The supplied finite-body coefficient, its primitive, and its composition with the existing R1 lens are internally consistent. Adding the same ordinary-GR 2PM correction to both models leaves the previously frozen screened excess delays nearly unchanged. This resolves the specific concern that omitting a common GR correction of comparable absolute size might numerically manufacture the reported model difference. It does **not** complete the scalar post-Newtonian field equations or finite-distance cosmological optics.

All calculations below retain the previously fixed body, source, redshifts, background maps, and three image branches. No parameter, source geometry, or image selection was changed after seeing the result.

| Relative delay | R1 minus control with 1PM base (s) | R1 minus control with shared GR2 base (s) | Change of model difference (s) |
|---|---:|---:|---:|
| Central image minus positive outer image | 25.9592581496 | 25.9593271397 | +0.0000689901640 |
| Negative outer image minus positive outer image | 8.64901937523 | 8.64900923806 | -0.0000101371711 |

The fractional changes are +2.657632 and -1.172060 parts per million of the respective R1 differences. The common GR correction itself changes the conventional delays by **8.2847277743 s** and **2.5524454784 s**. Those appreciable common shifts and the small changes of the differences are distinct results.

The corrected R1-minus-control image shifts, in microarcseconds along the image axis, are

\[
(6.12326934335,\;0.788107108546,\;-4.34858643117).
\]

Their changes from the original screened comparison are

\[
(-6.9055872\times10^{-6},\;-5.5760951\times10^{-8},\;5.0933702\times10^{-6})
\quad\text{microarcseconds}.
\]

## 1. The finite-body coefficient and its primitive

Let \(m=GM/c^2\), \(\epsilon=m/R\), \(y=b/R\), and \(t=\sqrt{1-y^2}\) inside the body. Here **\(b\) is the asymptotic invariant scattering impact**, not the turning areal radius. The supplied correction is

\[
\delta\widehat\alpha_2(b)=2\epsilon^2 I(y),\qquad
\delta\widehat\psi_2(b)=2\epsilon^2R J(y),\qquad J'=I,\quad J(1)=0.
\]

For \(0<y<1\),

\[
I(y)=\frac{15[\arcsin y-yt]}{4y^2}
+y\left(\frac{15}{2}t-4t^5-\frac{48}{7}t^7\right),
\]

\[
J(y)=\frac{15\pi}{8}-\frac{15}{4}\left(\frac{\arcsin y}{y}+t\right)
-\frac52t^3+\frac47t^7+\frac{16}{21}t^9.
\]

The primitive is correct analytically: the derivative of \(\arcsin(y)/y+t\) is \(-[\arcsin y-yt]/y^2\); differentiating its remaining powers supplies the rest of \(I\). Independent direct quadrature of \(I\) from 1 to ten representative radii, including the center, surface, and inherited image impacts, agrees with \(J\) within \(1.33\times10^{-15}\) in its dimensionless normalization. A separate finite-difference check of \(I'\) agrees within \(7.32\times10^{-10}\) relative.

The central expansion is

\[
I(y)=-\frac67y+31y^3-\frac{1065}{28}y^5
+\frac{385}{24}y^7-\frac{1285}{704}y^9+O(y^{11}).
\]

Consequently the extra deflection vanishes at the center, its radial and tangential Hessian limits agree and remain finite,

\[
\delta\widehat H_2(0)=-\frac{12\epsilon^2}{7R}\,\mathbf I,
\qquad J(0)=\frac{15\pi}{8}-\frac{26}{3}.
\]

At the surface, \(I(1)=15\pi/8\) and \(I'(1)=-15\pi/4\), matching the exterior expression

\[
\delta\widehat\alpha_2=\frac{15\pi}{4}\frac{m^2}{b^2},\qquad
\delta\widehat\psi_2=\frac{15\pi}{4}m^2\left(\frac1R-\frac1b\right).
\]

Thus potential, deflection, and Hessian join continuously. The coefficient is not infinitely differentiable at the surface: its interior-minus-exterior deflection coefficient starts as \(-10(1-y^2)^{3/2}\), so an unqualified claim of smoothness of all higher derivatives would be incorrect. The implemented lens equation and arrival derivatives require only the finite quantities above.

At the three original impacts the extra GR deflections are

\[
(1.84072751081\times10^{-11},\;-2.41589331705\times10^{-13},\;
2.97215800368\times10^{-11})\ \mathrm{rad},
\]

with radial signs quoted before multiplication by the image direction. A negative central correction is consistent with the independent exact static-metric check conducted by the coupling agent; it must not be replaced with the positive exterior point-mass coefficient.

## 2. Composition with the R1 calculation

The correct implementation is the one already supported by the new class:

```python
base_gr2 = TolmanVII2PMLens(original_TolmanVII_1PM_lens)
screened_gr2 = R1ScreenedLens(base_gr2, frozen_coefficients)
```

The conventional potential is \(\widehat\psi_1+\delta\widehat\psi_2\), and the screened potential is \(\widehat\psi_1+\delta\widehat\psi_2+\delta\widehat\psi_{\rm R1}\). The same additive relations must hold for deflection and Hessian. This composition preserves the original R1 force law and its frozen coefficient values; it does not replace the source mass profile by a mass inferred from a corrected image position.

For compatible distance matrices, the existing arrival functional is

\[
ct(\theta)=\frac12(\theta-\beta)^TT(\theta-\beta)
-(1+z_l)\widehat\psi(D_{ol}\theta),\qquad
TG=(1+z_l)D_{ol}^{T}.
\]

Therefore its stationary equation is exactly the implemented lens equation for **any** consistently differentiated additive encounter potential. The primitive must be included in the arrival calculation; changing deflection alone would violate this identity.

For a fixed-geometry perturbation, stationarity gives

\[
\delta(ct)_i=-(1+z_l)\delta\widehat\psi(D_{ol}\theta_i).
\]

The first-order common-GR changes and the re-solved changes differ by only \((+7.64518,-8.04016)\) microseconds, as expected for the nonlinear stationary-image correction. For R1 on the corrected base, the corresponding finite-minus-linear differences are \((+159.226,+56.0953)\) microseconds. These are finite-amplitude effects inside the defined optical model, not estimates of its physical truncation error.

There is a stronger independent check of the small four-case difference. If \(A\) is the extra GR contribution to the arrival action, \(B\) is the R1 contribution, and \(H\) is the unperturbed arrival Hessian, the mixed derivative of the on-shell action is

\[
\partial_{\eta}\partial_{\lambda}S_{\rm stationary}
=-(\nabla A)^TH^{-1}(\nabla B).
\]

Equivalently, its value for each image is

\[
-(1+z_l)\,\delta\widehat\alpha_{\rm R1}^{T}D_{ol}\delta\theta_{\rm GR2}.
\]

This predicts relative-delay changes of \((+68.9907790,-10.1375691)\) microseconds, agreeing with the four independently re-solved models to less than **0.62 nanoseconds**. The largest lens-equation residual among all four solutions is \(8.47\times10^{-22}\) radians. This confirms the numerical composition and explains why a substantial common correction changes the model difference only weakly.

The new mass and radius intrinsic derivatives also have the correct scaling. At fixed physical impact,

\[
\partial_{\log M}\delta\widehat\psi_2=2\delta\widehat\psi_2,
\quad \partial_{\log M}\delta\widehat\alpha_2=2\delta\widehat\alpha_2,
\]

\[
\partial_{\log R}\delta\widehat\psi_2=-\delta\widehat\psi_2-b\delta\widehat\alpha_2,
\quad
\partial_{\log R}\delta\widehat\alpha_2=-2\delta\widehat\alpha_2-b\delta\widehat\alpha_2'.
\]

Simply forwarding the old 1PM intrinsic derivatives would omit these contributions. The supplied class includes them.

## 3. What is and is not improved physically

The added coefficient is a genuine ordinary-GR 2PM isolated-scattering coefficient for the chosen finite body, including its central ray. Within the declared mapping of invariant impact to physical thin-plane impact, integrating that coefficient supplies a consistent local encounter generating potential. The additive constant fixed by \(J(1)=0\) cancels in relative delays.

However, identifying \(b\) with \(|D_{ol}\theta|\) remains the thin-plane convention. The metric turning point satisfies \(b=r_0e^{-\nu(r_0)}\); substituting \(r_0\) for \(b\) would alter the 2PM coefficient. A full physical calculation would also relate the observed local angle, the actual finite source and observer positions, the physical lens radius, and the matched cosmological clock. Those relations are not derived by integrating the asymptotic deflection alone. The same is true of nonlinear coupling to the smooth background.

More importantly, the corrected pair contains a common GR2 term plus the previously derived **leading local quasi-static** R1 potential. It does not contain a complete set of scalar–GR cross terms, higher-order scalar stress and gravitational slip, the associated hydrostatic readjustment at that order, or time-dependent and cosmological matching corrections. Re-solving rays in the additive potential includes some optical cross effects, such as the mixed stationary derivative above; it does not generate the missing field-equation contributions.

The result can therefore be called a **common ordinary-GR 2PM robustness control**. It cannot be called a complete 2PM R1 prediction, an exact cosmological ray-tracing result, or a quantitative error bar on the original scalar prediction. The very small numerical drift does not imply that every omitted physical correction is comparably small.

## 4. Observational equivalence and nuisance interpretation

The earlier positive-effective-halo construction remains valid at the retained leading quasi-static order:

\[
\delta M_{\rm eff}(r)=M\mu(r)F(r),\qquad
\delta\rho_{\rm eff}=\frac{M}{4\pi r^2}\frac{d[\mu F]}{dr}\ge0.
\]

Because \(\Phi=\Psi\) at that order, it reproduces both leading light deflection and leading nonrelativistic dynamics. The new control does not establish exact all-order equivalence. A conventional body plus the effective halo has its own GR post-Newtonian cross terms; an R1 body has its own scalar–metric cross terms. Neither complete set has been compared here. Conversely, the new control supplies no basis for claiming that the effective-halo degeneracy has been broken.

The existing nuisance-rank results remain explicitly results of their declared optical approximation and finite lens family. A small surviving template after restricting the ordinary density to Tolman-VII mass and radius is not evidence of a theory-specific measurable signal. If nuisance results are re-evaluated with GR2, use the corrected intrinsic derivatives and state that approximation consistently. Do not infer new identifiability solely from the common-GR robustness result.

## Recommended manuscript wording

> We add the finite-body ordinary-GR second post-Minkowskian deflection and its consistent encounter potential to both the control and scalar-on lenses, retaining the same body, source, background maps, and thin-plane impact convention. The common correction changes the two reference delays by 8.285 s and 2.552 s, while changing the scalar-minus-control differences by only +69.0 microseconds and -10.1 microseconds. This establishes numerical robustness of the conditional difference to the shared GR correction. It does not complete the scalar post-Newtonian or finite-distance cosmological calculation, and it does not remove the leading effective-mass-profile degeneracy.

## Review artifacts

- `check_common_gr_composition.py`: one fixed four-case composition and primitive check, plus its analytic mixed stationary derivative.
- `common_gr_composition_review_results.json`: numerical results and scope qualifications.
- The supplied and independently reviewed GR coefficient implementation is `code/strong_lensing_gr2.py`; the separate exact-static check is `scripts/r1_screened_review/check_gr_optical_control.py`.

No manuscript, main prediction module, or published repository was edited by this review.
