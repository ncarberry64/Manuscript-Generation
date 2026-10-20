# Conditional observed-lens model and withheld-light test

This note specifies an executable restricted comparison for HE 0435-1223. It preserves the already derived finite-body R1 scalar source, rather than adding that spherical scalar to an arbitrary elliptical matter potential. A rear quadratic optical plane supplies the asymmetric propagation map. Its optical response is conditional boundary information, not a derived global R1 environment. A poor fit of this restricted model is an empirical result; it must not be hidden by broadening astrometric errors after inspecting the held-out data.

The observational redshifts are versioned inputs: use the lens redshift `0.4546` and the source redshift `1.693` adopted by Bonvin et al. (2017), citing the latter measurement to Sluse et al. (2012). Some earlier analyses use `1.689`; that value must not be silently combined with maps labeled `1.693`. The rear-plane redshift is fixed to `0.8` by the frozen analysis protocol and is not an identification of an observed perturber.

## 1. Local body, scalar source and optical field

The local unknown is the spherical scalar section `pi(r)` for a declared nonnegative finite-mass density profile `rho(r;M,a)`. The conventional finite-body potential is `Phi_N`. At the measured lens redshift, evaluate the already adopted rolling reference to obtain `v`, `d`, `H`, and the original coupling `lambda`; do not fit a separate BHSM optical amplitude. In the notation of the existing local derivation,

\[
b_T=2v^2/\lambda,\qquad
C=1+4(d+2Hv)/\lambda-b_T^2/(2M_{\rm Pl}^2).
\]

The equations actually applied are

\[
C\Delta\pi-\frac2\lambda[(\Delta\pi)^2-\Pi:\Pi]
=-\frac{b_T\rho}{2M_{\rm Pl}^2},\qquad
2M_{\rm Pl}^2\Delta\Phi=\rho-b_T\Delta\pi,\quad\Psi=\Phi.
\]

Regularity, the weak-source branch and the isolated quasi-static matching prescription give

\[
Cr^2\pi'-\frac4\lambda r\pi'^2=-b_TGM(<r),\qquad
\pi'=-\frac{2b_TGM(<r)}{r^2[C+\sqrt{C^2+16b_TGM(<r)/(\lambda r^3)}]},
\]

\[
\Phi'=\frac{GM(<r)}{r^2}-\frac{b_T\pi'}{2M_{\rm Pl}^2}.
\]

Thus the same original `ProfileR1Lens` implementation can be evaluated at observed-lens parameters. The density is prescribed at this retained order; it is not a solved galaxy distribution function or a material equation of state. The potential, deflection and Hessian used in the image and delay likelihood must all come from the same projection. The conventional model uses `Phi_N`; the conditional R1 model uses the full `Phi` above.

For example a Plummer body has

\[
\rho=\frac{3Ma^2}{4\pi(r^2+a^2)^{5/2}},\quad
M(<r)=\frac{Mr^3}{(r^2+a^2)^{3/2}},\quad
\widehat\psi_N(b)=2m_g\ln\frac{b^2+a^2}{2a^2},\quad
\widehat\alpha_N(b)=\frac{4m_gb}{b^2+a^2}.
\]

The original quasi-static, weak-metric and thin-encounter restrictions remain. The unchanged `2e-10` numerical physical guard is not converted into an observational goodness-of-fit threshold or proof of physical admission.

## 2. An asymmetric propagation map with a separate optical plane

Let `l < p < s` denote the main body, a rear optical plane and the source. The scalar FLRW angular-diameter distances are `D_l,D_p,D_s,D_lp,D_ls,D_ps`. Use a symmetric physical rear-plane Hessian `U`, with units `Mpc^-1`, so that

\[
\widehat\psi_p(x_p)=\tfrac12x_p^TUx_p,\qquad
\widehat\alpha_p(x_p)=Ux_p.
\]

For the symmetric screen a convenient three-parameter normalization is

\[
U=\frac{D_s}{D_{ps}D_p}\Gamma,\qquad
\Gamma=\begin{pmatrix}\kappa+\gamma_1&\gamma_2\\\gamma_2&\kappa-\gamma_1\end{pmatrix}.
\]

No additional same-plane sheet or shear is then inserted in the main-body potential. In physical coordinates the two-plane ray maps are

\[
x_l=D_l\theta,\qquad
x_p=D_p\theta-D_{lp}\widehat\alpha_l(x_l),
\]

\[
x_s=D_s\theta-D_{ls}\widehat\alpha_l(x_l)-D_{ps}Ux_p.
\]

Eliminating `x_p` gives exactly

\[
x_s=\mathsf D_{os}\theta-\mathsf D_{ls}\widehat\alpha_l(D_l\theta),
\]

\[
\mathsf D_{ol}=D_lI,\quad
\mathsf D_{os}=D_sI-D_{ps}D_pU,\quad
\mathsf D_{ls}=D_{ls}I-D_{ps}D_{lp}U.                 \tag{1}
\]

These are the three maps to supply to the existing `LensGeometry`. They are compatible maps from a single explicit propagation problem. They must not be fitted independently. In terms of `q=D_ls/D_s` and `beta_lp=D_sD_lp/(D_pD_ls)`, they obey

\[
\mathsf D_{os}=D_s(I-\Gamma),\quad
\mathsf D_{ls}=D_{ls}(I-\beta_{lp}\Gamma),\quad
\mathsf G=q(I-\Gamma)^{-1}(I-\beta_{lp}\Gamma).
\]

Restrict the fit to the pre-conjugate branch of both maps. Requiring both eigenvalues of `I-Gamma` and `I-beta_lp Gamma` to be positive is an explicit, stronger sufficient condition than positive determinants alone. The chosen convergence/shear bounds are declared model priors, not results learned from the withheld light. Negative fitted screen convergence describes a tidal perturbation of the declared smooth reference; it is not assigned an independent nonnegative isolated surface-density interpretation.

## 3. Observer-time normalization follows from the two-plane action

Write `theta_p=x_p/D_p` and the true angular source position `beta=x_s/D_s`. At the retained paraxial order, the two-plane arrival length is

\[
ct=(1+z_l)\left[\frac{D_lD_p}{2D_{lp}}|\theta-\theta_p|^2-
\widehat\psi_l(D_l\theta)\right]
+(1+z_p)\left[\frac{D_pD_s}{2D_{ps}}|\theta_p-\beta|^2-
\tfrac12D_p^2\theta_p^TU\theta_p\right].             \tag{2}
\]

Define

\[
k_l=(1+z_l)D_lD_p/D_{lp},\quad
k_p=(1+z_p)D_pD_s/D_{ps},\quad
W=(1+z_p)D_p^2U,
\]

\[
H=(k_l+k_p)I-W,\qquad
\theta_p=H^{-1}(k_l\theta+k_p\beta).
\]

Substituting this stationary rear-plane section into (2) produces a Schur complement with quadratic matrix

\[
K=k_lI-k_l^2H^{-1},\qquad L=k_lk_pH^{-1}.
\]

The FLRW distance identity

\[
(1+z_l)D_lD_{ps}+(1+z_p)D_sD_{lp}
=(1+z_p)D_pD_{ls}
\]

then yields

\[
K=(1+z_l)D_l\mathsf D_{os}\mathsf D_{ls}^{-1}=\mathsf T,
\quad L=(1+z_l)D_lD_s\mathsf D_{ls}^{-1}.
\]

Therefore, up to a source-only constant that cancels between images,

\[
ct=\tfrac12(\theta-\beta_{\rm eff})^T\mathsf T
(\theta-\beta_{\rm eff})-(1+z_l)\widehat\psi_l(D_l\theta),
\quad \beta_{\rm eff}=\mathsf D_{os}^{-1}x_s.       \tag{3}
\]

All the matrices are commuting functions of the same symmetric `U`, so `T` is symmetric. Equation (3) is the existing `lens_application` arrival action with the maps (1). There is one observer clock, no independently fitted delay scale, and no duplicate tidal encounter. The two-plane action (2) supplies an independent check of image delay differences obtained from (3).

Four new exact SymPy checks were executed for this derivation: the curved-FLRW sine distance sum rule, the quadratic Schur coefficient, the source Schur coefficient, and their gradient/lens-equation identity. All simplified to exact zero. These are algebra checks of this new optical map, distinct from the twelve historical R1 field identities; observed-data residuals still have to be calculated by the actual fit.

## 4. What this conditional domain prescription does and does not fix

The main body's scalar Cauchy/matching data are the spherical branch above. The rear-plane `U` specifies total optical Weyl response outside that local body domain. It does not by itself specify remote three-dimensional density or its scalar section. In particular a two-dimensional measured tide does not imply that the environmental scalar Hessian vanishes in the body domain.

If an actual same-domain environment is to be included, the original equations instead require

\[
C\Delta\pi_L-\frac2\lambda[(\Delta\pi_L)^2-\Pi_L:\Pi_L]
-\frac4\lambda[(\Delta\pi_E)(\Delta\pi_L)-\Pi_E:\Pi_L]
=-\frac{b_T\rho_L}{2M_{\rm Pl}^2}.
\]

The consumed contraction is `[(tr Pi_E)delta_ij-Pi_E,ij] Pi_L,ij`, together with the environmental field's own scalar equation and Cauchy/matching data. The projected `U` does not determine that contraction. Choosing a remote plane and an isolated spherical body boundary is therefore a restricted optical conditional model, not a proof that the global environmental field of HE 0435-1223 has been solved. This is scientifically usable provided that restriction is explicit and a failed astrometric fit is retained as a failure of the model.

## 5. The conventional equivalence is a prediction, not an extra component

For the same body and state, define the effective added enclosed mass

\[
\delta M_{\rm eff}(r)=M(<r)\epsilon_F(r),\qquad
\epsilon_F=\frac{b_T^2}{M_{\rm Pl}^2[C+\sqrt{C^2+16b_TGM(<r)/(\lambda r^3)}]}.
\]

On the retained regular nonnegative-density branch the existing derivation gives a nonnegative added density. An ordinary model using `rho + delta rho_eff` has exactly the same `Phi=Psi`, projected potential, deflection and Hessian. With the same rear-plane map it consequently has identical image branches, magnifications and observer delays. This is true before any withheld flux is inspected.

Thus there are two distinct useful comparisons: (i) the restricted ordinary-body model versus the R1-body model, both fitted using the same training data and nuisance prescription; and (ii) R1 versus its fully matched ordinary effective-mass representation as a structural equality control. The latter equality prevents interpreting a better fit than the restricted ordinary-body family as unique evidence for BHSM. Different parameter priors can change Bayesian evidence even for equivalent likelihoods; equality of parameterized predictions is not an assertion about arbitrary priors.

## 6. A genuine withheld-light prediction

Use the historical 13-season HE 0435 data as training and the newly published seasons as withheld observations. Archive the exact row-level split before fitting. Published delays from an analysis that already includes those later seasons must not be used in training. Astrometry can be an external training constraint with its original covariance, not a published cosmological posterior. Train the conventional and conditional R1 macro parameters separately; hold fixed the background/calibration, redshifts, selected model family and rear-plane location. Report astrometric residuals and training residuals as well as held-out scores.

The common observation model may be written

\[
m_{ik}=c_i+s(t_{ik}-\tau_i)+\mu_i(t_{ik})+\epsilon_{ik},
\]

where `s` is the intrinsic quasar variability, `mu_i` is the slowly varying image-specific microlensing contribution, and the supplied measurement uncertainty enters `epsilon`. The macro delays `tau_i` come from (3), rather than being quietly refitted against the withheld images. Any empirical GP/spline law for `s` and `mu_i` is a stated observation model shared by the two hypotheses. It is not a derivation of a BHSM statistical-state kernel, a Wick replacement, or a prediction of intrinsic quasar fluctuations by the envelope equations.

The historical image offsets and extrinsic variability are held fixed for the comparison. Their empirical freedom means this test concerns temporal transfer conditional on photometric calibration and microlensing, rather than an independent prediction of macro magnification ratios. Point-image astrometry is also distinct from fitting the resolved Einstein-ring pixels. Neither unperformed test is credited to the executed likelihood.

A useful sharp transfer test freezes all parameters learned from historical data, then conditions on the *new-season image B observations only* as a measured driver and predicts new-season images A, C and D. These A/C/D observations remain withheld until scoring, as fixed by `preregistration/he0435_heldout_20261020.json`. This is a cross-image conditional prediction; it is not an all-channel forecast of future source variability. An additional all-channel chronological forecast is reported separately, with its broad intrinsic-source extrapolation uncertainty. Delayed-source input and microlensing are marginalized with the frozen training prescription; no per-image held-out polynomial or flux shift is fitted after looking at the answers.

Report uncertainty-aware held-out likelihood/coverage, standardized residuals and the numerical conventional/R1 delay difference. A zero-delay or incorrect-delay control can establish that the withheld data actually test the light-transfer operator. If both physical models make statistically indistinguishable predictions at the available cadence and precision, that is the result. Do not report a detected or uniquely mapped underlying envelope state.

## 7. Executed source activation and stationary-image audit before scoring

The separate `scripts/audit_he0435_macro.py` application consumes the frozen macro fit files and no held-out magnitudes. At the conventional fitted body, source, screen, geometry and clock, activating the original R1 scalar gives delay changes, ordered B/C/D minus A, of **19.29034420 s, 0.93856058 s and 29.03669445 s**. The original scalar and optical residuals are in `macro_forward_audit.json`.

The separately refitted models differ by 228.50610 s, 7.23433 s and 438.43556 s. That difference also includes nuisance-parameter reallocation and finite optimizer termination. It must not be described as the fixed-body R1 response. The audit performs no further fit.

The finite Plummer profile predicts a central stationary image in addition to the four quasar branches used in training. A search from 1060 Cartesian and polar seeds found five distinct roots in each model and verified every root using the original un-interpolated source projection. The central image in the conventional model lies at `(0.003736,-0.027160)` arcsec relative to the galaxy light center, with absolute macro flux about **0.06320 of image A**. The conditional R1 fit gives 0.06313. These ratios omit microlensing, extinction, finite source and lens-galaxy contamination. No observed central-image non-detection limit was introduced, so this is a model prediction requiring a separate image-level test, rather than a claim that it passes that test.

The screen changes absolute magnification. With `A=D_os_eff/D_s`, the actual source angular Jacobian is `A J_normalized`; hence the audit reports `1/det(A J_normalized)`. The common `det(A)` cancels only in relative image fluxes. The conventional central magnification is about 0.7001 and that of image A about 11.0774 in this convention.

For this positive circular mass family one can exclude images outside a finite radius. If `R=|theta|` exceeds both `|beta_eff|` and `|theta_center|`,

\[
(R-|\beta_{\rm eff}|)(R-|\theta_{\rm center}|)
\leq \frac{4\|\mathsf G\|m_{\rm total}}{D_l},
\]

where `m_total=m_g` conventionally and `m_g(1+epsilon_infinity)` for the derived positive R1 effective mass. This gives outer bounds 1.58183 and 1.58246 arcsec for the two fitted models. The multiseed search samples the interior of that rigorous bound, but is not an analytic proof that every interior basin was found. These results prevent promoting a four-branch timing model to a complete observed-galaxy reconstruction.

## Primary source trail

- Bonvin et al. (2017), *H0LiCOW V. New COSMOGRAIL time delays of HE 0435-1223*: historical 13-year photometry, source redshift version and time-delay method. https://doi.org/10.1093/mnras/stw3006
- Millon et al. (2020), *COSMOGRAIL XIX*: publicly released light curves and two additional HE 0435 monitoring seasons. https://doi.org/10.1051/0004-6361/202037740
- Courbin et al. (2011), *COSMOGRAIL IX*: original HST astrometry and the earlier source-redshift version. https://doi.org/10.1051/0004-6361/201015709
- Sluse et al. (2017), *H0LiCOW II*: independent spectroscopic environment, including `z_s=1.693 +/- .001`. https://doi.org/10.1093/mnras/stx1484
- Fleury, Larena and Uzan (2021), *Gravitational lenses in arbitrary space-times*: the compatible-map and arrival-time framework. https://arxiv.org/abs/2011.04440
- Existing repository derivations: `R1_LOCAL_COUPLING_DERIVATION_20261020.md`, `R1_REMAINING_FIELD_DERIVATION_20261020.md`, and `LENS_MODEL_ROBUSTNESS_20261020.md`. The two-plane Schur reduction above is derived explicitly here.
