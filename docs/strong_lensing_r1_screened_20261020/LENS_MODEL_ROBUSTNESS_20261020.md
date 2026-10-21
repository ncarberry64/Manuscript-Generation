# Does the state-conditioned transfer survive a change of lens profile?

## Controlled question

This bounded study tests one part of the paper's central argument: the light-transfer calculation remains well-defined when the ordinary lens profile changes, while its interpretation still depends on the admitted mass/source freedoms. It adds eight paired optical applications, not a new cosmological model or a fitted BHSM amplitude. All historical calculations remain unchanged.

Each pair uses the same frozen R1 coefficients, incident marker, source position, redshifts, observer clock and three physical distance maps. One member is the conventional response of the declared density; the other includes the action-owned scalar response of that same density. Ordinary input changes are declared before applying the response and are shared within each pair.

The design has two profiles and four settings: the reference, a 20% increase in total mass, a 20% increase in 3D half-mass radius, and a 20% increase in physical source displacement. Only one input changes at a time. These finite changes are controls, not posterior uncertainties or empirical priors.

## A common size convention

The reference remains the existing Tolman--VII density with total mass \(10^{11}M_\odot\) and material radius \(R=3\) kpc. Its half-mass point solves

\[
 \tfrac12(5x_h^3-3x_h^5)=\tfrac12,
 \qquad x_h=0.643138782163354,
 \qquad r_h=1.929416346490062\ {\rm kpc}.
\]

The alternative is the standard smooth Plummer sphere [Plummer (1911)](https://doi.org/10.1093/mnras/71.5.460), used as an independently established finite-mass spherical profile rather than a newly normalized ansatz:

\[
 \rho_P(r)=\frac{3M}{4\pi a^3}(1+r^2/a^2)^{-5/2},\qquad
 M_P(<r)=M\frac{r^3}{(r^2+a^2)^{3/2}}.
\]

Matching the same **total 3D mass and 3D half-mass radius** fixes

\[
 a=r_h\sqrt{2^{2/3}-1}=1.478745083254196\ {\rm kpc}.
\]

The Plummer scale is not a material boundary. Its density has a regular center and finite total mass without a hard cutoff. The study does not claim that either prescribed spherical profile fits an observed galaxy.

Direct projection of this density gives

\[
 \Sigma_P(b)=\frac{Ma^2}{\pi(b^2+a^2)^2},\qquad
 M_{P,2D}(<b)=M\frac{b^2}{b^2+a^2},
\]
\[
 \widehat\alpha_P(b)=\frac{4m_gb}{b^2+a^2},\qquad
 \widehat\psi_P(b)=2m_g\ln\frac{b^2+a^2}{2a^2},\quad m_g=GM/c^2.
\]

The potential zero at \(b=a\) is a common-arrival gauge. Its gradient and Hessian are derived from that same expression. No extra mass is inferred from a desired deflection.

## The same action and conservation law

For either enclosed-mass profile, the original retained scalar application is

\[
 Cr^2\pi_T'-\frac4\lambda r(\pi_T')^2=-b_TGM(<r),
\]

with the same regular root, frozen \(C,b_T,\lambda\), and metric constraint

\[
 \Phi'=\Psi'=\frac{GM(<r)}{r^2}-\frac{b_T}{2M_{\rm Pl}^2}\pi_T'.
\]

Replacing the enclosed mass is the only scalar-source change. The Plummer logarithmic mass slope is \(d\ln M(<r)/d\ln r=3a^2/(r^2+a^2)\), which also fixes the analytic deflection derivative. Both profiles have regular scalar sections at the origin and finite exterior charge. Their equal total reference mass gives the same \(r_V=212.0292249771\) kpc, although their optical transfers differ.

For prescribed density, the leading isotropic support follows the same conservation application,

\[
 p(r)/c^2=\int_r^{r_{\rm out}}\rho(s)\Phi'(s)\,ds,
\]

where \(r_{\rm out}=R\) for the finite Tolman body and infinity for Plummer. The fractional central pressure adjustments are \(3.56959\times10^{-6}\) and \(2.33994\times10^{-6}\), respectively, both positive. This checks leading prescribed-density support; it is not a shared fixed equation of state, an exact modified TOV solution or a collisionless distribution-function stability proof.

## Recomputed images and delays

The unchanged source has transverse displacement \(1.4576633674\) kpc. The reference-profile results are:

| Profile | Reference image directions (arcsec) | Reference delays (days) | R1 delay offsets (s) |
|---|---|---|---|
| Tolman VII | \(+0.835476034,-0.082161940,-0.657495525\) | \(21.716076566,14.271452387\) | \(25.95925815,8.64901938\) |
| Plummer, same \(M,r_h\) | \(+0.779552504,-0.042531153,-0.559040843\) | \(19.781760311,12.830475650\) | \(22.45096789,9.62003445\) |

Delays are central-minus-positive and negative-minus-positive, in that order. All eight applications retain three images. Their action-owned offsets are:

| Shared input control | Tolman VII offsets (s) | Plummer offsets (s) |
|---|---:|---:|
| Reference | \(25.95925815,8.64901938\) | \(22.45096789,9.62003445\) |
| \(M\) increased by 20% | \(32.69985084,9.89071878\) | \(28.67429799,10.78426324\) |
| \(r_h\) increased by 20% | \(24.26973517,9.94580400\) | \(20.27223884,10.19744776\) |
| Physical source offset increased by 20% | \(26.84524984,10.39452738\) | \(23.51104006,11.56337977\) |

At matched total mass and half-mass radius, changing the conventional profile changes the two reference delays by \(-167124.9245\) s and \(-124500.3901\) s, while changing the R1 offsets by \(-3.50829026\) s and \(+0.97101507\) s. The maximum ordinary-profile delay change is about 6438 times the largest R1 offset in these reference cases. The one-input 20% controls change conventional delays by thousands to hundreds of thousands of seconds.

These are deliberate sensitivity comparisons. They do not establish an uncertainty distribution for a real lens. They show why a quoted state-conditioned response and a measurement of that state are different questions: the transfer consumes the complete projected profile, and ordinary model freedoms can dominate its absolute arrival predictions.

Indeed the two reference profiles have exactly the same \(M(<r_h)=M/2\) and the same scalar force correction at \(r_h\). Their projected deflections and delays still differ because the ray samples the entire field. Matching a total mass and one size or force point does not match the light transfer.

## Effective-mass representation remains the same contribution

For each profile,

\[
 \delta M_{\rm eff}(r)=M(<r)F(r),\qquad
 \delta\rho_{\rm eff}(r)=\frac1{4\pi r^2}\frac{d\delta M_{\rm eff}}{dr}\ge0.
\]

The nonnegative-density proof uses \(n_M\ge0\) and is unchanged. The calculation integrates this density independently to verify its enclosed mass, then projects \(\delta M_{\rm eff}\) directly along \(z/b\), independently of the R1 angular-substitution projection and potential. It re-solves both reference image systems using the resulting conventional effective mass.

Across the two profiles, the direct projection agrees with the scalar projection to \(5.27\times10^{-13}\) relative. The matched effective-mass and scalar image directions agree within \(5.56\times10^{-17}\) arcsec and the two delays within \(6.14\times10^{-10}\) s. The three independently integrated effective enclosed-mass checks per profile agree within \(2.78\times10^{-16}\) relative. These figures characterize implementation agreement, not physical accuracy.

The effective mass is used **instead of** the equivalent explicit scalar contribution. It is not added to it. Thus the leading equivalence is not an artifact of the finite body's hard edge or of one selected radial profile.

## What can be inferred

The same source-position, total-mass, homologous-size, compensated mass-sheet and center response columns are evaluated separately on each profile. Both reference cases have one surviving scalar-template direction when the center is fixed and projected rank zero when both center coordinates are free. These are linearized one-source statements in the declared profile family; the numerical whitening scales remain illustrative, not observed precisions.

The larger positive-effective-mass equivalence holds independently of that local rank calculation. It explains why changing the radial profile does not manufacture a unique BHSM discriminator. Prepared light measures the retained envelope response, while physically justified restrictions on the ordinary lens model determine which source properties can be assigned separately.

## Which state coefficients enter the projected response?

This same transfer permits an explicit local coefficient sensitivity. At fixed mass, profile, impact and \(\lambda\), put \(S(r)=\sqrt{1+r_V^3f_M(r)/r^3}\). Differentiating the original force gives

\[
 \frac{\partial\ln\delta\Phi'}{\partial\ln b_T}
 =\frac32+\frac1{2S},\qquad
 \frac{\partial\ln\delta\Phi'}{\partial\ln C}=-\frac1S.
\]

The projected derivatives of the **extra** deflection are the positive line-of-sight weighted averages of these coefficients. At the three conventional reference impacts, in positive/central/negative-image order:

| Profile | \(\partial\ln\delta\widehat\alpha/\partial\ln b_T\) | \(\partial\ln\delta\widehat\alpha/\partial\ln C\) |
|---|---|---|
| Tolman VII | \(1.525461,1.516634,1.522208\) | \(-0.050923,-0.033268,-0.044417\) |
| Plummer | \(1.525605,1.515036,1.521978\) | \(-0.051211,-0.030071,-0.043956\) |

The deep-screening force suppresses sensitivity to \(C\). Its projected derivative remains nonzero because the ray also traverses the outer field. The exact coefficient identity \(2\partial_{\ln b_T}\ln\delta\widehat\alpha+\partial_{\ln C}\ln\delta\widehat\alpha=3\) is satisfied within \(1.6\times10^{-12}\). Finite coefficient variations independently check the analytic weighted derivatives.

These are coefficient partials, not independent newly selected physical states. An actual envelope-state derivative must compose them with its derived \(b_T,C\) Jacobian. The receipt also records the derivatives divided by the **total** deflection, which carry the additional small factor \(\delta\widehat\alpha/\widehat\alpha_{\rm total}\). No scalar amplitude, background or image solution is refitted by this supplementary derivative application.

## Verification and reproducibility

Across all eight recomputed cases, the maximum original scalar-flux relative residual is \(7.25\times10^{-16}\), and the maximum image-equation residual is \(2.55\times10^{-21}\) rad. Seven focused tests independently check density normalization and half-mass conventions, Plummer line projection, potential/gradient/Hessian consistency, intrinsic size derivatives, the scalar source flux, effective-mass reconstruction, the stated nuisance scope and projected coefficient sensitivities. The six initial tests and the additional coefficient-partial test were each run once; the latter passed after the supplemental derivative application was added.

```text
python scripts/run_strong_lensing_robustness.py
python -m pytest -q tests/test_strong_lensing_robustness.py
```

The result is `artifacts/strong_lensing_acceptance_20261020/robustness.json`. It preserves the exact frozen inputs and hashes. Existing historical outputs are not overwritten. No background shooting, new scalar calibration, common GR2 correction or environmental field solve is repeated. The approximations remain weak metric, frozen quasi-static scalar response, geometric optics and a thin optical plane. This study is not physical admission of a microscopic BHSM source state.
