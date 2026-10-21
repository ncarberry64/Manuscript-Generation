# Referee review: local R1 scalar prediction

> Scope reconciliation: the [current BHSM scope](../STRONG_LENSING_SCOPE_20261020.md)
> supersedes earlier project framing in this retained technical record. The
> governing question is surface contributions to envelope light transfer;
> conventional effective-mass equivalence is a constructive result at the
> retained order. Reference-body offsets are not discrepancies between fully
> matched representations. The full surface-to-R1 common-source derivation
> remains unfinished. Equations, numerical evidence and review chronology
> below keep their original scope.

> **Review disposition.** This is an AI-assisted technical review, not an external human journal referee report. It records the accuracy issue identified before the common GR2 calculation. That requested common conventional correction has now been implemented and checked, including the central ray; see `FINITE_BODY_GR_2PM_DERIVATION_20261020.md` and `gr2_control.json`. Statements below that this common correction has not yet been calculated describe the earlier review state. Missing scalar higher-order terms, global matching and other physical scope limitations remain as stated.

Date: 2026-10-20. Reviewed the completed local force/projection module, replay script, numerical result, and common-action derivation. This review does not modify the project. Its one additional numerical check is the analytic ordinary-GR second-post-Minkowskian scale at the two exterior images; no lens solver or background fit was rerun.

## Recommendation

The local calculation can be retained as a reproducible **matched leading-order benchmark of the supplied cubic R1 effective action**. I find no factor, sign, or unit correction required in the implemented source coefficient, constrained stiffness, regular scalar root, Vainshtein radius, Weyl-force projection, or leading hydrostatic correction. The major corrections are to accuracy and novelty claims. In particular, numerical convergence must not be equated to absolute microarcsecond physical accuracy.

## 1. Must quantify the omitted ordinary-GR scale

`code/strong_lensing_local_family.py:159–202` supplies the leading weak-field, thin-plane lens potential. Its separate exact GR fluid solution does not make this optical potential exact. `code/strong_lensing_r1_screened.py` adds the leading quasi-static scalar correction to that same leading optical model.

For either ray remaining outside the finite body, the asymptotically flat Schwarzschild deflection in terms of invariant impact parameter is

\[
\hat\alpha(b)=\frac{4m}{b}+\frac{15\pi}{4}\left(\frac{m}{b}\right)^2
+O\big((m/b)^3\big),\qquad m=GM/c^2.
\]

Thus the first omitted ordinary-GR term relative to the retained deflection is

\[
\epsilon_{2\rm PM}=\frac{15\pi}{16}\frac{m}{b}.
\]

This follows, for example, from Iyer and Petters, *Light's Bending Angle due to Black Holes: From the Photon Sphere to Infinity*, arXiv:gr-qc/0611086v2, Eq. (23), using their \(1-b'=3\sqrt3m/b\): <https://arxiv.org/pdf/gr-qc/0611086>.

Using the unchanged mass and stored image impacts gives:

| Image | Impact (kpc) | Omitted 2PM / retained GR deflection | Illustrative 2PM image displacement | Reported R1 displacement |
|---|---:|---:|---:|---:|
| Positive outer | 3.82838408 | 3.68150 ppm | +1.35457 microarcsec | +6.12328 microarcsec |
| Negative outer | 3.01282778 | 4.67807 ppm | −1.72124 microarcsec | −4.34859 microarcsec |

For this scale estimate only, the image response is evaluated through the existing leading radial Jacobian:

\[
\delta\theta_{2\rm PM}
=\operatorname{sgn}(\theta)
\frac{(D_{ls}/D_s)(15\pi/4)(m/b)^2}
{1+(D_{ls}/D_s)D_l\,4m/b^2}.
\]

The two rays have \(b>R\); this estimate was not applied to the central ray traversing the fluid. It is also not a complete finite-distance, closed-cosmology second-order lens calculation. Its purpose is to demonstrate the size of an omitted known term. The values use `results.json` keys `body.mass_Msun`, `image_impact_radii_Mpc[0,2]`, `physical_maps_Mpc`, and `conventional.images_rad[0,2]`.

**Interpretation:** this does not by itself invalidate the *difference* between consistently matched leading-order models. Common GR 2PM terms can cancel in that difference to leading order; a scalar-dependent mixed-order remainder is expected to be smaller. The implementation has not established that cancellation through a full higher-order calculation. Therefore the 2PM values must not be added to the scalar shifts as if they were additional scalar effects, or presented as an error bar on the differential prediction.

**Suggested main-text replacement:**

> The quoted shifts are differences between matched leading weak-field, thin-lens models. They do not imply absolute astrometric accuracy at the displayed microarcsecond scale. For the two exterior rays, the first omitted ordinary Schwarzschild term is 3.68 and 4.68 ppm of the leading deflection, corresponding illustratively to 1.35 and 1.72 microarcseconds through the reference image Jacobian. Common higher-order GR terms may cancel in the matched difference, but that cancellation and the remaining scalar-dependent corrections have not been calculated here.

Round headline values to about 6.12, 0.788, and −4.35 microarcsec; 26.0 and 8.65 seconds; and 6.82 ppm. Retain full precision in the reproducibility data.

## 2. The matching-tail check is not a global error bound

`R1ScreenedLens._integral` extends a frozen local spherical profile to infinity. Its `checks` method removes the radial contribution beyond 20 Mpc and finds fractions of order \(10^{-6}\) **of the scalar deflection**, not of the total deflection. This correctly measures an integration tail within that profile. It does not bound a physically different compensating environment, external scalar Hessian, time-dependent solution, or global closed-spatial matching prescription.

The local solution is admissible as an overlap-domain construction with \(r_V\ll r\ll H^{-1},R_{\rm curvature}\). An isolated positive excess on compact \(S^3\) is not by itself a completed global perturbation; a global source/background or compensation prescription remains necessary. The stored \(H_l r_V/c\simeq5.64\times10^{-5}\) and curvature scale support the local approximation but are not themselves certified observable error bounds.

**Suggested wording:**

> Extending the frozen isolated profile beyond 20 Mpc changes its calculated scalar deflection by less than two parts in a million. This is a tail-convergence test for the specified profile, not a bound on cosmological matching or environmental uncertainty.

For a future nonzero environmental scalar perturbation, nonlinear superposition must also be checked. If \(\pi=\pi_{\rm loc}+\pi_{\rm env}\), the retained cubic equation includes

\[
-\frac4\lambda\left[
(\Delta\pi_{\rm loc})(\Delta\pi_{\rm env})
-\partial_i\partial_j\pi_{\rm loc}\,
\partial_i\partial_j\pi_{\rm env}\right].
\]

This does not alter the present isolated-body calculation, but it forbids claiming that an arbitrary scalar environment can be added only through an independent optical distance map with no change to the local force.

## 3. Density, pressure, and equilibrium claims

The replay script computes

\[
\delta p(r)=\int_r^R\rho(s)\frac{GM(<s)}{s^2}\epsilon_F(s)\,ds.
\]

Its units and Tolman-VII mass normalization are correct. It is a leading hydrostatic correction at prescribed density. Fixing that density while changing the required pressure also changes the inferred barotropic relation; it is not a comparison of two equilibria of an independently fixed common material equation of state. The exact GR TOV family must remain separately labeled.

**Suggested wording:**

> We hold the finite density profile fixed and calculate the pressure change needed by leading hydrostatic balance. We do not impose a common material equation of state across the two descriptions or construct an exact modified-gravity TOV family.

As for the optical shift, the computed pressure *difference* should not be mistaken for the absolute accuracy of the Newtonian pressure: ordinary relativistic pressure corrections have not been included in that leading-pressure table.

## 4. Principal positivity is a limited result, and the superluminality belongs in the text

The algebra establishing \(K_r,K_\Omega,K_t>0\) on the retained regular positive-density branch is sound. It establishes local static ellipticity and a positive leading scalar time coefficient. It is not a full coupled-mode ghost/gradient theorem, a global stability analysis, or a proof of causal ultraviolet completion.

The same retained principal calculation gives a radial scalar characteristic speed squared of approximately 1.47 near 100 kpc. This result should appear beside the positivity statement, not only in a supplementary file. Time-radial connection terms and higher weak-field terms were omitted, so the many stored digits are not exact cone measurements. The order-unity excess over the metric light speed is nevertheless too large to be dismissed as numerical roundoff. Photons still follow metric null geodesics.

**Suggested wording:**

> The retained local scalar principal coefficients are positive on the selected branch. The same approximation admits radial scalar characteristics faster than the metric null cone, with a maximum squared speed of about 1.47 for the benchmark body. We therefore make no claim of a complete stability or ultraviolet causality certificate.

## 5. Novelty and field identification

The Poisson/scalar equations and screened spherical branch reproduce the standard cubic-Galileon quasi-static reduction. Barreira et al., arXiv:1306.3219v3, Eqs. (11)–(15), give an explicit comparison; their metric signature differs, but their background energy density fixes the matching \(c_2=1\), \(c_3/\mathcal M^3=1/\lambda\): <https://arxiv.org/html/1306.3219v3>. Do not claim discovery of a new screening mechanism, a first Galileon lensing signal, or an independently derived fundamental BHSM coupling.

The supported contribution is narrower: an explicitly reproducible local finite-body realization on the previously frozen R1 background, propagated through the same strong-lens comparison and accompanied by its identifiability/effective-halo limit. The exact effective-halo equivalence at retained \(\Phi=\Psi\) order must remain prominent; independent dynamical information measuring the same potential alone does not remove it.

The distinctions between geometric \(\psi=-\Psi\), independent action field \(T\), target phase \(\eta\), and cosmological time displacement \(\pi_{\rm time}\) are appropriate and should be maintained. The pure cubic specialization supplies a conditional result without proving that the broader reduced fourth-order coefficient is zero.

## Final assessment

No numerical result needs to be discarded on the basis of this targeted review. The calculation supports a **conditional leading-order forward-model result**, with explicit scalar-cone and matching limitations and no halo-independent detection claim. The quantitative ordinary-GR 2PM scale is the most important addition needed to prevent the reported tiny shifts from overstating the calculation's physical accuracy.
