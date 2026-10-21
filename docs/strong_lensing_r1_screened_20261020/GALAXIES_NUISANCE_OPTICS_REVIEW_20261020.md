# Research review: optical comparison and identifiability for the Galaxies revision

> Scope reconciliation: the [current BHSM scope](../STRONG_LENSING_SCOPE_20261020.md)
> supersedes earlier project framing in this retained technical record. The
> governing question is surface contributions to envelope light transfer;
> conventional effective-mass equivalence is a constructive result at the
> retained order. Reference-body offsets are not discrepancies between fully
> matched representations. The full surface-to-R1 common-source derivation
> remains unfinished. Equations, numerical evidence and review chronology
> below keep their original scope.

> This is an AI-assisted technical review of the leading-order calculation. The common GR2 control was added afterward; it does not complete scalar higher-order terms or lift the leading mass-profile equivalence.

## Scope and overall assessment

Reviewed the supplied `49a1dab` main file and local-family section, the supplied earlier `strong_lensing_20261020.tex` body, and the completed frozen-R1 screened calculation and three-source review results. The `49a1dab` review snapshot does not contain its full `strong_lensing_20261020.tex`; the requested earlier body is therefore treated as text to reconcile, not as proof of a missing file in the actual repository. No manuscript or prediction code was edited for this review. No new numerical calculation was performed.

**Recommendation:** make the sourced R1 cubic-scalar encounter and its mass-profile degeneracy the paper's main argument. The earlier transport, SIS and constant-target calculations are useful controls, but their current prominence obscures the result that now answers a concrete physical question. The strongest manuscript is a conditional finite-body optical prediction with a rigorous account of what it can and cannot distinguish. It is not a claim that the constant target phase acquired a source, that the full BHSM closure has been completed, or that a new effect has been detected.

The optics and nuisance algebra reviewed so far are internally consistent within their declared domains. The corrections below concern several genuinely ambiguous equations, stale statements, and the interpretation of the new results. I found no new sign or normalization error in the physical lens equation, stationary arrival derivative, finite-body deflection, or screened projection.

## 1. Required corrections before a research submission

### 1.1 Update the abstract, scope and conclusion to describe the calculation now performed

The present abstract foregrounds the six conditional axial-state columns, rank-zero projection, the 17.094532-day SIS control, and the source-free target-family Hessian. It does not contain the new matter-sourced R1 local prediction. The supplied earlier body also says that an action-derived local encounter is unresolved, that the executed intrinsic encounter derivative is zero, and that the next calculation must supply a local source. Those statements describe the earlier checkpoint, not the new restricted R1 result.

Replace them with a three-way distinction:

1. The constant `S^6` target family still has exactly zero additional target stress and optical response.
2. The separately specified rolling scalar `T` of the R1 cubic compatibility action supplies a nonzero, regular, screened local response under a frozen-background, quasi-static approximation.
3. Neither result completes the general target/topographic/cap constitutive closure or a globally matched BHSM galaxy solution.

Do not describe the new R1 scalar as a nonzero continuation of the old constant-target phase. Its source and action are different. Also do not describe the conventional comparison as a complete independent GR cosmology: the local conventional encounter and screened encounter share the same prescribed external R1 propagation geometry.

The specific SIS delay should leave the abstract. It is an analytic validation datum, not the main physical outcome. The abstract should instead identify the shared finite density, screening scale, approximate conditional delay shifts, and the unresolved mass-profile degeneracy.

### 1.2 Define the linearized data and covariance correctly

At `eq:sl-fisher`, the text calls `y` the observable vector but writes `y=BX+N eta+epsilon`. Raw image positions and delays are nonzero at the reference state. Write either

\[
d=d_0+B\,\delta p+N\,\delta\eta+\epsilon,
\]

or define `y=d-d0` explicitly. State `C` is fixed and positive definite for this calculation before defining a square whitening map. The implemented routine uses Cholesky whitening, not a generic singular-covariance prescription.

The projection

\[
P=I-(WN)(WN)^+,
\qquad F=(WB)^TP(WB)
\]

is correct for the local Gaussian/quadratic problem with unrestricted nuisance directions. Call it the nuisance-projected local information or profile information. A proper marginal likelihood with flat priors on an overcomplete nuisance parameterization is a separate normalization question, which the projector does not need to resolve.

Keep the existing warning that null Fisher directions are not zero-uncertainty estimates. Explicitly state that numerical rank is not a detection significance, and that the illustrative weighting is not an observed likelihood. Exact algebraic rank is invariant under invertible changes of units, but thresholded singular values require the declared state normalization.

### 1.3 Distinguish a predicted difference vector from a physical parameter Jacobian

For the screened test, the single column fed to the projector is the fixed conventional-to-screened optical difference at first order. The model parameters and background are already fixed. A surviving rank of one means **one nonzero nuisance-orthogonal template direction**. It does not identify a newly fitted fundamental scalar parameter or establish one measured BHSM degree of freedom.

If a diagnostic amplitude multiplying that template is introduced, label it as a diagnostic amplitude. A Fisher forecast for an actual coupling requires differentiating the physical model with respect to that coupling, including the induced background and screening changes.

### 1.4 Repair the constant-background/constant-perturbation ambiguity

The earlier envelope section says that a constant texture has zero self entry. Read literally as a statement about the background target, this conflicts with the local-family section's positive phase gradient stiffness. The correct distinction is:

> A constant reference target has zero derivative stress and zero mixed metric/phase response. A nonconstant phase perturbation on that background has a positive quadratic gradient energy. Only a spatially constant phase perturbation has a vanishing derivative-sector self entry, and it is admissible only if the boundary conditions permit it.

The Dirichlet annulus mode in the new local-family section is nonconstant, so its positive self entry is consistent. Do not suggest that its positivity requires adding another action term to an otherwise zero derivative Hessian.

### 1.5 Disambiguate the harmonic amplitudes in the axial transport equations

At `eq:sl-metric`, `Phi` and `Psi` denote physical spacetime potentials. The axial equations subsequently write `P_g=1+2 Phi q` and `M_g=1-2 Psi q`, where `q` is already the spatial harmonic. Without an explicit change of notation this appears to multiply the harmonic twice.

Define

\[
\Phi(N,\mathbf x)=\phi(N)Q(\mathbf x),
\qquad \Psi(N,\mathbf x)=\psi(N)Q(\mathbf x),
\]

and use the lower-case coefficient functions in the axial transport equations. The code uses coefficient rows and performs the harmonic multiplication once; the manuscript should say so.

The symbol `q` is also used for a background action parameter, a target meridian angle, an angular separation vector, and a radial/profile abbreviation. A notation table or systematic renaming is needed. In the new core sections use a distinct symbol for the target angle and reserve `q_s` only for `D_ls/D_os` if that quantity is retained.

### 1.6 Do not present scale estimates as a proved optical error bound

The finite-body section says that the fractional omitted scales are `GM/(Rc^2)` and `R/D_ab`. The first is the natural weak-field expansion parameter, but the second is a geometric thin-plane scale, not a universally derived leading fractional error coefficient. Symmetry and the matching prescription affect the leading thickness correction.

Suggested replacement:

> The calculation retains the leading weak-field thin-plane projection. Compactness and the ratio of the effective encounter scale to the propagation distances are expansion parameters; this work does not derive a complete bound on all finite-thickness and global matching corrections.

For the screened field, the effective encounter is not confined to the material radius `R=3 kpc`. Its screening radius is approximately `212 kpc`, and its exterior tail is explicitly integrated. Use `H_l r_V/c` and the finite-tail check when discussing frozen local coefficients. The material-radius background-curvature estimate does not bound a rolling-scalar or cap response.

Keep exact GR TOV matching separate from the new approximation: the screened encounter is not an exact modified-gravity TOV solution or an exact Schwarzschild exterior. The same prescribed density is retained; any pressure adjustment is only the leading hydrostatic response already described in the R1 derivation.

## 2. Main identifiability result: show the structural factorization

The rank-zero smooth-propagation result should be explained by a short analytic statement, not primarily by singular values near machine zero.

For the finite radial profile with scalar maps, define

\[
\beta=s/D_{os},\quad
\theta_R=R/D_{ol},\quad
\theta_E^2=4mD_{ls}/(D_{ol}D_{os}),\quad
T=(1+z_l)D_{ol}D_{os}/D_{ls}.
\]

With `f_ab=delta ln D_ab`, the data response lies in the source, mass, radius and compensated-MST nuisance span. In the axial source-x convention,

\[
B=N_{\min}\begin{pmatrix}
-\beta_{\rm as}f_{os}\\
f_{ls}-f_{os}-f_{ol}\\
-f_{ol}\\
f_{ol}+f_{os}-f_{ls}
\end{pmatrix},
\qquad
N_{\min}=(N_{s_x},N_{\ln M},N_{\ln R},N_{\rm MST}).
\]

It follows directly that `PWB=0` for every positive-definite covariance under this nuisance model. This establishes the result independently of the particular six cosmological columns. For the SIS, source-x, log strength and compensated MST already span the three active response rows; the larger nuisance model additionally saturates all five observable rows.

State the scope precisely: this is a scalar-map, specified-profile, retained-order identifiability result. It is not an exclusion of all cosmological lensing information, all non-axial distortions, all multi-source configurations, or all higher-order effects.

A compensated mass sheet is an admitted extra model freedom, not a member of the finite Tolman-VII `(M,R)` family. It must be identified as such. Also distinguish one physical sheet shared by several source planes from arbitrary independent per-source sheets; they are not interchangeable nuisance models.

## 3. Put the new conditional optical prediction in the main text

The shared finite-body comparison is now the meaningful numerical result. A compact table can show:

| Quantity | Conventional local encounter | Screened R1 local encounter / difference |
|---|---:|---:|
| Central relative delay | 21.716076566 d | increase 25.959258 s |
| Negative outer relative delay | 14.271452387 d | increase 8.649019 s |
| Three image shifts | reference | +6.123276, +0.788107, -4.348592 microarcsec |
| Central/outer delay ratio | 1.5216444674 | 1.5216548468 |
| Ratio fractional change | reference | 6.82122e-6 |
| Screening radius | not applicable | 212.029 kpc |
| Unscreened force enhancement | not applicable | 0.001507886 |

Use fewer digits in the narrative. The machine-readable output can retain full precision. Identify the two image pairs defining every delay and ratio; the earlier softened-isothermal example used the inverse ratio ordering, so copying its numerical convention would confuse readers.

State the source of the prediction: the already frozen R1 compatibility action and background, the prescribed finite density, and the regular screened branch. The result does not use a tuned extra lensing amplitude. It also does not establish that the full BHSM microscopic theory selects this effective action or body.

The optical projection deserves an explicit equation,

\[
\delta\hat\alpha(b)=4b\int_0^\infty
\frac{\delta\Phi'(\sqrt{b^2+z^2})}{\sqrt{b^2+z^2}}\,dz,
\]

with the potential specified by a finite radial primitive or a difference relative to a reference impact. An unregularized absolute line integral of a `1/r` potential has a logarithmic constant divergence; the code correctly fixes an irrelevant common potential constant at `b=R` instead. Do not describe `R` as a cutoff of the exterior gravitational field.

One sentence can summarize validation: direct-z projection independently agrees with the angular-substitution implementation at maximum relative error `4.43e-13`, and independent image roots agree with the stationary first-order response at the expected second-order level. The long residual census belongs in the supplement.

## 4. The effective-halo equivalence is a central scientific result

At the retained order the new scalar has `Phi=Psi`, and

\[
\delta\Phi'(r)=\frac{GM}{c^2r^2}\,\mu(r)F(r),
\qquad \mu(r)=M(<r)/M.
\]

Therefore one conventional spherical halo with

\[
\delta M_{\rm eff}(r)=M\mu(r)F(r),\qquad
\delta\rho_{\rm eff}(r)=\frac{M}{4\pi r^2}\frac{d[\mu(r)F(r)]}{dr}
\]

reproduces the full extra optical potential up to its irrelevant additive constant. This is stronger than merely finding a nuisance fit for one lens system.

For the retained solution, put `s=sqrt(1+r_V^3 mu/r^3)`. Then

\[
\frac{d\ln(\mu F)}{d\ln r}
=\frac12(1+s^{-1})\frac{d\ln\mu}{d\ln r}
 +\frac32(1-s^{-1})\geq0.
\]

The effective halo has nonnegative density and total additional mass

\[
\delta M_{\rm eff}(\infty)=M F_\infty
=1.50788624\times10^8 M_\odot
\]

for this body. It is one common halo for every source plane, not a separately adjusted mass distribution per source.

Because both metric potentials agree, the same construction reproduces leading nonrelativistic dynamics as well as lensing. Kinematic information can constrain a mass-sheet nuisance, but it does not by itself distinguish this scalar solution from an unrestricted conventional halo with the same effective gravitational mass. A claim that kinematics automatically provides a theory discriminator would therefore be too strong.

This result should appear next to the new prediction, not only in a limitations footnote. It defines the correct interpretation: a fixed-matter conditional prediction, not a lensing-only identification of modified gravity against arbitrary matter models.

## 5. How to present the single-source and three-source projections

### Single source

The old source-free target null is exact throughout its family. The new screened difference is nonzero at fixed physical inputs, but its first-order observable vector is absorbed when source, mass, radius, compensated MST and center are allowed to vary.

The fixed-center result is a small residual direction, not the full roughly 26-second shift. Once the center is free, the local projected rank is zero. A conventional nuisance change involving a center displacement of only `1.77e-9 arcsec` reproduces the first-order vector. Report the structural conclusion; the detailed compensation coefficients can be supplemental.

Do not report zero rank for log mass/radius columns as evidence that mass and radius are unmeasurable: those same columns were also inserted into the nuisance matrix. That zero is built into the chosen question.

### Frozen three-source design

The three-source calculation is useful because its geometry was fixed before testing the new local response. It retains source positions, common mass/radius, common center, common physical convergence/shear, and independent delay-scale columns. Adding independent source-convergence columns removes the cosmographic distance response but leaves a small screened radial-shape template.

The raw whitened screened norm is `0.01515279989`. It falls to `2.02264543e-5` with the shared nuisance set and to `1.82168870e-5` after adding the independent source-convergence columns: only `0.12022%` remains in the latter calculation. With the illustrative weighting, the largest fitted residuals are about `6 nanoarcsec` and `0.086 s`. These are model residual scales, not measured precisions or a forecast.

This result can receive one main-text paragraph and a compact row in an identifiability table, with the geometry, all columns and the nuisance-orthogonal statistic in the supplement. It must immediately be limited to the shared Tolman-VII profile family. The positive effective-halo construction removes the distinction once a flexible common radial halo is admitted.

Also state that the shared convergence/shear are phenomenological weak tidal fields in the optical model. The scalar field remains the isolated spherical solution. The calculation does not solve the nonlinear cubic scalar in a three-dimensional external tidal environment, nor is each independent convergence column a complete multi-plane matter realization. It is an identifiability design for the specified forward model.

## 6. Recommended main-text versus supplementary allocation

| Material | Recommended position | Reason |
|---|---|---|
| R1 action, source/mixing constraint, regular branch and screening function | Main text | Identifies the actual physical mechanism and fixes the prediction. |
| Shared finite density, fixed physical source, redshift/clock convention | Main text | Makes the comparison reproducible and fair. |
| Compact Jacobi/lens/arrival equations and complete intrinsic derivative | Main text | Explains how the metric produces the observables. |
| Screened image/delay changes and force-profile figure | Main text | Concrete conditional outcome. |
| Effective positive-halo equivalence | Main text | Determines what the prediction can distinguish. |
| Structural scalar-map factorization and single-source local degeneracy | Main text, concise | Explains why response is not the same as inference. |
| Frozen three-source shape result | Brief main-text paragraph; full supplement | Shows the difference between finite-family and unrestricted profile identification. |
| Full six-anchor cosmological reduction/transport and raw arrays | Supplement | Valid machinery, but secondary to the new local prediction. |
| SIS 17.094532-day benchmark and arbitrary zeta=1e-3 trial | Supplement | Validation and normalization controls, not a selected physical anomaly. |
| Softened-isothermal central-image ratio example | Supplement or omit | Repeats the finite-profile lesson with another assumed galaxy. |
| Constant-target derivation | Short main-text recovery control; full appendix | Preserves the earlier null without competing with the sourced scalar. |
| Two-mode annulus Hessian, all numerical entries, ADM auxiliary inverse checks | Appendix | Necessary provenance/control, not the mechanism of the sourced R1 effect. |
| Exact GR TOV residuals and induced EOS details | Appendix | Supports the control; leading optics only uses the prescribed density. |
| Synthetic prior-assisted rank and illustrative Fisher singular values | Supplement | Not data, not a forecast, and basis/prior dependent. |
| Three-source angular crossratio versus the declared closed-LambdaCDM comparator | Separate supplement or omit | A distinct background-cosmography question; not the local screened effect. |
| Hashes, solver tolerances, test counts, replay receipts | Reproducibility supplement | Essential audit trail, distracting in the narrative. |

## 7. Suggested paper structure and framing

Suggested title direction: **Finite-Body Strong Lensing in a Screened Cubic Scalar Model: Matched Predictions and Mass-Profile Degeneracies.** A BHSM connection can be stated accurately in the introduction and model provenance without implying that the complete fundamental framework has been established by this calculation.

1. Physical question and precise effective-model scope.
2. Frozen R1 action/background and the sourced local screened solution.
3. Shared finite-body optical experiment and approximation hierarchy.
4. Conditional image, delay and ratio predictions, with concise validation.
5. Nuisance projection, the three-source design, and the effective-halo equivalence.
6. Conclusions: what has been predicted, what is not identifiable, and which independent physical information would be required next.

The paper's contribution should be the reproducible application, internally consistent matching conventions, specified-model prediction, and explicit degeneracy analysis. Do not claim novelty for Vainshtein screening, mass-sheet invariance, generic Fisher projection, or scalar-to-effective-mass rewriting without a literature comparison. Use the primary screened-scalar, strong-lensing degeneracy and multi-plane references selected by the literature review.

Remove project history from the scientific argument: the previous muon preparation, inclusive pair-production calculation, internal checkpoint language, and accepted-contribution/submission scheduling do not explain the lensing result. Replace version-specific claims such as '51 focused tests pass' in the main text with a compact reproducibility statement tied to the final archived release.

## Recommended conclusion sentence

> For the declared R1 effective action and fixed finite matter distribution, the regular screened solution produces nonzero, reproducible image and time-delay changes. Their interpretation remains conditional: source and lens freedoms absorb the single-source difference, a small radial-shape component survives the stated three-source finite-profile model, and a positive shared effective halo reproduces the full leading-order optical modification when general radial matter profiles are admitted.
