# Strong lensing scope and assumptions

The first comparison retains a conventional galaxy-potential encounter and evaluates the available Berger–Hopf geometric response with the same incoming optical preparation. It tests image positions and relative arrival times, including degeneracies. A local envelope replacement of the galaxy encounter is not yet derived. The submission target is **25 February 2027**.

## Accepted invitation scope

> Thank you for the invitation to contribute to “Strong Gravitational Lens Modeling.” I am pleased to accept and am preparing a theoretical and computational study of strong-lensing time delays in a spatially modulated closed cosmology within the Berger–Hopf framework. The planned manuscript will compare a conventional galaxy-potential benchmark with a proposed geometric envelope response, using the same incident-light preparation to examine observable differences. Particular attention will be given to time-delay predictions, mass-sheet and line-of-sight degeneracies, and reproducible numerical calculations, with assumptions and derived results clearly distinguished. I aim to submit the manuscript by 25 February 2027.

## Governing physical comparison

For a common incident optical preparation and a specified measurement protocol, determine which image-position and relative-arrival-time differences follow from the conventional galaxy encounter and the available envelope equations, then test which survive source, galaxy-model, mass-sheet, line-of-sight and cosmological nuisance freedoms.

In the BHSM interpretation, gravitational lensing is the optical response of an envelope event traversed by environmental variables, with light providing the probe here. The calculation relates prepared incoming light to outgoing observables, conditional on the envelope state and surrounding geometry. Light probes state dependence, effective geometry and any response actually derived from the equations. Traversal alone establishes no encapsulation, de-encapsulation, reflection, absorption or new photon interaction law.

The executed first realization treats the conventional galaxy potential as the main local thin encounter. Available smooth BHSM metric perturbations modify reference propagation. A local envelope contribution enters separately only after its action projection and source normalization are supplied. No gravitational or line-of-sight mean contribution is counted twice.

## Matched preparation and conventions

- Use one short test pulse with common local emission spectrum, polarization and physical transverse source coordinate; neglect its backreaction.
- Use orthonormal observer/lens/source frames, parallel Sachs screens, radians for angles and physical lengths for angle-to-separation maps.
- Normalize the affine frequency magnitude to one at the observer and define positive observed redshift by the frequency ratio.
- Fix measured redshifts at the illustrative lens/source values 0.3 and 1.0. Include endpoint motion and moving lens-to-source initialization in any response claimed to hold these redshifts fixed.
- Fix the transverse physical source coordinate in the declared screens. Recompute inferred angular source coordinates and image positions when maps change. Do not freeze all three independently.
- Let longitudinal position and the source emission event follow the fixed observed-redshift surface in each realization. This transverse matching does not additionally prescribe a complete three-dimensional source world-line across varied geometries.
- Keep the main encounter out of all three reference maps if represented as a separate potential. Initialize the lens-to-source map at the lens; do not subtract observer distances.

The incoming optical pulse, the coupled cosmological radiation-fluid perturbations and muon/reference-state preparation data are separate objects. The current muon preparation is not an admitted lensing calibration and is not required by the optical equations executed here.

## Inherited state and domain

The current action-native coupled provider is primary. For one selected normalized n = 2 spatial harmonic, the physical anchor is `(zeta, Pi_2, delta_m_ADM, delta_r_ADM, v_m_ADM, v_r_ADM)` at redshift 2.1. `Pi_2` is the reference-normalized scalar momentum and is mapped to the full coupled canonical state; it is not silently identified with `p_zeta`. Nine spatial n = 2 components form the full eigenspace; selecting one profile leaves six initial amplitudes. The historical two-column seed is a restricted slice.

Use the frozen R1 background with `Omega_m0 = 0.31`, `Omega_r0 = 9e-5`, `Omega_k0 = -0.018`, `H0 = 67.4 km/s/Mpc`, `lambda = 1`, `q = 3`, and `V0/(Mpl^2 H0^2) = 2.3986073449708067`. Reuse the supported redshift range 0 through 2.1 and the existing documented background/stability work.

For a normalized observer-dipole profile on the reference three-sphere, `Q = sqrt(3)/(pi R_H^(3/2)) sin(2 psi) (n dot p)` has integral `Q^2 dV = 1`. On the unit sphere omit the radius factor. The chosen axis is conditional. It is not selected by isotropy or established by past supernova fits. Lensing uses the conventional `Phi + Psi` potential combination; legacy provider/sign conventions require an adapter before reuse.

## Claim boundaries and dependencies

| Object | Status at the first comparison |
| --- | --- |
| Conventional galaxy model and physical source coordinate | Declared illustrative inputs, not an observational fit |
| Frozen R1 background and coupled action evolution | Existing numerical owner reused |
| Three compatible Jacobi maps, relative arrival functional and complete first variation | Mathematical transport construction and executed controls |
| Selected coherent cosmological anchor amplitudes | Conditional state inputs; unit sensitivities are not fitted physical amplitudes |
| Static spatial phase mixed source with zero spatial current or constant local texture | Exactly zero in this channel |
| Local physical metric/topographic Hessian entries | Require the normalized common-domain action projection; lapse/shift auxiliary blocks are not substitutes |
| Static seam export | Inherited radial profile; common-action normalization/orientation remains unresolved |
| Moving seam time-derivative source | Proposed kinematic continuation, not used as an action-selected law |
| Local galaxy replacement by an envelope | Uncomputed until the encounter metric/source is derived |
| Stochastic optical covariance | Requires source statistics, including correlations; no Gaussian/Wick closure inferred |

The local inverse response requires `a,b,c,u,v,c_phi`, basis/profile norms, the determinant or constrained invertibility condition, response residual, Schur curvature and common boundary/domain status. These are action projections, not fitted parameters. Preserve the explicit mixed source and do not infer its absence from an unavailable numerical array.

## Observable information and nuisance assumptions

Whiten the declared covariance, project physical response columns away from the admitted nuisance span, and report surviving singular values, rank and parameter-space null directions. An unconstrained direction has no finite likelihood-only uncertainty; a Fisher pseudoinverse must not manufacture a zero error bar. Identify information contributed by priors separately.

An independent unrestricted mass-sheet scale per lens absorbs a common rescaling of that lens's relative delays. Additional sky directions do not by themselves remove it. A pure-dipole Born-order shear null does not imply zero redshift, focusing or total lens response. A low-rank covariance result applies only to the stated coherent amplitudes and fixed response.

The familiar inverse Hubble-constant delay scaling assumes fixed dimensionless background and reduced angular galaxy/source conventions. It does not automatically apply to a variation holding the physical transverse source position fixed.

## Reproducibility and source continuity

The focused manuscript is `manuscript/strong_lensing_main.tex`, with `manuscript/sections/strong_lensing_20261020.tex`. The existing Universe manuscript and bibliography remain intact. The inspected source revision is `3427c5297fac71faa501685bd38e5d420619c761`, on `codex/strong-lensing-first-comparison-20261020`; the older `5fd65ea18a6055dbc6028839c2918415e3c00f71` is provenance only. No reset is implied.

The supplied request identifies `04ajb_coupled_environmental_state.tex`; that path is absent at the inspected revision. The existing action-native implementation and the explicitly supplied six-anchor contract provide the corresponding state convention; no missing file is fabricated.

Primary external checks: [Fleury, Larena and Uzan, arbitrary-spacetime lensing, equations 2–3 and 11–14](https://arxiv.org/html/2011.04440v4), [their dominant-lens treatment of line-of-sight effects](https://arxiv.org/abs/2104.08883v5), and [Teodori et al., mass-sheet and differential line-of-sight degeneracies](https://arxiv.org/abs/2201.05111v3). These references support the conventional optical/degeneracy framework; they do not derive BHSM envelope coefficients.
