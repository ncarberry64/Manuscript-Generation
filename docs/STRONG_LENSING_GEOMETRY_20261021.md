# HE 0435 geometry correction and first unused JWST application

The old Gaia failure was principally a restricted lens model failure, already
visible in the fitted NIC2 data. A positive explicit nearby G1 gives a useful
geometry correction. It improves the first unused JWST diagnostic substantially,
but **neither conventional nor R1 is validated**. The instrumental response and
exposure-consistency checks fail. No errors were enlarged to obtain a pass.

## What was diagnosed and corrected

The old conventional prediction has Gaia coordinate RMS 3.25560 mas, versus
3.09131 mas on its own fitted NIC2 positions. Gaia and NIC2 differ by only
1.24691 mas RMS. Exact gnomonic mapping changes at most 0.003212 mas; allowing
global rotation/scale changes Gaia RMS only to 3.23931 mas. Labels, parity,
RA-cosine convention, shared-reference covariance and stationary rays were checked.
Refitting all nine old parameters against Gaia alone still gives a local
quadratic residual about 1559 and RMS 2.986 mas. This is a bounded/local diagnosis,
not a theorem excluding every circular profile.

For the circular deflector and constant symmetric effective tide, in its
eigenframe, beta_j=u_j[1-g_j k(r)] implies

    (g1-g2)u1 u2 + g2 beta1 u2 - g1 beta2 u1 = 0.

The omitted conventional freedom tested here is spatial variation of the nearby
G1/G22 tide. G1 is observed at (-2.585,-3.637) arcsec east/north relative to the
main light center, z=0.7821. The old west-positive catalogue coordinate is flipped
once. The old paper's main-lens alias G1 is a different object.

The rear potential is a positive Plummer perturber plus a **separately specified**
quadratic remainder, rather than adding G1 to the old complete environment:

    psi_p(x) = x^T U_rem x/2 + 2 GM_G/c^2 log[(|x-c_G|^2+a_G^2)/(2a_G^2)]
    Sigma_G(R) = M_G a_G^2/[pi(R^2+a_G^2)^2] >= 0,  a_G=0.2 kpc.

The scale is a prescribed regularization of a distant compact optical perturber,
not a measured internal galaxy profile. The fixed remainder convergence anchors
the mass-sheet coordinate; it is not measured total environmental convergence.

## Actual equations and applications

The local unknown radial fields are the regular R1 scalar gradient p and metric
gradient Phi', sourced by the **total** outer+nuclear Plummer mass once. In the
original code's hatted dimensionless variables the evaluated equations are

    rho(r) = sum_j 3 M_j a_j^2/[4 pi(r^2+a_j^2)^(5/2)]
    Menc(r) = sum_j M_j r^3/(r^2+a_j^2)^(3/2)
    C x^2 p - 4 x p^2/lambda_hat = -b_hat m_total(x)
    2 x^2 Phi' + b_hat x^2 p = 2 m_total(x).

The original regular origin/isolated radial branch and relative projected
potential are retained. The compact component is 1e8 Msun, scale 0.03 kpc.
This conditional quasistatic boundary application is not a completed global
Cauchy evolution or microscopic preparation solve.

For each observed-angle unknown theta the actual domain maps and residuals are

    xp(theta) = D_op theta - D_lp alpha_l(D_ol theta)
    c_G = xp(theta_G)
    R(theta) = theta - D_ls/D_os alpha_l
               - D_ps/D_os [U_rem xp + alpha_G(xp-c_G)] - beta0 = 0
    J = I - D_ls D_ol/D_os H_l
          - D_ps/D_os (U_rem+H_G)(D_op I-D_lp D_ol H_l).

Products are ordered; the Jacobian is not artificially symmetrized. Observer
angles map to the rear plane and the same source coordinate. G1's physical center
is traced through each matching foreground field. Both independent stationary
variations of the full two-plane arrival action are checked:

    ct = k_l |theta-theta_p|^2/2 + k_p |theta_p-beta0|^2/2
         -(1+z_l)psi_l(D_ol theta) -(1+z_p)psi_p(D_op theta_p)
    theta_p=xp/D_op
    k_l=(1+z_l)D_ol D_op/D_lp,  k_p=(1+z_p)D_op D_os/D_ps.

The development optimization consumes the exact original light likelihood, the
full moving-origin gnomonic Gaia covariance, NIC2 light center and prescribed
3 mas mass-center prior. Nine parameters versus six relative Gaia coordinates
remain underidentified. No new BHSM response amplitude is introduced.

At r=0.03 kpc the executed original local R1 section is:

| Quantity | Computed value |
|---|---:|
| rho/(Mpl^2 H0^2) | 3.73605305235e9 |
| reduced J_eff/(Mpl H0^2) | -1.71279915113e8 |
| scalar gradient/(Mpl H0) | -3.59919198068e-5 |
| scalar integrated lhs and source | -3.49557028640e-17 |
| metric Gauss lhs and source | 7.62473289969e-16 |

Across the saved sections the scalar relative residual is at most 6.67e-16 and
Gauss at most 3.05e-16. Potential-derivative checks are separately about 3.6e-9;
they are not described as machine precision. Bright-ray residuals are below
1.6e-21 rad and arrival stationarity below 1.5e-17 Mpc. Finite projection tails,
material support and global matching remain conditional.

Development Gaia RMS becomes 0.01055/0.01051 mas (conventional/R1), using unchanged
catalogue errors. This is fitting, not validation. Both finite censuses contain
seven stationary images; central/A ratios are 0.00165054/0.00165038. The consumed
HST F814W test permits them under the same conditional prescription, with upper
scale multipliers about 7.46.

A frozen-body counterfactual replaces only G1 by its value/gradient/Hessian
Taylor polynomial at the rear origin. It preserves the local body and quadratic
remainder. Removing the spatial variation changes relative images by up to
22.56 mas and delays by 12,266 s; solved counterfactual rays/action gradients
still satisfy their original equations. Its consumed-Gaia RMS is 14.08 mas.
This isolates a response within the frozen model, not unique physical causality.

## Unused-photon application and calibration chronology

The physical predictions were frozen and published at **fc3a449 before downloading**
four archival 2022-10-12 JWST/NIRCam F115W NRCB3 native CAL exposures, program 1198 obs007,
1803.776 s. Their SHA identities are retained. Quasar values were excluded by a
10 arcsec guard during field work. The first field ePSF method failed too few
training stars; its separate serialization failure and recovery remain visible.

Independent LMC STDPSF calibration was frozen at **29dcff9 before library array
use**. All four field-star fits failed known-ERR adequacy (chi2/dof 43.44-52.77).
Correct sampling/origin/orientation/spatial interpolation were independently
verified. Native wings differ by 10-23% near radius 4-5 pixels. Additional
positive jitter/halo responses and empirical dilation remain inadequate; none
was transferred as a validated response. This is a mean-response problem under
the tested likelihood, not a source of permission to inflate ERR.

Before opening quasar values, **287b284** froze an explicitly diagnostic amendment
with the original fixed independent PSFs, field-only association translations,
quadratic smooth backgrounds justified from consumed HST and synthetic data,
and unchanged native errors. Every validation flag is forced false. The original
precision-validation experiment remains blocked; its failure is not erased.

The once-only target marker precedes the actual numerical application. All four
exposures and sixteen point fits are retained, with exact native SCI/ERR/DQ stamps
and portable official nonlinear GWCS. No target physical parameter, scale, shear,
PSF shape or error floor is fitted.

| Frozen prediction | Raw five-coordinate RMS (mas) | Formal Q5 | Validated |
|---|---:|---:|---|
| Old conventional | 3.313314 | 18667.775 | No |
| Old R1 | 3.313156 | 18667.225 | No |
| G1 conventional | 0.779400 | 958.434 | No |
| G1 R1 | 0.779454 | 958.428 | No |

The raw reduction is about 76.48%, without wider errors. All sixteen pixel fits
nevertheless fail (chi2/dof 24.72-79.97). Exposure consistency gives Q=432.262
for 15 coordinates, conditional p=1.13e-82. A common predicted geometry cannot
remove this between-exposure scatter: weighted total squared residual equals
scatter about the weighted mean plus the mean-to-prediction squared residual.
Formal Q values under this inadequate extraction are diagnostic, not calibrated
physical-model rejection probabilities. No BHSM preference is established.

## Precise remaining dependency and preserved scope

For source i, exposure e and native pixel p, the instrument mean consumes the
spectral integral of the pupil intensity convolved with image-motion/detector
response and integrated over the native pixel. It depends on the incident SED,
throughput, wavefront, detector response and illumination history. A fixed LMC
ePSF supplies one empirical instance; this field application demonstrates that
it does not supply an adequate response here. Four failed stellar wings do not
uniquely determine the wavefront/SED/detector decomposition or its transfer to
quasar+host light. An adequate mean, diffuse-light model, distortion treatment and
noise/calibration uncertainty are needed before precision geometry validation.
The tested finite response corrections do not solve that problem.

The local R1 effective source is implemented and applied here; positive G1 is a
prescribed total optical potential. Its BHSM scalar environmental state has not
been derived. Material EOS/support, global matching and the cap-to-canonical-field
microscopic source/state embedding remain explicit. In particular the original
bilinear source j_c=mu/2 sum_s[U_s^dagger(iG_c)U_s-V_s^dagger(iG_c)V_s]
is not evaluated by substituting these classical densities for its sections.
No dressed-state factorization or vacuum is silently selected. The exclusive
chronological two-point, inclusive contour four-point and complete matched
occupation tasks retain their distinctions; no later four-point is imposed as
a prerequisite for this conditional classical optical experiment.

The physical guard stays **2e-10**, and **physical_admission=false**. Old retarded
arrays/native heat functions do not become statistical state kernels. Original
background, source, calibration and once-only local/native accounting are retained.
Old LCO coverage belongs to its old frozen geometry and is not assigned to G1.
Other bands remain unopened; no primary switch or target-refit success is claimed.

## Evidence and reproduction

The executable geometry is `code/strong_lensing_geometry_revision.py`; the
development application is `scripts/fit_he0435_geometry_revision.py`. Frozen
protocols are in `preregistration/he0435*20261021.json`. Exact receipts, failed
calibrations, source sections and the target result are in
`artifacts/strong_lensing_geometry_20261021/`. Public native download identities
are retained, while small native stamps and GWCS permit portable inspection.

The focused suite passes **123 tests**; the final diagnostic guard amendment is
separately checked by its 15-test file. Internal independent agent reviews cover
the physical equations, conventions, preopening diagnostic and postopening
scores. These are internal computational checks, not external human peer review.
The inheritance audit confirms original science/data/template bytes are retained;
two inherited README checkout line endings are recorded and not repaired.
