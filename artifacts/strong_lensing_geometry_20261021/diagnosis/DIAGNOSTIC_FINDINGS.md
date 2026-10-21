# Consumed-observation geometry diagnosis

Every numerical observation used in this diagnosis was already exposed in the
20261020 work. None is a new validation target. No uncertainty floor is fitted
or added to Gaia. All original files and numerical arrays are preserved.

The original conditional model misses its own NIC2 development coordinates by
`[2.603, .712, 1.915, -4.171, -5.223, -1.311]` mas (NIC2 minus model).
Gaia minus NIC2 is only `[-1.980, .990, -.00047, -.856, -.349, 1.891]` mas.
A global scale and rotation fit reduces the old-model Gaia coordinate RMS from
3.25560 to 3.23931 mas; its scale is 1.000206 and rotation .005348 degrees.
Thus a coordinate rotation/scale error is not the dominant cause.

With all original nine parameters free, no delays or centroid priors, and the
actual six-coordinate Gaia covariance, the old positive radial-source plus
constant quadratic-screen family reaches Q=1559.083 and RMS=2.98605 mas.
Its outer Plummer scale approaches the .2 kpc lower bound. Fixing the mass
center and screen convergence gives essentially the same result (Q=1559.065).
Adding the historical light center, center prior and a one-day development
delay control does not remove this obstruction. These are numerical local
development fits over the declared bounded family, not a theorem excluding
every possible radial matter profile.

In the eigenframe of a constant symmetric optical map G, a circular deflector
obeys beta_j=u_j[1-g_j k(r)], u=theta-c. Eliminating k gives

`(g1-g2)u1u2 + g2 beta1 u2 - g1 beta2 u1 = 0`.

All circular radial profiles retain this rectangular-hyperbola relation; their
radial equations impose further restrictions. Near the best old-family fit,
the six Gaia rows have five strong angular directions and one tiny singular
direction. This explains why merely tightening the optimizer does not remove
the discrepancy in that family. Four relative point positions alone do not
establish that the unique real missing field is a particular perturber.

The new correction uses the independently documented nearby galaxy G1/G22 at
observed east/north position (-2.585,-3.637) arcsec and redshift .7821. Its
positive regular Plummer potential supplies spatially varying tidal terms
missing from a quadratic environmental potential. The measured center is
traced through the foreground lens:

`x_p = D_op theta - D_lp alpha_l(D_ol theta)`

`c_G = D_op theta_G - D_lp alpha_l(D_ol theta_G)`.

The original full two-plane residual is

`r = theta - (D_ls/D_os)alpha_l - (D_ps/D_os)[U x_p + alpha_G(x_p-c_G)] - beta0`.

Its ordered Jacobian is

`I - (D_ls/D_os) H_l D_ol - (D_ps/D_os)(U+H_G)[D_op I-D_lp H_l D_ol]`.

The action is

`ct = .5 k_l |theta-theta_p|^2 + .5 k_p |theta_p-beta0|^2 - (1+z_l)psi_l - (1+z_p)[.5 x_p^T U x_p+psi_G]`,

where `k_l=(1+z_l)D_ol D_op/D_lp`, `k_p=(1+z_p)D_op D_os/D_ps`,
and `theta_p=x_p/D_op`. Both independent action gradients and the unsymmetrized
ordered ray Jacobian are tested. A zero-mass perturber at the old .8 screen
redshift reproduces the original conventional and original R1 applications.

The added source is positive:

`Sigma_G(R)=M_G a_G^2/[pi(R^2+a_G^2)^2]`,

`psi_G=2 G M_G/c^2 log[(R^2+a_G^2)/(2a_G^2)]`.

The scale .2 kpc is a fixed regular point-mass surrogate, not an inferred
galactic profile. The fitted quadratic remainder is a different jointly
specified environmental term; the old total optical potential is not inserted
again on top of the new positive source. Its kappa is anchored to the inherited
value to fix a mass-sheet degree. This is a declared conditional choice, not a
measurement that the new total convergence equals the previous one.

The final development calculation uses six Gaia coordinates, the old two
lens-light coordinates, the prescribed 3 mas mass-center prior, and the exact
old development-light likelihood. There are nine free parameters: eight of
the old local/source/screen/center parameters plus G1 log mass, with kappa
fixed. It is not a validation using Gaia and not a complete posterior.

The first exploratory G1 fit used the previous small-field Gaia coordinate
map. It reached G1 mass 2.4157771e11 solar masses, Q_Gaia=.01038592 and
RMS=.00661347 mas for the conventional model, and Q=.01271245 and
RMS=.00751532 mas for R1. Those exploratory figures are preserved under
`diagnosis/`; they are not the final frozen prediction.

The final fit uses exact gnomonic coordinates centered on A and the analytic
six-by-eight covariance pushforward from all eight catalogue alpha-star and
declination perturbations, including the uncertain reference origin. Gaia's
RA errors already describe alpha-star; each catalogue row's own local tangent
basis is used. The largest coordinate correction is .00321176 mas, far too
small to explain the original approximately 3.26 mas failure. The old frozen
Gaia arrays and results remain unchanged.

The final `positive_g1_development_fit.json` and matching image censuses give:

| Model | G1 mass, solar masses | Gaia Q | Coordinate RMS, mas | Central/A |
| --- | ---: | ---: | ---: | ---: |
| Conventional | 235148467074 | 0.026511601827 | 0.0105512092357 | 0.00165054058905 |
| R1 | 235146069817 | 0.0262799707485 | 0.0105059067542 | 0.00165037739633 |

The final conventional six-coordinate whitened derivative has six strong
directions (smallest singular value 3.77491 versus .002 in the exploratory
old-family trial). Both finite searches find seven stationary images, with
all three additional-image contributions individually retained.

The precision rerun also reached a better weak-nuisance optimum: the
conventional exact development objective changed from -19200.5925683 in the
exploratory fit to -19201.1919554 in the final fit. The parameter and delay
changes cannot be attributed physically to the tiny coordinate correction;
these fits contain weak nuisance directions and are not a converged complete
posterior. The final R1 objective is -19201.1918397. The tiny conventional/R1
objective difference does not establish a state preference.

The physical point predictions were frozen at commit `fc3a449` before the
new native-file acquisition and target application. The unusually small
in-sample Gaia Q is not evidence of predictive success. No uncertainty floor
or wider-error validation is used.

The local R1 calculation retains the original nonlinear scalar equation on
the TOTAL local two-Plummer density once. Its section residuals remain near
floating-point precision. The external G1 is a prescribed total optical
potential on another plane; no independent nonlinear scalar fields are summed
in a shared domain, and no empirical BHSM interface-response amplitude is
introduced. A globally matched scalar environment and its physical state are
not derived by this fit. Physical admission remains false and the physical
guard remains 2e-10.
