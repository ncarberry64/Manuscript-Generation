# 2MRS luminosity-weighted foreground transfer

Classification: **not_established**

- catalog: official 2MRS J/ApJS/199/26/table3
- photometric field: Kcmag -> analysis field Kmag
- foreground coefficient:
  $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.primary.beta_foreground_mag_per_score_sigma) +/- 0.007744918616606419 mag / score sigma
- t: $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.primary.t_foreground)
- Delta chi2: $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.primary.delta_chi2)
- blocked permutation p: $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.primary.blocked_permutation_p)
- leave-one-shell-out sign stable:
  $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.primary.sign_stable_all_shell_dropouts)
- held-out Delta chi2:
  $(@{baseline=; claim_boundary=System.Object[]; classification=not_established; data=; foreground=; k_magnitude_column=Kmag; low_to_high_transfer=; primary=; schema=r1-2mrs-luminosity-transfer-result-v1}.low_to_high_transfer.validation_delta_chi2)

## Interpretation

The K-band luminosity-weighted foreground coefficient has a stable positive
sign under every leave-one-shell-out check, but the preregistered blocked
permutation threshold is not met and the low-redshift coefficient fails to
transfer to the higher-redshift validation set.

Accordingly, the simple cumulative luminosity-weighted foreground-well proxy
is **not established** as the explanation of the Pantheon+ residual structure.

This result does not test an SMBH-specific or action-normalized BHSM
de-encapsulation field. It tests a stellar-mass/gravitational-well proxy with
a single frozen path weighting.

## Claim boundary

- No weighting exponent was fit.
- No angular scale was rescanned.
- No shell-specific coefficient was chosen after seeing the result.
- Conventional peculiar velocity, density, lensing and survey effects remain
  alternative explanations for the weak same-sign tendency.
