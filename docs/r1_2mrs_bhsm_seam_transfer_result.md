# BHSM seam-charge-weighted 2MRS foreground transfer

Classification: **not_established**

- BHSM source scaling: V_Sigma^2 proportional to sqrt(M_b)
- fixed K-band mass proxy: M_b proportional to L_K
- object weight: 10^(-0.2 m_K) d_L
- foreground coefficient:
  $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.primary.beta_foreground_mag_per_score_sigma) +/- 0.007804864269566229 mag / score sigma
- t: $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.primary.t_foreground)
- Delta chi2: $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.primary.delta_chi2)
- blocked permutation p: $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.primary.blocked_permutation_p)
- leave-one-shell-out sign stable:
  $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.primary.sign_stable_all_shell_dropouts)
- held-out Delta chi2:
  $(@{claim_boundary=System.Object[]; classification=not_established; data=; foreground=; low_to_high_transfer=; primary=; schema=r1-2mrs-bhsm-seam-transfer-result-v1; source_chain=}.low_to_high_transfer.validation_delta_chi2)

This is the first foreground test whose object weighting is taken directly
from the retained spherical BHSM seam-radius plus seam-charge scaling rather
than from raw counts or linear luminosity. It is still only a proxy because
the K-band luminosity-to-baryonic-mass relation and the export of seam charge
into the dynamical cosmological topographic field are not action-normalized.
