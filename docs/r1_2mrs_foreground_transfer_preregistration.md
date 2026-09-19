# R1 / BHSM 2MRS foreground-transfer test

This test is the first real-data replacement for the illustrative `z_c=0.03`
coherence cutoff in the LOS-topography candidate.

For SN `i`, define shell density contrast `delta_s(n_i)` from a fixed
10-degree Gaussian-smoothed 2MRS map. The shell contributes only to the
fraction of its comoving path lying in front of the SN. The single frozen
foreground score is

    F_raw,i
      = [C_K(chi_i)/S_K(chi_i)]
        sum_s DeltaChi_visible(i,s) delta_s(n_i).

`F_raw` is standardized without using SN residuals.

The primary regression adds this one score to the fixed-axis dipole and the
available host/survey nuisance controls. No shell-specific coefficient and no
new angular-scale scan are allowed.

Primary significance is the improvement in weighted chi-square relative to
the same model without the foreground score, evaluated against redshift-blocked
permutations of the foreground score.

A secondary transfer check fits the one foreground amplitude at
`0.01 <= z <= 0.03`, freezes it, and carries it to `0.03 < z <= 0.15`.

Interpretation remains deliberately narrow: 2MRS density is a proxy for
foreground structure. Even a positive result is not by itself evidence that
the source is BHSM de-encapsulation rather than conventional peculiar velocity,
lensing, density, or selection effects.
