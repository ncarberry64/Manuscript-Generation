# Final methods audit

No frozen BHSM parameter, direction, amplitude, transfer function, cosmology, original sample cut, or environment density split was optimized during this work.

Corrections addressed data parsing or mathematical consistency:

- Parse SNANA VARLIST headers with and without a colon; use exact SNID before an IAUC alias and the released LIST file where duplicate versions occur.
- Apply the released SNANA 0.27 magnitude normalization. Use unstarred named observer passbands rather than overwriting them with auxiliary starred curves. Honor the DES LOWZ FILTERMAP and slash-suffixed band codes; retain DES nominal passband-proxy limitations.
- Preserve rank-failed jackknifes without discarding the nominal fit. Do not assign numerical noise in an exhausted residual subspace to a physical phase component. The separate orthogonal-decomposition table is authoritative for additive residual powers; the original fit-summary band/phase diagnostic powers are separate projections and are not additive.
- Transform foreground-cell covariance perturbations through the same training/test operator as the residual vector. A kernel on untransformed rows cannot be applied directly to transformed residuals.
- Keep each permuted SN's exact redshift transfer; exchange only angular-template labels inside fixed survey/redshift/sector strata and retain duplicate blocks.
- Remove the external PV difference from fitted nuisance features. Its fixed mean enters only the independently labeled coefficient-one sensitivity. The superseded development output is retained in an explicitly labeled audit archive, not used in the final report or likelihood tables.

The observed 2MRS galaxy-count supplement was added after the primary field design and does not modify any tested nuisance feature, threshold, prediction or likelihood. Its incompleteness and redshift-space coordinates are explicit.

The entire current analysis remains conditional on measured source/environment features. It is not a replacement for a jointly propagated SNANA/calibration/BBC and external-field likelihood.
