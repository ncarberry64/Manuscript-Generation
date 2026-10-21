# Converged conditional static baseline

The full existing mean has now been optimized jointly, rather than stopped after one continuation step. It retains the supplied 318-mode positive field response, 36 nonnegative shared sky coefficients, independent point positions and positive fluxes, detector background planes, and all 24 field-point nuisance parameters. The observed SCI/ERR arrays, priors, finite ePSF domain, and source family are unchanged. Neither conventional nor BHSM lens predictions enter this fit.

Convergence requires both a dimensionless projected KKT gradient of the half-objective no greater than `1e-3` and a final relative full-objective change no greater than `1e-10`. This threshold concerns numerical optimization; it changes neither detector ERR nor the `2e-10` physical guard. The first all-data stop satisfied relative objective convergence but missed the KKT threshold, so it was continued with stricter numerical termination. That first termination receipt remains available.

| Application | Raw target chi-square | Raw field chi-square | Projected KKT maximum |
|---|---:|---:|---:|
| All four target exposures | 33639.0078831 | 9995.6544215 | 0.000048959 |
| Target exposures 1–3 | 23295.0188837 | 10593.4672608 | 0.000353515 |
| Held exposure 4, strict independent point fit | 11018.9817616 | — | 0.000166121 |

The all-data objective is `44420.2979866`; the training objective is `34631.5436856`. The full linearized joint ranks are 438/438 and 423/423 respectively, including the declared priors. Five shared sky coefficients lie at zero in each fit. Covariance is explicitly unconstrained local curvature, not a truncated posterior or evidence of physical identification.

Convergence does **not** repair the held mean. Its raw chi-square is slightly worse than the earlier one-step application (`10906.4832918`). Per-image held chi-square is A 2738.4, B 3638.7, C 2517.1, D 2124.7. Only about 15–28% of that residual power lies within three native pixels of the point centers; adjacent standardized residual correlations are approximately 0.08–0.28. This is broad spatial error in the current conditional mean. A larger covariance alone cannot count as repairing it.

`static_diagnostics/` retains native target and field maps, including SCI, ERR, DQ, masks, fitted means, raw residuals and standardized residuals. The full 32-coordinate train–held covariance includes the shared response/sky cross blocks; its minimum eigenvalue is `2.1308573e-6 pixel²`. Raw and finite-grid moment conventions remain separate. Exposure 4 targets never update training response or sky. All four separate field-star exposures remain calibration data. This is consumed F115W development, not unused validation.

Exact commands from the worktree root:

```powershell
.venv-geometry/Scripts/python.exe scripts/run_he0435_response_completion.py --label all_data --max-nfev 400
.venv-geometry/Scripts/python.exe scripts/run_he0435_response_completion.py --label all_data --max-nfev 400 --resume
.venv-geometry/Scripts/python.exe scripts/run_he0435_response_completion.py --label training --max-nfev 400
.venv-geometry/Scripts/python.exe scripts/finalize_he0435_static_response.py
.venv-geometry/Scripts/python.exe -m pytest tests/test_strong_lensing_response_completion.py -q
```

Three focused synthetic tests check the full coupled analytic Jacobian, bound signs, and converged nonlinear recovery. They pass. Existing provisional results are preserved; new results live only in this completion namespace.

This static result supplies an initializer for the actual read-history calculation. It does not supply the missing read intercepts, actual segment weights, detector counting covariance, or microscopic envelope statistical kernel. Those must come from the matching native instrument equations and recovered histories. Any parameter-dependent Gaussian detector likelihood must include its covariance determinant and derivatives, with immutable raw native-ERR residual metrics retained separately. No physical admission or reserved-data test is claimed.
