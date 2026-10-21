# Coupled native-history application

This executes the existing positive 318-coefficient response, 36 shared positive
sky coefficients, free point/field centroids and fluxes, exposure planes, and a
shared nonnegative collection coefficient. All F115W observations here were
already consumed. No lens prediction enters the objective, no original ERR is
rescaled, and no reserved band is opened.

For effective optical source rate `f P_eta + S_eta a`, the native-reported SCI
conversion is multiplied by `TGROUP/(5 TFRAME)` before the physical read-time
application. The electronic rate plane is kept separate. The actual mean is:

1. `lambda_pre = gain * G * (lambda_effective/gain)`, using the exact CRDS IPC
   deconvolution and a fixed positive harmonic exterior on the 65-pixel domain.
   Trusted source squares are 33 pixels for targets and 17 for FIELD. The
   induced initial native-ROI mean changes are recorded, not silently discarded.
2. Original leading collection mean with conditional translated ceilings:
   `H_p(t)=integral_0^t min(lambda_p s,C_p) ds`,
   `D_p=sum_q (lambda_p+lambda_q)[H_q-H_p]`,
   `mu_p(t)=min(lambda_p t+kappa D_p(t),C_p)`.
   This is a first-order storage-ceiling hypothesis. Lost/overflow charge,
   physical well capacities and an exact nonlinear stochastic completion are
   not inferred by imposing the SAT reference cutoff.
3. Fit the electronic intercept beta from actual G0 through
   `mean_r L^-1(beta + mu(t_0r)/gain + dark_0r) = G ADC_G0`.
   Beta is differentiated implicitly, not held as an independent noiseless
   reset value.
4. Apply `L_native[H mean_r L^-1(beta+mu/gain+dark)] - native_group_dark`,
   then the original native branch-dependent weighted OLS equations and the
   flat/PHOTMJSR map. `H` is the order8 positive Neumann approximation to `G^-1`,
   with omitted l1 mass bounded by 2.07573e-15. No extra optical PSF convolution
   is added.
   ZEROFRAME follows its own actual branch and native dark convention.

The optimizer uses the continuous double extension of the exact native OLS
source equations. Official float32 OLS_C values are replayed at the final saved
mean. Their small numerical difference is retained explicitly. The full reverse
derivative includes both harmonic maps, translated caps, implicit beta, IPC,
linearity, ZEROFRAME and actual segment/median weighting. An actual target and
FIELD directional check has relative L2 errors below 7e-10.

| Application | Target raw Q with original ERR | Field raw Q | kappa (1/electron) | Scaled projected KKT |
| --- | ---: | ---: | ---: | ---: |
| All data, continuous mean | 33330.090404 | 9916.532669 | 1.594530921e-7 | 1.32149e-4 |
| E1--3 targets, all four FIELD | 23005.301878 | 10532.883972 | 1.654003367e-7 | 1.64766e-4 |
| E4 targets, frozen trained shared state | 10974.065189 | Not fitted again | Fixed from training | 1.91219e-4 |

The all-data joint rank is 439/439; training rank is 424/424 under the declared
priors. Positive sky bounds
remain active. Covariances are full local UNCONSTRAINED curvature, not a
truncated posterior or a coverage result. Held propagation retains trained
sky/eta/kappa crossblocks and their correlations with trained centroids; the
full32 covariance minimum eigenvalue is 2.115898519e-6 pixel squared.

The strict held official-native replay is Q=10974.065420, compared with the
converged static held Q=11018.981762. The 0.41% improvement is not a repaired
profile. Un-JUMPed held pixels alone contribute Q=9867.836391; about79% of the
held Q lies outside3pixels. The same frozen optical body with kappa=0 and its
fitted collection response retain most observed near-wing temporal curvature.
Raw native pixel means, errors and residuals are in native_replay_maps.npz;
spatial and temporal diagnostic sidecars retain the unresolved components.

The bounded background-partition counterfactual assigned the entire fitted
plane to photon carriers before IPC/collection and removed its late electronic
addition. All24 applications stayed in the positive trusted domain. Held Q
changed from10974.065420 to10972.215765, only0.017%; near-wing curvature changed
by roughly0.02--0.05%. This extreme does not repair the failure. The actual
native adaptive estimator was reapplied and its derivative changes retained;
the electronic intercept was re-solved from the same original G0 equation.

The RAPID-dark stage check keeps its nonlocal support explicit: direct E4 ROIs
have zero reference dark mean, but their IPC halos contain isolated nonzero
entries. Removing the nominal preIPC dark-mean assignment while retaining the
actual native dark subtraction and re-solving G0 changes the continuous native
ROI mean by at most2.46e-9SCI; all original float32 native outputs are identical.
Held Q is unchanged. This measures that specific stage sensitivity; it does not
derive a general physical dark/overflow law or its stochastic covariance.

The first failed optimization, which trusted the compact-template edge merely
because an initial rate was nonnegative, is preserved under
all_data/attempt_1_compact_edge_mask. It is not a scientific milestone or the
final fit. The subsequent fixed-domain fit has zero domain rejections and
meets both declared convergence gates: scaled projected KKT <=1e-3 and last
relative objective change <=1e-10. Kappa bound proximity is tested in scaled
coordinates, avoiding a false zero-bound classification at small positive kappa.

Remaining dependencies are explicit: the photon/electronic background
partition; actual post-full-well/overflow response; deposited-charge/segment
state for native JUMP histories; and full nonlocal G0-anchored noise, including
the SAT-cutoff boundary's raw/post-SIRS G0 cross-covariance. A fixed native-ERR
estimating metric is not a converged parameter-dependent Gaussian likelihood.
The separate measured dark covariance cannot simply be added to native ERR or
reapplied through IPC. None of these classical detector calculations supplies
the BHSM U/V bilinear source or its optical state kernel. Physical admission
and validation remain false; the 2e-10 physical guard is unchanged.

Reproduce with the existing environments, without reinstalling dependencies:

```powershell
.venv-native-history/Scripts/python.exe scripts/run_he0435_native_completion.py --label all_data
.venv-native-history/Scripts/python.exe scripts/run_he0435_native_completion.py --label training
.venv-native-history/Scripts/python.exe scripts/finalize_he0435_native_held.py
.venv-geometry/Scripts/python.exe scripts/check_he0435_native_completion_direction.py --label all_data
.venv-geometry/Scripts/python.exe scripts/diagnose_he0435_native_spatial_residuals.py --label held_target
.venv-geometry/Scripts/python.exe scripts/diagnose_he0435_native_temporal_residuals.py --label all_data
.venv-geometry/Scripts/python.exe scripts/diagnose_he0435_held_native_temporal.py
.venv-native-history/Scripts/python.exe scripts/compare_he0435_native_background_partition.py
.venv-native-history/Scripts/python.exe scripts/diagnose_he0435_native_dark_stage.py
.venv-geometry/Scripts/python.exe -m pytest tests/test_strong_lensing_native_completion.py -q
```
