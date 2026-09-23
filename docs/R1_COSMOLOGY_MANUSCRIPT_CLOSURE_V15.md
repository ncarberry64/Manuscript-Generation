# R1 Cosmology Manuscript Correction / Closure Ledger â€” V15

**Status:** manuscript consolidation freeze
**Branch:** `theory/cosmology-universe-mdpi-submission-pass`
**Starting commit:** `f6fb600662d5ab838d72a50d92a830eba5637b5c`

This ledger freezes the scientific interpretation reached after the canonical-basis
and BAO observable-class audits. It is a manuscript-status document, not a new fit
and not a new cosmological parameter choice.

## 1. Canonical state basis â€” CLOSED

All current action-native common-state statements use

\[
X_2=(q_2,\Pi_2),
\qquad
q_2=\zeta_{\rm anchor},
\qquad
\Pi_2=2a_{\rm anchor}^3\mathcal G_{S,2}\dot\zeta_{\rm anchor}.
\]

The Pantheon high-redshift state is

\[
X_2=(-2.646722099677e-01,\,+2.300532232547e-02).
\]

Legacy/raw response rows written in \((\zeta,d\zeta/dN)\) coordinates must not
be dotted directly into this canonical state.

## 2. Background geometry â€” CLOSED

The frozen R1 background remains

\[
H_0=67.4\ {\rm km\,s^{-1}\,Mpc^{-1}},
\qquad
\Omega_{k0}=-0.018,
\]

with

\[
R_H=33.153131317\ {\rm Gpc}.
\]

No background parameter was changed by this reconciliation.

## 3. Low-redshift supernova amplitude â€” CLOSED AS EMPIRICAL CALIBRATION

At \(z_\star=0.02\),

\[
\frac{\delta D_L}{D_L}
=-1.897330116627e-02\pm4.144653167389e-03,
\]

corresponding kinematically to

\[
|\Delta H_0|
=1.278800499\pm0.279349623
\ {\rm km\,s^{-1}\,Mpc^{-1}}.
\]

This is the empirical local calibration. It is not replaced by the homogeneous
high-redshift state extrapolation.

## 4. Global-state versus local response â€” SEPARATED

Propagating the high-redshift canonical state through the canonical action-native
luminosity-distance row gives at \(z=0.02\)

\[
|\Delta H_0|_{\rm global}
=0.010172329492\ {\rm km\,s^{-1}\,Mpc^{-1}},
\]

whereas the empirical local calibration is \(1.278800499\). The central-value
ratio is approximately \(125.714\).

Therefore the present homogeneous/global realization is **not** used as a
substitute for the local calibration. A local/topographic contribution must be
modeled and tested separately if the framework is to predict that amplitude.

This is a limitation of the minimal global-only realization, not a license to
retune the canonical state.

## 5. BAO observable status

### 5.1 Transverse BAO distance response â€” CLOSED

On the retained late-time branch with \(d\ln r_d=0\),

\[
\boxed{\mathcal F_{D_M/r_d}(z)=\mathcal F_{D_L}(z)}.
\]

This is a distance/light-cone response. At \(z=0.93\), the current canonical
high-redshift state gives

\[
\delta\ln(D_M/r_d)=+7.849412258190e-03.
\]

### 5.2 Historical curvature formula â€” RETAINED, REINTERPRETED

The historical relation

\[
\frac{\Delta r_{\rm BAO}}{r_{\rm BAO}}
\simeq -\frac{K\chi^2(z)}{6}
\]

is retained as a **background curvature distance-projection** expression. It is
not interpreted as a physical stretching of the 150-Mpc sound horizon.

For reference at \(z=0.93\),

\[
\delta_{\rm BAO,bg}=-1.561712887951e-03.
\]

By contrast, curvature acting intrinsically across a 150-Mpc ruler is only

\[
-\frac16\left(\frac{150\ {\rm Mpc}}{R_H}\right)^2
=-3.411788963831e-06.
\]

The preprint-sized \(10^{-3}\) effect is therefore a Gpc-baseline projection
scale, not a \(10^{-3}\) deformation of the sound horizon itself.

### 5.3 Historical \(R_{BS}\sim700\text{--}800\) â€” HISTORICAL ONLY

The old direct BAO--SN proportionality is **not action-native closed**. The modern
R1 realization has a two-component canonical state and observable-specific
kernels, so a universal coefficient cannot be asserted merely by eliminating one
scalar amplitude.

The 700--800 number remains in the historical scientific record as a
phenomenological forecast/representation. It must not be presented as the final
modern R1 prediction.

### 5.4 Radial BAO â€” OPEN

The exact radial relation is

\[
\mathcal F_{D_H/r_d}(z)=-\mathcal F_H(z)
\]

on the same late-time branch, but \(\mathcal F_H\) has not been action-native
closed. Numerical derivatives of the transverse distance response are diagnostics
only and are not promoted to a final \(H(z)\) kernel.

### 5.5 Survey projection â€” OPEN DOWNSTREAM OPERATOR

DESI footprint, pair-window, estimator, and covariance effects enter **after**
the cosmic observable kernels are fixed. They do not define the underlying
cosmic amplitude.

## 6. Other action-native observables

The current canonical matter/growth and Weyl/lensing rows remain available as
cross-observable predictions from the same canonical state. They should be
presented with their existing retrospective/prospective caveats and not
reinterpreted through the superseded mixed-basis bridge.

## 7. Superseded intermediate interpretations

The following are explicitly superseded and must not be used for manuscript
claims:

- canonical Pantheon \(X_2\) dotted into legacy/raw \((\zeta,\zeta_N)\)
  luminosity-distance rows;
- BAO window target bands derived from that mixed-basis construction;
- raw dust-transfer rows dotted directly into canonical \(X_2\);
- identifying the gauge-invariant long-mode spatial dilation itself with the
  BAO peak shift.

## 8. What is allowed to remain open at manuscript consolidation

The manuscript does **not** require these objects to be closed before the present
R1 paper can be finalized honestly:

1. exact \(\mathcal F_H(z)\) / radial \(D_H/r_d\) kernel;
2. DESI pair-window and covariance projection;
3. a completed local/topographic-to-global decomposition that predicts the
   empirical low-z amplitude from first principles;
4. independent future-data validation.

These are stated as future tests/extensions rather than silently filled with
surrogate quantities.

## 9. Manuscript policy

The paper should describe the current R1 result as a **testable action-native
observable framework**, not as proof that every historical phenomenological
coefficient survived unchanged.

No frozen cosmological parameter was retuned in reaching this status.
