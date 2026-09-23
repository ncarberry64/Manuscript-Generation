# BHSM / R1 coupled environmental-state audit

## Result

I independently reconstructed the frozen R1 six-state action-native system and changed coordinates at the anchor \(z=2.1\) from the canonical state

\[
(\zeta,D_m,D_r,p_\zeta,p_{D_m},p_{D_r})
\]

to the physically interpretable coordinates

\[
(q_2,\Pi_2,\delta_m,\delta_r,v_m,v_r),
\qquad
q_2=\zeta,\quad
\Pi_2=2a^3G_{S,2}\dot\zeta .
\]

The coordinate inversion reproduces the repository's original two-column topographic seed to \(3.1\times10^{-13}\). An independent integration reproduces the stored full canonical propagator at \(z=1.5\) to \(3.1\times10^{-9}\).

With the initial topographic pair held at zero and the four environmental coordinates varied independently, the later projected topographic pair \((\zeta,\dot\zeta)\) receives a nonzero \(2\times4\) response.

## Temporal rank

| z | s1 | s2 | condition s1/s2 | det matter-only [delta_m,v_m] |
|---:|---:|---:|---:|---:|
| 1.50 | 9.65325193 | 0.00168033957 | 5744.82 | 0.0162206969 |
| 1.00 | 7.21390025 | 0.00988894461 | 729.491 | 0.0713376389 |
| 0.80 | 5.80181931 | 0.0153657497 | 377.581 | 0.0891490049 |
| 0.50 | 4.13861138 | 0.0228942773 | 180.771 | 0.0947501431 |
| 0.30 | 3.5300629 | 0.0235622417 | 149.819 | 0.0831758119 |
| 0.10 | 3.30044344 | 0.0191582562 | 172.273 | 0.0632303811 |
| 0.02 | 3.28211299 | 0.0165162724 | 198.72 | 0.0542079304 |
| 0.00 | 3.28191222 | 0.0158257248 | 207.378 | 0.0519383028 |

At every reported epoch after the anchor, the environmental-to-topographic transfer has **rank 2**. The two matter inputs \((\delta_m,v_m)\) alone already have nonzero determinant, so the two-dimensional temporal topographic response does not require radiation to become full rank.

This is a rank statement, not an amplitude prior. The singular-value ratio depends on the normalization chosen for the four environmental coordinates.

## Spatial reduction

The full scalar \(n=2\) harmonic sector on \(S^3\) has dimension 9. On an isotropic R1 background, each harmonic component obeys the same temporal system. Writing the environmental initial coefficients as a \(4\times9\) matrix \(E_i\) and the topographic coefficients as a \(2\times9\) matrix \(X\),

\[
X(z)=U_{XX}(z)X_i+U_{XE}(z)E_i .
\]

Therefore:

- the environment can carry and transmit an \(n=2\) spatial orientation;
- isotropic R1 dynamics do **not** generate a preferred axis from nothing;
- the traditional one-profile reduction \(X=x\,e^T\) requires \(\mathrm{rank}_{spatial}(X)\le1\);
- generic multi-pattern environmental data can instead give spatial rank 2 because \(\mathrm{rank}(U_{XE})=2\).

So the precise remaining spatial requirement is not “find one isolated eigenmode.” It is:

\[
oxed{	ext{specify/derive the environmental }n=2	ext{ spatial state and justify the reduction used for observables}.}
\]

If the environmental \(n=2\) density/velocity coefficients are all aligned with one spatial profile, the one-profile picture is preserved automatically. If they are not aligned, the full coupled realization contains more than one spatial \(n=2\) profile and the observable treatment must retain that fact.

## Claim boundary

This calculation begins at the existing R1 anchor \(z=2.1\). It therefore **does not derive the observed/inferred anchor \(X_2\) from earlier cosmological initial conditions**. It establishes that, once the coupled initial state is specified at the anchor, environmental density and velocity data can control both temporal components of the later topographic response.

The remaining predictive inputs are:

1. physical/environmental initial conditions;
2. their full \(n=2\) spatial coefficients;
3. a justified spatial reduction or population rule;
4. the microscopic BHSM-to-R1 normalization bridge.

No observational parameter was used to choose the environmental state in this audit.
