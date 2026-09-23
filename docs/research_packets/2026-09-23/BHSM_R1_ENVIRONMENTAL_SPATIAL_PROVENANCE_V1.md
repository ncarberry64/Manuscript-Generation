# BHSM / R1 environmental spatial provenance audit v1

## Executive result

The coupled R1 calculation already shows that environmental density/velocity data can span both temporal topographic degrees of freedom. The remaining spatial question can be stated exactly.

Let the four environmental fields have degree-two \(S^3\) coefficient vectors collected into

\[
E_i\in\mathbb R^{4\times9},
\]

and let the two topographic temporal components across the same nine spatial coefficients be

\[
X(z)\in\mathbb R^{2\times9}.
\]

On the isotropic R1 background,

\[
X(z)=U_{XX}(z)X_i+U_{XE}(z)E_i .
\]

Therefore a one-profile reduction is **not fundamental**. It is the special condition

\[
\operatorname{rank}_{spatial}X\le1.
\]

At a fixed epoch the environmental contribution has rank at most two because \(U_{XE}\) is \(2\times4\). Generic environmental data need not be rank one.

## Multi-epoch result

Using the independently reconstructed environmental transfer matrices at \(z=1.5,1.0,0.5,0\), the stacked map has singular values

\[
13.157526737,\quad 0.0710417139623,\quad
2.39053722135e-05,\quad 2.57783336054e-08.
\]

It is full column rank in double precision, although the fourth direction is extremely weak. The matter-only stacked map \((\delta_m,v_m)\) has singular values

\[
13.1575173113,\quad 0.0710415006827
\]

and is robustly rank two.

Hence a redshift-persistent single spatial profile requires, at minimum for the matter sector, aligned \(n=2\) density and velocity coefficient vectors unless some additional physical reduction or cancellation is derived.

## Exact degree-two tensor projection

For a scalar field \(F(X)\) on the unit \(S^3\subset\mathbb R^4\), define \(X_AX_A=1\). The exact symmetric trace-free degree-two coefficient tensor is

\[
\boxed{
A_{AB}[F]
=
\frac{6}{\pi^2}
\int_{S^3}
F(X)
\left(X_AX_B-\frac14\delta_{AB}\right)d\Omega_3
}.
\]

The reconstructed harmonic is

\[
F_{n=2}(X)=A_{AB}X_AX_B,
\]

with

\[
\int_{S^3}F_{n=2}^2d\Omega_3
=
\frac{\pi^2}{6}A_{AB}A_{AB}.
\]

The pure observer-dipole sector is \(A_{44}=0\), trace-free \(S_{ij}=0\), leaving the three \(A_{4i}\) coefficients.

This gives the exact object that should be extracted from independently specified environmental density and velocity fields.

## Existing independent low-z diagnostic — not an initial-condition derivation

The retained project already contains independently constructed low-redshift environmental directions:

- 2MRS cumulative K-flux density proxy: RA=165.454 deg, Dec=-4.096 deg.
- 2MRS NN velocity-field descriptive dipole: RA=224.739 deg, Dec=-33.381 deg.

Their angular separation is

\[
\boxed{62.308^\circ}.
\]

Treating only their directions as equal-norm vectors, the two-vector matrix has singular values

\[
1.210256,\quad 0.731628,
\]

so the best rank-one approximation captures only 73.24% of the squared norm, leaving 26.76% in the second spatial direction.

The five individual 2MRS K-flux shell directions are also not rank one: their principal squared-norm fractions are approximately 55.5%, 35.9%, and 8.6%.

These facts are **diagnostic only**. They do not represent the \(z=2.1\) environmental initial state and cannot be used to choose the theory. They do show that a generic multi-profile environment is not an artificial possibility: the existing local environmental proxies are visibly multi-directional.

For posthoc context only, the velocity-field dipole is 23.86 deg from the frozen SN axis, while the K-flux density dipole is 46.28 deg away. Those comparisons do not select either input.

## Resulting interpretation

The cleanest current architecture is

\[
\boxed{
\text{specified coupled environmental state}
\rightarrow
\text{physical }n=2\text{ spatial projection}
\rightarrow
\text{frozen coupled R1 propagation}
\rightarrow
\text{observable state}
}.
\]

Amplitude, temporal phase, and orientation may all belong to the coupled state. A uniquely isolated scalar eigenline is unnecessary.

What remains is **provenance**, not another attractor ansatz:

1. specify or derive the environmental initial state;
2. project its density and velocity fields onto the complete \(n=2\) tensor basis;
3. measure the spatial rank rather than assuming one profile;
4. justify any reduction to the pure dipole or one-axis sector;
5. retain the separate microscopic BHSM-to-R1 normalization bridge.

No observational amplitude or axis was used to construct the coupled-state theory in this audit.
