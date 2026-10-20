# Common finite-body GR optical correction through second post-Minkowskian order

> Scope reconciliation: the [current BHSM scope](../STRONG_LENSING_SCOPE_20261020.md)
> supersedes earlier project framing in this retained technical record. The
> governing question is surface contributions to envelope light transfer;
> conventional effective-mass equivalence is a constructive result at the
> retained order. Reference-body offsets are not discrepancies between fully
> matched representations. The full surface-to-R1 common-source derivation
> remains unfinished. Equations, numerical evidence and review chronology
> below keep their original scope.

## Result and scope

The same finite Tolman-VII density used in the R1 comparison admits an explicit ordinary-GR deflection correction through second post-Minkowskian (2PM) order. It includes the central ray through the fluid. A matching analytic primitive supplies a consistent physical thin-plane potential and Hessian.

Add this **same** GR correction to both the conventional lens and the conventional lens plus the previously computed leading R1 scalar response. Re-solve both image systems and compute both sets of stationary arrivals with the same potential convention. This improves the common ordinary-GR optical control; it does not derive scalar-dependent 2PM terms or an exact finite-distance cosmological time-transfer function.

The independent impact coordinate \(b\) is the invariant asymptotic impact of an isolated static scattering problem. The comparison identifies it with its existing physical thin-plane coordinate. It must not be replaced by the closest areal radius \(r_0\); their distinction contributes at the very order being computed.

## 1. Exact isolated scattering integral without subtracting two nearly equal angles

Write the exact GR fluid and exterior metric as

\[
ds^2=-e^{2\nu(r)}dt^2+\frac{dr^2}{f(r)}+r^2d\Omega^2,
\qquad f(r)=1-\frac{2m(r)}r.
\]

Here \(m(r)=GM(<r)/c^2\) and the lapse is normalized at spatial infinity. The existing `TolmanVII` class solves the TOV equations inside the body and matches to Schwarzschild outside. Conservation of photon energy and angular momentum gives

\[
\frac{d\varphi}{dr}
=\frac{b e^{\nu(r)}}{r^2\sqrt{f(r)}
\sqrt{1-b^2e^{2\nu(r)}/r^2}}.
\]

Introduce the optical circumference radius

\[
\rho=r e^{-\nu(r)},\qquad
\frac{d\rho}{dr}=e^{-\nu}(1-r\nu').
\]

The turning point now has \(\rho=b\). On the weakly compact branch, \(1-r\nu'>0\), so this change of variable is monotone. With \(\rho=b/\sin\theta\), the exact deflection becomes

\[
\hat\alpha_{\rm exact}(b)
=2\int_0^{\pi/2}\left[
\frac{1}{\sqrt{f(r)}(1-r\nu'(r))}-1
\right]d\theta,
\qquad r e^{-\nu(r)}=\frac b{\sin\theta}.
\tag{1}
\]

This has neither a moving integration endpoint nor a turning-point square-root singularity. The flat subtraction is performed inside the integral. The implementation evaluates its integrand as

\[
\operatorname{expm1}\left[-\tfrac12\operatorname{log1p}(-2m/r)
-\operatorname{log1p}(-r\nu')\right].
\]

It inverts the optical-radius relation by Newton iteration and splits the quadrature at the optical image of the material surface. This provides an independent check of the analytic expansion below.

## 2. Expansion with the same prescribed finite density

Set \(\epsilon=m/R\), \(x=r/R\), and

\[
\mu(x)=\frac{m(r)}m=
\begin{cases}
\tfrac12x^3(5-3x^2),&x<1,\\
1,&x\ge1.
\end{cases}
\]

Use \(u(x)=\mu(x)/x\), and expand

\[
\nu(r)=\epsilon n_1(x)+O(\epsilon^2),\qquad
p_g(r)=\frac{m^2}{R^4}\big[p_0(x)+O(\epsilon)\big],
\]

where \(p_g=Gp_{\rm physical}/c^4\). The Newtonian TOV quantities, with the correct exterior lapse normalization, are

\[
n_1(x)=
\begin{cases}
-\frac{15}{8}+\frac54x^2-\frac38x^4,&x<1,\\
-1/x,&x\ge1,
\end{cases}
\]

\[
p_0(x)=\frac{15}{32\pi}(2-5x^2+4x^4-x^6),\quad x<1,
\qquad p_0=0\quad x\ge1.
\]

The exact radial Einstein/TOV row gives

\[
r\nu'=\epsilon u+
\epsilon^2\big[2u^2+4\pi x^2p_0\big]+O(\epsilon^3).
\]

At fixed optical radius \(s=\rho/R\), the areal coordinate is
\(x=s[1+\epsilon n_1(s)]+O(\epsilon^2)\). Expanding the integrand of (1) therefore gives

\[
\frac{1}{\sqrt f(1-r\nu')}-1
=2\epsilon u(s)+\epsilon^2 F_2(s)+O(\epsilon^3),
\]

\[
F_2(x)=\frac{11}{2}u^2+4\pi x^2p_0+2x n_1u'.
\tag{2}
\]

The final term is the optical-to-areal coordinate conversion. Omitting it would fail the exterior Schwarzschild check. The pressure term is also required by the common GR solution.

For the specified density this simplifies to

\[
F_2(x)=
\begin{cases}
-15x^2+60x^4-\frac{105}{2}x^6+15x^8,&x<1,\\
\frac{15}{2x^2},&x\ge1.
\end{cases}
\tag{3}
\]

The \(O(\epsilon)\) term in (1) reproduces the existing projected-density lens. The additional deflection is

\[
\delta\hat\alpha_{\rm GR,2PM}(b)
=2\epsilon^2 I(y),\qquad
I(y)=\int_0^{\pi/2}F_2(y/\sin\theta)d\theta,
\quad y=b/R.
\tag{4}
\]

## 3. Closed deflection, derivative, and potential

For \(0<y<1\), let \(t=\sqrt{1-y^2}\). The exact coefficient of the 2PM term is

\[
I(y)=\frac{15}{4y^2}\big[\arcsin y-yt\big]
+y\left[\frac{15}{2}t-4t^5-\frac{48}{7}t^7\right].
\tag{5}
\]

Outside the body,

\[
I(y)=\frac{15\pi}{8y^2},\qquad y\ge1,
\]

so (4) recovers exactly
\(\delta\hat\alpha_{\rm GR,2PM}=15\pi(m/b)^2/4\), the standard invariant-impact Schwarzschild coefficient. See Iyer and Petters, Eq. (23): <https://arxiv.org/pdf/gr-qc/0611086>.

The derivative needed for the lens Hessian is

\[
I'(y)=-\frac{15}{2y^3}(\arcsin y-yt)
+15t+20t^3+24t^5-\frac{384}{7}t^7,
\quad y<1,
\tag{6}
\]

and \(I'=-15\pi/(4y^3)\) outside. The primitive normalized by \(J(1)=0\) is

\[
J(y)=\frac{15\pi}{8}
-\frac{15}{4}\left[\frac{\arcsin y}{y}+t\right]
-\frac52t^3+\frac47t^7+\frac{16}{21}t^9,
\quad y<1,
\tag{7}
\]

\[
J(y)=\frac{15\pi}{8}(1-1/y),\quad y\ge1.
\]

Consequently,

\[
\delta\hat\psi_{\rm GR,2PM}(b)=2\epsilon^2R J(b/R),
\qquad
\frac{d\delta\hat\psi}{db}=\delta\hat\alpha,
\qquad
\frac{d\delta\hat\alpha}{db}=\frac{2\epsilon^2}{R}I'(b/R).
\tag{8}
\]

The physical projected potential has units of length, the deflection is dimensionless, and its derivative has units of inverse length. The common potential constant at \(b=R\) drops out of delay differences.

The central limits are regular:

\[
I(y)=-\frac67y+31y^3-\frac{1065}{28}y^5+O(y^7),
\qquad J(0)=\frac{15\pi}{8}-\frac{26}{3}.
\]

Thus the additional central slope is \(-12\epsilon^2/(7R)\). A small **negative** central 2PM deflection correction is expected for this fixed areal density profile and asymptotic impact convention; the full GR deflection remains positive. At the surface, \(I\) and \(I'\) match their exterior values. The correction potential, deflection, and lens Hessian are continuous; higher derivatives need not be smooth because of the finite density cutoff.

The implementation uses a convergent small-argument series for \(\arcsin y-y\sqrt{1-y^2}\) and a stable rearrangement of (7) near the surface. A separate AI-assisted computational review differentiated (7) and confirmed both central and surface limits.

## 4. Focused independent check against exact GR geodesics

The physical compactness is

\[
\epsilon=1.595138613862614\times10^{-6}.
\]

At the three unchanged reference impacts:

| Impact (kpc) | Additional radial deflection, 2PM (rad) | Exact static GR minus 1PM minus 2PM (rad) |
|---|---:|---:|
| 3.82838407985 | \(+1.84072751081140\times10^{-11}\) | \(+8.33312033332\times10^{-17}\) |
| 0.376488911045, inside fluid | \(-2.41589331705463\times10^{-13}\) | \(-2.11955918329\times10^{-18}\) |
| 3.01282777544 | \(+2.97215800368347\times10^{-11}\) | \(+1.70972779711\times10^{-16}\) |

These are radial correction functions; the vector direction is supplied by the image's position relative to the lens center. The remainder includes ordinary-GR orders beyond 2PM and numerical error. Its smallness is a check of the ordinary-GR expansion for this body, not a bound on missing scalar or cosmological effects.

For an independent coefficient test of the interior formula, the exact TOV/scattering calculation was repeated at the same central impact ratio \(y=0.12549630368161785\) and compactness values \(5\times10^{-4}\), \(2.5\times10^{-4}\), and \(1.25\times10^{-4}\). Extrapolating \((\alpha_{\rm exact}-\alpha_{\rm 1PM})/\epsilon^2\) quadratically to zero compactness gives

\[
-0.09494692325461365,
\]

compared with the analytic \(2I(y)=-0.09494692324937426\). The absolute difference is \(5.24\times10^{-12}\). Both exterior coefficient checks agree exactly at displayed double precision. The full compactness sequence, quadrature estimates, physical rows, and source hashes are in `gr_optical_control_results.json`.

## 5. Callable interface and consistent comparison

The standalone module is `code/strong_lensing_gr2.py`. It exports:

- `TolmanVII2PMCorrection(base)`: `extra_alpha(b)`, `extra_alpha_prime(b)`, `extra_potential(b)`, and `extra_application(x)`.
- `TolmanVII2PMLens(base)`: the original 1PM body plus that 2PM correction, with a consistent `application(x)` returning potential, gradient, and Hessian. Mass/radius response derivatives are included.
- `exact_static_deflection(body,b)`: the independent flat-subtracted exact static-scattering check using a solved `TolmanVII` body.

The intended shared composition is:

```python
base_1pm = TolmanVIILens()
common_gr = TolmanVII2PMLens(base_1pm)
gr_solution = solve_images(geometry, common_gr, source)
scalar_on = R1ScreenedLens(base=common_gr, coefficients=unchanged_coefficients)
r1_solution = solve_axis_images(
    geometry, scalar_on, source, gr_solution['images_rad'])
```

Compute both stationary arrivals with these same applications; then compare the new R1-minus-GR difference with its old 1PM counterpart. Adding only an exterior angle estimate, omitting the central ray, or changing deflection without its potential would not provide this consistency check.

This is a controlled improvement in common ordinary-GR optics. The leading R1 response is deliberately held fixed; mixed scalar/GR post-Minkowskian terms, exact finite-distance optical matching, environmental response, and full modified-equilibrium changes remain outside this calculation.
