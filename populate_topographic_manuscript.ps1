$ErrorActionPreference = 'Stop'

if (-not (Test-Path '.git')) { throw 'Run this script from the root of Manuscript-Generation.' }

New-Item -ItemType Directory -Force -Path 'manuscript' | Out-Null

$content = @'
\documentclass[11pt]{article}
\usepackage[margin=1in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern}
\usepackage{amsmath,amssymb,mathtools,booktabs}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage[numbers,sort&compress]{natbib}
\setlength{\emergencystretch}{2em}

\title{Boundary-Mode Topographic Cosmology:\\
A Reduced-Response Bridge from Berger--Hopf Geometry to Late-Time Structure Growth}
\author{Norman P. Carberry\\Independent Researcher}
\date{Working draft --- September 2026}

\newcommand{\Mpl}{M_{\rm Pl}}
\newcommand{\OmK}{\Omega_k}
\newcommand{\RH}{R_H}
\newcommand{\Kphys}{\mathcal K_{\rm phys}}

\begin{document}
\maketitle

\begin{abstract}
We develop a late-time cosmological bridge between a scalar topographic effective field theory on compact spatial sections \(S^3(R_H)\) and the constraint-reduced response architecture of the Berger--Hopf program. The construction keeps a covariant cubic-Galileon background with second-order time evolution and introduces the higher-spatial-gradient response as the low-energy Schur/Feshbach reduction of a stable physical complement rather than as a fundamental higher-time-derivative interaction. On \(S^3\), scalar harmonics obey \(-\nabla^2Y_n=n(n+2)Y_n/R_H^2\). The reduced response can therefore soften the lowest nonconstant \(n=1\) branch while stiffening \(n\ge2\), providing a mechanism for a large low-redshift dipolar susceptibility without a comparable quadrupole or large modification of ordinary redshift-space-distortion modes. We define a reference branch R1, specify an auditable preregistration protocol, and derive prospective growth and lensing targets. The low-redshift supernova amplitude is used only as an observational normalization of the \(n=1\) response; the curvature forecast and structure-growth predictions are kept logically separate. The framework is explicitly falsifiable through spatial curvature, harmonic hierarchy, fixed-axis anisotropy, growth, and lensing tests.
\end{abstract}

\input{sections/01_introduction}
\input{sections/02_hyperspherical_parent}
\input{sections/03_reduced_response}
\input{sections/04_scalar_metric_bridge}
\input{sections/05_harmonic_selection}
\input{sections/06_observables}
\input{sections/07_predictions_falsification}
\input{sections/08_numerical_protocol}
\input{sections/09_results_R1}
\input{sections/10_discussion}
\input{sections/11_conclusions}

\bibliographystyle{unsrtnat}
\bibliography{references}
\end{document}
'@

Set-Content -Path 'manuscript\main.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Introduction}

Late-time cosmological tensions and directional residuals motivate tests that are capable of separating local structure, survey geometry, and genuinely new large-scale dynamics. The scalar-topographic program introduced a compact hyperspherical parent \(S^3(R_H)\) and a long-wavelength scalar response whose largest observational effect occurs at low redshift \citep{Carberry2026Topographic}. A later covariant formulation placed the topographic scalar in the cubic-Galileon sector of Horndeski theory, preserving second-order time evolution while allowing nontrivial late-time scalar--metric braiding.

The present work asks a narrower question: can the long-wavelength topographic response be represented as the large-scale branch of the same constraint-reduced response logic used in the Berger--Hopf program, without assuming that the full microscopic BHSM construction has already been compactified into cosmology? We answer this at the level of a controlled bridge model.

The central mathematical chain is
\[
S[\Phi]
\longrightarrow
\widehat{\mathcal H}_{\rm phys}
\longrightarrow
\mathcal K_{\rm topo}
\longrightarrow
\lambda_{\rm topo}
\longrightarrow
T_n
\longrightarrow
\delta g_{\mu\nu}
\longrightarrow
\{D_L,H,f\sigma_8,\Phi+\Psi\}.
\]

The distinction between established input and bridge extension is essential. The \(S^3\) harmonic spectrum and cubic-Galileon background are treated as the cosmological starting point. The Schur/Feshbach reduction is an exact mathematical operation once a physical complement is specified. The identification of its low-energy coefficients with the topographic response is a bridge assumption to be tested, not a declaration that the complete microscopic BHSM action has already been reduced to four-dimensional cosmology.

The paper makes four main contributions. First, it standardizes the \(S^3\) harmonic convention and shows how the observed sky dipole arises from the four-dimensional \(n=1\) eigenspace rather than by identifying a radial mode with an angular dipole. Second, it derives a reduced kernel in which a positive higher-spatial-gradient term can emerge after integrating out a stable complement. Third, it identifies a parameter interval in which the \(n=1\) response is selectively softened while all \(n\ge2\) harmonics are stiffened. Fourth, it gives an executable preregistration protocol for curvature, supernova anisotropy, growth, and lensing tests.

The claim boundary is therefore:
\begin{quote}
This work does not assume that the scalar-topographic cosmological field has already been derived from the complete microscopic Berger--Hopf action. Exact reductions are distinguished from cross-scale identifications and exploratory reference-branch assumptions. Prospective observables are frozen before comparison with designated independent data.
\end{quote}
'@

Set-Content -Path 'manuscript\sections\01_introduction.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Hyperspherical Parent and Harmonic Convention}

We use a closed FLRW geometry
\[
ds^2=-dt^2+a^2(t)R_H^2\left[d\chi^2+\sin^2\chi\,d\Omega_2^2\right],
\]
with \(a_0=1\). Thus \(R_H\) is the present-day comoving curvature radius and the physical curvature radius is \(aR_H\).

Scalar harmonics satisfy
\[
-\nabla_\gamma^2Y_n=\frac{n(n+2)}{R_H^2}Y_n,
\qquad
n=0,1,2,\ldots,
\]
with degeneracy \((n+1)^2\). The physical spatial eigenvalue entering the time-dependent perturbation equations is
\[
p_n(a)=\frac{n(n+2)}{a^2R_H^2}.
\]

The first nonconstant eigenspace is therefore
\[
p_1=\frac{3}{a^2R_H^2}.
\]
This convention is used throughout. Older scalar-topographic notes employed more than one normalization for the lowest mode; no numerical forecast is carried forward here without first converting it to the convention above.

\subsection{Exact light-cone form of the \(n=1\) mode}

Embed the unit three-sphere in \(\mathbb R^4\). The \(n=1\) harmonics are the four coordinate functions restricted to \(S^3\). A general real \(n=1\) mode can therefore be written
\[
T_1(\chi,\hat{\mathbf n})
=
A\left[
\cos\alpha\cos\chi
+
\sin\alpha\sin\chi\,
(\hat{\mathbf n}\cdot\hat{\mathbf p})
\right].
\]
At fixed redshift this contains an exact sky monopole plus an exact sky dipole and no higher angular multipoles. The dipole part is
\[
T_{1,{\rm dip}}(z,\hat{\mathbf n})
=
A\sin\alpha\,
\sin\!\left[\frac{\chi(z)}{R_H}\right]
(\hat{\mathbf n}\cdot\hat{\mathbf p}),
\]
where \(\chi(z)\) denotes the comoving radial distance in dimensional units when written inside the sine as \(\chi/R_H\).

This removes an ambiguity present in simplified radial descriptions: the zonal term \(\cos\chi\) alone is isotropic for an observer at its pole. A sky dipole is supplied by a generic orientation of the full \(n=1\) eigenspace relative to the observer.
'@

Set-Content -Path 'manuscript\sections\02_hyperspherical_parent.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Constraint-Reduced Response and the Topographic Kernel}

Let \(\widehat{\mathcal H}\) be the physical Hessian after constraints and gauge directions have been removed. For a selected topographic branch \(P\) and physical complement \(Q\), the exact Schur/Feshbach reduced operator is
\[
\mathcal K
=
P\left[
\widehat{\mathcal H}
-
\widehat{\mathcal H}Q
(Q\widehat{\mathcal H}Q)^{-1}
Q\widehat{\mathcal H}
\right]P.
\]

For one scalar harmonic, write the quadratic block as
\[
\widehat{\mathcal H}_n
=
\begin{pmatrix}
A(p) & C(p)\\
C^\dagger(p) & D(p)
\end{pmatrix}.
\]
The exact reduced stiffness is
\[
K_{\rm red}(p)=A(p)-C(p)D^{-1}(p)C^\dagger(p).
\]

A minimal long-wavelength realization is
\[
A(p)=Z_0p,\qquad
D(p)=M_c^2+c_cp,\qquad
C(p)=g\sqrt p.
\]
Then
\[
K_{\rm red}(p)
=
Z_0p-\frac{g^2p}{M_c^2+c_cp}.
\]
For \(c_cp\ll M_c^2\),
\[
K_{\rm red}(p)
=
\left(Z_0-\frac{g^2}{M_c^2}\right)p
+
\frac{g^2c_c}{M_c^4}p^2
+O(p^3).
\]

Thus a positive \(p^2\) stiffness can arise from a coupled second-order parent system after a stable complement is integrated out. In this interpretation the old static operator
\[
-\nabla^2+B\nabla^4
\]
is a low-energy spatial derivative expansion, not evidence for a fundamental fourth-order time equation.

Define
\[
\chi\equiv\frac{g^2}{M_c^2},
\qquad
L_c^2\equiv\frac{c_c}{M_c^2}.
\]
After canonical normalization of the leading spatial term, the reduced kernel may be written schematically as
\[
K_{\rm topo}(p)
=
m_T^2+
(c_s^2-\chi)p
+
\chi L_c^2p^2+\cdots.
\]
The signs of the temporal kinetic term and of the full physical gradient matrix must be checked independently; a positive spatial \(p^2\) coefficient by itself is not a complete ghost analysis.
'@

Set-Content -Path 'manuscript\sections\03_reduced_response.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Cubic-Galileon Scalar--Metric Bridge}

The cosmological parent action is taken to be
\[
S=
\int d^4x\sqrt{-g}
\left[
\frac{M_{\rm Pl}^2}{2}R
-\frac12\nabla_\mu T\nabla^\mu T
+\frac{1}{\Lambda^3}\Box T(\nabla T)^2
-V(T)
+\mathcal L_m
\right].
\]
Using
\[
X\equiv-\frac12\nabla_\mu T\nabla^\mu T,
\]
this corresponds to
\[
G_2=X-V(T),\qquad
G_3=\frac{2X}{\Lambda^3},\qquad
G_4=\frac{M_{\rm Pl}^2}{2},\qquad
G_5=0
\]
in the convention \(\mathcal L_3=-G_3\Box T\).

The model belongs to the kinetic-gravity-braiding subclass of Horndeski theory \citep{BelliniSawicki2014}. Since \(G_4\) is constant and \(G_5=0\),
\[
M_*^2=M_{\rm Pl}^2,\qquad
\alpha_M=0,\qquad
\alpha_T=0.
\]
The scalar--metric mixing is controlled by the braiding function
\[
\alpha_B
=
\frac{2\dot{\bar T}^{\,3}}
{HM_{\rm Pl}^2\Lambda^3},
\]
and
\[
\alpha_K
=
\frac{\dot{\bar T}^{\,2}}
{H^2M_{\rm Pl}^2}
+
6\alpha_B.
\]
The scalar no-ghost combination is
\[
D_{\rm kin}
=
\alpha_K+\frac32\alpha_B^2.
\]

For numerical work we softly break shift symmetry with
\[
V(T)=V_0-\mu^3T.
\]
The scalar energy density and pressure are
\[
\rho_T
=
\frac12\dot T^2+V(T)
+\frac{6H\dot T^3}{\Lambda^3},
\]
\[
P_T
=
\frac12\dot T^2-V(T)
-\frac{2\dot T^2\ddot T}{\Lambda^3}.
\]
The shift current
\[
J=\dot T+\frac{6H\dot T^2}{\Lambda^3}
\]
obeys
\[
\dot J+3HJ=\mu^3.
\]

The closed-FLRW background equation is
\[
3M_{\rm Pl}^2
\left(
H^2+\frac{1}{a^2R_H^2}
\right)
=
\rho_m+\rho_r+\rho_T.
\]

The bridge does not require the bare Horndeski sound speed to become negative. Instead, the BHSM-inspired reduced response can lower the effective long-wavelength stiffness while the positive \(p^2\) term stabilizes shorter wavelengths.
'@

Set-Content -Path 'manuscript\sections\04_scalar_metric_bridge.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Selective Harmonic Response}

Let the complement correlation length scale with the parent radius,
\[
L_c(a)=aR_H\,\ell,
\]
where \(\ell\) is dimensionless. Since
\[
p_n=\frac{n(n+2)}{a^2R_H^2},
\]
we have
\[
L_c^2p_n=\ell^2n(n+2).
\]

In the effectively massless topographic sector, define
\[
K_n=p_nC_n
\]
with
\[
\boxed{
C_n(a)
=
c_s^2(a)
-
\chi(a)\left[1-\ell^2n(n+2)\right].
}
\]

If
\[
\boxed{
\frac18<\ell^2<\frac13,
}
\]
then
\[
1-3\ell^2>0
\]
for \(n=1\), whereas
\[
1-\ell^2n(n+2)<0
\]
for every \(n\ge2\). Hence increasing \(\chi\) softens the \(n=1\) response but stiffens all higher harmonics:
\[
C_1<c_s^2,\qquad
C_n>c_s^2\quad(n\ge2).
\]

This produces a stable near-critical alternative to a true gradient instability. If
\[
0<C_1\ll1,
\]
then a sourced mode
\[
T_1\simeq\frac{S_1}{p_1C_1}
\]
can be strongly enhanced while all \(n\ge2\) responses remain regular.

For the quadrupole,
\[
\frac{T_2}{T_1}
=
\frac{S_2}{S_1}
\frac{3C_1}{8C_2}.
\]
Therefore the response per unit source satisfies
\[
\left|
\frac{T_2/S_2}{T_1/S_1}
\right|
=
\frac{3C_1}{8C_2}
<
\frac38
\]
on the stable selective-filter branch, with much stronger suppression as \(C_1/C_2\to0\).

This is the primary harmonic null prediction of the model.
'@

Set-Content -Path 'manuscript\sections\05_harmonic_selection.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Observable Projection}

The observable distance perturbation is not identified directly with the scalar amplitude. At linear order it is a light-cone projection of metric, velocity, and scalar perturbations:
\[
\frac{\delta D_L}{\bar D_L}
=
\int d\chi\,
\left[
W_\Phi\Phi+
W_\Psi\Psi+
W_vv_\parallel+\cdots
\right].
\]
For an \(n=1\) topographic source, the angular dependence factors into
\[
\hat{\mathbf n}\cdot\hat{\mathbf p},
\]
so that
\[
\frac{\delta D_L}{\bar D_L}
=
\epsilon_D(z)\,
\hat{\mathbf n}\cdot\hat{\mathbf p}.
\]

For supernovae,
\[
\mu=5\log_{10}D_L+\mathrm{const},
\]
and therefore a small magnitude dipole \(A_\mu\) corresponds to
\[
\frac{\delta D_L}{D_L}
=
\frac{\ln10}{5}A_\mu\cos\theta.
\]

For structure growth we define
\[
k^2\Psi
=
-4\pi Ga^2\mu(k,a)\rho_m\delta_m,
\qquad
\Phi=\gamma(k,a)\Psi,
\]
and
\[
\Sigma(k,a)
=
\frac{\mu(k,a)}2[1+\gamma(k,a)].
\]
The matter growth factor obeys
\[
D''+
\left(2+\frac{H'}H\right)D'
-\frac32\Omega_m(a)\mu(k,a)D=0,
\]
with prime denoting \(d/d\ln a\). The observable is
\[
f\sigma_8(k,z)
=
\frac{D'}{D}D\,\sigma_{8,\rm ini}
=
D'(k,z)\sigma_{8,\rm ini}.
\]

The quasi-static limit is used only on scales inside the scalar sound horizon \citep{SawickiBellini2015}. Horizon-scale \(n=1\) observables are treated through the full long-wavelength response rather than by extrapolating the quasi-static formulas.
'@

Set-Content -Path 'manuscript\sections\06_observables.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Predictions and Falsification Criteria}

\subsection{Low-redshift calibration}

The existing Pantheon+SH0ES tomographic analysis provides the observational normalization
\[
A_\mu^\star=-0.0412\pm0.009
\]
in the bin
\[
0.01<z<0.03,
\]
with best-fit direction
\[
({\rm RA},{\rm Dec})
=
(211.48^\circ,-12.81^\circ).
\]
This corresponds to
\[
\frac{\delta D_L}{D_L}
=
(-0.01895\pm0.00414)\cos\theta,
\]
and at sufficiently low redshift,
\[
\frac{\delta H}{H}
=
(+0.01895\pm0.00414)\cos\theta.
\]

This datum normalizes only the observed \(n=1\) response,
\[
A_\mu
=
\frac{\mathcal K_{\rm SN}S_1}{p_1C_1},
\]
and does not independently determine \(S_1\), \(\mathcal K_{\rm SN}\), or \(C_1\).

\subsection{Curvature}

The pre-existing prospective curvature forecast is retained as
\[
\boxed{\Omega_k=-0.018\pm0.004}.
\]
At \(H_0=67.4\ {\rm km\,s^{-1}\,Mpc^{-1}}\),
\[
R_H
=
\frac{c}{H_0\sqrt{|\Omega_k|}}
\simeq33.15\ {\rm Gpc}.
\]
The sign prediction is \(\Omega_k<0\) in the usual FLRW curvature convention.

\subsection{Harmonic hierarchy}

The bridge predicts selective \(n=1\) softening with no analogous \(n\ge2\) enhancement. A future source-normalized \(n=2\) susceptibility comparable to or larger than the \(n=1\) susceptibility would reject the minimal harmonic-filter branch.

\subsection{Fixed-axis test}

A prospective independent supernova sample should first test the frozen axis
\[
(211.48^\circ,-12.81^\circ)
\]
rather than scanning freely over the sky. A free-axis search may be reported secondarily.

\subsection{High-\(n\) recovery}

At large \(n\),
\[
C_n\simeq\chi\ell^2n^2,
\]
so the extra scalar susceptibility falls as \(n^{-2}\). The model therefore predicts a strong separation between the horizon-scale dipole and conventional clustering modes.

\subsection{Failure conditions}

The minimal bridge is disfavored if one or more of the following are established with adequate precision:
\begin{enumerate}
\item independent curvature measurements exclude the frozen negative-\(\Omega_k\) forecast;
\item the source-normalized \(n=2\) response is comparable to or exceeds \(n=1\);
\item an independent low-redshift sample excludes the frozen dipole axis while finding a stable incompatible physical axis;
\item a coherent dipole of comparable amplitude persists into a redshift regime in which the model predicts loss of enhanced susceptibility;
\item ordinary RSD modes show a large scale-independent modification incompatible with the high-\(n\) suppression;
\item lensing requires a metric response incompatible with the same scalar propagator that controls matter clustering.
\end{enumerate}
'@

Set-Content -Path 'manuscript\sections\07_predictions_falsification.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Numerical Implementation and Preregistration Protocol}

The numerical pipeline is separated into
\[
\boxed{
\text{background}
\rightarrow
\text{physical scalar response}
\rightarrow
\text{observables}
\rightarrow
\text{comparison}.
}
\]
Prediction generation and data comparison are implemented in separate scripts.

\subsection{Background system}

Let \(v=\dot T\). Expanding the scalar current equation gives
\[
\left(1+\frac{12Hv}{\Lambda^3}\right)\dot v
+
\frac{6v^2}{\Lambda^3}\dot H
=
\mu^3-3H\left(v+\frac{6Hv^2}{\Lambda^3}\right).
\]
The Raychaudhuri equation gives
\[
\frac{2v^2}{\Lambda^3}\dot v
-
2M_{\rm Pl}^2\dot H
=
\rho_m+\frac43\rho_r+v^2+\frac{6Hv^3}{\Lambda^3}
-\frac{2M_{\rm Pl}^2}{a^2R_H^2}.
\]
At each integration step these equations form a \(2\times2\) linear system for \((\dot v,\dot H)\). The remaining equations are
\[
\dot T=v,\qquad
\dot a=aH.
\]

The Friedmann equation is retained as an independent residual check:
\[
3M_{\rm Pl}^2
\left(H^2+\frac{1}{a^2R_H^2}\right)
-
(\rho_m+\rho_r+\rho_T)
=0.
\]

\subsection{Dimensionless variables}

We integrate in \(N=\ln a\) with
\[
E=\frac{H}{H_0},
\quad
\phi=\frac{T}{M_{\rm Pl}},
\quad
\lambda=\frac{\Lambda^3}{M_{\rm Pl}H_0^2},
\quad
q=\frac{\mu^3}{M_{\rm Pl}H_0^2}.
\]

The exploratory reference branch R1 uses
\[
\Omega_{m0}=0.31,\quad
\Omega_{r0}=9.0\times10^{-5},\quad
\Omega_{k0}=-0.018,\quad
\lambda=1,\quad q=3.
\]
The constant \(V_0\) is determined only by the shooting condition \(E(a=1)=1\).

\subsection{Initial conditions}

The numerical integration begins at
\[
z_i=100.
\]
The decaying homogeneous current solution is removed using the matter-era particular solution
\[
J_i\simeq\frac{2\mu^3}{9H_i}.
\]
The initial scalar velocity is the positive small-velocity root of
\[
v_i+\frac{6H_iv_i^2}{\Lambda^3}=J_i.
\]

\subsection{Growth grid}

Before prospective comparison, predictions are generated on the fixed grid
\[
k=\{0.02,0.05,0.10,0.15,0.20\}\ {\rm Mpc}^{-1}
\]
and
\[
z=\{0.5,0.8,1.0,1.5,2.1\}.
\]

A deterministic JSON manifest containing the action version, parameters, numerical tolerances, calibration values, and predicted observables is hashed with SHA-256. Any later scientific or numerical revision creates a new manifest version rather than overwriting the previous freeze.
'@

Set-Content -Path 'manuscript\sections\08_numerical_protocol.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Results: Reference Branch R1}

R1 is an exploratory reference branch rather than a uniquely derived BHSM cosmology. Its purpose is to demonstrate that the bridge can simultaneously realize late-time scalar activation, a selectively enhanced \(n=1\) response, and near-standard high-\(n\) structure growth.

For
\[
\Omega_{m0}=0.31,\quad
\Omega_{r0}=9.0\times10^{-5},\quad
\Omega_{k0}=-0.018,\quad
\lambda=1,\quad q=3,
\]
the shooting condition \(E(1)=1\) gives
\[
\boxed{
\frac{V_0}{M_{\rm Pl}^2H_0^2}=2.39861.
}
\]

The scalar braiding and no-ghost combination evolve as shown in Table~\ref{tab:r1background}.

\begin{table}[ht]
\centering
\begin{tabular}{rrrr}
\toprule
\(z\) & \(\alpha_B\) & \(D_{\rm kin}\) & \(\Omega_T\)\\
\midrule
10  & \(2.07\times10^{-7}\) & \(1.64\times10^{-6}\) & 0.00194\\
3   & \(8.60\times10^{-5}\) & \(6.80\times10^{-4}\) & 0.03915\\
2.1 & \(3.72\times10^{-4}\) & 0.00294 & 0.08060\\
1.5 & 0.00122 & 0.00959 & 0.14321\\
1.0 & 0.00380 & 0.02988 & 0.24551\\
0.8 & 0.00622 & 0.04889 & 0.30787\\
0.5 & 0.01346 & 0.10536 & 0.43200\\
0.3 & 0.02271 & 0.17740 & 0.53538\\
0   & 0.04940 & 0.38485 & 0.70791\\
\bottomrule
\end{tabular}
\caption{Reference-branch background evolution.}
\label{tab:r1background}
\end{table}

The scalar equation of state evolves from nearly vacuum-like behavior at high redshift to
\[
w_T(z=0)\simeq-0.906,
\]
and the background acceleration transition occurs at
\[
z_{\rm acc}\simeq0.683.
\]

Relative to a curved \(\Lambda\)CDM model with identical \(H_0,\Omega_{m0},\Omega_{r0},\Omega_{k0}\), the R1 expansion rate differs by only order one percent over the late-time interval. The maximum difference in the tabulated range is about \(1.7\%\) near \(z\sim0.5\).

\subsection{Prospective growth target}

In the minimal high-\(n\) realization, the direct topographic fifth-force term is negligible on ordinary RSD modes and the dominant growth shift is produced by the altered background history. With identical matter-era normalization,
\[
R_{f\sigma_8}(z)
=
\frac{f\sigma_{8,\rm R1}}
{f\sigma_{8,\Lambda{\rm CDM}}}
\]
is

\begin{table}[ht]
\centering
\begin{tabular}{rrr}
\toprule
\(z\) & \(R_{f\sigma_8}\) & suppression\\
\midrule
2.1 & 0.99346 & 0.65\%\\
1.5 & 0.98869 & 1.13\%\\
1.0 & 0.98164 & 1.84\%\\
0.8 & 0.97795 & 2.20\%\\
0.5 & 0.97271 & 2.73\%\\
0.3 & 0.97168 & 2.83\%\\
0.0 & 0.98482 & 1.52\%\\
\bottomrule
\end{tabular}
\caption{R1 growth forecast relative to matched curved \(\Lambda\)CDM. These values are reference-branch predictions and become formally preregistered only when reproduced by the canonical code and written to the hashed manifest.}
\label{tab:r1growth}
\end{table}

The central phenomenological result is therefore not a large generic modification of gravity. It is the combination
\[
\boxed{
\text{enhanced horizon-scale }n=1\text{ susceptibility}
+
\text{rapid high-}n\text{ recovery}
+
\text{few-percent late-time growth shift}.
}
\]
'@

Set-Content -Path 'manuscript\sections\09_results_R1.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Discussion}

The bridge developed here differs from several common late-time alternatives in a specific way. A local-flow explanation can generate a low-redshift dipole but does not by itself predict a compact \(S^3\) harmonic hierarchy, a frozen curvature sign, and a linked high-\(n\) response. Ordinary quintessence has a canonical smooth gradient sector and does not naturally single out one discrete global spatial harmonic. A pure cubic Galileon supplies kinetic braiding but, without an additional scale-selective completion, a negative bare gradient coefficient would not isolate only the deep infrared. The reduced-response term introduced here instead permits the bare Horndeski branch to remain kinetically healthy while the physical \(n=1\) susceptibility is selectively enhanced.

The compact spectrum is central. With \(R_H\simeq33\) Gpc, ordinary RSD wavenumbers correspond to harmonics hundreds or thousands of levels above \(n=1\). This makes it possible for a large global dipolar response to coexist with nearly standard galaxy clustering. The model is therefore not equivalent to a scale-independent fifth force.

The same mechanism also clarifies the role of BHSM. The manuscript does not require that a cosmological scalar be identified directly with one microscopic Berger--Hopf field. What is imported is the reduced-response logic: a selected physical branch is evaluated only after the coupled complement has adjusted. In a derivative expansion, this produces a negative long-wave Schur correction together with a positive higher-spatial-gradient stiffness. That is the cross-scale structural claim made here.

Several limitations remain. The source law \(S_1\) and the light-cone projection kernel must be derived in a complete perturbation treatment before the Pantheon normalization can be inverted into a unique microscopic susceptibility. The minimal no-slip closure used for some exploratory lensing estimates is not assumed to be a final theorem. The R1 values \(q=3\) and \(\lambda=1\) are reference-branch choices, not yet action-derived constants. The curvature forecast should be tested with directional-aware likelihoods and not retuned in response to a discrepant isotropic estimator.

These limitations are scientifically useful because each corresponds to a concrete next calculation rather than an unrestricted fitting freedom.
'@

Set-Content -Path 'manuscript\sections\10_discussion.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript\sections' | Out-Null

$content = @'
\section{Conclusions}

We have constructed a reduced-response bridge between hyperspherical scalar-topographic cosmology and the Berger--Hopf response architecture. The model retains a covariant cubic-Galileon background and interprets the higher-spatial-gradient topographic stiffness as the low-energy effect of integrating out a stable physical complement.

On \(S^3(R_H)\), the exact scalar spectrum
\[
p_n=\frac{n(n+2)}{a^2R_H^2}
\]
allows a simple harmonic filter. For
\[
\frac18<\ell^2<\frac13,
\]
the same response that softens \(n=1\) stiffens every \(n\ge2\) mode. This yields a large dipolar susceptibility without requiring a comparable quadrupole and without forcing a large modification of ordinary structure-growth modes.

The low-redshift Pantheon+SH0ES dipole is used only to normalize the observed \(n=1\) response. The pre-existing curvature forecast
\[
\Omega_k=-0.018\pm0.004
\]
is kept separate. An exploratory reference branch R1 gives a late-time braiding component that is negligible at high redshift and predicts only a few-percent suppression of \(f\sigma_8\) relative to matched curved \(\Lambda\)CDM over the main survey range.

The strongest prediction is therefore a pattern rather than a single number:
\[
\boxed{
\text{strong }n=1\text{ response}
\;\oplus\;
\text{suppressed }n\ge2\text{ response}
\;\oplus\;
\text{near-GR high-}k\text{ growth}.
}
\]
The accompanying preregistration protocol is designed so that curvature, harmonic hierarchy, fixed-axis anisotropy, growth, and lensing can falsify the bridge without post-comparison retuning.
'@

Set-Content -Path 'manuscript\sections\11_conclusions.tex' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'manuscript' | Out-Null

$content = @'
@article{Carberry2026Topographic,
  author = {Norman P. Carberry},
  title = {A Hyperspherical Scalar--Topographic Framework for Late-Time Cosmological Anomalies},
  year = {2026},
  note = {Preprints.org, posted 20 January 2026, doi:10.20944/preprints202601.1427.v1}
}

@article{BelliniSawicki2014,
  author = {Emilio Bellini and Ignacy Sawicki},
  title = {Maximal freedom at minimum cost: linear large-scale structure in general modifications of gravity},
  journal = {JCAP},
  year = {2014},
  volume = {07},
  pages = {050},
  eprint = {1404.3713},
  archivePrefix = {arXiv}
}

@article{SawickiBellini2015,
  author = {Ignacy Sawicki and Emilio Bellini},
  title = {Limits of quasistatic approximation in modified-gravity cosmologies},
  journal = {Phys. Rev. D},
  year = {2015},
  volume = {92},
  pages = {084061},
  eprint = {1503.06831},
  archivePrefix = {arXiv}
}

@article{Scolnic2022,
  author = {D. Scolnic and others},
  title = {The Pantheon+ Analysis: The Full Dataset and Light-Curve Release},
  journal = {Astrophys. J.},
  year = {2022},
  volume = {938},
  pages = {113},
  eprint = {2112.03863},
  archivePrefix = {arXiv}
}

@article{Brout2022,
  author = {D. Brout and others},
  title = {The Pantheon+ Analysis: Cosmological Constraints},
  journal = {Astrophys. J.},
  year = {2022},
  volume = {938},
  pages = {110},
  eprint = {2202.04077},
  archivePrefix = {arXiv}
}
'@

Set-Content -Path 'manuscript\references.bib' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'docs' | Out-Null

$content = @'
# Claim status — BHSM/topographic cosmology manuscript

## Established source inputs
- Closed hyperspherical parent `S^3(R_H)`.
- Cubic-Galileon/Horndeski cosmological action used in the later Topographic Dark Energy formulation.
- Standardized scalar spectrum: `-ΔY_n = n(n+2)Y_n/R_H^2`.
- Existing low-z Pantheon+SH0ES tomographic calibration: `A_mu = -0.0412 ± 0.009` in `0.01 < z < 0.03`, axis `(RA,Dec)=(211.48°, -12.81°)`.
- Existing prospective curvature forecast: `Omega_k = -0.018 ± 0.004`.

## Exact bridge mathematics
- Physical Schur/Feshbach reduction: `K_red = A - C D^{-1} C†`.
- Long-wave expansion of a stable complement can generate a negative `p` correction and positive `p^2` stiffness.
- Exact `n=1` light-cone decomposition on `S^3`: generic mode = monopole + sky dipole.
- For `1/8 < ell^2 < 1/3`, the filter softens `n=1` and stiffens every `n>=2`.

## Exploratory bridge assumptions
- `L_c(a)=a R_H ell`.
- Minimal source-response form `T_n=S_n/[p_n C_n]`.
- Reference Branch R1: `Omega_m0=0.31`, `Omega_r0=9e-5`, `Omega_k0=-0.018`, `lambda=1`, `q=3`.
- Any no-slip closure `gamma=1` is R1-specific unless separately derived.

## New derived reference-branch quantities
- R1 shooting solution: `V0/(M_Pl^2 H0^2) ≈ 2.39861`.
- R1 growth ratios are generated by `code/growth.py`; they are not framework-wide predictions until the manifest is frozen.

## Not claimed
- Full BHSM microscopic-to-cosmological compactification theorem.
- Unique derivation of `q`, `lambda`, `ell`, or the source normalization.
- Unique microscopic origin of the observed supernova dipole axis.
- Proof that the entire Hubble tension is explained by the measured ~1.9% low-z distance dipole.
'@

Set-Content -Path 'docs\claim_status.md' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'docs' | Out-Null

$content = @'
# Derivation ledger

## D1 — S3 spectrum
Input:
`ds_3^2 = a^2 R_H^2[dchi^2 + sin^2(chi)dOmega_2^2]`.

Result:
`p_n(a)=n(n+2)/(a^2 R_H^2)`.

Status: exact geometry.

## D2 — n=1 observer-sky projection
Using the R4 embedding of S3:
`T_1=A[u_4 cos chi + u_vec·n_hat sin chi]`.

Result: fixed-redshift sky contains monopole + dipole only.

Status: exact geometry.

## D3 — Schur-generated higher-spatial stiffness
Parent block:
`A=Z0 p`, `D=M_c^2+c_c p`, `C=g sqrt(p)`.

Exact:
`K_red=Z0 p - g^2 p/(M_c^2+c_c p)`.

Long-wave:
`K_red=(Z0-g^2/M_c^2)p + (g^2 c_c/M_c^4)p^2 + ...`.

Status: exact algebra given the declared block.

## D4 — harmonic filter
Set `L_c=a R_H ell`.

Then:
`C_n=c_s^2-chi[1-ell^2 n(n+2)]`.

If `1/8 < ell^2 < 1/3`:
- n=1 softened;
- n>=2 stiffened.

Status: exact bridge result.

## D5 — R1 background
Action:
cubic Galileon + linear soft-breaking potential.

Numerical branch:
`Omega_m0=0.31`, `Omega_r0=9e-5`, `Omega_k0=-0.018`, `lambda=1`, `q=3`.

Shooting:
`E(a=1)=1`.

Status: exploratory numerical reference branch.

## D6 — R1 growth
Growth uses identical early normalization and matched curved LCDM reference.

Output:
`R_fsigma8 = D'_R1 / D'_LCDM`.

Status: prospective reference-branch forecast; freeze only after canonical regeneration and hashing.
'@

Set-Content -Path 'docs\derivation_ledger.md' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'code' | Out-Null

$content = @'
"""Reference Branch R1 background integrator.

Units: H0 = M_Pl = 1.
The script solves the cubic-Galileon + linear-potential background on closed FLRW,
shoots V0 so E(a=1)=1, and prints a deterministic JSON summary.
"""
from __future__ import annotations
import json
import math
from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OMEGA_M0 = 0.31
OMEGA_R0 = 9.0e-5
OMEGA_K0 = -0.018
KCURV = -OMEGA_K0
LAMBDA = 1.0
Q = 3.0
Z_INITIAL = 100.0
A_INITIAL = 1.0 / (1.0 + Z_INITIAL)
N_INITIAL = math.log(A_INITIAL)

@dataclass(frozen=True)
class Params:
    omega_m0: float = OMEGA_M0
    omega_r0: float = OMEGA_R0
    omega_k0: float = OMEGA_K0
    lam: float = LAMBDA
    q: float = Q
    z_initial: float = Z_INITIAL

def _initial_state(v0: float):
    phi = 0.0
    E = math.sqrt(
        OMEGA_M0*A_INITIAL**-3
        + OMEGA_R0*A_INITIAL**-4
        + OMEGA_K0*A_INITIAL**-2
        + v0/3.0
    )
    for _ in range(100):
        J = 2.0*Q/(9.0*E)
        A = 6.0*E/LAMBDA
        v = (-1.0 + math.sqrt(1.0 + 4.0*A*J))/(2.0*A)
        rho_t = 0.5*v*v + v0 + 6.0*E*v**3/LAMBDA
        E_new = math.sqrt(
            OMEGA_M0*A_INITIAL**-3
            + OMEGA_R0*A_INITIAL**-4
            + rho_t/3.0
            - KCURV*A_INITIAL**-2
        )
        if abs(E_new-E) < 1e-15*E:
            E = E_new
            break
        E = E_new
    return np.array([phi, v, E], dtype=float)

def rhs(N: float, y: np.ndarray) -> np.ndarray:
    phi, v, E = y
    a = math.exp(N)
    rho_m = 3.0*OMEGA_M0*a**-3
    rho_r = 3.0*OMEGA_R0*a**-4
    J = v + 6.0*E*v*v/LAMBDA

    A11 = 1.0 + 12.0*E*v/LAMBDA
    A12 = 6.0*v*v/LAMBDA
    b1 = Q - 3.0*E*J

    A21 = 2.0*v*v/LAMBDA
    A22 = -2.0
    b2 = (
        rho_m + (4.0/3.0)*rho_r + v*v
        + 6.0*E*v**3/LAMBDA
        - 2.0*KCURV*a**-2
    )
    det = A11*A22 - A12*A21
    if abs(det) < 1e-14:
        raise RuntimeError("Background derivative matrix became singular.")

    dot_v = (b1*A22 - A12*b2)/det
    dot_E = (A11*b2 - b1*A21)/det
    return np.array([v/E, dot_v/E, dot_E/E], dtype=float)

def integrate(v0: float):
    return solve_ivp(
        rhs, (N_INITIAL, 0.0), _initial_state(v0),
        method="DOP853", rtol=2e-13, atol=1e-15,
        dense_output=True, max_step=0.005,
    )

def shoot_v0() -> float:
    def target(v0):
        return float(integrate(v0).y[2, -1] - 1.0)
    return brentq(target, 2.3, 2.5, xtol=1e-13)

def diagnostics(sol, v0: float, z: float):
    N = math.log(1.0/(1.0+z))
    phi, v, E = [float(x) for x in sol.sol(N)]
    a = math.exp(N)
    dphi_dN, dv_dN, dE_dN = rhs(N, np.array([phi, v, E]))
    dot_v = dv_dN*E
    alpha_b = 2.0*v**3/(E*LAMBDA)
    x = v/E
    alpha_k = x*x + 6.0*alpha_b
    dkin = alpha_k + 1.5*alpha_b**2
    V = v0 - Q*phi
    rho_t = 0.5*v*v + V + 6.0*E*v**3/LAMBDA
    p_t = 0.5*v*v - V - 2.0*v*v*dot_v/LAMBDA
    omega_t = rho_t/(3.0*E*E)
    w_t = p_t/rho_t

    rho_m = 3.0*OMEGA_M0*a**-3
    rho_r = 3.0*OMEGA_R0*a**-4
    friedmann = (E*E + KCURV*a**-2) - (rho_m+rho_r+rho_t)/3.0
    residual = abs(friedmann)/(3.0*E*E)

    return {
        "z": z, "E": E, "alpha_B": alpha_b, "D_kin": dkin,
        "Omega_T": omega_t, "w_T": w_t, "friedmann_residual": residual,
    }

def main():
    v0 = shoot_v0()
    sol = integrate(v0)
    zs = [10.0, 3.0, 2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]
    rows = [diagnostics(sol, v0, z) for z in zs]
    payload = {
        "model": "R1",
        "params": asdict(Params()),
        "V0_over_Mpl2H02": v0,
        "rows": rows,
    }
    print(json.dumps(payload, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
'@

Set-Content -Path 'code\background.py' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'code' | Out-Null

$content = @'
"""Growth forecast for Reference Branch R1 relative to matched curved LCDM."""
from __future__ import annotations
import json, math
import numpy as np
from scipy.integrate import solve_ivp
import background

V0 = background.shoot_v0()
R1 = background.integrate(V0)
OM = background.OMEGA_M0
OR = background.OMEGA_R0
OK = background.OMEGA_K0
ODE = 1.0 - OM - OR - OK
NI = background.N_INITIAL
AI = background.A_INITIAL

def r1_dlnH_dN(N):
    y = R1.sol(N)
    return float(background.rhs(N, y)[2] / y[2])

def e_lcdm(N):
    a = math.exp(N)
    return math.sqrt(OM*a**-3 + OR*a**-4 + OK*a**-2 + ODE)

def dlnE_lcdm(N):
    a = math.exp(N)
    e2 = OM*a**-3 + OR*a**-4 + OK*a**-2 + ODE
    de2 = -3*OM*a**-3 - 4*OR*a**-4 - 2*OK*a**-2
    return 0.5*de2/e2

def growth_rhs_r1(N, y):
    D, Dp = y
    E = float(R1.sol(N)[2])
    om = OM*math.exp(-3*N)/E**2
    return [Dp, -(2.0+r1_dlnH_dN(N))*Dp + 1.5*om*D]

def growth_rhs_lcdm(N, y):
    D, Dp = y
    E = e_lcdm(N)
    om = OM*math.exp(-3*N)/E**2
    return [Dp, -(2.0+dlnE_lcdm(N))*Dp + 1.5*om*D]

def integrate_growth(fun):
    return solve_ivp(
        fun, (NI, 0.0), [AI, AI],
        method="DOP853", rtol=1e-12, atol=1e-14,
        dense_output=True, max_step=0.01,
    )

def main():
    gr1 = integrate_growth(growth_rhs_r1)
    gl = integrate_growth(growth_rhs_lcdm)
    zs = [2.1, 1.5, 1.0, 0.8, 0.5, 0.3, 0.0]
    rows = []
    for z in zs:
        N = math.log(1.0/(1.0+z))
        Dr, Dpr = gr1.sol(N)
        Dl, Dpl = gl.sol(N)
        rows.append({
            "z": z,
            "R_fsigma8": float(Dpr/Dpl),
            "D_ratio": float(Dr/Dl),
        })
    print(json.dumps({"model":"R1","rows":rows}, indent=2, sort_keys=True))

if __name__ == "__main__":
    main()
'@

Set-Content -Path 'code\growth.py' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'code' | Out-Null

$content = @'
"""Minimal harmonic-filter utilities for the BHSM/topographic bridge."""
from __future__ import annotations
import math

def nu(n: int) -> int:
    if n < 0:
        raise ValueError("n must be non-negative")
    return n*(n+2)

def p_n(n: int, a: float, R_H: float) -> float:
    return nu(n)/(a*a*R_H*R_H)

def C_n(n: int, c_s2: float, chi: float, ell: float) -> float:
    return c_s2 - chi*(1.0 - ell*ell*nu(n))

def selective_filter_condition(ell: float) -> bool:
    e2 = ell*ell
    return (1.0/8.0) < e2 < (1.0/3.0)

def response_per_source(n: int, a: float, R_H: float, c_s2: float, chi: float, ell: float) -> float:
    p = p_n(n, a, R_H)
    c = C_n(n, c_s2, chi, ell)
    if p == 0.0:
        raise ZeroDivisionError("n=0 has no inverse-gradient response in this massless form")
    return 1.0/(p*c)

def quadrupole_to_dipole_response_ratio(c1: float, c2: float) -> float:
    return 3.0*c1/(8.0*c2)
'@

Set-Content -Path 'code\perturbations.py' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'code' | Out-Null

$content = @'
"""Generate deterministic preregistration manifest and SHA-256 digest."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import background, growth
from scipy.integrate import solve_ivp
import math

MANIFEST = Path(__file__).resolve().parents[1] / "preregistration" / "prediction_manifest.json"

def canonical_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")

def build_manifest():
    v0 = background.shoot_v0()
    r1 = background.integrate(v0)

    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)
    zgrid = [0.5, 0.8, 1.0, 1.5, 2.1]
    grow = {}
    for z in zgrid:
        N = math.log(1/(1+z))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        grow[str(z)] = float(dpr/dpl)

    return {
        "protocol": "BHSM_TOPO_COSMOLOGY_PREREG_V1",
        "status": "REFERENCE_BRANCH_R1_NOT_YET_PROMOTED_TO_FRAMEWORK_WIDE_FREEZE",
        "geometry": {
            "Omega_k": -0.018,
            "Omega_k_sigma": 0.004,
            "H0_km_s_Mpc": 67.4,
            "laplacian": "-Delta Y_n = n(n+2)/R_H^2 Y_n",
        },
        "background_R1": {
            "Omega_m0": background.OMEGA_M0,
            "Omega_r0": background.OMEGA_R0,
            "Omega_k0": background.OMEGA_K0,
            "lambda": background.LAMBDA,
            "q": background.Q,
            "z_initial": background.Z_INITIAL,
            "V0_over_Mpl2H02": v0,
        },
        "topographic_calibration": {
            "redshift_range": [0.01, 0.03],
            "A_mu": -0.041156,
            "A_mu_sigma": 0.009,
            "axis_RA_deg": 211.48,
            "axis_Dec_deg": -12.81,
        },
        "harmonic_filter": {
            "condition": "1/8 < ell^2 < 1/3",
            "softened_mode": 1,
            "higher_modes": "stiffened",
        },
        "growth": {
            "k_Mpc_inverse": [0.02,0.05,0.10,0.15,0.20],
            "redshift": zgrid,
            "R_fsigma8": grow,
            "comparison": "matched curved LCDM",
            "parameter_refit_after_freeze": False,
        },
    }

def main():
    manifest = build_manifest()
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    digest = hashlib.sha256(canonical_bytes(manifest)).hexdigest()
    print(f"manifest={MANIFEST}")
    print(f"sha256={digest}")

if __name__ == "__main__":
    main()
'@

Set-Content -Path 'code\preregister.py' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'preregistration' | Out-Null

$content = @'
# Preregistration

`prediction_manifest.json` is generated by `code/preregister.py`.

The current protocol is a **reference-branch freeze candidate**, not yet a framework-wide final freeze. Before comparison with a designated independent dataset:

1. run `python code/background.py`;
2. run `python code/growth.py`;
3. run `python code/preregister.py`;
4. record the printed SHA-256 digest in the manuscript/release notes;
5. do not alter upstream parameters after inspecting the comparison data.

Any scientific or numerical correction creates a new protocol version and preserves the old manifest in Git history.
'@

Set-Content -Path 'preregistration\README.md' -Value $content -Encoding UTF8

New-Item -ItemType Directory -Force -Path 'preregistration' | Out-Null

$content = @'
{
  "protocol": "BHSM_TOPO_COSMOLOGY_PREREG_V1",
  "status": "GENERATE_WITH_code/preregister.py_BEFORE_PROSPECTIVE_COMPARISON"
}
'@

Set-Content -Path 'preregistration\prediction_manifest.json' -Value $content -Encoding UTF8

$content = @'
numpy>=1.26
scipy>=1.11
'@

Set-Content -Path 'requirements.txt' -Value $content -Encoding UTF8

Write-Host 'Populated manuscript and reference implementation.'

python -c "import numpy, scipy" 2>$null
if ($LASTEXITCODE -ne 0) {
    python -m pip install -r requirements.txt
}

python code\background.py

python code\growth.py

python code\preregister.py

git add .

git commit -m "Populate topographic cosmology manuscript and R1 reference implementation"

git push

git status
