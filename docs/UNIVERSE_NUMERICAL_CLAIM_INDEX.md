# Numerical claim/source index

Inventory for review, not an assertion of new numerical validation. Status and provenance families are in UNIVERSE_CLAIM_PROVENANCE.md.

## manuscript/main.tex

```text
1: \documentclass[11pt]{article}
2: \usepackage[margin=1in]{geometry}
3: \usepackage[T1]{fontenc}
9: \setlength{\emergencystretch}{2em}
13: \date{September 2026}
24: We develop a conditional cosmological response chain on a closed $S^3$
27: $n=1$ sector; the first ordinary physical scalar carrier is $n=2$.
30: $\sin(2\psi)\hat{\mathbf n}\cdot\hat{\mathbf p}$.
39: each at $p<0.01$ under a preregistered survey design.
41: for 3,521 supernova rows and adds independent foreground information.
42: With $\Delta\chi^2=\chi^2_{\rm frozen}-\chi^2_{\rm null}$, conditional held-out
43: Pantheon+ source/environment adjustment gives $-0.883$, whereas DES transfer
44: gives $+1.321$. Fixing an external peculiar-velocity mean at coefficient one
45: gives $-4.909$ for Pantheon+, but incomplete joint uncertainty propagation
56: \input{sections/01_introduction}
57: \input{sections/01a_geometry_before_fields}
58: \input{sections/01c_reproducibility_and_ai_workflow}
59: \input{sections/01b_claim_status_and_falsifiability}
60: \input{sections/02_hyperspherical_parent}
61: \input{sections/03_reduced_response}
62: \input{sections/04aa_bhsm_native_n2_transport}
63: \input{sections/04aj_r1_n2_reduced_propagator}
64: \input{sections/04ak_r1_n2_growth_lensing_transfer}
65: \input{sections/04al_r1_n2_luminosity_distance_kernel}
66: \input{sections/04ao_bhsm_optical_transfer}
67: \input{sections/04as_bhsm_optical_low_rank_covariance}
68: \input{sections/04at_effective_representation_audit}
69: \input{sections/09_results_R1}
70: \input{sections/09d_sn_sightline_tomography}
71: \input{sections/09c_negative_results}
72: \input{sections/07_predictions_falsification}
73: \input{sections/08b_make_or_break_prediction}
74: \input{sections/08a_prospective_test_matrix}
75: \input{sections/10_discussion}
76: \input{sections/09_geometry_first_conclusion}
77: \input{sections/09b_universe_submission_back_matter}
80: \input{sections/04_scalar_metric_bridge}
81: \input{sections/04a_closed_horndeski_n2}
82: \input{sections/04b_matter_sourced_n2}
83: \input{sections/04c_bianchi_dictionary_n2}
84: \input{sections/05_harmonic_selection}
85: \input{sections/06_observables}
86: \input{sections/04am_bhsm_los_topographic_transfer}
87: \input{sections/04an_2mrs_foreground_transfer}
88: \input{sections/04ap_bhsm_optical_constitutive}
89: \input{sections/04aq_bhsm_optical_hessian_projection}
90: \input{sections/04ar_bhsm_optical_current_gate}
91: \input{sections/08_numerical_protocol}
92: \input{sections/09a_reference_figures}
```

## manuscript/sections/01_introduction.tex

```text
15: Our contribution is their organization into an explicit closed-$S^3$
21: The regular one-scalar closed-FLRW quotient removes the exceptional $n=1$
22: scalar sector. The first ordinary scalar carrier is $n=2$ in our convention.
23: One selected normalized harmonic carries $X_2=(q_2,\Pi_2)^T$; the entire
24: $n=2$ eigenspace is not two dimensional. Its profile and axis must be fixed
29: $\sin(2\psi)\hat{\mathbf n}\cdot\hat{\mathbf p}$ without inserting a new
40: with separate prospective $p<0.01$ rejection gates.
43: Preprints.org study \citep{Carberry2026Topographic}. The physical carrier,
```

## manuscript/sections/01a_geometry_before_fields.tex

```text
19: $P_{\rm phys}D^2\mathcal L_{\rm KKT}P_{\rm phys}$.
21: $\mathcal K_{\rm phys}=A-BK^{-1}B^\dagger$.
25: \mathcal A_g h+\mathcal B_{g2}X_2+J_g&=0,\\
26: \mathcal R_{g2}&=-\mathcal A_g^{-1}\mathcal B_{g2},\\
27: h_2&=h_{\rm src}+\mathcal R_{g2}X_2,\qquad
28: h_{\rm src}=-\mathcal A_g^{-1}J_g.
```

## manuscript/sections/01c_reproducibility_and_ai_workflow.tex

```text

```

## manuscript/sections/01b_claim_status_and_falsifiability.tex

```text
5: \begin{tabular}{@{}p{0.22\linewidth}p{0.73\linewidth}@{}}
11: and Jacobi/Sachs identities within linear geometric optics.\\[0.5em]
14: channel count.\\[0.5em]
16: $n=2$ profile; R1 coefficients, angular axis, SN amplitude and
18: undetermined.\\[0.5em]
21: joint uncertainty and are sensitivity only.\\[0.5em]
23: 2MRS number-density, K-band luminosity and seam-charge foreground transfer
24: candidates.\\[0.5em]
```

## manuscript/sections/02_hyperspherical_parent.tex

```text
5: ds^2=-dt^2+a^2(t)R_H^2\left[d\psi^2+\sin^2\psi\,d\Omega_2^2\right],
7: with \(a_0=1\). Define
9: D_C(z)=c\int_0^z\frac{dz'}{H(z')},
15: -\nabla_\gamma^2Y_n=\frac{n(n+2)}{R_H^2}Y_n,
16: \qquad n=0,1,2,\ldots,
20: p_n(a)=\frac{n(n+2)}{a^2R_H^2}.
22: A common closed-universe convention uses \(k^2=N^2-1\), hence
24: \boxed{N=n+1}.
28: For our \(n=1\), \(k^2=3\) and \(N=2\). The degree-one scalar harmonics on the unit three-sphere obey
30: D_iD_jQ^{(1)}=-\gamma_{ij}Q^{(1)}.
34: \left(D_iD_j+\frac{k^2}{3}\gamma_{ij}\right)Q^{(1)}=0.
36: The scalar-derived trace-free metric harmonic is absent. For the regular generally covariant one-scalar sector, the remaining scalar gauge transformations and lapse/shift constraints leave no ordinary propagating scalar degree of freedom in this exceptional sector. Modern gauge-invariant closed-universe analyses likewise take the propagating scalar spectrum to begin at alternative-convention \(N\ge3\) \citep{KieferVardanyan2022,SpecognaEtAl2026}.
38: Accordingly, the minimal bridge no longer uses our \(n=1\) as the physical scalar--metric carrier. An additional boundary degree of freedom could evade this conclusion, but it is not required for the construction below.
43: \boxed{n_\star=2},\qquad \boxed{N_\star=3},
44: \qquad p_2=\frac{8}{a^2R_H^2}.
46: Let the unit \(S^3\subset\mathbb R^4\) embedding coordinates be
48: X^4=\cos\psi,\qquad X^i=\sin\psi\,\hat n^i.
52: T_2=A_{AB}X^AX^B,
53: \qquad A_{AB}=A_{BA},\qquad \delta^{AB}A_{AB}=0.
55: Writing \(c=\cos\psi\), \(s=\sin\psi\), \(a_i=A_{4i}\), and
57: A_{ij}=S_{ij}-\frac{A_{44}}{3}\delta_{ij},\qquad S_i{}^i=0,
62: T_2=A_{44}\left(c^2-\frac{s^2}{3}\right)
63: +2cs\,\mathbf a\!\cdot\!\hat{\mathbf n}
64: +s^2S_{ij}\hat n^i\hat n^j .
69: A_{44}=0,\qquad S_{ij}=0,\qquad a_i\ne0
74: T_{2,\mathrm{dip}}=A\sin(2\psi)\,(\hat{\mathbf n}\!\cdot\!\hat{\mathbf p}).
```

## manuscript/sections/03_reduced_response.tex

```text
11: (Q\widehat{\mathcal H}Q)^{-1}
27: K_{\rm red}(p)=A(p)-C(p)D^{-1}(p)C^\dagger(p).
32: A(p)=Z_0p,\qquad
33: D(p)=M_c^2+c_cp,\qquad
40: Z_0p-\frac{g^2p}{M_c^2+c_cp}.
42: For \(c_cp\ll M_c^2\),
46: \left(Z_0-\frac{g^2}{M_c^2}\right)p
48: \frac{g^2c_c}{M_c^4}p^2
49: +O(p^3).
52: Thus a positive \(p^2\) stiffness can arise from a coupled second-order parent system after a stable complement is integrated out. In this interpretation the old static operator
54: -\nabla^2+B\nabla^4
60: \chi\equiv\frac{g^2}{M_c^2},
62: L_c^2\equiv\frac{c_c}{M_c^2}.
68: m_T^2+
69: (c_s^2-\chi)p
71: \chi L_c^2p^2+\cdots.
73: The signs of the temporal kinetic term and of the full physical gradient matrix must be checked independently; a positive spatial \(p^2\) coefficient by itself is not a complete ghost analysis.
```

## manuscript/sections/04aa_bhsm_native_n2_transport.tex

```text
31: D^2\mathcal L_{\rm KKT}(z_\star)
35: If $P_2$ denotes the projector onto the first admissible physical
39: \mathcal K_2
41: P_2
45: \widehat{\mathcal H}Q_2
47: Q_2\widehat{\mathcal H}Q_2
48: \right)^{-1}
49: Q_2\widehat{\mathcal H}
51: P_2 ,
52: \label{eq:bhsm-K2}
55: Here $Q_2=1-P_2$ on the already gauge-reduced tangent space.
59: \paragraph{Proposition 1 (action-owned $n=2$ envelope transport).}
64:       $\mathcal E_{\rm cos}\simeq S^3(R_H)$;
68: \item the homogeneous $n=0$ component is absorbed into the background and
69:       the exceptional scalar $n=1$ sector is removed by the regular
71: \item the reduced $n=2$ kinetic form is nondegenerate and positive on the
75:       retained $n=2$ sector.
79: $n=2$ envelope mode.  It does not require the introduction of an additional
85: Let $e_2(x)$ be a normalized element of
86: $P_2T_{\rm phys}$ and write the retained perturbation locally as
88: \delta\Phi_{\rm topo}(t,x)=q_2(t)e_2(x).
93: S^{(2)}_2
95: \frac12
96: \int dt\,a^3
98: \mathcal G_2(t)\,
99: (D_t q_2)^2
101: \mathcal M_2^2(t)\,q_2^2
103: 2\,\mathcal J_2(t)\,q_2
105: \label{eq:bhsm-q2-action}
107: The coefficients are the $e_2$ matrix elements of the reduced operator
108: \eqref{eq:bhsm-K2}; explicitly,
109: $\mathcal G_2$ is the kinetic quadratic form,
110: $\mathcal M_2^2$ the reduced spatial/restoring form, and
111: $\mathcal J_2$ the projection of the action-owned source onto $e_2$.
114: \Pi_2
116: a^3\mathcal G_2 D_tq_2 ,
119: \eqref{eq:bhsm-q2-action} yields
121: D_t\!\left(a^3\mathcal G_2D_tq_2\right)
123: a^3\mathcal M_2^2 q_2
125: a^3\mathcal J_2 .
126: \label{eq:bhsm-q2-eom}
128: Thus the dynamical phase-space pair $(q_2,\Pi_2)$ is the reduced description
135: Sec.~\ref{sec:closed-horndeski-n2} is used here as a late-time effective
138: convention used in that section, the exact physical-$n=2$ kinetic coefficient
141: \mathcal G_{S,2}
143: M_{\rm Pl}^2
145: 5D_{\rm kin}
149: 10\left(1-\alpha_B^{\rm BS}/2\right)^2
153: For $D_{\rm kin}>0$ and nonsingular braiding denominator this coefficient is
155: effective kinetic matching condition for assumption 4 of Proposition~1 on the
159: treated as a coordinate representative of the retained response of $q_2$.
163: \subsection{Exact observer dipole inside the physical $n=2$ sector}
165: A degree-two harmonic on $S^3$ may be written
167: T_2=A_{AB}X^AX^B,\qquad A_A{}^A=0 .
171: T_2
173: A_{44}\left(c^2-\frac{s^2}{3}\right)
175: 2cs\,\mathbf a\!\cdot\!\hat{\mathbf n}
177: s^2 S_{ij}\hat n^i\hat n^j .
179: The pure $A_{4i}$ subspace is therefore
182: T_{2,\rm dip}
184: A\sin(2\psi)
189: \label{eq:bhsm-n2-exact-dipole}
201: At the present stage, the literal curvature-aware $00$, $0i$, spatial-trace,
210: Eq.~\eqref{eq:bhsm-q2-eom}: propagation is generated by the action-reduced
211: physical $n=2$ operator.  A future exact identification of the curved EFT
212: transfer map with $\mathcal K_2$ must reproduce this reduced phase space
217: Proposition~1 states the consequences of assumed mode ownership and the form of its reduced
221: $\mathcal M_2^2$, $\mathcal J_2$, or the observable transfer coefficient.
225: Consequently, the exact $n=2$ dipole geometry, the physical-harmonic
```

## manuscript/sections/04aj_r1_n2_reduced_propagator.tex

```text
1: \section{First numerical R1 physical-$n=2$ reduced propagator}
2: \label{sec:r1-n2-reduced-propagator}
4: The conditional reference response is evaluated at the reduced gravity--scalar level. For the physical $n=2$ harmonic, $-(\mathcal D^2+3K)=5K$.
10: +\left(3+\frac{d\ln H}{dN}+\frac{d\ln\mathcal G_{S,2}}{dN}\right)\zeta_{,N}
11: +\frac{5K}{a^2H^2}\frac{\mathcal F_{S,2}}{\mathcal G_{S,2}}\zeta=0.
15: The fundamental matrix is normalized to the identity at $z=2.1$ in the state $(\zeta,d\zeta/dN)^T$.
18: U(0,2.1)=
20: -0.186499277 & 0.019664322\\
21: 0.013019584 & -0.005651284
25: Across $0\le z\le2.1$ the numerical scan gives $\min\mathcal G_{S,2}=1.471470\times10^{-3}$ and $\min\mathcal F_{S,2}=5.037267\times10^{-1}$.
26: The maximum absolute residual of the exact lapse/shift constraints over both fundamental solutions is $2.720\times10^{-14}$, and the Liouville determinant identity is reproduced with relative residual $4.171\times10^{-14}$.
29: This is the first numerical phase-space transfer kernel for the physical $n=2$ R1 mode.  It uses the exact reduced closed-Horndeski $n=2$ kinetic coefficient and constraints.
30: The current $\mathcal F_{S,2}$ implementation is explicitly the gravity--scalar expression; the complete matter+radiation scalar stability/response system is not yet included.  Consequently this result is not yet the directional $f\sigma_8$, lensing, or luminosity-distance observable kernel.
```

## manuscript/sections/04ak_r1_n2_growth_lensing_transfer.tex

```text
1: \section{Direct $n=2$ matter-growth and Weyl transfer}
2: \label{sec:r1-n2-growth-lensing-transfer}
4: The reduced physical $n=2$ propagator can be coupled directly to minimally coupled matter through covariant energy--momentum conservation, without invoking the curvature-defective Newtonian-gauge DAE.
9: \dot\delta_m+3\frac{d}{dt}(Hv_m)&=-\left(3\dot\zeta+\frac{k_2^2}{a^2}\chi\right)+\frac{k_2^2}{a^2}v_m,
11: with $k_2^2/a^2=8K/a^2$ for the physical manuscript $n=2$ harmonic.
19: The response is integrated for both basis columns of the physical state $(\zeta,d\zeta/dN)^T$, with zero matter response at the $z=2.1$ anchor.  Thus the result is a row-vector kernel on the same two-dimensional physical state space.
23: \mathcal F^{\rm state}_{f\sigma_8}(0)&=(2.172303247,\,3.071170118\times10^{-2}),\\
24: \mathcal F^{\rm state}_{\rm Weyl}(0)&=(6.167895805\times10^{-2},\,-2.498933170\times10^{-2}).
27: These are not two new amplitudes.  They are the two components of linear functionals acting on the same initial $n=2$ state.  One amplitude gives unique predictions only after the normalized state direction is independently fixed; the distance row alone has rank one.
30: The $f\sigma_8$ object above is a state-space fractional response kernel, using the R1 baseline $dD/dN$ normalization.  It is not yet the final single-number directional prediction because one $D_L$ calibration row cannot select the complete initial-state vector.
```

## manuscript/sections/04al_r1_n2_luminosity_distance_kernel.tex

```text
2: \label{sec:r1-n2-DL-kernel}
4: The remaining light-cone observable is evaluated with the gauge-invariant geometric luminosity-distance construction.\citep{YooScaccabarozzi2016}
5: For positive spatial curvature the background transverse distance and lensing efficiency use $S_K(\chi)=\sin(\sqrt K\chi)/\sqrt K$, with the corresponding curved-FLRW line-of-sight kernel.\citep{PyneBirkinshaw2004}
7: For the pure $A_{4i}$ physical $n=2$ mode, per unit $\mu=\hat{\mathbf n}\cdot\hat{\mathbf p}$,
9: Q_{2,\rm dip}(\chi)=\sin(2\sqrt K\,\chi),
10: \qquad \widehat\nabla^2 Q_{2,\rm dip}=-2Q_{2,\rm dip}.
20: with $C_K(\chi)=\cos(\sqrt K\chi)$ and the redshift, radial and convergence distortions evaluated along the same R1 physical-$n=2$ metric response. Here $\delta z$ denotes the dimensionless perturbation $\delta\ln(1+z)$, not an unnormalized change in the measured redshift.
23: At the reference bin center $z_\star=0.02$, the resulting luminosity-distance row on the $(\zeta,d\zeta/dN)$ state anchored at $z=2.1$ is
26: \mathcal F_{D_L}(z_\star)=(1.842472667\times10^{-1},\,-3.520884651\times10^{-3}).
32: \mathcal F_{D_L}(z_\star)\cdot X_{2,\rm anchor}=-1.897330117\times10^{-2}\pm4.144653167\times10^{-3}.
37: n_{\rm cal}=(1.910607424\times10^{-2},\,9.998174623\times10^{-1}),\qquad \mathcal F_{D_L}\cdot n_{\rm cal}=0.
41: This numerical light-cone calculation exposes a logical condition implicit in the single-calibration argument: one amplitude calibration determines the model only after the normalized physical mode trajectory $\widehat X_2$ has been fixed independently.
47: \text{normalized state direction }\widehat X_2.
```

## manuscript/sections/04ao_bhsm_optical_transfer.tex

```text
5: $S^3(R_H)$, the physical BHSM $n=2$ carrier, and local line-of-sight
17: k^\nu\nabla_\nu k^\mu=0,
19: k^\mu k_\mu=0,
21: and let $e_A^\mu$, $A=1,2$, be a parallel-transported Sachs screen basis.
24: \frac{d^2\mathcal D}{d\lambda^2}
63: 2H_K\mathcal A'
79: \gamma_1 & \gamma_2\\
80: \gamma_2 & -\gamma_1
87: -\frac12\,{\rm Tr}\,\delta\mathcal T,
89: \mathcal G_1
91: -\frac12(\delta\mathcal T_{11}-\delta\mathcal T_{22}),
93: \mathcal G_2
95: -\delta\mathcal T_{12}.
99: \kappa''+2H_K\kappa'&=\mathcal F,\\
100: \gamma_1''+2H_K\gamma_1'&=\mathcal G_1,\\
101: \gamma_2''+2H_K\gamma_2'&=\mathcal G_2.
106: \zeta_\gamma=\delta\ln(1+z),
114: (\zeta_\gamma,\kappa,p_\kappa,\gamma_1,p_1,\gamma_2,p_2)^T.
129: S_{\rm opt}=(\mathcal Z,\mathcal F,\mathcal G_1,\mathcal G_2)^T,
135: 0&0&0&0&0&0&0\\
136: 0&0&1&0&0&0&0\\
137: 0&0&-2H_K&0&0&0&0\\
138: 0&0&0&0&1&0&0\\
139: 0&0&0&0&-2H_K&0&0\\
140: 0&0&0&0&0&0&1\\
141: 0&0&0&0&0&0&-2H_K
149: 1&0&0&0\\
150: 0&0&0&0\\
151: 0&1&0&0\\
152: 0&0&0&0\\
153: 0&0&1&0\\
154: 0&0&0&0\\
155: 0&0&0&1
161: The physical global mode is $X_2=(q_2,\Pi_2)^T$ and the existing
164: h_2=\mathcal R_{g2}X_2,
166: \mathcal R_{g2}=-\mathcal A_g^{-1}\mathcal B_{g2}.
176: -\mathcal A_g^{-1}\mathcal B_{g\psi},
194: $K_T=-\nabla^2+B\nabla^4$ and
198: V_\Sigma^2
200: \frac{c^2}{8\pi}
208: $(\mathcal Z,\mathcal F,\mathcal G_1,\mathcal G_2)$ and
209: $\mathcal G_\psi=\mathcal L_\psi^{-1}$, then
219: \mathcal R_{g2}X_2
247: \langle\xi_{\rm opt}\rangle=0,
281: complexity is carried by a $7\times7$ path covariance rather than by an
299: D_L=(1+z)^2D_A
305: 2\zeta_\gamma-\kappa.
308: state is passed through the existing closed-$S^3$ source remapping,
321: The exact closed-$S^3$ optical propagator, the seven-state reduction and the
```

## manuscript/sections/04as_bhsm_optical_low_rank_covariance.tex

```text
31: \sigma_q^2 m_Om_O^\dagger,
36: C_{ii}C_{jj}-|C_{ij}|^2=0.
41: R_1=\frac{\lambda_1}{\sum_i\lambda_i}=1,
43: \rho_{21}=\frac{\lambda_2}{\lambda_1}=0.
```

## manuscript/sections/04at_effective_representation_audit.tex

```text
6: the development of the R1 cosmological bridge.  A conventional
26: X_2
28: \mathcal R_{g2}X_2
```

## manuscript/sections/09_results_R1.tex

```text
1: \section{Results: Reference Branch R1}
3: R1 is an exploratory reference branch rather than a uniquely derived BHSM cosmology. Its purpose is to demonstrate that the bridge can simultaneously realize late-time scalar activation, a selected physical \(n=2\) response, and near-standard high-\(n\) structure growth.
7: \Omega_{m0}=0.31,\quad
8: \Omega_{r0}=9.0\times10^{-5},\quad
9: \Omega_{k0}=-0.018,\quad
10: \lambda=1,\quad q=3,
12: the shooting condition \(E(1)=1\) gives
15: \frac{V_0}{M_{\rm Pl}^2H_0^2}=2.39861.
19: The scalar braiding and no-ghost combination evolve as shown in Table~\ref{tab:r1background}.
27: 10  & \(2.07\times10^{-7}\) & \(1.64\times10^{-6}\) & 0.00194\\
28: 3   & \(8.60\times10^{-5}\) & \(6.80\times10^{-4}\) & 0.03915\\
29: 2.1 & \(3.72\times10^{-4}\) & 0.00294 & 0.08060\\
30: 1.5 & 0.00122 & 0.00959 & 0.14321\\
31: 1.0 & 0.00380 & 0.02988 & 0.24551\\
32: 0.8 & 0.00622 & 0.04889 & 0.30787\\
33: 0.5 & 0.01346 & 0.10536 & 0.43200\\
34: 0.3 & 0.02271 & 0.17740 & 0.53538\\
35: 0   & 0.04940 & 0.38485 & 0.70791\\
39: \label{tab:r1background}
44: w_T(z=0)\simeq-0.906,
48: z_{\rm acc}\simeq0.683.
51: Relative to a curved \(\Lambda\)CDM model with identical \(H_0,\Omega_{m0},\Omega_{r0},\Omega_{k0}\), the R1 expansion rate differs by only order one percent over the late-time interval. The maximum difference in the tabulated range is about \(1.7\%\) near \(z\sim0.5\).
57: R_{f\sigma_8}(z)
59: \frac{f\sigma_{8,\rm R1}}
60: {f\sigma_{8,\Lambda{\rm CDM}}}
68: \(z\) & \(R_{f\sigma_8}\) & suppression\\
70: 2.1 & 0.99346 & 0.65\%\\
71: 1.5 & 0.98869 & 1.13\%\\
72: 1.0 & 0.98164 & 1.84\%\\
73: 0.8 & 0.97795 & 2.20\%\\
74: 0.5 & 0.97271 & 2.73\%\\
75: 0.3 & 0.97168 & 2.83\%\\
76: 0.0 & 0.98482 & 1.52\%\\
79: \caption{R1 growth forecast relative to matched curved \(\Lambda\)CDM. These values are reference-branch predictions and become formally preregistered only when reproduced by the canonical code and written to the hashed manifest.}
80: \label{tab:r1growth}
87: \text{enhanced horizon-scale }n=2\text{ susceptibility}\\
```

## manuscript/sections/09d_sn_sightline_tomography.tex

```text
6: multi-observable gates. The axis $(211.48^\circ,-12.81^\circ)$, amplitude
7: $A=-0.04321623436$ mag, $z_{\rm ref}=0.02$, $z_c=0.03$, local redshift
8: transfer, R1 background and original sample cuts remain fixed.
12: Individual epoch/band photometry for all 3,521 released Pantheon+ and DES
14: \citep{Scolnic2022,Sanchez2024DES,Abbott2024DES}.
15: There are 3,480 converged solutions and 41 non-converged local fits, all
25: the published standardized distance modulus minus the frozen R1 distance.
32: 2MRS neural reconstruction, trained on simulation mocks
33: \citep{LilowEtAl2024Fields}. Its valid radius is
34: $200\,h^{-1}\mathrm{Mpc}$; distant foregrounds remain unmeasured.
35: Observed 2MRS counts \citep{Huchra2012} are archived separately from the
58: \begin{tabular}{@{}p{0.37\linewidth}rrrp{0.18\linewidth}@{}}
60: Test & Rows & $\Delta\chi^2$ & Mock $p$ & Status\\
62: Pantheon, released distances & 1208 & $+2.985$ & 0.247 & PRIMARY\\
63: Pantheon, source/environment & 1208 & $-0.883$ & 0.0725 & CONDITIONAL\\
64: Pantheon, fixed external PV sensitivity & 1208 & $-4.909$ & 0.0110 & SENSITIVITY ONLY\\
65: DES, released distances & 1467 & $+0.240$ & 0.771 & PRIMARY\\
66: DES, Pantheon nuisance transfer & 1467 & $+1.321$ & 0.590 & CONDITIONAL\\
70: $\Delta\chi^2=\chi^2_{\rm frozen}-\chi^2_{\rm null}$ favors the prediction.
71: One-sided score tails use 2,000 full-covariance Gaussian mocks,
74: fixed at coefficient 1; complete joint uncertainty propagation unavailable.
75: All $\Delta\chi^2$ and $p$ values are dimensionless.}
87: their nominal 0.05 tail occurs in 0.09 of the Pantheon null calibration
92: of distance covariance \citep{HuiGreene2006,Smith2014Lensing}.
93: The released PV/redshift corrections \citep{Peterson2022PV,Carr2022Redshift}
94: and the 2M++ reconstruction \citep{Carrick2015} motivate comparison with
97: correction onto $+0.033953716$ mag, compared with $+0.009368488$ mag for
102: $A_{\rm fit}=-0.032870\pm0.012398$ mag with direction and redshift dependence
103: fixed. The source/environment nuisance space retains 79.0\% of the
104: bin-profiled template information, inflating its uncertainty by 1.125.
108: $10\,h^{-1}\mathrm{Mpc}$ Cartesian cells within the unchanged radial shells,
109: using $2\,h^{-1}\mathrm{Mpc}$ midpoint steps. The local matrix has 3,931
110: sampled columns, with numerical ranks 1,268 and 238 for the Pantheon and
111: DES row sets at relative singular tolerance $10^{-7}$.
114: smaller than 0.328 (Pantheon) and 0.531 (DES). A 0.0195 tail in one of 22
122: \texttt{analyses/sn\_sightline\_tomography\_v1} in the repository.
123: The full data-bearing local archive is identified there by SHA-256;
134: mean fixed at coefficient 1; complete joint uncertainty propagation
```

## manuscript/sections/09c_negative_results.tex

```text
12: Using the 2MRS release \citep{Huchra2012}, three frozen foreground proxies also fail held-out transfer.
18: Foreground proxy & Permutation $p$ & Held-out $\Delta\chi^2$\\
20: 2MRS number density & 0.076 & $+4.107958$\\
21: K-band luminosity & 0.092 & $+4.911833$\\
22: BHSM seam-charge proxy & 0.091 & $+4.968593$\\
26: Here $\Delta\chi^2=\chi^2_{\rm fixed\ foreground}-\chi^2_{\rm null}$;
31: Each uses 277 supernovae; the number-density catalog has 38600 galaxies,
32: the photometric variants 38601. The latter have 83 calibration and 194
```

## manuscript/sections/07_predictions_falsification.tex

```text
6: The Pantheon+ data release and cosmological analysis provide the supernova basis for the carried-forward tomographic calibration used here \citep{Scolnic2022,Brout2022}. The earlier scalar-topographic analysis reported
8: A_\mu^\star=-0.0412\pm0.009
12: 0.01<z<0.03,
18: (211.48^\circ,-12.81^\circ).
24: (-0.01895\pm0.00414)\cos\theta,
30: (+0.01895\pm0.00414)\cos\theta.
36: $A=-0.04321623436$ mag. The two calibrations must not be conflated, and
38: normalize only the selected physical \(n=2\) response,
42: \frac{\mathcal K_{\rm SN}S_2}{p_2C_2},
44: and do not independently determine \(S_2\), \(\mathcal K_{\rm SN}\), or \(C_2\). The physical distance row is now explicit; one scalar calibration leaves a state-space null direction.
50: \boxed{\Omega_k=-0.018\pm0.004}.
52: At \(H_0=67.4\ {\rm km\,s^{-1}\,Mpc^{-1}}\),
56: \frac{c}{H_0\sqrt{|\Omega_k|}}
57: \simeq33.15\ {\rm Gpc}.
59: The sign prediction is \(\Omega_k<0\) in the usual FLRW curvature convention.
61: This value is a frozen model forecast, not a current empirical best fit. DESI DR2 BAO results are reported as well described by flat \(\Lambda\)CDM, even while combinations with other probes motivate tests of extensions such as dynamical dark energy \citep{DESIDR2ResultsII2025}. The purpose of retaining \(-0.018\pm0.004\) is therefore prospective falsification: it must be compared against independent curvature inference without being moved toward the observed result.
65: The bridge predicts selective \(n=2\) softening with no analogous \(n\ge3\) enhancement. A future source-normalized \(S^3\) \(n=3\) susceptibility comparable to or larger than the physical \(n=2\) susceptibility would reject the minimal harmonic-filter branch. This statement is made at the compact-space harmonic level; comparison to a measured sky quadrupole requires an explicit observer-shell projection.
71: (211.48^\circ,-12.81^\circ)
79: C_n\simeq\chi\ell^2n^2,
81: so the extra scalar susceptibility falls as \(n^{-2}\). The model therefore predicts a strong separation between the first physical compact-space level and conventional clustering modes.
88: \item the source-normalized \(S^3\) \(n=3\) response is comparable to or exceeds physical \(n=2\);
93: \item the selected physical \(n=2\) response cannot satisfy the common-domain constraint or stability requirements.
```

## manuscript/sections/08b_make_or_break_prediction.tex

```text
4: \path{preregistration/bhsm_geometry_first_make_or_break_v1.json}.
7: hashed design before unblinding. Existing Pantheon+ and 2MRS analyses are
13: d=FX_2+\epsilon,\qquad
14: \chi^2_\perp=\min_X(d-FX)^T\Sigma^{-1}(d-FX).
17: With rank-$2$ $F$, known positive-definite Gaussian $\Sigma$ and no fitted
18: nuisances, $\chi^2_\perp\sim\chi^2_{N-2}$. Thus $N>2$ overidentifies the state.
19: Cholesky whitening and least squares avoid forming $\Sigma^{-1}$.
20: Reject the specific two-state deterministic realization at $p_A<0.01$.
26: data-selected axes and nonlinear kernels invalidate the simple $N-2$
32: $(211.48^\circ,-12.81^\circ)$; its frame and selection treatment belong in
39: \operatorname{rank}C_{\rm topo}\le1.
43: eigenvalues, $\rho_{21}=\lambda_2/\lambda_1$ when $\lambda_1>0$, and
44: $C_{ii}C_{jj}-|C_{ij}|^2$.
51: realization at $p_B<0.01$. Ratios/minors are diagnostics, not extra
56: fixed total null covariance and 9999 Monte Carlo replicates:
57: $p=(1+\#\{T_{\rm sim}\ge T_{\rm obs}\})/10000$.
68: decisive target inspection and cannot rescue this $m=1$ realization.
72: when data overlap. Each uses a $1\%$ threshold. Under both nulls, the
73: probability of at least one false rejection is bounded by $2\%$, not $1\%$.
```

## manuscript/sections/08a_prospective_test_matrix.tex

```text
11: \begin{tabular}{@{}p{0.17\linewidth} p{0.22\linewidth} p{0.28\linewidth} p{0.20\linewidth}@{}}
14: SN / standard candles & $F_{D_L}(z)X_2$ &
16: Amplitude or redshift evolution inconsistent with the independently selected state\\[0.5em]
18: Growth / RSD & $F_{f\sigma_8}(z)X_2$ &
20: Requires a second independent state not present in the reduction\\[0.5em]
22: Weak lensing / Weyl & $F_W(z)X_2$ and optical shear &
24: Axis or transfer incompatible with the distance response\\[0.5em]
28: Observed covariance structure incompatible with propagated optical state\\[0.5em]
40: C_{ii}C_{jj}-|C_{ij}|^2=0,
42: for every $2\times2$ minor of an isolated rank-one topographic covariance.
```

## manuscript/sections/10_discussion.tex

```text
4: not constrained-perturbation mathematics. The regular exclusion of $n=1$
9: The degree-two scalar eigenspace on $S^3$ has dimension nine; its dipole
13: The R1 coordinates $(\zeta,\zeta_{,N})$ are response coordinates; numerical
14: identification with canonical BHSM $(q_2,\Pi_2)$ requires a normalization
23: The positive $p^2$ term in a long-wavelength derivative expansion cannot
```

## manuscript/sections/09_geometry_first_conclusion.tex

```text
4: a selected $n=2$ response, induced metric/matter kernels, fixed-observed-redshift
11: stochastic texture amplitude. Rejection at $p<0.01$ eliminates the relevant
```

## manuscript/sections/09b_universe_submission_back_matter.tex

```text
9: \url{https://github.com/ncarberry64/Manuscript-Generation}.
13: archive is identified by SHA-256. A public permanent archive remains unavailable. Missing historical foreground-input checksums are not
28: Preprints.org 202601.1427, cited under its verified published
29: title \citep{Carberry2026Topographic}. That unreviewed preprint is not
```

## manuscript/sections/04_scalar_metric_bridge.tex

```text
4: The conditional R1 effective background action, used as a compatibility representation rather than a completed BHSM parent action, is taken to be
7: \int d^4x\sqrt{-g}
9: \frac{M_{\rm Pl}^2}{2}R
10: -\frac12\nabla_\mu T\nabla^\mu T
11: +\frac{1}{\Lambda^3}\Box T(\nabla T)^2
18: X\equiv-\frac12\nabla_\mu T\nabla^\mu T,
22: G_2=X-V(T),\qquad
23: G_3=\frac{2X}{\Lambda^3},\qquad
24: G_4=\frac{M_{\rm Pl}^2}{2},\qquad
25: G_5=0
27: in the convention \(\mathcal L_3=-G_3\Box T\).
29: The model belongs to the kinetic-gravity-braiding subclass of Horndeski theory \citep{DeffayetKGB2010,BelliniSawicki2014}. Since \(G_4\) is constant and \(G_5=0\),
31: M_*^2=M_{\rm Pl}^2,\qquad
32: \alpha_M=0,\qquad
33: \alpha_T=0.
39: \frac{2\dot{\bar T}^{\,3}}
40: {HM_{\rm Pl}^2\Lambda^3},
46: \frac{\dot{\bar T}^{\,2}}
47: {H^2M_{\rm Pl}^2}
49: 6\alpha_B.
55: \alpha_K+\frac32\alpha_B^2.
57: Their perturbation equations are formulated for a spatially flat background \citep{BelliniSawicki2014}. We therefore use \(D_{\rm kin}>0\) below as an inherited kinetic diagnostic for the cubic-Galileon sector, not as a complete proof of perturbative stability for the present closed background plus the additional reduced spatial response. A full stability proof requires the constrained gauge-invariant kinetic and gradient matrices, including the exceptional lowest scalar sector identified above.
61: V(T)=V_0-\mu^3T.
67: \frac12\dot T^2+V(T)
68: +\frac{6H\dot T^3}{\Lambda^3},
73: \frac12\dot T^2-V(T)
74: -\frac{2\dot T^2\ddot T}{\Lambda^3}.
78: J=\dot T+\frac{6H\dot T^2}{\Lambda^3}
82: \dot J+3HJ=\mu^3.
84: These background expressions are the standard kinetic-gravity-braiding relations specialized to the present \(G_2,G_3\) choice \citep{DeffayetKGB2010}.
88: 3M_{\rm Pl}^2
90: H^2+\frac{1}{a^2R_H^2}
96: The bridge does not require the bare Horndeski sound speed to become negative. Instead, the BHSM-inspired reduced response is used to lower the effective long-wavelength spatial stiffness while a positive \(p^2\) term stiffens shorter wavelengths. This statement concerns the reduced spatial kernel; it is not, by itself, a complete closed-FLRW stability theorem.
```

## manuscript/sections/04a_closed_horndeski_n2.tex

```text
2: \label{sec:closed-horndeski-n2}
4: The previous section uses the Bellini--Sawicki functions to organize the cubic-Galileon background. For the curvature-scale mode, however, the flat-FLRW perturbation formulas are not sufficient. We therefore specialize the non-flat Horndeski quadratic action derived by Akama and Kobayashi \citep{AkamaKobayashi2019}.
8: G_2=X-V(T),\qquad
9: G_3=\frac{2X}{\Lambda^3},\qquad
10: G_4=\frac{M_{\rm Pl}^2}{2},\qquad
11: G_5=0,
15: \mathcal F_T=\mathcal G_T=M_{\rm Pl}^2.
17: Because \(G_4\) is constant and \(G_5=0\), the explicit curvature contributions to their \(\Theta\) and \(\Sigma\) vanish. The remaining coefficients are
22: HM_{\rm Pl}^2-\frac{\dot T^3}{\Lambda^3}
30: \frac12\dot T^2
31: +\frac{12H\dot T^3}{\Lambda^3}
32: -3M_{\rm Pl}^2H^2.
39: HM_{\rm Pl}^2
40: \left(1-\frac{\alpha_B}{2}\right).
47: \frac{M_{\rm Pl}^2\Sigma}{\Theta^2},
49: \kappa(t)\equiv\frac{K}{a^2}.
51: The first physical scalar harmonic has alternative-convention index \(N=3\), so
53: k^2=K(N^2-1)=8K.
58: +8\Theta\kappa\,\chi
59: +3\Theta\dot\zeta
60: +5M_{\rm Pl}^2\kappa\,\zeta=0,
64: -M_{\rm Pl}^2\dot\zeta
65: -M_{\rm Pl}^2\kappa\,\chi=0.
72: \frac{5M_{\rm Pl}^2}{\Theta(r+8)}
75: -\frac{M_{\rm Pl}^2\kappa}{\Theta}\zeta
84: -\frac{r+3}{\kappa(r+8)}\dot\zeta
85: -\frac{5M_{\rm Pl}^2}{\Theta(r+8)}\zeta.
89: The reduced kinetic coefficient of this \(n=2\) mode is therefore
92: \mathcal G_{S,2}
94: M_{\rm Pl}^2
95: \frac{5(r+3)}{r+8}.
100: r+3
103: {2(1-\alpha_B/2)^2},
108: \mathcal G_{S,2}
110: M_{\rm Pl}^2
111: \frac{5D_{\rm kin}}
113: +10(1-\alpha_B/2)^2}.
116: Thus, provided \(\Theta\neq0\), positivity of the previously computed \(D_{\rm kin}\) is sufficient for positivity of the pure gravity--scalar kinetic coefficient of the first physical closed-universe scalar harmonic.
120: \mathcal F_{S,2}
122: \frac1a\frac{d}{dt}
124: \frac{5}{r+8}
126: \frac{aM_{\rm Pl}^4}{\Theta}
129: -M_{\rm Pl}^2
131: \frac{5}{r+8}
132: \frac{M_{\rm Pl}^6}{\Theta^2}
133: \frac{K}{a^2}.
135: The repository evaluates this expression as a gravity--scalar diagnostic only. The late-time R1 background contains matter and radiation, whose perturbations introduce additional dynamical variables and source terms. Consequently, \(\mathcal F_{S,2}\) is not promoted here to a complete matter-era stability or light-cone transfer result.
```

## manuscript/sections/04b_matter_sourced_n2.tex

```text
1: \section{Matter-Sourced Constraint Transfer at n=2}
3: The reduced kernel describes the pure gravity--scalar constrained dynamics of the first physical closed-universe harmonic. Here we introduce minimally coupled matter sources without importing a flat-space or quasi-static approximation. We use the curvature-aware Newtonian-gauge EFT equations written by Gambino and Pace \citep{GambinoPace2024}. Their metric convention is
5: ds^2
7: -(1+2\Phi)dt^2
9: a^2(1-2\Psi)\gamma_{ij}dx^idx^j,
13: \delta T^0{}_0=-\delta\rho,
15: \delta T^0{}_i=(\rho+p)D_i v.
20: \alpha_M=\alpha_T=\alpha_H=0,
22: M^2=M_{\rm Pl}^2.
29: -\frac12\alpha_B^{\rm BS}
38: \alpha_K+6(\alpha_B^{\rm GP})^2
40: \alpha_K+\frac32(\alpha_B^{\rm BS})^2
49: On the physical \(n=2\) harmonic,
51: D^2Q=-8KQ,
53: (D^2+3K)Q=-5KQ.
57: \kappa=\frac{K}{a^2},
59: the sourced \(00\) equation becomes
61: 6(1+\alpha_B^{\rm GP})H\dot\Phi
63: (\alpha_K-6\alpha_B^{\rm GP})H^2\dot\pi
66: (6-\alpha_K+12\alpha_B^{\rm GP})H^2
68: 10\kappa
73: 6H
75: \frac{\rho+p}{2M_{\rm Pl}^2}
77: \dot H(1+\alpha_B^{\rm GP})
79: \kappa\left(1-\frac83\alpha_B^{\rm GP}\right)
82: -\frac{\delta\rho}{M_{\rm Pl}^2}.
85: The \(0i\) equation is
87: 2\dot\Phi
89: 2H\alpha_B^{\rm GP}\dot\pi
91: 2(1+\alpha_B^{\rm GP})H\Phi
94: 2\dot H+\frac{\rho+p}{M_{\rm Pl}^2}
97: -\frac{q}{M_{\rm Pl}^2},
105: \frac43\rho_r v_r,
112: These equations form a \(2\times2\) local constraint-transfer system for
121: -2H^2
123: \alpha_K+6(\alpha_B^{\rm GP})^2
126: \boxed{-2H^2D_{\rm kin}}.
135: \frac{k^2}{a^2}v_m
137: 3\dot\Phi,
145: \frac43\frac{k^2}{a^2}v_r
147: 4\dot\Phi,
151: -\Phi-\frac14\delta_r,
155: \frac{k^2}{a^2}=8\kappa.
158: The repository implements these equations as a local linear transfer matrix with the additional topographic response explicitly switched off. This is deliberately not yet promoted to a globally integrated observable solution. The curvature-aware EFT equations and the reduced ADM kernel use different perturbation variables and normalizations; the remaining evolution equations must first be matched analytically and then subjected to a Bianchi/constraint-preservation audit. Only after that compatibility gate is passed will we integrate the full matter-era \(n=2\) system and add the independent reduced topographic operator.
```

## manuscript/sections/04c_bianchi_dictionary_n2.tex

```text
9: N=1+\delta n,\qquad
11: g_{ij}=a^2e^{2\zeta}\gamma_{ij},
17: and give the curvature-aware field equations in Newtonian gauge \citep{AkamaKobayashi2019,GambinoPace2024}.
41: The first implementation attempt treated the \(00\) and \(0i\) equations as a complete local evolution system and then demanded that the spatial-trace and scalar equations vanish for an arbitrary pointwise state. That test is too strong and uses the constraint equations in the wrong role.
47: while the \(00\) and \(0i\) equations define constraints on the initial data. Covariant consistency requires those constraints to be preserved by the evolution.
51: \alpha_M=\alpha_T=\alpha_H=0,\qquad
60: 2 & -2H\alpha_B\\
61: 6H\alpha_B & H^2\alpha_K
70: 2H^2\left(\alpha_K+6\alpha_B^2\right)
72: 2H^2D_{\rm kin}>0.
77: The repository now integrates the trace/scalar equations together with ideal dust and radiation conservation, starting from data that satisfy the sourced \(00\) and \(0i\) constraints. The Bianchi audit is then the preservation of those two constraints over the interval
79: 2.1\ge z\ge0.
```

## manuscript/sections/05_harmonic_selection.tex

```text
6: \qquad L_c^2p_n=\ell^2n(n+2),
12: C_n=c_s^2-\chi\left[1-\ell^2n(n+2)\right].
14: The physical scalar spectrum begins at \(n=2\). Selective softening of \(n=2\) with stiffening of every propagating \(n\ge3\) requires
16: 1-8\ell^2>0,
18: 1-15\ell^2<0,
22: \boxed{\frac1{15}<\ell^2<\frac18}.
26: C_2<c_s^2,
28: C_n>c_s^2\quad(n\ge3).
30: The formal \(n=1\) coefficient also lies on the soft side, but \(n=1\) is not a propagating scalar degree of freedom in the regular one-field closed-FLRW sector.
32: If \(0<C_2\ll1\),
34: T_2\simeq\frac{S_2}{p_2C_2}
38: \frac{T_3}{T_2}=\frac{S_3}{S_2}\frac{8C_2}{15C_3},
43: \left|\frac{T_3/S_3}{T_2/S_2}\right|
44: =\frac{8C_2}{15C_3}<\frac8{15}.
47: For the pure-dipole \(A_{4i}\) subspace of the selected \(n=2\) harmonic, the observer-sky pattern is exactly dipolar with geometric radial envelope \(\sin2\psi\).
```

## manuscript/sections/06_observables.tex

```text
15: If the selected \(n=2\) curvature-scale sector admits a physical gauge-invariant metric response, its exact scalar angular dependence can factor as
30: \mu=5\log_{10}D_L+\mathrm{const},
36: \frac{\ln10}{5}A_\mu\cos\theta.
41: k^2\Psi
43: -4\pi Ga^2\mu(k,a)\rho_m\delta_m,
51: \frac{\mu(k,a)}2[1+\gamma(k,a)].
56: \left(2+\frac{H'}H\right)D'
57: -\frac32\Omega_m(a)\mu(k,a)D=0,
61: f\sigma_8(k,z)
63: \frac{D'}{D}D\,\sigma_{8,\rm ini}
65: D'(k,z)\sigma_{8,\rm ini}.
68: The quasi-static limit is used only on scales inside the scalar sound horizon \citep{SawickiBellini2015}. Horizon-scale \(n=2\) observables require the full constrained long-wavelength response and must not be obtained by extrapolating quasi-static formulas into the first physical curvature-scale mode.
```

## manuscript/sections/04am_bhsm_los_topographic_transfer.tex

```text
4: The smooth $S^3(R_H)$ background and the physical $n=2$ response do not exhaust the geometry sampled by a photon. We therefore separate the global parent from a local deformation field $\psi(\mathbf x,t)$:
6: \gamma_{ij}=a^2 e^{2\psi}\bar\gamma^{S^3}_{ij}.
10: \delta{}^{(3)}R=-\frac{4}{a^2}\left(\bar\nabla^2+\frac{3}{R_H^2}\right)\psi,
12: so the physical manuscript $n=2$ harmonic carries the shifted factor $8-3=5$.  The local volume rate is
17: The retained BHSM topographic flux candidate is $F_T=\nabla T-B\nabla(\nabla^2T)$ with $\nabla\cdot F_T=S$.  The existing seam-charge candidate is
19: V_\Sigma^2=\frac{c^2}{8\pi}\int_{\Sigma_\star}H_\Sigma\Delta K_{uu}\,dA,
21: giving $V_\Sigma^2=GM/r_\star$ at leading spherical weak-field order.  The exact action-level export from this boundary charge into the dynamical $\psi$ source remains open; no normalization is invented here.
24: For a minimal coherent local-source proxy the accumulated radial perturbation grows only to a coherence distance $\chi_c$ and then saturates.  The closed-$S^3$ radial distance term therefore has the normalized shape
28: The value $z_c=0.03$ originated as an illustrative local-coherence choice
32: The table below retains the historical $-0.0412$ mag normalization;
39: 0.02 & -0.041200\\
40: 0.10 & -0.012598\\
41: 0.30 & -0.004435\\
42: 0.50 & -0.002810\\
43: 0.80 & -0.001904\\
44: 1.00 & -0.001605\\
45: 1.50 & -0.001212\\
46: 2.10 & -0.000991\\
51: The carried-forward low-redshift SN scalar amplitude normalizes the local topographic response in this exploratory decomposition.  It cannot simultaneously determine an independent global $n=2$ amplitude.  The global/local split must ultimately be fixed by the common BHSM action or an independent preregistered observable.
```

## manuscript/sections/04an_2mrs_foreground_transfer.tex

```text
1: \section{Frozen 2MRS foreground-transfer test}
2: \label{sec:2mrs-foreground-transfer}
6: field constructed from the 2MRS galaxy distribution.  For supernova \(i\),
14: where \(\delta_s\) is the fixed \(10^\circ\) Gaussian-smoothed density contrast
19: The frozen analysis used 277 Pantheon+ Hubble-flow
20: supernovae and 38600 2MRS galaxies with
21: \(z<0.05\).  Relative to the fixed-axis plus nuisance-control model, adding
24: \Delta\chi^2_{\rm fitted-null} = -1.619359,
26: \beta_F = 0.015286
27: \pm 0.007846\ {\rm mag},
31: p_{\rm perm}=0.076.
37: \(0.01\le z\le0.03\) and then held fixed for \(0.03<z\le0.15\).  The held-out
40: \Delta\chi^2_{\rm heldout}
42: +4.107958.
46: The 2MRS field is a matter-density proxy, not yet an action-normalized BHSM
```

## manuscript/sections/04ap_bhsm_optical_constitutive.tex

```text
15: \mathbf F_T=\nabla T-B\nabla(\nabla^2T),
25: Q S_{\rm enc}(r)=\frac{Q_\Sigma}{r^2},
27: Q_\Sigma=V_\Sigma^2,
34: \mathcal E_\Sigma^{(0)}[Q_\Sigma]
38: \frac{\Theta(r-r_\star)}{r^2}.
42: The $r^{-2}$ shape is inherited from the conserved-throughput neighborhood
49: \mathcal E_\Sigma^{(1)}
53: \dot Q_\Sigma\frac{\Theta(r-r_\star)}{r^2}
56: \frac{\delta(r-r_\star)}{r_\star^2}
68: Here the conventional spatial-potential sign gives Weyl response $(\Phi+\Psi)/2$. The growth section instead uses the positive spatial-curvature perturbation, giving $\Psi-\Phi$; convert the spatial-potential sign before composing maps. Choose a two-coordinate physical scalar metric basis
84: For $\Delta_g=ac-b^2\ne0$,
90: (\mathcal A_g^{\rm loc})^{-1}
93: \frac1{\Delta_g}
107: {2(ac-b^2)}.
129: With $\mathcal G_\psi=\mathcal L_\psi^{-1}$,
```

## manuscript/sections/04aq_bhsm_optical_hessian_projection.tex

```text
31: 2F''(X)(j\!\cdot\!d\varphi)A_{\mu\nu}
33: 2F'(X)j_{(\mu}\partial_{\nu)}\varphi
70: (\mathcal A_g^{\rm loc})^{-1}
81: together with the basis/profile normalizations, $\Delta_g=ac-b^2$, the
93: \mathcal R_{g\psi}\to S^{-1}\mathcal R_{g\psi},
```

## manuscript/sections/04ar_bhsm_optical_current_gate.tex

```text
12: 2F''(X)(j^k\nabla_k\varphi)A_{ij}
14: 2F'(X)j_{(i}\nabla_{j)}\varphi
22: The static conditions $D_0\eta=0$ and $j_0=0$ remove the static shift/phase
28: j_i=0
30: \nabla_i\varphi=0
32: u=v=0.
```

## manuscript/sections/08_numerical_protocol.tex

```text
21: \left(1+\frac{12Hv}{\Lambda^3}\right)\dot v
23: \frac{6v^2}{\Lambda^3}\dot H
25: \mu^3-3H\left(v+\frac{6Hv^2}{\Lambda^3}\right).
29: \frac{2v^2}{\Lambda^3}\dot v
31: 2M_{\rm Pl}^2\dot H
33: \rho_m+\frac43\rho_r+v^2+\frac{6Hv^3}{\Lambda^3}
34: -\frac{2M_{\rm Pl}^2}{a^2R_H^2}.
36: At each integration step these equations form a \(2\times2\) linear system for \((\dot v,\dot H)\). The remaining equations are
44: 3M_{\rm Pl}^2
45: \left(H^2+\frac{1}{a^2R_H^2}\right)
48: =0.
55: E=\frac{H}{H_0},
59: \lambda=\frac{\Lambda^3}{M_{\rm Pl}H_0^2},
61: q=\frac{\mu^3}{M_{\rm Pl}H_0^2}.
64: The exploratory reference branch R1 uses
66: \Omega_{m0}=0.31,\quad
67: \Omega_{r0}=9.0\times10^{-5},\quad
68: \Omega_{k0}=-0.018,\quad
69: \lambda=1,\quad q=3.
71: The constant \(V_0\) is determined only by the shooting condition \(E(a=1)=1\).
77: z_i=100.
81: J_i\simeq\frac{2\mu^3}{9H_i}.
85: v_i+\frac{6H_iv_i^2}{\Lambda^3}=J_i.
92: k=\{0.02,0.05,0.10,0.15,0.20\}\ {\rm Mpc}^{-1}
96: z=\{0.5,0.8,1.0,1.5,2.1\}.
99: A deterministic JSON manifest containing the action version, parameters, numerical tolerances, calibration values, and predicted observables is hashed with SHA-256. Any later scientific or numerical revision creates a new manifest version rather than overwriting the previous freeze.
```

## manuscript/sections/09a_reference_figures.tex

```text
6: \includegraphics[width=0.68\linewidth]{figures/r1_alphaB.png}
7: \caption{Conditional reference calculation, not an observational fit: kinetic-braiding function \(\alpha_B(z)\) for Reference Branch R1. The logarithmic vertical scale emphasizes the suppression of the braiding sector toward high redshift.}
8: \label{fig:r1alphab}
13: \includegraphics[width=0.68\linewidth]{figures/r1_omegaT.png}
14: \caption{Conditional reference calculation, not an observational fit: scalar energy fraction \(\Omega_T(z)\) for R1. The scalar component is subdominant at high redshift and becomes dynamically important only at late times.}
15: \label{fig:r1omegat}
21: \includegraphics[width=0.68\linewidth]{figures/harmonic_filter.png}
22: \caption{Sign of the Schur-response factor \(1-\ell^2 n(n+2)\) for the illustrative filtering value \(\ell=0.32\), matching the canonical code. The physical \(n=2\) branch is softened and \(n\ge3\) stiffened; the displayed \(n=1\) coefficient is exceptional and excluded from propagation.}
28: \includegraphics[width=0.68\linewidth]{figures/r1_growth_ratio.png}
29: \caption{Prospective R1 growth ratio relative to a matched curved-\(\Lambda\)CDM background with the same early normalization.}
30: \label{fig:r1growth}
```
