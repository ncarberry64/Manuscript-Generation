
# Derivation ledger

## D1 Ã¢â‚¬â€ S3 spectrum
Input:
`ds_3^2 = a^2 R_H^2[dpsi^2 + sin^2(psi)dOmega_2^2]`.

Result:
`p_n(a)=n(n+2)/(a^2 R_H^2)`.

Status: exact geometry.

## D2 Ã¢â‚¬â€ n=1 scalar observer-sky projection
Using the R4 embedding of S3:
`T_1=A[u_4 cos psi + u_vecÃ‚Â·n_hat sin psi]`.

Result: a fixed-radius shell of the scalar function contains sky monopole + dipole only.

Status: exact scalar geometry.

## D2G Ã¢â‚¬â€ exceptional closed-FLRW metric sector
Common literature convention:
`k^2=N^2-1`, `N=1,2,...`, with `N=n+1` relative to this manuscript.

Our `n=1` therefore maps to literature `N=2`. Standard gauge-invariant closed-FLRW scalar Hamiltonians contain `1/(N^2-4K)` and the usual scalar-derived tensor-harmonic construction is exceptional for the lowest modes when `K=+1`.

Result:
the scalar geometry in D2 survives, but its physical metric/light-cone sourcing is not yet derived.

Status: open bridge calculation; all n=1 metric observables are conditional on closure.

## D3 Ã¢â‚¬â€ Schur-generated higher-spatial stiffness
Illustrative parent block:
`A=Z0 p`, `D=M_c^2+c_c p`, `C=g sqrt(p)`.

Exact:
`K_red=Z0 p - g^2 p/(M_c^2+c_c p)`.

Long-wave:
`K_red=(Z0-g^2/M_c^2)p + (g^2 c_c/M_c^4)p^2 + ...`.

Status: exact algebra given the declared illustrative block; the block itself is exploratory.

## D4 Ã¢â‚¬â€ harmonic filter
Set `L_c=a R_H ell`.

Then:
`C_n=c_s^2-chi[1-ell^2 n(n+2)]`.

If `1/8 < ell^2 < 1/3`:
- n=1 softened;
- n>=2 stiffened.

Status: exact reduced-spatial bridge result.

Important:
the n=2/n=1 source-normalized ratio is an S3 harmonic-level ratio, not automatically a sky quadrupole/dipole ratio.

## D5 Ã¢â‚¬â€ cubic-Galileon kinetic diagnostic
For the declared KGB/Horndeski model:
`D_kin=alpha_K+3 alpha_B^2/2`.

Status:
verified inherited Horndeski kinetic diagnostic. It is not by itself a complete stability theorem for the closed background plus reduced spatial kernel; Bellini-Sawicki formulate the standard perturbation equations on spatially flat FRW.

## D6 Ã¢â‚¬â€ R1 background
Numerical branch:
`Omega_m0=0.31`, `Omega_r0=9e-5`, `Omega_k0=-0.018`, `lambda=1`, `q=3`.

Shooting:
`E(a=1)=1`.

Status: exploratory numerical reference branch.

## D7 Ã¢â‚¬â€ R1 growth
Growth uses identical early normalization and matched curved LCDM reference.

Output:
`R_fsigma8 = D'_R1 / D'_LCDM`.

Status:
prospective R1 forecast under the stated high-n/background-only approximation.

## D8 Ã¢â‚¬â€ carried-forward observational calibration
Pantheon+SH0ES low-z amplitude and direction are carried forward from the earlier tomographic analysis.

Status:
calibration input, not independently reproduced by the current repository.

## D9 Ã¢â‚¬â€ frozen curvature target
`Omega_k=-0.018 Ã‚Â± 0.004`.

Status:
pre-existing prospective model forecast. Current DESI DR2 BAO results are reported as well described by flat LambdaCDM, so the frozen value must not be presented as a current empirical best fit or retuned after comparison.


## D10 â€” exceptional n=1 no-go in regular one-field closed FLRW
For manuscript n=1: `k^2=3` and `D_iD_jQ=-gamma_ijQ`, hence `(D_iD_j+k^2 gamma_ij/3)Q=0`. The scalar-derived trace-free metric harmonic is absent; under the regular one-scalar diffeomorphism constraints this sector is nondynamical/pure gauge.

## D11 â€” first physical n=2 pure-dipole subspace
`T_2=A_AB X^A X^B`, `A_A^A=0`, decomposes on an observer shell as `A44(c^2-s^2/3)+2cs a_i n^i+s^2 S_ij n^i n^j`. The pure cross term is `T_2,dip=A sin(2 psi) p_hat.n_hat`.

## D12 â€” corrected physical harmonic filter
`n=2` softened and `n>=3` stiffened gives `1/15 < ell^2 < 1/8`, with next physical response `|(T3/S3)/(T2/S2)|=8 C2/(15 C3)`.
## D13 â€” exact closed-Horndeski n=2 constraint kernel

Source: Akama & Kobayashi (2019), non-flat Horndeski quadratic scalar action.

Specialization:
`G2=X-V`, `G3=2X/Lambda^3`, `G4=M_Pl^2/2`, `G5=0`.

Exact:
`F_T=G_T=M_Pl^2`,
`Theta=H M_Pl^2-dot(T)^3/Lambda^3`,
`Sigma=dot(T)^2/2+12H dot(T)^3/Lambda^3-3M_Pl^2 H^2`.

For manuscript n=2 / literature N=3: `k^2=8K`.

Define `r=M_Pl^2 Sigma/Theta^2`.

Then:
`G_S,2=M_Pl^2*5(r+3)/(r+8)`.

Using the alpha dictionary:
`r+3=D_kin/[2(1-alpha_B/2)^2]`.

Hence:
`G_S,2=M_Pl^2*5D_kin/[D_kin+10(1-alpha_B/2)^2]`.

Status: exact gravity-scalar bridge result for the first physical closed
harmonic. Matter/radiation perturbations and the additional topographic
response remain to be coupled before claiming a complete late-time observable
transfer matrix.
