# Derivation ledger

## D1 â€” S3 spectrum
Input:
`ds_3^2 = a^2 R_H^2[dchi^2 + sin^2(chi)dOmega_2^2]`.

Result:
`p_n(a)=n(n+2)/(a^2 R_H^2)`.

Status: exact geometry.

## D2 â€” n=1 observer-sky projection
Using the R4 embedding of S3:
`T_1=A[u_4 cos chi + u_vecÂ·n_hat sin chi]`.

Result: fixed-redshift sky contains monopole + dipole only.

Status: exact geometry.

## D3 â€” Schur-generated higher-spatial stiffness
Parent block:
`A=Z0 p`, `D=M_c^2+c_c p`, `C=g sqrt(p)`.

Exact:
`K_red=Z0 p - g^2 p/(M_c^2+c_c p)`.

Long-wave:
`K_red=(Z0-g^2/M_c^2)p + (g^2 c_c/M_c^4)p^2 + ...`.

Status: exact algebra given the declared block.

## D4 â€” harmonic filter
Set `L_c=a R_H ell`.

Then:
`C_n=c_s^2-chi[1-ell^2 n(n+2)]`.

If `1/8 < ell^2 < 1/3`:
- n=1 softened;
- n>=2 stiffened.

Status: exact bridge result.

## D5 â€” R1 background
Action:
cubic Galileon + linear soft-breaking potential.

Numerical branch:
`Omega_m0=0.31`, `Omega_r0=9e-5`, `Omega_k0=-0.018`, `lambda=1`, `q=3`.

Shooting:
`E(a=1)=1`.

Status: exploratory numerical reference branch.

## D6 â€” R1 growth
Growth uses identical early normalization and matched curved LCDM reference.

Output:
`R_fsigma8 = D'_R1 / D'_LCDM`.

Status: prospective reference-branch forecast; freeze only after canonical regeneration and hashing.
