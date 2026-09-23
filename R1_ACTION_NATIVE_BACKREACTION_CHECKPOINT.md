# R1 action-native backreaction checkpoint

**STOPPED FOR USAGE PRESERVATION — candidate action, not certified closure.**

Branch: `theory/cosmology-universe-mdpi-submission-pass`  
Exact current/base commit at capture: `8d0d4877a0de8a9145ed4175d25cba7111315ebe`. New action/checkpoint files are durable working-tree files and are not committed; this commit therefore remains HEAD.

## Resume here

- Full equations, assumptions and derivation: `docs/R1_ACTION_NATIVE_MATTER_BACKREACTION_V1_DERIVATION.md`.
- Deterministic module: `code/action_native_matter_n2.py`.
- Complete 401-point scan and propagator: `artifacts/R1_ACTION_NATIVE_MATTER_BACKREACTION_V1.json`.
- Machine-readable handoff: `R1_ACTION_NATIVE_BACKREACTION_CHECKPOINT.json`.

Variables: `x=(zeta,D_m,D_r)`, `u=(alpha,chi,v_m,v_r)`, canonical state `(x,p)`.
`D_I=delta rho_I/(rho_I+p_I)+3zeta`. Schutz–Sorkin convention:
`S_I=-integral[sqrt(-g)rho(n)+J^mu partial_mu ell]`, `delta ell=mu v`.
Dust w=cs²=0 exactly; radiation w=cs²=1/3. Frozen R1 densities and V0 are read unchanged.

## Matrices already derived

Write `b=8K/a²`, `s=5K/a²`, `h_I=rho_I+p_I`, `h=h_m+h_r`.

```text
A = [ 2Sigma  2bTheta       0       0 ]
    [ 2bTheta -2b*kappa    b*h_m   b*h_r ]
    [ 0        b*h_m     -b*h_m    0 ]
    [ 0        b*h_r       0     -b*h_r ]
B = [6Theta,0,0; -2b,0,0; 0,h_m,0; 0,0,h_r]
D = [2s+3h,-h_m,-h_r; 0,0,0; 0,0,0; 0,0,0]
K0=diag(-6,0,0)
E=diag(2s,0,0)-sum h_I*w_I*d_I*d_I^T
  d_m=(-3,1,0), d_r=(-3,0,1)
```

Full reciprocal Hessian in `(u,xdot,x)` is `[A B D; B^T K0 0; D^T 0 E]`.
After velocity elimination the lapse/shift matrix is
`[2Sigma,2bTheta; 2bTheta,b(h-2kappa)]`.
`J=B xdot+D x`; `u*=-A^-1 J`; reduced action is
`a³/2 [xdot^T K0 xdot+x^T E x-J^T A^-1 J]`.
Equivalently `M=K0-B^T A^-1 B`, `L=-B^T A^-1 D`, `V=E-D^T A^-1 D`;
`p=a³(M xdot+Lx)`, `xdot=M^-1(p/a³-Lx)`, `pdot=a³(L^T xdot+Vx)`.
No new phenomenological topographic/fluid force and no numerical projection.

## Checks completed, gates not yet certified

- A: quick vacuum constraint and `M=2G_S,2` checks passed at z=2.1,1,0; full vacuum evolution identity still pending.
- B: fluid conservation derived; native residual maximum `4.55202e-10` after state scaling.
- C: Hessian symmetry quick checks pass; independent polynomial differentiation test pending.
- D/E: all sampled auxiliary blocks invertible; all three physical kinetic eigenvalues positive across 401 points. Minimum `0.00248499050267`. These are sampled candidate checks, not certified interval/gradient stability.
- F: no independently fitted or added matter/topographic coupling.
- G: reconstructed action-constraint residual `2.85249e-11`; symplectic residual `4.53021e-09`; no-slip row residual `9.89644e-10`. No stepwise projection. Independent Noether certification remains pending.
- **H FAILS:** GP local-constraint residual/state scale is `363.624097949`. This is a large unresolved mismatch, not a rounding error. Do not force agreement or patch GP Eq.12.

The full six-state propagator and two unfitted topographic seed columns already exist.
Seeds set ADM fluid density/velocity to zero at z=2.1, using
`Pi2=2 a_anchor³ G_S,2 zdot_anchor`. They do not assert unique primordial fluid conditions.
No new automated test suite has been run; no observational return is authorized by these incomplete gates.

## Unresolved algebra / exact next step

Run from `C:\Users\carbe\Manuscript-Generation`:

```powershell
python scripts/probe_action_native_checkpoint.py
```

The probe isolates the failed GP comparison column by column without changing equations.
Audit signs in `v_Newtonian=v_ADM-chi`, density time shift, `pi=-chi`, source normalization,
the audit's use of no-slip to obtain Phi_dot, and total-background tadpole cancellation.
Do not assume which side is wrong. Then independently check the complete polynomial Hessian,
sequential Schur reduction and the vacuum reduced equation. Add small focused tests only.
Do not rerun the broad suite, GP localization campaign, observations or LaTeX.

User-reported prior GP localization is preserved verbatim in the numerical artifact;
its later detailed artifacts were not found in this checkout. No observational parameter,
frozen manifest, original owner, GP equation or manuscript was changed. Pre-existing
untracked SN inputs remain untouched. No unresolved free coefficient has been introduced;
the unresolved item is mathematical certification and the compatibility mismatch.
