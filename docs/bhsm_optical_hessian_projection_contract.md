# BHSM local optical Hessian extraction contract

The remaining local optical coefficients are not fit parameters.

Use one common closed-FLRW physical domain and evaluate:

    a = <e_A,   H_gg e_A>
    b = <e_A,   H_gg e_psi>
    c = <e_psi, H_gg e_psi>

and, for the normalized local topographic profile phi_loc,

    u = B_gphi[e_A,   phi_loc]
    v = B_gphi[e_psi, phi_loc].

The mixed bilinear is the retained action variation

    integral dmu w gamma^{mu nu} [
        2 F''(X) (j.d phi) A_{mu nu}
      + 2 F'(X) j_(mu partial_nu) phi
      - F'(X) g_{mu nu} (j.d phi)
    ].

The optional topographic self-entry

    c_phi = <phi_loc, H_phiphi phi_loc>

must be evaluated on the same domain so that the local Schur curvature can be
checked.

Required output:
- a,b,c,u,v,c_phi
- metric-basis normalization
- topographic-profile normalization
- Delta_g = ac-b^2
- ||A_g R_gpsi + B_gpsi||
- local Schur curvature
- Legendre-gate status
- common-domain status

Do not use the algebra-witness coefficients in older helper tests as physical
numbers, and do not calibrate these entries to SN data.
