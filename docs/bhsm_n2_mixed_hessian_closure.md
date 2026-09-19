# Physical n=2 mixed-Hessian closure

The finite observable-kernel inventory exposed the actual upstream blocker.

BHSM v5.11 already supplies the formal scalar/topographic operator

diag(-partial_rho^2+5, -Delta_B+5)

with off-diagonal trace/extension maps -I_mix and -I_mix^dagger.

On the physical n=2 S3 harmonic, -Delta_B -> 8, so the boundary diagonal is
13/L^2.  If lambda_rho is the selected collar radial eigenvalue and I_2 the
projected trace/extension matrix element, then

C_2,ST = (1/L^2) [[lambda_rho+5, -I_2],
                  [-I_2*,        13]]

and this block is positive iff

|I_2|^2 < 13 (lambda_rho+5).

The remaining common action data are:

- A_g^(2): physical scalar metric Hessian;
- B_g2: metric/topographic mixed block;
- I_2: nonhomogeneous collar-to-boundary trace/extension matrix element;
- lambda_rho: collar radial eigenvalue on the selected domain;
- D_phys^(2): common self-adjoint physical domain.

Once these are derived from the retained BHSM action, one Schur/Feshbach
reduction fixes K_2 and R_g2, covariant conservation fixes R_m2, and all five
observable kernels follow.

The existing growth.py R_fsigma8 quantity is an isotropic/reference
R1-to-LCDM growth ratio and is not the directional BHSM-native F_fsigma8
kernel.
