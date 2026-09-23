# BHSM optical spatial-current support gate

A static zero-momentum eta branch has zero ADM momentum density and zero
static shift/phase mixed block, but this does not imply zero spatial
metric/topographic mixing.

For a phase texture phi, the retained spatial mixed source is

    S_ij =
        2 F'' (j.grad phi) A_ij
      + 2 F'  j_(i grad_j) phi
      - F' g_ij (j.grad phi).

Therefore

    j_i = 0  OR  grad_i phi = 0
       => u = v = 0.

A nonzero spatial target Noether current and a nonconstant texture are
necessary, though not sufficient, for local metric sourcing.

If the texture has one fixed shape and one stochastic amplitude q,

    B_gpsi = q b_loc,

then its mixed-source covariance is rank one:

    C_B = Var(q) b_loc b_loc^T.

With m independent texture channels, the topographic contribution has rank at
most m before adding ordinary-gravity or other independent source channels.

This is the next useful data reduction: reconstruct a small number of
action-owned texture/current channels rather than an arbitrary optical
covariance or an object-by-object foreground ledger.
