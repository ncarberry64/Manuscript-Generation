# BHSM-native cosmology bridge

This note records the repository lineage used to reformulate the cosmology
manuscript around the native BHSM mode ontology.

## BHSM source lineage

- PR #84: modes are physical boundary-local fields acted on by sector
  projectors; standalone 4D promotion is not automatic.
- PR #132: conditional unified action includes a scalar/topographic sector,
  symbolic variational equations, boundary terms, and quadratic/Hessian
  extraction.
- PR #134: scalar/topographic vacuum condition is derived conditionally from
  an action-reduced mode functional.
- PR #135: normalized scalar/topographic profile-boundary BVP is closed
  conditionally.
- PR #177: a boundary/interface mode is not promoted to a positive-norm 4D
  field until the relevant kinetic/operator closure is established; a
  constraint-determined transition trace can exist without such promotion.

## Cosmology consequence

The physical n=2 topographic carrier is therefore represented first as an
action-owned envelope/boundary mode on the physical tangent space.  A
four-dimensional Horndeski/Stueckelberg variable may represent its effective
late-time response, but is not counted as an additional ontological scalar.

The reduced evolution theorem is obtained by projecting the full Hessian onto
the physical n=2 subspace and Schur/Feshbach eliminating constrained
directions.  The resulting amplitude q2 and its action-generated conjugate
momentum form the retained reduced phase-space pair.

## Claim boundary

This reformulation does not manufacture the still-open numerical transfer
normalization between the BHSM reduced operator and every curved-EFT
observable.  It changes what is foundational: the BHSM reduced operator is
the transport law; the curved EFT DAE is a compatibility representation that
must reproduce the same physical mode count.
