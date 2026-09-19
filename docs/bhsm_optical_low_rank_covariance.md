# BHSM optical low-rank covariance prediction

The compact full-preimage spatial-current background needed to evaluate the
local mixed coefficients numerically is not yet available in the retained
BHSM record. The existing degree-one surrogate must not be promoted to that
physical background.

A normalization-independent prediction nevertheless survives.

If m physical topographic texture channels source an observable vector O,

    delta O = M_O q,

then

    C_O^topo = M_O C_q M_O^T

and therefore

    rank(C_O^topo) <= m.

For one channel,

    C_O^topo = sigma_q^2 m m^T,

so all 2x2 minors vanish and the isolated covariance has one nonzero
eigenvalue.

The observational test must subtract independently modeled conventional and
measurement covariance first. It must not tune bins, observable selection or
channel count to produce a low-rank answer.
