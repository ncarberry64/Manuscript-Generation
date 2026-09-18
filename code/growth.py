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
