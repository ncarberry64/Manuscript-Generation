"""Generate deterministic preregistration manifest and SHA-256 digest."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import background, growth
from scipy.integrate import solve_ivp
import math

MANIFEST = Path(__file__).resolve().parents[1] / "preregistration" / "prediction_manifest.json"

def canonical_bytes(obj) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")

def build_manifest():
    v0 = background.shoot_v0()
    r1 = background.integrate(v0)

    gr1 = growth.integrate_growth(growth.growth_rhs_r1)
    gl = growth.integrate_growth(growth.growth_rhs_lcdm)
    zgrid = [0.5, 0.8, 1.0, 1.5, 2.1]
    grow = {}
    for z in zgrid:
        N = math.log(1/(1+z))
        _, dpr = gr1.sol(N)
        _, dpl = gl.sol(N)
        grow[str(z)] = float(dpr/dpl)

    return {
        "protocol": "BHSM_TOPO_COSMOLOGY_PREREG_V1",
        "status": "REFERENCE_BRANCH_R1_NOT_YET_PROMOTED_TO_FRAMEWORK_WIDE_FREEZE",
        "geometry": {
            "Omega_k": -0.018,
            "Omega_k_sigma": 0.004,
            "H0_km_s_Mpc": 67.4,
            "laplacian": "-Delta Y_n = n(n+2)/R_H^2 Y_n",
        },
        "background_R1": {
            "Omega_m0": background.OMEGA_M0,
            "Omega_r0": background.OMEGA_R0,
            "Omega_k0": background.OMEGA_K0,
            "lambda": background.LAMBDA,
            "q": background.Q,
            "z_initial": background.Z_INITIAL,
            "V0_over_Mpl2H02": v0,
        },
        "topographic_calibration": {
            "redshift_range": [0.01, 0.03],
            "A_mu": -0.041156,
            "A_mu_sigma": 0.009,
            "axis_RA_deg": 211.48,
            "axis_Dec_deg": -12.81,
        },
        "harmonic_filter": {
            "condition": "1/8 < ell^2 < 1/3",
            "softened_mode": 1,
            "higher_modes": "stiffened",
        },
        "growth": {
            "k_Mpc_inverse": [0.02,0.05,0.10,0.15,0.20],
            "redshift": zgrid,
            "R_fsigma8": grow,
            "comparison": "matched curved LCDM",
            "parameter_refit_after_freeze": False,
        },
    }

def main():
    manifest = build_manifest()
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    digest = hashlib.sha256(canonical_bytes(manifest)).hexdigest()
    print(f"manifest={MANIFEST}")
    print(f"sha256={digest}")

if __name__ == "__main__":
    main()
