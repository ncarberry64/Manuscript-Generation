#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, math, sys
from pathlib import Path
import numpy as np

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo",default=str(Path.home()/"Manuscript-Generation"))
    ap.add_argument("--output",default=str(Path.home()/"Downloads"/"BHSM_R1_COUPLED_ENVIRONMENTAL_STATE_AUDIT_LOCAL.json"))
    args=ap.parse_args()
    repo=Path(args.repo).expanduser().resolve()
    sys.path.insert(0,str(repo/"code"))
    import action_native_matter_n2 as an
    import closed_horndeski_n2 as ch

    bg=an.frozen_background()
    prop=an.integrate_propagator(bg)
    N0=an.N_START

    def physical_map(N):
        y=bg.sol(N)
        vel,aux=an.canonical_maps(N,y)
        a=math.exp(N)
        G=ch.gs_n2(y[1],y[2])
        R=np.zeros((6,6))
        R[0,0]=1.0
        R[1,:]=2*a**3*G*vel[0,:]
        R[2,0],R[2,1]=-3.0,1.0
        R[3,0],R[3,2]=-4.0,4.0/3.0
        R[4,:]=aux[2,:]
        R[5,:]=aux[3,:]
        return R

    R0=physical_map(N0)
    S=np.linalg.inv(R0)
    seed0=an.topographic_seed(bg)
    seed_error=float(np.max(np.abs(S[:,:2]-seed0)))

    zs=[1.5,1.0,0.8,0.5,0.3,0.1,0.02,0.0]
    rows=[]
    for z in zs:
        N=-math.log1p(z)
        U=np.asarray(prop.sol(N),float).reshape(6,6)
        y=bg.sol(N)
        vel,_=an.canonical_maps(N,y)
        O=np.vstack([np.eye(6)[0],vel[0]])  # zeta,zeta_dot
        T=O@U@S
        E=T[:,2:6]
        sv=np.linalg.svd(E,compute_uv=False)
        dm_vm=E[:,[0,2]]
        rows.append({
            "z":z,
            "U_XE_zeta_zdot":E.tolist(),
            "singular_values":sv.tolist(),
            "rank_rel_1e-8":int(np.sum(sv > 1e-8*sv[0])),
            "matter_delta_v_determinant":float(np.linalg.det(dm_vm)),
        })

    payload={
        "schema":"BHSM-R1-coupled-environmental-state-audit-local-v1",
        "repo":str(repo),
        "physical_anchor_coordinates":["q2","Pi2","delta_m","delta_r","v_m","v_r"],
        "seed_inverse_residual":float(np.max(np.abs(R0@S-np.eye(6)))),
        "topographic_seed_recovery_max_abs_difference":seed_error,
        "environment_inputs":["delta_m","delta_r","v_m","v_r"],
        "output_pair":["zeta","zeta_dot"],
        "rows":rows,
        "retuning":False,
    }
    out=Path(args.output).expanduser().resolve()
    out.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))

if __name__=="__main__":
    main()
