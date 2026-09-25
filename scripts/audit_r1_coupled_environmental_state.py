#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, math, sys
from pathlib import Path
import numpy as np

def build_audit(repo):
    repo=Path(repo).expanduser().resolve()
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
    stored=json.loads((repo/"artifacts/R1_ACTION_NATIVE_MATTER_BACKREACTION_V1.json").read_text())
    stored_z1p5=next(row for row in stored["samples"] if row["z"]==1.5)
    reproduction=float(np.max(np.abs(prop.sol(-math.log1p(1.5)).reshape(6,6)
                                    -np.asarray(stored_z1p5["full_canonical_propagator"]))))
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
        physical=physical_map(N)@U@S
        physical_E=physical[:2,2:]
        factor=2*math.exp(3*N)*ch.gs_n2(y[1],y[2])
        physical_sv=np.linalg.svd(physical_E,compute_uv=False)
        rows.append({
            "z":z,
            "U_XE_zeta_zdot":E.tolist(),
            "singular_values":sv.tolist(),
            "rank_rel_1e-8":int(np.sum(sv > 1e-8*sv[0])),
            "matter_delta_v_determinant":float(np.linalg.det(dm_vm)),
            "physical_propagator":physical.tolist(),
            "U_XE_q2_Pi2":physical_E.tolist(),
            "q2_Pi2_singular_values":physical_sv.tolist(),
            "q2_Pi2_rank_rel_1e-8":int(np.sum(physical_sv>1e-8*physical_sv[0])),
            "q2_Pi2_matter_determinant":float(np.linalg.det(physical_E[:,[0,2]])),
            "Pi2_over_zeta_dot":float(factor),
            "output_basis_conversion_error":float(np.max(np.abs(physical_E-np.diag([1.,factor])@E))),
        })

    payload={
        "schema":"BHSM-R1-coupled-environmental-state-audit-local-v2",
        "repo":str(repo),
        "physical_anchor_coordinates":["q2","Pi2","delta_m","delta_r","v_m","v_r"],
        "seed_inverse_residual":float(np.max(np.abs(R0@S-np.eye(6)))),
        "topographic_seed_recovery_max_abs_difference":seed_error,
        "stored_full_propagator_z1p5_max_abs_difference":reproduction,
        "environment_inputs":["delta_m","delta_r","v_m","v_r"],
        "output_pair":["zeta","zeta_dot"],
        "rows":rows,
        "observational_environment_selection":False,
        "scope":"Linear basis responses only; no physical E_i inferred or fitted. Pi2 is the reference-normalized coordinate, not the full coupled canonical p_zeta.",
        "input_SHA256":{p.relative_to(repo).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [repo/"preregistration/prediction_manifest.json",repo/"code/action_native_matter_n2.py",
                      repo/"code/closed_horndeski_n2.py",repo/"code/background.py",
                      repo/"artifacts/R1_ACTION_NATIVE_MATTER_BACKREACTION_V1.json"]},
        "retuning":False,
    }
    return payload

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--repo",default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--output",default=str(Path.home()/"Downloads"/"BHSM_R1_COUPLED_ENVIRONMENTAL_STATE_AUDIT_LOCAL.json"))
    args=ap.parse_args()
    payload=build_audit(args.repo)
    out=Path(args.output).expanduser().resolve()
    out.write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(payload,indent=2))

if __name__=="__main__":
    main()
