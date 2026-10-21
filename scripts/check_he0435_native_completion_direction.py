"""No-fit actual target/FIELD full-chain derivative check on the saved mean."""
from pathlib import Path
import argparse,hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'code'))
import run_he0435_native_completion as run


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training'],default='all_data');args=parser.parse_args()
    dest=run.OUT/args.label;solution=dest/'solution.npz';x=np.load(solution)['parameters']
    problem,_,_,_,_,_=run.build(args.label)
    problem.contexts=[problem.contexts[0],problem.contexts[-1]]
    r,j=problem.evaluate(x,True);rng=np.random.default_rng(81569)
    direction=rng.normal(size=len(x))*problem.scales*.02
    # Positive sky bound coefficients are intentionally held fixed; a centered
    # numerical derivative must not cross their physical parameter boundary.
    direction[problem.target.sky_start:problem.target.bg_start]=0.
    direction[-1]=1e-8
    analytic=j@direction;rows=[]
    for step in [2e-4,1e-4]:
        plus=problem.evaluate(x+step*direction)[0];minus=problem.evaluate(x-step*direction)[0]
        finite=(plus-minus)/(2*step);delta=finite-analytic
        rows.append({'dimensionless_step':step,'relative_L2_error':float(np.linalg.norm(delta)/np.linalg.norm(analytic)),
            'maximum_absolute_whitened_direction_error':float(abs(delta).max()),'L2_analytic':float(np.linalg.norm(analytic))})
    run.prior.dump(dest/'actual_direction_check.json',{'scope':'Saved all-native coupled mean, one actual target and one actual FIELD patch. No optimizer, native read processing, unused data or change of saved parameters.',
        'records':[ctx['name'] for ctx in problem.contexts],'direction_includes':'centroids, positive log-flux, shared318eta, electronic planes and kappa; bound sky coefficients held fixed.',
        'checks':rows,'pass_relative_L2_below_1e-5':all(row['relative_L2_error']<1e-5 for row in rows),
        'solution_sha256':hashlib.sha256(solution.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'physical_validation_pass':False,'physical_admission':False})
    print(rows,flush=True)


if __name__=='__main__':main()
