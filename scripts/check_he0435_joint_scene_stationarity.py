"""No-fit necessary target-nuisance KKT checks at saved provisional means."""
from pathlib import Path
import sys,json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_joint_scene as run
js=run.js


def check(problem,p,eta):
    run.set_eta(problem.records,eta)
    r,j,_,_=problem.evaluate(p);rp,jp=problem.prior(p)
    gradient=j.T@r+jp.T@rp
    lo,hi=problem.bounds();atlo=np.isfinite(lo)&np.isclose(p,lo,rtol=0,atol=1e-8)
    athi=np.isfinite(hi)&np.isclose(p,hi,rtol=0,atol=1e-8)
    projected=gradient.copy()
    projected[atlo]=np.minimum(projected[atlo],0.)
    projected[athi]=np.maximum(projected[athi],0.)
    scaled=gradient*problem.scales;pg=projected*problem.scales
    index={'point_coordinate':np.array([3*i+d for i in range(problem.n) for d in [0,1]]),
           'log_flux':np.array([3*i+2 for i in range(problem.n)]),
           'detector_plane':np.arange(problem.bg_start,problem.size),
           'shared_sky':np.arange(problem.sky_start,problem.bg_start)}
    groups={}
    for name,take in index.items():
        groups[name]={'parameters':len(take),'lower_bound_count':int(atlo[take].sum()),'upper_bound_count':int(athi[take].sum()),
                     'unconstrained_scaled_gradient_max_abs':float(np.max(abs(scaled[take]))) if len(take) else 0.,
                     'projected_KKT_scaled_gradient_max_abs':float(np.max(abs(pg[take]))) if len(take) else 0.}
    return {'target_raw_chi2':float(r@r),'sky_prior_penalty':float(rp@rp),
            'groups':groups,'all_target_nuisance_projected_KKT_scaled_gradient_max_abs':float(np.max(abs(pg))),
            'parameter_gradient_of_half_objective':gradient,'scaled_projected_KKT_gradient':pg,
            'active_lower_indices':np.flatnonzero(atlo),'active_upper_indices':np.flatnonzero(athi)}


def main():
    records,fields,eta0,wcss,origin=run.load_inputs();allproblem=js.SceneProblem(records);answers={}
    unchanged=[run.OUT/'all_data.json',run.OUT/'held_exposure_training.json',run.OUT/'all_data_linearized_arrays.npz',run.OUT/'held_exposure_training_linearized_arrays.npz']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    before={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in unchanged}
    for label,problem in [('all_data',allproblem),('held_exposure_training',js.SceneProblem([r for r in records if r['exposure']!=allproblem.exposures[-1]]))]:
        a=np.load(run.OUT/f'{label}_linearized_arrays.npz')
        answers[label]=check(problem,a['parameters_updated'],a['eta_updated'])
    after={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in unchanged}
    if before!=after:raise AssertionError('Saved original result bytes changed')
    result={'scope':'No-fit saved-mean necessary stationarity check of target point/flux/plane/shared-sky parameters at fixed saved eta. Actual original target residuals and declared sky prior are differentiated. No field nuisance refit, joint response-gradient inference, optimizer or unused observation is used.',
            'gradient_convention':'Gradient of1/2(target rawchi2+sky prior penalty), multiplied by the same declared dimensionless parameter scales as the joint SVD. For a lower bound, positive gradient meets the KKT sign and is set tozero in the projected violation; for an upper bound, negative gradient meets it. Active tolerance1e-8 in native parameter units.',
            'applications':answers,'original_saved_result_sha256':before,'original_results_byte_identical':before==after,
            'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [Path(__file__).resolve(),ROOT/'code/strong_lensing_joint_scene.py',ROOT/'scripts/run_he0435_joint_scene.py']},
            'interpretation':'A nonzero free point/flux/plane or projected sky gradient is a necessary-condition failure for target-nuisance stationarity even with response fixed. This does not provide the missing total eta gradient because the field response rows are not refitted here. It is not a coverage or adequacy result.',
            'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10}
    run.dump(run.OUT/'target_nuisance_stationarity.json',result)
    for name,a in answers.items():print(name,json.dumps(a['groups']),a['all_target_nuisance_projected_KKT_scaled_gradient_max_abs'])


if __name__=='__main__':main()
