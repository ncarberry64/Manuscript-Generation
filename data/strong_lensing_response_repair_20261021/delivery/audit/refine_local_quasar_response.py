"""Converge the same local-response objective with an exact dense trust solve.

LSMR reached its outer-iteration cap. No response family, prior, error, window,
held exposure or physical constraint changes. Initializers are each fold's own
training-only iterates, never iterates trained on the held exposure.
"""
from pathlib import Path
import sys,json,time
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from local_quasar_response import get_targets,read_primary_fits,ART,FineBasis,fit,serial,profile_information,residuals

def main():
    records=get_targets();names=sorted(set(r['exposure'] for r in records))
    grid,_=read_primary_fits(ART/'standard_psf/jw01198007001_02101_00001_nrcb3_FIELD_standard_epsf.fits')
    basis=FineBasis(grid+max(0.,-float(grid.min()))*1.1+1e-8);basis.iterative_solver=False
    initial=json.loads((HERE/'local_quasar_results.json').read_text())
    report={k:v for k,v in initial.items() if k not in ['folds','all_image_fit','all_rank','all_pixel_results','held_exposure_chi2_total','held_exposure_nominal_dof_conditional_on_shape']}
    report['optimizer_refinement']='Exact dense trust-region continuation from each same-fold LSMR iterate; unchanged objective.';report['folds']=[]
    for hold,name in enumerate(names):
        tr=[r for r in records if r['exposure']!=name];te=[r for r in records if r['exposure']==name]
        before=np.load(HERE/f'local_quasar_fold{hold+1}.npz');start=time.time()
        train=fit(tr,basis,initial=before['parameters'],max_nfev=25)
        test=fit(te,basis,a_fixed=train['shape'],max_nfev=70)
        row={'held_exposure':name,'train':serial(train),'test':serial(test),'rank':profile_information(tr,basis,train),
             'held_pixel_results':residuals(te,basis,test,free_shape=False),'validation_pass':False}
        report['folds'].append(row)
        np.savez(HERE/f'local_quasar_fold{hold+1}_refined.npz',shape=train['shape'],parameters=train['parameters'],test_parameters=test['parameters'])
        (HERE/'local_quasar_refined_results.json').write_text(json.dumps(report,indent=2))
        print(json.dumps({'fold':hold+1,'seconds':time.time()-start,'training_success':train['success'],'training_nfev':train['nfev'],
              'held_chi2_per_dof':[p['chi2']/p['dof_conditional_on_shape'] for p in test['fits']]}),flush=True)
    before=np.load(HERE/'local_quasar_shape.npz');allfit=fit(records,basis,initial=before['parameters'],max_nfev=30)
    report['all_image_fit']=serial(allfit);report['all_rank']=profile_information(records,basis,allfit)
    report['all_pixel_results']=residuals(records,basis,allfit)
    report['held_exposure_chi2_total']=sum(r['chi2'] for f in report['folds'] for r in f['test']['fits'])
    report['held_exposure_nominal_dof_conditional_on_shape']=sum(r['dof_conditional_on_shape'] for f in report['folds'] for r in f['test']['fits'])
    np.savez(HERE/'local_quasar_shape_refined.npz',shape=allfit['shape'],parameters=allfit['parameters'],nodes=basis.nodes,gauge=basis.gauge,width=basis.width)
    (HERE/'local_quasar_refined_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'held_total':report['held_exposure_chi2_total'],'nominal_dof':report['held_exposure_nominal_dof_conditional_on_shape'],
                     'all_success':allfit['success'],'all_nfev':allfit['nfev'],'all_rank':report['all_rank']}),flush=True)

if __name__=='__main__':main()
