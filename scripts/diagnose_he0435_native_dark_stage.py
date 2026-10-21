"""Bounded frozen-body sensitivity to the nominal preIPC RAPID-dark assignment."""
from pathlib import Path
import hashlib,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_native_completion as run
from apply_he0435_bounded_native_history import native_apply
from stdatamodels.jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from stcal.ramp_fitting import ramp_fit


def main():
    dest=run.OUT/'dark_stage_support';dest.mkdir(exist_ok=True)
    path=run.OUT/'held_target/solution.npz';x=np.load(path)['full_parameters']
    problem,_,_,_,_,_=run.build('all_data');e4=problem.exposure_order[-1]
    problem.contexts=[c for c in problem.contexts if c['original']['exposure']==e4]
    problem.evaluate(x);rows=[];arrays={}
    with datamodels.RampModel(run.BASE/'native_cache'/e4.replace('_cal.fits','_uncal.fits')) as source:
        template=ramp_fit.create_ramp_fit_class(source,'OLS_C',dqflags.pixel,False)
    for ctx,(state,offset,_) in zip(problem.contexts,problem.last_states):
        h=ctx['history'];p=ctx['patch'];cal=ctx['calibration'];sel=ctx['selected'];rate=problem.optical(x,ctx)[0]
        before_cont=problem.native_ols(ctx,state)[0]*ctx['native_SCI_factor']+offset
        old=cal.dark.copy()
        try:
            # This is a declared mean-stage sensitivity, not the assertion that
            # true physical dark charge is zero. Actual native dark subtraction
            # and its source reference are retained in both applications.
            cal.dark=np.zeros_like(old);alternative=cal.evaluate(rate,x[-1])
        finally:cal.dark=old
        after_cont=problem.native_ols(ctx,alternative)[0]*ctx['native_SCI_factor']+offset
        before=native_apply(template,h,p,state['groups_DN'],state['zeroframe_DN']).ravel()[sel]*ctx['native_SCI_factor']+offset
        after=native_apply(template,h,p,alternative['groups_DN'],alternative['zeroframe_DN']).ravel()[sel]*ctx['native_SCI_factor']+offset
        data=ctx['original']['data'];err=ctx['original']['error'];shape=(65,65)
        reference=h['dark_RAPID42_DN'][p];direct=np.zeros(shape,bool);direct.ravel()[sel]=True
        halo=np.zeros(shape,bool);halo[16:49,16:49]=True
        row={'name':ctx['name'],'field':ctx['field'],
            'nonzero_RAPID_reference_pixels_direct_ROI':int(np.any(reference!=0,axis=0)[direct].sum()),
            'nonzero_RAPID_reference_pixels_IPC_support33':int(np.any(reference!=0,axis=0)[halo].sum()),
            'maximum_abs_RAPID_reference_DN_in_support33':float(np.max(abs(reference[:,halo]))),
            'max_abs_continuous_native_mean_change_SCI':float(np.max(abs(after_cont-before_cont))),
            'max_abs_original_float32_native_mean_change_SCI':float(np.max(abs(after-before))),
            'nominal_original_ERR_Q':float(np.sum(((before-data)/err)**2)),
            'zero_preIPC_dark_assignment_original_ERR_Q':float(np.sum(((after-data)/err)**2)),
            'native_dark_subtraction_unchanged':True,'G0_intercept_recomputed':True}
        rows.append(row);arrays[ctx['name']+'_native_group_mean_change_DN']=alternative['groups_DN']-state['groups_DN']
        arrays[ctx['name']+'_continuous_ROI_mean_change_SCI']=after_cont-before_cont;arrays[ctx['name']+'_original_native_ROI_mean_change_SCI']=after-before
    np.savez_compressed(dest/'sensitivity_maps.npz',**arrays)
    targets=[r for r in rows if not r['field']]
    result={'scope':'Frozen trained shared state and independently measured E4 point state; nominal RAPID reference assigned before inverseL/IPC versus no preIPC dark-mean assignment. Native post-linearization dark subtraction retained exactly. A bounded stage-sensitivity diagnostic, not a claim of zero physical dark charge or a calibrated alternative detector law.',
        'equation':'Compare Lnative[H mean Linv(beta+mu/gain+dark_read)]-dark_native with Lnative[H mean Linv(beta+mu/gain)]-dark_native, re-solving the same G0 beta equation in each case.',
        'support_scope':'Direct17ROI has zero RAPID dark reference, but actual33IPC-support halo contains isolated nonzero entries in A/C/FIELD. Their propagation is explicitly evaluated, not excluded by a local assertion.',
        'rows':rows,'held_target_Q_nominal':sum(r['nominal_original_ERR_Q'] for r in targets),'held_target_Q_alternative':sum(r['zero_preIPC_dark_assignment_original_ERR_Q'] for r in targets),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'held_solution_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'remaining_scope':'Does not identify a general physical dark/overflow/counting law or dark stochastic covariance. The nominal dark-stage ambiguity is not a demonstrated cure of the broad held profile failure.',
        'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,'error_rescaling':False}
    run.prior.dump(dest/'result.json',result);print(rows,flush=True)
    print('HELD_Q',result['held_target_Q_nominal'],result['held_target_Q_alternative'],flush=True)


if __name__=='__main__':main()
