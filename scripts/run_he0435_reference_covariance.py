"""Execute native-reference relation replay, conditional J and unit noise Gram.

No detector step or fit is rerun. This is NOT a calibrated raw noise estimate:
READNOISE0282 was produced after reference/bias/IPC correction. Its marginal
CDS values do not supply the missing raw spatiotemporal covariance.
"""
from pathlib import Path
import sys,json,hashlib,time
import numpy as np
from stdatamodels.jwst import datamodels
from jwst.refpix.optimized_convolution import make_kernels
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_reference_covariance import conditional_sirs,predicted_native_post,export_sparse
BASE=ROOT/'data/strong_lensing_response_completion_20261021';OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/reference_covariance'
DATA=BASE/'reference_covariance'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)
    kernelpath=BASE/'native_cache/jwst_nircam_sirskernel_0002.asdf'
    with datamodels.SIRSKernelModel(kernelpath) as model:kernels=make_kernels(model,'NRCB3',1.,30)
    rows=[]
    for exposure in range(1,5):
        started=time.time();stem=f'jw01198007001_02101_0000{exposure}_nrcb3';path=BASE/'histories'/f'{stem}_histories.npz';dpath=BASE/'histories'/f'{stem}_native_diagnostics_v2.npz'
        output=DATA/f'{stem}_reference_operator.npz';receipt=OUT/f'{stem}.json';earlyout=DATA/f'{stem}_early_reference_sensitivity.npz'
        if output.exists() and receipt.exists():rows.append(json.loads(receipt.read_text()));continue
        a=np.load(path);d=np.load(dpath);controls=a['pre_refpix_reference_control_DN'];dq=d['reference_control_pre_sirs_pixeldq']
        xy=np.c_[a['x_pixel'].ravel(),a['y_pixel'].ravel()];arrays={'labels':a['labels'],'bounds_xyxy':a['bounds_xyxy'],'kernel_left':np.asarray(kernels[0]),'kernel_right':np.asarray(kernels[1]),'control_xy_native':a['reference_control_native_xy'],'control_pre_sirs_DQ':dq}
        group_records=[];early=[];corexy=np.c_[a['x_pixel'][:,24:41,24:41].ravel(),a['y_pixel'][:,24:41,24:41].ravel()]
        dt=float(np.diff(a['native_read_times_s'].mean(axis=1))[0]);branch=d['linearity_DQ_mapped']
        valid=((a['groupdq'][:,:2,24:41,24:41]&7)==0).all(axis=1)&((a['pixeldq'][:,24:41,24:41]&1)==0)&((d['flat_DQ_mapped'][:,24:41,24:41]&1)==0)
        for group in range(9):
            # NIRDataset.group starts float32. SIRS returns float64 and the same
            # work buffer is reused for subsequent groups, before final SCI cast.
            work=np.float32 if group==0 else np.float64
            o=conditional_sirs(controls[group],a['reference_control_native_xy'],dq,xy,kernels,np.ones(controls.shape[1]),native_work_dtype=work)
            pred=predicted_native_post(a['pre_refpix_group_DN'][:,group].ravel(),xy,o);actual=a['post_refpix_group_DN'][:,group].ravel()
            difference=pred-actual
            prefix=f'group{group}';export_sparse(arrays,prefix+'_T',o['T']);export_sparse(arrays,prefix+'_D',o['D']);arrays[prefix+'_U']=o['U']
            arrays[prefix+'_unit_variance_feature_covariance']=o['covariance'];arrays[prefix+'_unit_variance_feature_factor']=o['factor']
            arrays[prefix+'_side_good_left']=o['side_good_left'];arrays[prefix+'_side_good_right']=o['side_good_right'];arrays[prefix+'_top_bottom_clip_counts']=o['top_bottom_clip_counts']
            if group==0:
                arrays['features_detector_row_amp_parity']=o['features'];arrays['pixel_feature_index_full65']=o['pixel_feature_index'].reshape(5,65,65)
            record={'group':group,'native_pre_post_replay_max_abs_DN':float(np.max(abs(difference))),'native_pre_post_bitwise_equal':bool(np.array_equal(pred,actual)),
                    'native_work_dtype':o['native_work_dtype'],
                    'float64_Jr_vs_native_float32_correction_max_abs_DN':float(np.max(abs(o['linear_operator_reconstruction']-o['native_correction']))),
                    'features':len(o['features']),'unit_variance_factor_rank':o['factor'].shape[1],'minimum_Gram_eigenvalue':float(o['eigenvalues'][0]),
                    'factor_replay_max_abs':float(np.max(abs(o['factor']@o['factor'].T-o['covariance']))),'roundoff_eigenvalue_tolerance':float(o['eigenvalue_roundoff_tolerance'])}
            group_records.append(record)
            if group<2:
                mapping=o['pixel_feature_index'].reshape(5,65,65)[:,24:41,24:41].ravel()
                lp=a['linearity_polynomial_first_derivative'][:,group,24:41,24:41].copy()
                no_lin=(branch[:,24:41,24:41]&1048576)!=0
                lp[no_lin]=1.
                scale=lp.ravel()/dt*float(a['photmjsr'])/a['flat'][:,24:41,24:41].ravel()
                # Corrected G has derivative -Lprime*Jref with respect to raw
                # controls; early rate is (G1-G0)/dt.
                early.append((1 if group==0 else -1)*scale[:,None]*o['factor'][mapping])
                if group==1:
                    factor=np.column_stack(early)
                    if not earlyout.exists():np.savez_compressed(earlyout,unit_white_reference_sensitivity_factor_SCI=factor,labels=a['labels'],native_xy=corexy,valid_early17=valid,interval_s=dt,PHOTMJSR=a['photmjsr'])
                    print('EARLY_UNIT_REFERENCE_SENSITIVITY_READY',earlyout.name,factor.shape,flush=True)
        selected=controls[:,(dq&1)==0].astype(float);centered=selected-selected.mean(axis=0,keepdims=True);centered-=centered.mean(axis=1,keepdims=True)
        arrays['reference_only_centered_temporal_Gram_DN2']=centered@centered.T/max(centered.shape[1]-1,1)
        arrays['reference_readnoise0282_CDS_DN']=d['reference_control_readnoise_CDS_DN']
        np.savez_compressed(output,**arrays)
        row={'status':'EXECUTED_EXACT_CONDITIONAL_SIRS_OPERATOR_AND_UNIT_VARIANCE_SENSITIVITY','source_complete':False,'ready_for_native_ERR_addition':False,
             'operator_file':str(output.relative_to(ROOT)),'operator_sha256':sha(output),'early_sensitivity_file':str(earlyout.relative_to(ROOT)),'early_sensitivity_sha256':sha(earlyout),
             'input_history_sha256':sha(path),'mapped_reference_dq_sha256':sha(dpath),'SIRS_reference_sha256':sha(kernelpath),
             'source_sha256':{p:sha(ROOT/p) for p in ['code/strong_lensing_reference_covariance.py','scripts/run_he0435_reference_covariance.py']},
             'relation':'corrected_prelinearity_science = raw_post_superbias_science - Jref*r; Jref=D+U*T, conditional on actual native rejection sets',
             'factor':'Jref*I*Jref^T with unit variance per independent reference control. This is a conditional operator sensitivity, NOT identified raw ADC noise.',
             'early_map':'Native Lprime, first-two-group difference / exact53.68385s, /flat*PHOTMJSR; all5patch cross terms retained in shared factor columns.',
             'groups':group_records,'all_native_mean_replays_bitwise_equal':all(x['native_pre_post_bitwise_equal'] for x in group_records),
             'reference_only_temporal_Gram':'Descriptive PSD9x9 statistic across nonDNU referencepixels after per-pixel temporal and per-group pixel centering; rank<=8 by construction. Not used as a fully identified covariance or a replacement state kernel.',
             'unresolved_base_covariance':'READNOISE0282 HISTORY states CDS estimated AFTER UofA ncdhas reference/bias/IPC correction, from disjoint pairs. Its marginal variances do not identify rawSigma or justify adding Jref Cref JrefT to nativeERR. Calibration operator/full correlations or raw calibration darks are required.',
             'other_limits':['Rejection masks and native float32 branches conditioned on observed controls; mask uncertainty and quantization derivative omitted.','Kernel/reference calibration covariance is not supplied by reference ERR marginals.','The same reference controls create cross-stamp and row correlations; DΣT cross terms retained.','No target residual was used to fit noise or an error scale.'],
             'physical_admission':False,'native_processing_rerun':False,'elapsed_s':time.time()-started}
        receipt.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n');rows.append(row)
        print('SIRS_OPERATOR_READY',stem,row['all_native_mean_replays_bitwise_equal'],output.stat().st_size,flush=True)
    (OUT/'results.json').write_text(json.dumps({'records':rows,'source_complete':False,'ready_for_native_ERR_addition':False,'physical_admission':False},indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
