"""Competing actual-trajectory contractions: common drift versus collection."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'));sys.path.insert(0,str(ROOT/'scripts'))
from strong_lensing_joint_scene import conservative_charge_laplacian,svd_fit
from run_he0435_joint_scene import dump
BASE=ROOT/'artifacts/strong_lensing_response_completion_20261021'


def main():
    arrays=[];inputs={};out=BASE/'temporal_operators';out.mkdir(exist_ok=True)
    paths=sorted((ROOT/'data/strong_lensing_response_completion_20261021/histories').glob('*_histories.npz'))
    for exposure,hpath in enumerate(paths):
        h=np.load(hpath);stem=hpath.name.replace('_histories.npz','')
        ppath=BASE/'common_early/pilot_means'/f'training_plus_held_{stem}_pilot.npz';pilot=np.load(ppath)
        cpath=BASE/'temporal_curvature'/f'{stem}_curvature.npz';curves=np.load(cpath)
        for pp in [hpath,ppath,cpath]:inputs[str(pp.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(pp.read_bytes()).hexdigest()
        for k,label in enumerate(h['labels']):
            key=f'{label}_minus_data_anchored_Jensen';observed=curves[key+'_coefficients_beta_rate_curvature'][2];variance=curves[key+'_curvature_variance'];stencil=curves[f'{label}_all9valid_neighbour_stencil']
            I=pilot['positive_point_plus_coherent_sky_SCI'][k];gain=h['gain_e_per_DN'][k];conversion=gain*h['flat'][k]/float(h['photmjsr']);rate=I*conversion
            dy,dx=np.gradient(I);diffusion=conservative_charge_laplacian(I)
            basis=np.stack([.5*conservative_charge_laplacian(rate**2)/gain,
                            .5*conversion*dx/gain,.5*conversion*dy/gain,.5*conversion*diffusion/gain],axis=-1)
            x=h['x_pixel'][k];y=h['y_pixel'][k];center=h['prior_guess_native_xy'][k];r=np.hypot(x-center[0],y-center[1])
            selected=stencil&(r>=3)&(r<8.6)&(abs(x-round(center[0]))<=8)&(abs(y-round(center[1]))<=8)&((h['archived_cal_DQ'][k]&1)==0)&np.isfinite(observed)&np.isfinite(variance)&(variance>0)&np.isfinite(basis).all(axis=-1)
            a=np.zeros((selected.sum(),13));a[:,0]=basis[selected,0];a[:,1+3*exposure:4+3*exposure]=basis[selected,1:]
            for i,(z,v,row) in enumerate(zip(observed[selected],variance[selected],a)):
                arrays.append((exposure,str(label),z,v,row))
    exp=np.array([a[0] for a in arrays]);labels=np.array([a[1] for a in arrays]);y=np.array([a[2] for a in arrays]);variance=np.array([a[3] for a in arrays]);a=np.array([a[4] for a in arrays]);train=(exp<3)|(labels=='FIELD');held=~train
    result={};zero={'train_Q':float(np.sum(y[train]**2/variance[train])),'held_Q':float(np.sum(y[held]**2/variance[held])),'held_raw_curvature_RMS_DN_s2':float(np.sqrt(np.mean(y[held]**2)))}
    for name,columns in [('collection_only',[0]),('common_linear_drift',list(range(1,13))),('joint_collection_and_drift',list(range(13)))]:
        columns=np.array(columns);weighted=a[train][:,columns]/np.sqrt(variance[train,None]);rhs=y[train]/np.sqrt(variance[train]);norm=np.linalg.norm(weighted,axis=0);scale=np.where(norm>0,1/norm,1.)
        lower=np.full(len(columns),-np.inf);upper=np.full(len(columns),np.inf)
        if 0 in columns:lower[np.flatnonzero(columns==0)[0]]=0
        fit=lsq_linear(weighted*scale,rhs,bounds=(lower,upper),tol=1e-12,lsmr_tol=1e-12,max_iter=1000)
        theta=np.zeros(13);theta[columns]=fit.x*scale;prediction=a@theta;res=y-prediction;info=svd_fit(weighted,scales=scale)
        result[name]={'parameters_order':['shared_kappa_e_inverse']+[f'E{i+1}_{s}' for i in range(4) for s in ['vx_pixel_s','vy_pixel_s','signed_diffusion_pixel2_s']],
             'parameters':theta,'rank':info['rank'],'columns':len(columns),'scaled_singular_values':info['singular_values'],'success':bool(fit.success),'optimality':float(fit.optimality),
             'train_Q':float(np.sum(res[train]**2/variance[train])),'held_Q':float(np.sum(res[held]**2/variance[held])),
             'held_raw_curvature_RMS_DN_s2':float(np.sqrt(np.mean(res[held]**2))),'held_raw_curvature_RMS_fraction_of_zero':float(np.sqrt(np.mean(res[held]**2))/zero['held_raw_curvature_RMS_DN_s2']),
             'train_residual_RMS_DN_s2':float(np.sqrt(np.mean(res[train]**2))),'positive_kappa_bound_active':bool(0 in columns and theta[0]<1e-18)}
        np.savez_compressed(out/f'{name}.npz',prediction_curvature_DN_s2=prediction,residual_curvature_DN_s2=res,observed_curvature_DN_s2=y,variance_DN2_s4=variance,train=train,exposure=exp,labels=labels,basis=a)
    dump(out/'result.json',{'scope':'Actual9group temporalcurvature contrasted against physically motivated commontranslation/linear-diffusion versus originalbrightness-squared neighbourcontraction. E1–3targets+all4FIELD train; E4targets held. This is consumeddevelopment andconditionalcurvaturemetric, not fullnativeforward or calibratedphysicalvalidation.',
         'equations':'a_DN=κ/(2g)L[(c I)^2]+(c/2g)[vx d_x I+vy d_y I+D L I], c=gain*flat/PHOTMJSR. Detector conversion c stays fixed under optical translation/diffusion; it is not differentiated as sky. One sharedpositiveκ plus3signedopticalrate derivatives perexposure. Only all9calibrated neighbourstencils3–8.6pixel used.',
         'parameter_scope':'SignedD is a local commonwidth derivative, not necessarily positive diffusion. The linear-in-time optical expansion and first-order collection expansion require separate positivity/small-change tests before full-intervalforward use. Staticfieldshape maydifferbysource; no newPSFfamily orlensparameters fitted.',
         'missing_noise_relation':'Frozenconditionalκ0photon/CDSmarginal curvaturecovariance omits unresolvedrawscience/referencecrosscovariance; Q and coefficienterrors are not calibrated significance. Read/reference/linearity temporal drift can mimic some contraction; no uniquephysicalcause claimed.',
         'pixels_training':int(train.sum()),'pixels_held':int(held.sum()),'zero_curvature_predictor':zero,'models':result,'input_sha256':inputs,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    print(json.dumps({'zero':zero,'models':{k:{q:v[q] for q in ['parameters','train_Q','held_Q','held_raw_curvature_RMS_fraction_of_zero','rank']} for k,v in result.items()}},default=lambda x:x.tolist(),indent=2))


if __name__=='__main__':main()
