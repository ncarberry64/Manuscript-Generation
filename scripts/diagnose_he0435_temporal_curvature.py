"""Actual native-group temporal curvature versus the original neighbour law.

Individual later reads were not downlinked. Only the nine observed two-read
averages are fitted; their eighteen known time supports enter the design/cov.
"""
from pathlib import Path
import hashlib,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'));sys.path.insert(0,str(ROOT/'scripts'))
from strong_lensing_joint_scene import conservative_charge_laplacian
from run_he0435_joint_scene import dump
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021/temporal_curvature'


def fit_curves(data,valid,rate,gain,rnoise,times):
    n,h,w=data.shape;tf=times.mean(axis=1);t2=(times**2).mean(axis=1);scale=times.max()
    design=np.c_[np.ones(n),tf/scale,t2/scale**2]
    mincov=np.minimum(times[:,None,:,None],times[None,:,None,:]).mean(axis=(2,3))
    coef=np.full((3,h,w),np.nan);cov=np.full((3,3,h,w),np.nan);q=np.full((h,w),np.nan);qlin=q.copy()
    bits=np.sum(valid*(1<<np.arange(n))[:,None,None],axis=0)
    for code in np.unique(bits):
        groups=np.flatnonzero(code&(1<<np.arange(n)))
        if len(groups)<3:continue
        select=(bits==code)&np.isfinite(rate)&(rate>0)&(gain>0)&np.isfinite(rnoise)
        if not select.any():continue
        c=(rate[select]/gain[select])[:,None,None]*mincov[np.ix_(groups,groups)][None]
        c+=np.eye(len(groups))[None]*rnoise[select,None,None]**2/4
        x=design[groups];y=data[groups][:,select].T
        invx=np.linalg.solve(c,np.broadcast_to(x,(len(y),*x.shape)));invy=np.linalg.solve(c,y[:,:,None])[:,:,0]
        normal=np.einsum('gi,ngj->nij',x,invx);v=np.linalg.inv(normal);p=np.einsum('nij,gj,ng->ni',v,x,invy)
        residual=y-p@x.T;qq=np.einsum('ni,ni->n',residual,np.linalg.solve(c,residual[:,:,None])[:,:,0])
        xl=x[:,:2];vl=np.linalg.inv(normal[:,:2,:2]);pl=np.einsum('nij,gj,ng->ni',vl,xl,invy)
        rl=y-pl@xl.T;qqlinear=np.einsum('ni,ni->n',rl,np.linalg.solve(c,rl[:,:,None])[:,:,0])
        conversion=np.array([1,scale,scale**2]);coef[:,select]=(p/conversion).T
        cov[:,:,select]=(v/conversion[None,:,None]/conversion[None,None,:]).transpose(1,2,0)
        q[select]=qq;qlin[select]=qqlinear
    return coef,cov,q,qlin


def main():
    OUT.mkdir(exist_ok=True);rows=[];inputs={}
    for path in sorted((ROOT/'data/strong_lensing_response_completion_20261021/histories').glob('*_histories.npz')):
        a=np.load(path);inputs[str(path.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(path.read_bytes()).hexdigest()
        biaspath=ROOT/'artifacts/strong_lensing_response_completion_20261021/group_average'/path.name.replace('.npz','_group_bias.npz')
        inputs[str(biaspath.relative_to(ROOT)).replace('\\','/')]=hashlib.sha256(biaspath.read_bytes()).hexdigest()
        bias=np.load(biaspath)['correction_DN'];products={}
        for k,label in enumerate(a['labels']):
            data=a['corrected_group_DN'][k].astype(float);valid=a['valid_group_DNU_SAT_only'][k].copy()
            # JUMP delimits a new intercept. A single shared-beta trajectory is
            # inapplicable across that event, so reject the whole pixel here.
            valid[:,np.any((a['groupdq'][k]&4)!=0,axis=0)]=False
            native=a['native_rate_DN_s'][k].astype(float);gain=a['gain_e_per_DN'][k].astype(float);rn=a['readnoise_CDS_DN'][k].astype(float)
            t=a['native_read_times_s'];rate=native*float(a['group_time_header_s'])/(t[1].mean()-t[0].mean())
            x=a['x_pixel'][k];y=a['y_pixel'][k];center=a['prior_guess_native_xy'][k];radius=np.hypot(x-center[0],y-center[1])
            roi=(abs(x-round(center[0]))<=8)&(abs(y-round(center[1]))<=8)&((a['archived_cal_DQ'][k]&1)==0)
            # The original unsaturated neighbour contraction, in DN/s² per κ_e:
            # (1/2g) L[(g*lambda_DN)^2]. It is tested only where the full local
            # stencil remains calibrated through all nine groups.
            eigenrate=gain*rate;basis=.5*conservative_charge_laplacian(eigenrate**2)/gain
            full=valid.all(axis=0)
            stencil=full.copy()
            stencil[1:] &=full[:-1];stencil[:-1]&=full[1:];stencil[:,1:]&=full[:,:-1];stencil[:,:-1]&=full[:,1:]
            for variant,values in [('native_corrected',data),('minus_data_anchored_Jensen',data-bias[k])]:
                p,v,q,ql=fit_curves(values,valid,rate,gain,rn,t);products[f'{label}_{variant}_coefficients_beta_rate_curvature']=p
                products[f'{label}_{variant}_curvature_variance']=v[2,2]
                products[f'{label}_{variant}_Q_quadratic']=q;products[f'{label}_{variant}_Q_linear']=ql
                radial=[]
                for lower,upper in [(0,3),(3,6),(6,8.6)]:
                    region=roi&(radius>=lower)&(radius<upper)&np.isfinite(p[2])
                    safe=region&stencil&np.isfinite(v[2,2])&(v[2,2]>0)
                    if safe.any():
                        weight=1/v[2,2,safe];b=basis[safe];curv=p[2,safe];den=float(np.sum(weight*b*b))
                        estimate=float(np.sum(weight*b*curv)/den) if den>0 else None
                        correlation=float(np.corrcoef(b,curv)[0,1]) if len(b)>2 and np.std(b)>0 and np.std(curv)>0 else None
                        uncertainty=float(1/np.sqrt(den)) if den>0 else None
                    else:estimate=correlation=uncertainty=None
                    fixed=region&full
                    # Remove each pixel's fitted static intercept+linear rate;
                    # the remaining temporal shape cannot be a static PSF error.
                    trajectory=[]
                    for g in range(9):
                        observed=values[g]-p[0]-p[1]*t[g].mean()
                        trajectory.append(float(np.mean(observed[fixed])) if fixed.any() else None)
                    radial.append({'radius_pixel':[lower,upper],'pixels_with_fit':int(region.sum()),'fixed_all9valid_pixels':int(fixed.sum()),
                         'median_quadratic_DN_s2':float(np.nanmedian(p[2,region])) if region.any() else None,
                         'raw_deltaQ_linear_to_quadratic':float(np.nansum(ql[region]-q[region])),
                         'fully_unsaturated_neighbour_stencil_pixels':int(safe.sum()),'conditional_leading_kappa_e_inverse':estimate,
                         'conditional_kappa_standard_error':uncertainty,'curvature_vs_original_neighbour_basis_correlation':correlation,
                         'mean_group_trajectory_after_static_intercept_and_rate_DN':trajectory})
                rows.append({'exposure':path.name,'label':str(label),'variant':variant,'radial_results':radial})
            products[f'{label}_original_neighbour_basis_DN_s2_per_kappa']=basis;products[f'{label}_all9valid_neighbour_stencil']=stencil
        np.savez_compressed(OUT/path.name.replace('_histories.npz','_curvature.npz'),**products)
    dump(OUT/'result.json',{'scope':'Actual9downlinked two-read-average trajectories on consumed4target+field exposures; no invented individual readvalues or optical/refit. Staticper-pixel intercept/rate profiled out. Actual18time supports enter quadraticmean and conditionalκ0Poisson/CDSmarginal metric.',
       'equations':'y_g=beta+lambda*mean(t_gr)+a*mean(t_gr^2); original neighbour low-charge law a=κ_e/(2gain) L[(gain*lambda_DN)^2]. The neighbour test is restricted to completeunsaturated9group centre+fourneighbours; JUMPpixels excluded.','noise_scope':'Independent post-correction CDS-marginal+conditionalκ0Poisson comparator only, no SIRS addition/no sigmafit. Errorbars/deltaQ are not source-complete significance or validation.',
       'Jensen_scope':'Nativecorrected and subtraction ofroot data-anchored knownnonlinearaveraging diagnostic are shown separately. This is not a joint latentβ/readlikelihood solve.',
       'lambda_basis_scope':'The neighbour basis uses the actual native-rate pilot converted to physicalTFRAME timing, not the independently profiled initial-rate coefficient. Diagnosticκ is not a joint incidentλ/κ optimum or uniquely identified collection coefficient. Staticoptical shape errors alone are linear-in-time and are profiled out, but coefficient-calibration/reference temporal drift can also produce curvature.',
       'rows':rows,'input_sha256':inputs,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    for row in rows:
        if row['variant']=='minus_data_anchored_Jensen':print(row['exposure'][21:26],row['label'],[(r['radius_pixel'],r['conditional_leading_kappa_e_inverse'],r['curvature_vs_original_neighbour_basis_correlation']) for r in row['radial_results']])


if __name__=='__main__':main()
