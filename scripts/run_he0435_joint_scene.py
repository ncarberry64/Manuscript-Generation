"""Execute consumed F115W coherent-sky and calibration-uncertainty development."""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,hashlib,time
import numpy as np
from scipy.optimize import least_squares
from astropy.io import fits
import asdf

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import strong_lensing_joint_scene as js
from strong_lensing_jwst_rate_mask import jwst_native_rate_mask
from strong_lensing_pixel_geometry import tangent_plane_arcsec,gwcs_centroid,intrinsic_geometry,intrinsic_jacobian

OUT=ROOT/'artifacts/strong_lensing_response_repair_20261021/joint_scene'
DELIVERY=ROOT/'data/strong_lensing_response_repair_20261021/delivery'
OLD=ROOT/'artifacts/strong_lensing_geometry_20261021'


def dump(path,value):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(value,indent=2,default=convert,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def make_record(h,guess,response,name,exposure,patch,flux,radius=8):
    cx,cy=guess;ix,iy=np.rint(guess).astype(int)
    x0,x1=ix-radius,ix+radius+1;y0,y1=iy-radius,iy+radius+1
    data,error,dq=[np.array(h[k].data[y0:y1,x0:x1]) for k in ['SCI','ERR','DQ']]
    valid,receipt=jwst_native_rate_mask(data,error,dq,primary_header=h[0].header)
    yy,xx=np.mgrid[y0:y1,x0:x1]
    return {'name':name,'exposure':exposure,'patch':patch,'guess':np.asarray(guess,float),
            'x':xx[valid].astype(float),'y':yy[valid].astype(float),
            'data':data[valid].astype(float),'error':error[valid].astype(float),
            'response':response,'bounds':[x0,y0,x1,y1],
            'flux_guess':float(flux*response.normalization),'mask_receipt':receipt}


def load_inputs():
    saved=np.load(DELIVERY/'audit/positive_field_shape.npz')
    eta=np.array(saved['shape']);nodes=np.array(saved['nodes']);null=np.array(saved['null']);width=float(saved['width'])
    fieldrows=json.loads((OLD/'standard_psf_convention_diagnostic.json').read_text())['records']
    fieldmap={r['exposure']:r for r in fieldrows}
    records=[];fields=[];wcss={};origin=None;anchors=[]
    for path in sorted((OLD/'jwst_target_diagnostic/exposures').glob('*.json')):
        item=json.loads(path.read_text());native=Path(item['file']);name=native.name
        with asdf.open(OLD/'jwst_target_diagnostic/wcs'/f'{native.stem}_wcs.asdf',lazy_load=False) as a:
            wcss[name]=a['wcs']
        w=wcss[name]
        if origin is None:origin=np.asarray(w(*item['points'][0]['centroid_xy_pixel'])[:2],float)
        field=fieldmap[name]
        grid=fits.getdata(OLD/'standard_psf'/name.replace('_cal.fits','_FIELD_standard_epsf.fits'))
        fr=js.PositiveResponse(grid,nodes,null,width,eta)
        with fits.open(native,memmap=False) as h:
            fields.append(make_record(h,field['reproduced_global_centroid'],fr,name,name,0,field['flux_SCI_pixel_units']))
            for k,p in enumerate(item['points']):
                grid=fits.getdata(ROOT/p['field_psf_path']);response=js.PositiveResponse(grid,nodes,null,width,eta)
                rec=make_record(h,p['centroid_xy_pixel'],response,f'{native.stem}_{p["label"]}',name,k,p['flux_in_SCI_pixel_units'])
                records.append(rec)
                if len(anchors)<4:anchors.append(tangent_plane_arcsec(np.asarray(w(*p['centroid_xy_pixel'])[:2]),origin))
    # Fixed sky knots use only exposure1's prior measurement coordinates;
    # exposure4 targets do not select the held-out experiment's sky basis.
    for r in records:
        w=wcss[r['exposure']]
        def sky(x,y,w=w):
            ra,dec=w(x,y)[:2]
            return tangent_plane_arcsec(np.stack([ra,dec],axis=-1),origin)
        js.build_scene_projection(r,sky,np.asarray(anchors[r['patch']]))
    for e in wcss:
        center=np.mean([r['guess'] for r in records if r['exposure']==e],axis=0)
        for r in records:
            if r['exposure']==e:r['detector_background']=np.c_[np.ones(len(r['x'])),(r['x']-center[0])/100,(r['y']-center[1])/100]
    return records,fields,eta,wcss,origin


def field_information(fields,eta):
    blocks=[];residuals=[];rows=[]
    for r in fields:
        resp=r['response'];resp.update(eta);x,y=r['x'],r['y'];origin=r['guess']
        def application(p):
            q,dx,dy=resp.values(x-p[0],y-p[1],True);f=np.exp(p[2])
            background=np.c_[np.ones(len(x)),x-origin[0],y-origin[1]]
            mean=f*q+background@p[3:]
            j=np.c_[-f*dx,-f*dy,f*q,background]/r['error'][:,None]
            return (mean-r['data'])/r['error'],j
        p0=np.r_[origin,np.log(r['flux_guess']),np.median(r['data']),0.,0.]
        lo=np.r_[origin-.7,-5.,[-np.inf]*3];hi=np.r_[origin+.7,20.,[np.inf]*3]
        fit=least_squares(lambda p:application(p)[0],p0,jac=lambda p:application(p)[1],bounds=(lo,hi),x_scale='jac',max_nfev=70,ftol=2e-9,xtol=2e-9,gtol=1e-7)
        resid,j=application(fit.x)
        cal=np.exp(fit.x[2])*resp.cal_values(x-fit.x[0],y-fit.x[1])/r['error'][:,None]
        free=fit.active_mask==0
        prof,rr=js.profile_columns(j[:,free],cal,resid,scales=np.asarray([1,1,1,10,10,10])[free])
        blocks.append(prof);residuals.append(rr)
        rows.append({'name':r['name'],'chi2':float(resid@resid),'pixels':len(resid),'success':bool(fit.success),
                     'centroid':fit.x[:2].tolist(),'active_bound_mask':fit.active_mask.tolist(),
                     'optimality':float(fit.optimality),'unconstrained_nuisance_gradient_max':float(np.max(abs(j.T@resid))),
                     'profiled_residual_norm_squared':float(rr@rr)})
    a=np.vstack(blocks);r=np.concatenate(residuals)
    full=np.vstack([a,np.eye(len(eta))/.2]);fr=np.r_[r,eta/.2]
    return {'matrix':full,'residual':fr,'covariance':js.svd_fit(full,scales=np.full(len(eta),.2))['covariance'],
            'profiled_data_matrix':a,'records':rows,'raw_field_chi2':sum(x['chi2'] for x in rows),
            'profiled_field_chi2':float(r@r),'nuisance_active_bound_count':sum(np.count_nonzero(x['active_bound_mask']) for x in rows),
            'field_eta_gradient_scaled_max':float(np.max(abs(full.T@fr)*.2)),
            'data_rank':js.svd_fit(a,scales=np.full(len(eta),.2))['rank']}


def set_eta(records,eta):
    for r in records:r['response'].update(eta);js.refresh_scene_projection(r)


def active_svd_step(j,r,p,scales,lo,hi):
    free=np.ones(j.shape[1],bool);step=np.zeros(j.shape[1])
    for iteration in range(8):
        solved=js.svd_fit(j[:,free],-r-j[:,~free]@step[~free],scales=scales[free])
        step[free]=solved['solution']
        violation=free&((p+step<lo)|(p+step>hi))
        if not violation.any():break
        step[violation]=np.clip(p[violation]+step[violation],lo[violation],hi[violation])-p[violation]
        free[violation]=False
    return step,free


def make_joint(problem,p,eta,field):
    r,a,b,rows=problem.evaluate(p,True);rp,ap=problem.prior(p);n=problem.size;d=len(eta)
    j=np.vstack([np.c_[a,b],np.c_[ap,np.zeros((len(rp),d))],np.c_[np.zeros((len(field['residual']),n)),field['matrix']]])
    residual=np.r_[r,rp,field['residual']]
    return j,residual,r,a,b,rows


def measurement_covariances(problem,p,eta,field):
    r,a,b,rows=problem.evaluate(p,True);rp,ap=problem.prior(p)
    conditional,cut,sensitivity,info=js.propagated_cut_covariance(a,b,field['covariance'],ap,problem.scales)
    take=np.array([3*i+j for i in range(problem.n) for j in [0,1]])
    moments=np.vstack([rec['response'].moment_derivative for rec in problem.records])
    physical_sensitivity=sensitivity[take]+moments
    physical_cut=conditional[np.ix_(take,take)]+physical_sensitivity@field['covariance']@physical_sensitivity.T
    return {'conditional_parameter_covariance':conditional,'cut_parameter_covariance':cut,
            'cut_raw_centroid_covariance_pixel2':cut[np.ix_(take,take)],
            'cut_moment_centroid_covariance_pixel2':physical_cut,
            'conditional_centroid_covariance_pixel2':conditional[np.ix_(take,take)],
            'response_sensitivity_centroid_pixel_per_eta':sensitivity[take],
            'parameter_rank':info['rank'],'response_gradient_scaled_max':float(np.max(abs(b.T@r)*.2))}


def geometry_from_centroids(problem,p,cov,wcss,origin,mode='moment'):
    n=problem.n;positions=[];j=np.zeros((2*n,2*n))
    for i,rec in enumerate(problem.records):
        shift=rec['response'].moment if mode=='moment' else (rec['response'].moment-rec['response'].base_moment if mode=='baseline_relative_moment' else 0.)
        xy=p[3*i:3*i+2]+shift
        sky,_,jac=gwcs_centroid(wcss[rec['exposure']],xy,np.eye(2),origin)
        positions.append(sky);j[2*i:2*i+2,2*i:2*i+2]=jac
    skycov=j@cov@j.T;positions=np.array(positions)
    q=[];t=np.zeros((5*len(problem.exposures),2*n))
    for k,e in enumerate(problem.exposures):
        take=[i for i,r in enumerate(problem.records) if r['exposure']==e]
        pts=positions[take];q.extend(intrinsic_geometry(pts))
        cols=np.array([2*i+d for i in take for d in [0,1]])
        t[np.ix_(np.arange(5*k,5*k+5),cols)]=intrinsic_jacobian(pts)
    q=np.array(q);c=t@skycov@t.T
    if len(problem.exposures)>1:
        weight=np.linalg.inv(c);b=np.tile(np.eye(5),(len(problem.exposures),1));mean_cov=np.linalg.inv(b.T@weight@b)
        mean=mean_cov@b.T@weight@q;res=q-b@mean
        scatter={'Q':float(res@weight@res),'dof':len(q)-5,'scope':'Formal conditional correlated-centroid diagnostic only; inadequate mean/no physical validation.'}
    else:mean=q;mean_cov=c;scatter=None
    return {'tangent_positions_arcsec':positions,'intrinsic_by_exposure_arcsec':q.reshape(-1,5),
            'centroid_convention':{'moment':'Raw fitted center plus FULL finite-grid response moment; a declared finite-grid convention, not established physical centroid truth.',
                                   'raw':'Raw point-response nuisance center.',
                                   'baseline_relative_moment':'Raw fitted center plus response finite-grid moment minus the positive-base finite-grid moment, matching the delivered relative-moment convention.'}[mode],
            'FULL_intrinsic_covariance_arcsec2':c,'mean_intrinsic_arcsec':mean,'mean_covariance_arcsec2':mean_cov,
            'between_exposure_scatter':scatter,'response_and_sky_correlations_retained':True,'physical_validation_pass':False}


def apply(problem,fields,eta0,wcss,origin,label):
    started=time.perf_counter();set_eta(problem.records,eta0)
    field=field_information(fields,eta0)
    print(label,'fielddata_rank',field['data_rank'],flush=True)
    fit=problem.fit(max_nfev=90);p=fit.x.copy()
    print(label,'baselinefit',fit.nfev,float(fit.fun@fit.fun),flush=True)
    cut=measurement_covariances(problem,p,eta0,field)
    cutgeo=geometry_from_centroids(problem,p,cut['cut_moment_centroid_covariance_pixel2'],wcss,origin)
    cutgeoraw=geometry_from_centroids(problem,p,cut['cut_raw_centroid_covariance_pixel2'],wcss,origin,'raw')
    cutgeorelative=geometry_from_centroids(problem,p,cut['cut_moment_centroid_covariance_pixel2'],wcss,origin,'baseline_relative_moment')
    j,r,raw,a,b,rows=make_joint(problem,p,eta0,field)
    scales=np.r_[problem.scales,np.full(len(eta0),.2)]
    joint0=js.svd_fit(j,scales=scales)
    lo,hi=problem.bounds();allp=np.r_[p,eta0]
    step,free=active_svd_step(j,r,allp,scales,np.r_[lo,np.full(len(eta0),-np.inf)],np.r_[hi,np.full(len(eta0),np.inf)])
    # The profiled residual is used in the GN rows. Compare actual nonlinear
    # line-search trials against the SAME RAW field likelihood, not its
    # projected residual (equal only at an exact interior nuisance optimum).
    rp=problem.prior(p)[0]
    base_objective=float(raw@raw+rp@rp+field['raw_field_chi2']+np.sum((eta0/.2)**2))
    trials=[];accepted=0.;newfield=field
    for factor in [1.,.5,.25,.125,.0625,.03125]:
        trial=allp+factor*step;pp=trial[:problem.size];eta=trial[problem.size:]
        try:
            set_eta(problem.records,eta);fnew=field_information(fields,eta)
            rr=problem.evaluate(pp)[0];rp=problem.prior(pp)[0]
            objective=float(rr@rr+rp@rp+fnew['raw_field_chi2']+np.sum((eta/.2)**2))
        except (ValueError,FloatingPointError):objective=1e100
        trials.append({'factor':factor,'original_nonlinear_objective':objective})
        if objective<=base_objective:
            accepted=factor;pnew=pp;etanew=eta;newfield=fnew;break
    if not accepted:
        pnew=p;etanew=eta0.copy();set_eta(problem.records,eta0);newfield=field_information(fields,eta0)
    jj,rr,rawnew,aa,bb,newrows=make_joint(problem,pnew,etanew,newfield)
    joint=js.svd_fit(jj,scales=scales)
    n=problem.size;take=np.array([3*i+d for i in range(problem.n) for d in [0,1]])
    transform=np.zeros((2*problem.n,n+len(eta0)))
    transform[np.arange(len(take)),take]=1.
    transform[:,n:]=np.vstack([rec['response'].moment_derivative for rec in problem.records])
    cov=transform@joint['covariance']@transform.T
    geom=geometry_from_centroids(problem,pnew,cov,wcss,origin)
    rawcov=joint['covariance'][np.ix_(take,take)]
    null_observable=(transform*scales)@joint['null_vectors_dimensionless'].T
    lower,upper=problem.bounds();near_lower=np.isfinite(lower)&np.isclose(pnew,lower,rtol=0,atol=1e-8)
    near_upper=np.isfinite(upper)&np.isclose(pnew,upper,rtol=0,atol=1e-8)
    actual_updated_objective=float(rawnew@rawnew+problem.prior(pnew)[0]@problem.prior(pnew)[0]+newfield['raw_field_chi2']+np.sum((etanew/.2)**2))
    result={'scope':'Consumed development only; fixed delivered response baseline plus ONE actually executed coupled Gauss-Newton continuation, not a converged nonlinear posterior or validation.',
       'baseline_mean_fit_success':bool(fit.success),'baseline_mean_fit_nfev':int(fit.nfev),
       'baseline_mean_fit_optimality':float(fit.optimality),'baseline_active_bound_mask':fit.active_mask,
       'baseline_known_ERR_fit':js.failure_summary(raw,a),'baseline_records':rows,
       'cut_covariance_scope':'Field-only/cut uncertainty propagation at unchanged delivered eta0. Includes all shared response modes plus profiled coherent sky/point flux/background. No target-conditioned response information used in this baseline covariance.',
       'cut_covariance':cut,'cut_geometry':cutgeo,
       'cut_geometry_raw_center':cutgeoraw,
       'cut_geometry_baseline_relative_moment':cutgeorelative,
       'field_calibration_data_rank':field['data_rank'],'field_response_coefficients':len(eta0),
       'field_prior_sd':.2,'field_eta_gradient_at_delivered_mean_scaled_max':field['field_eta_gradient_scaled_max'],
       'field_at_delivered_mean':{k:field[k] for k in ['records','raw_field_chi2','profiled_field_chi2','nuisance_active_bound_count']},
       'coupled_GN':{'accepted_factor':accepted,'line_search':trials,'baseline_joint_objective':base_objective,
          'updated_joint_objective':actual_updated_objective,'updated_GN_profiled_row_norm_squared':float(rr@rr),
          'delta_eta_prior_norm':float(np.linalg.norm((etanew-eta0)/.2)),
          'gradient_dimensionless_max_after_step':float(np.max(abs(jj.T@rr)*scales)),
          'original_nonlinear_known_ERR_fit':js.failure_summary(rawnew,aa),'records':newrows,
          'joint_linearized_rank':joint['rank'],'joint_parameters':len(scales),'svd_rcond':js.SVD_RCOND,
          'joint_full_column_rank':joint['full_column_rank'],'centroid_observable_null_projection_max':float(np.max(abs(null_observable))) if null_observable.size else 0.,
          'active_bound_count':int(near_lower.sum()+near_upper.sum()),'active_bound_parameter_indices':np.flatnonzero(near_lower|near_upper),
          'GN_step_clamped_parameter_indices':np.flatnonzero(~free),'field_nuisance_diagnostics':newfield['records'],
          'column_scales':scales,'covariance_scope':'UNCONSTRAINED local curvature covariance of the specified linearized joint field+target+priors problem at the actual provisional iterate. Positive-coefficient active bounds are recorded but this is neither a truncated posterior nor a tangent-cone covariance. Convergence, mean adequacy and coverage are not established.',
          'FULL_moment_centroid_covariance_pixel2':cov,'FULL_raw_centroid_covariance_pixel2':rawcov,'geometry':geom,
          'geometry_raw_center':geometry_from_centroids(problem,pnew,rawcov,wcss,origin,'raw'),
          'geometry_baseline_relative_moment':geometry_from_centroids(problem,pnew,cov,wcss,origin,'baseline_relative_moment')},
       'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10,'seconds':time.perf_counter()-started}
    dump(OUT/f'{label}.json',result)
    np.savez(OUT/f'{label}_linearized_arrays.npz',eta_initial=eta0,eta_updated=etanew,parameters_initial=p,parameters_updated=pnew,
             joint_covariance=joint['covariance'],joint_column_scales=scales,
             baseline_cut_centroid_covariance=cut['cut_moment_centroid_covariance_pixel2'],joint_centroid_covariance=cov)
    print(label,'done',result['coupled_GN']['original_nonlinear_known_ERR_fit']['chi2'],accepted,flush=True)
    return pnew,etanew,joint['covariance'],result


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    marker=OUT/'execution_started.json'
    if marker.exists():
        if '--resume-development' not in sys.argv:raise RuntimeError('Prior development attempt exists; explicit --resume-development required')
        if not (OUT/'runtime_attempt_1.json').exists():
            dump(OUT/'runtime_attempt_1.json',{'status':'FAILED_BEFORE_FIRST_FIT','reason':'Native unsigned DQ BZERO scaling cannot use Astropy memmap=True; changed only file access to memmap=False. No response/scene parameters executed or changed.','resumed_utc':datetime.now(timezone.utc).isoformat()})
    else:
        with marker.open('x',encoding='utf-8') as stream:json.dump({'started_utc':datetime.now(timezone.utc).isoformat(),'no_reserved_data':True},stream)
    records,fields,eta0,wcss,origin=load_inputs()
    allproblem=js.SceneProblem(records)
    apply(allproblem,fields,eta0,wcss,origin,'all_data')
    train=[r for r in records if r['exposure']!=allproblem.exposures[-1]]
    held=[r for r in records if r['exposure']==allproblem.exposures[-1]]
    training=js.SceneProblem(train)
    p,eta,vcov,result=apply(training,fields,eta0,wcss,origin,'held_exposure_training')
    set_eta(held,eta);testing=js.SceneProblem(held,scene_fixed=training.sky(p));heldfit=testing.fit(max_nfev=70)
    rr,a,b,rows=testing.evaluate(heldfit.x,True)
    # Propagate both trained sky and calibration modes, including their cross block.
    scene_j=np.zeros((len(rr),36));start=0
    for rec in held:
        count=len(rec['x']);scene_j[start:start+count,9*rec['patch']:9*rec['patch']+9]=rec['scene_matrix']/rec['error'][:,None];start+=count
    shared=np.r_[np.arange(training.size,training.size+len(eta)),np.arange(training.sky_start,training.bg_start)]
    vshared=vcov[np.ix_(shared,shared)]
    conditional,marginal,sensitivity,info=js.propagated_cut_covariance(a,np.c_[b,scene_j],vshared,scales=testing.scales)
    take=np.array([3*i+d for i in range(testing.n) for d in [0,1]])
    moment=np.c_[np.vstack([r['response'].moment_derivative for r in held]),np.zeros((2*len(held),36))]
    transfer=sensitivity[take]+moment
    physical=conditional[np.ix_(take,take)]+transfer@vshared@transfer.T
    dump(OUT/'held_exposure_application.json',{'scope':'One held whole target exposure. The trained coherent sky and response are fixed before this fit; only independent held point positions/fluxes and one detector background plane are fitted. This remains consumed development data.',
      'held_exposure':allproblem.exposures[-1],'fit_success':bool(heldfit.success),'known_ERR_fit':js.failure_summary(rr,a),'records':rows,
      'FULL_centroid_covariance_pixel2_propagated_training_sky_and_response':physical,
      'geometry':geometry_from_centroids(testing,heldfit.x,physical,wcss,origin),
      'geometry_raw_center':geometry_from_centroids(testing,heldfit.x,marginal[np.ix_(take,take)],wcss,origin,'raw'),
      'geometry_baseline_relative_moment':geometry_from_centroids(testing,heldfit.x,physical,wcss,origin,'baseline_relative_moment'),
      'fit_optimality':float(heldfit.optimality),'active_bound_mask':heldfit.active_mask,
      'shared_training_covariance_scope':'UNCONSTRAINED local linearized training joint curvature at a provisional one-step iterate; not a converged calibrated posterior or active-bound/truncated posterior.',
      'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    summary={'status':'EXECUTED_PROVISIONAL_COHERENT_SKY_AND_SHARED_CALIBRATION_UNCERTAINTY_DEVELOPMENT',
       'scene':'Four disjoint fixed observed-sky patches, each3x3positive Gaussian coefficients shared across four exposures; width.10arcsec/spacing.18arcsec, proper coefficient priorSD10SCI. No unique source/galaxy decomposition or full host inferred.',
       'operator':'Sky basis sampled through exact archived GWCS on1/4nativepixel grid, convolved once with already pixel-integrated normalized positive field response; no extra pixel integration/PAM/ERR scaling.',
       'response':'Delivered318coefficient coarse positive field family only: spline(positive base)(u)*exp(B(u)eta), divided by finite-grid normalization (free flux reparameterization), with analytic normalization/moment derivatives and supplied first-order3gauge removal; finite25pixel support declaredzero outside, no spline extrapolation. Canonical17pixel likelihood and padded sky integration.',
       'calibration':'All4actualfieldstars constrain shared response precision+declared.2coefficientprior. Cut propagation and one re-evaluated nonlinear coupled GN application are distinct. Full covariance retains image/exposure correlations and finite-grid centroid-moment response.',
       'history':'history_rate accepts actualreadmeans andpixel-specificweights. Target histories are unavailable in current compact inputs; no historyarrays or transferredchargecoefficient fabricated. ChargeexchangeLaplacian helper conserves per-read electron charge (or common-gain units); raw DN with varying gain requires conversion first. It is not fitted or claimed identified.',
       'no_lens_prediction_loaded':True,'no_reserved_band_opened':True,'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10}
    summary['source_sha256']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['code/strong_lensing_joint_scene.py','scripts/run_he0435_joint_scene.py']}
    paths=[DELIVERY/'audit/positive_field_shape.npz',DELIVERY/'audit/positive_field_response.py',OLD/'standard_psf_convention_diagnostic.json',ROOT/'code/strong_lensing_jwst_rate_mask.py']
    native=[]
    for p in sorted((OLD/'jwst_target_diagnostic/exposures').glob('*.json')):
        item=json.loads(p.read_text());native.append({'file':item['file'],'sha256':item['sha256'],'hash_scope':'Inherited already-consumed native-file hash from this exposure receipt; not recomputed by this application.'})
        paths.extend([p,OLD/'jwst_target_diagnostic/wcs'/f'{Path(item["file"]).stem}_wcs.asdf',OLD/'standard_psf'/Path(item['file']).name.replace('_cal.fits','_FIELD_standard_epsf.fits')])
        paths.extend(ROOT/r['field_psf_path'] for r in item['points'])
    summary['input_sha256']={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    summary['native_input_receipt_hashes']=native
    summary['moment_convention']='Full finite-grid moment, raw center, and delivered baseline-relative finite-grid moment are saved separately; none establishes physical centroid accuracy.'
    summary['covariance_bound_scope']='Local UNCONSTRAINED curvature with active bound counts/indices reported. Full-column rank and observable null projections are reported; pseudoinverse zeros are not interpreted as known null directions.'
    dump(OUT/'execution_summary.json',summary)


if __name__=='__main__':main()
