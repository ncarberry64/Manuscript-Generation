"""Conditional direct-group charge redistribution diagnostic; never physical admission.

Rates are DN/s, group values DN, kappa electron^-1. The graph Laplacian has
off-diagonal +1 and diagonal minus degree, with no exterior edges.
"""
from __future__ import annotations

import numpy as np
from scipy import linalg, sparse
from scipy.optimize import least_squares
from scipy.interpolate import RectBivariateSpline
from scipy.stats import chi2


def group_time_moments(read_times):
    t = np.asarray(read_times, float)
    if t.ndim != 2 or np.any(t <= 0):
        raise ValueError("positive group-by-read times required")
    return t.mean(axis=1), np.square(t).mean(axis=1)


def temporal_covariance(read_times, rate, gain, readnoise_cds):
    """Ideal independent-pixel Poisson+CDS covariance, after unit linearity gain.

    Shared SIRS, dark-reference, and native linearity Jacobians are NOT included.
    The caller must freeze rate from training data, not from the held group.
    """
    t = np.asarray(read_times, float)
    n = t.shape[1]
    kernel = np.minimum(t[:, None, :, None], t[None, :, None, :]).mean(axis=(2, 3))
    return max(float(rate), 0.) / float(gain) * kernel + np.eye(len(t)) * float(readnoise_cds)**2 / (2*n)


def zero_flux_laplacian(shape):
    ny, nx = shape
    rows, cols, vals = [], [], []
    for y in range(ny):
        for x in range(nx):
            p = y*nx+x
            degree = 0
            for yy, xx in ((y-1,x),(y+1,x),(y,x-1),(y,x+1)):
                if 0 <= yy < ny and 0 <= xx < nx:
                    rows.append(p); cols.append(yy*nx+xx); vals.append(1.)
                    degree += 1
            rows.append(p); cols.append(p); vals.append(-float(degree))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(ny*nx, ny*nx))


def electron_conserving_operator(gain):
    """Convert DN rates to electrons, exchange charge, then return DN units.

    B(f_DN^2)=diag(1/g) L diag(g^2) f_DN^2. Thus sum(g*B)=0;
    sum(B) generally does not vanish when adjacent electronic gains differ.
    """
    g=np.asarray(gain,float)
    if not np.all(np.isfinite(g)) or np.any(g<=0):
        raise ValueError("finite positive gains required on the full charge graph")
    return sparse.diags(1/g.ravel())@zero_flux_laplacian(g.shape)@sparse.diags(g.ravel()**2)


def selected_history(groupdq, pixeldq=0, stop=None):
    """Observed-mask conditional history: DNU/SAT excluded; JUMP starts intercept.

    A gap also starts a fresh intercept. A JUMP group itself is retained.
    """
    dq = np.asarray(groupdq)
    if int(pixeldq) & 1:
        return np.array([], int), np.array([], int)
    valid = (dq & 3) == 0
    if stop is not None:
        valid[np.arange(len(dq)) >= stop] = False
    idx = np.flatnonzero(valid)
    seg = []
    k = -1
    last = -2
    for j in idx:
        if j != last+1 or (int(dq[j]) & 4):
            k += 1
        seg.append(k); last = j
    return idx, np.asarray(seg, int)


def _design(segments):
    return np.eye(int(np.max(segments))+1)[segments]


def pilot_rates(groups, groupdq, pixeldq, mean_times, train_stop=8):
    """Training-only unweighted segmented-line pilot; no held data consumed."""
    ng, ny, nx = groups.shape
    out = np.zeros((ny, nx), float)
    rank_fail = 0
    for p in range(ny*nx):
        y, x = divmod(p, nx)
        idx, seg = selected_history(groupdq[:,y,x], pixeldq[y,x], train_stop)
        idx_finite = np.isfinite(groups[idx,y,x])
        if not np.all(idx_finite):
            rank_fail += 1; continue
        if len(idx) <= 1 or len(idx)-len(np.unique(seg)) < 1:
            rank_fail += 1; continue
        h = _design(seg)
        design = np.column_stack((h, mean_times[idx]))
        solution, _, rank, _ = np.linalg.lstsq(design, groups[idx,y,x], rcond=1e-12)
        if rank == design.shape[1]:
            out[y,x] = solution[-1]
        else:
            rank_fail += 1
    return out, {"training_only": True, "rank_deficient_or_unusable_pixels": rank_fail,
                 "negative_pilot_rates": int(np.sum(out < 0)), "covariance_rate_clipped_at_zero": True}


def conditional_innovation(ytrain, yhold, mutrain, muhold, ctt, cht, chh):
    r = linalg.solve(ctt, cht.T, assume_a="pos").T
    return yhold-muhold-r@(ytrain-mutrain), chh-r@cht.T


def fit_history(groups, groupdq, pixeldq, gain, readnoise_cds, read_times,
                roi=(25,42,25,42), train_stop=8, redistribution=True, max_nfev=80,
                optical=None):
    """Joint positive-rate fit with profiled segment intercepts and a full graph.

    Covariance is frozen once from a segmented TRAIN-only pilot. Padded rates
    are fixed at that nonnegative pilot. All ROI rates and kappa are coupled.
    The kappa upper bound guarantees nonnegative accumulated source mean for
    EVERY rate in the predeclared rate box, including zero-flux graph edges.
    """
    g = np.asarray(groups, float)
    ng, ny, nx = g.shape
    m, v = group_time_moments(read_times)
    pilot, pilot_note = pilot_rates(g, groupdq, pixeldq, m, train_stop)
    physical_pilot = np.maximum(pilot, 0.)
    y0,y1,x0,x1 = roi
    roi_flat = np.arange(ny*nx).reshape(ny,nx)[y0:y1,x0:x1].ravel()
    lap = electron_conserving_operator(gain).tocsr()
    full_gain=np.asarray(gain,float).ravel()
    # All selected ROI pixels must carry an estimable slope after intercepts.
    infos = []
    active = []
    for flat in roi_flat:
        yy,xx = divmod(int(flat),nx)
        idx, seg = selected_history(groupdq[:,yy,xx], pixeldq[yy,xx], train_stop)
        if len(idx)-len(np.unique(seg)) < 1 or not np.all(np.isfinite(g[idx,yy,xx])):
            continue
        if not np.isfinite(gain[yy,xx]) or gain[yy,xx] <= 0 or not np.isfinite(readnoise_cds[yy,xx]) or readnoise_cds[yy,xx] <= 0:
            continue
        h = _design(seg)
        c = temporal_covariance(read_times, physical_pilot[yy,xx],gain[yy,xx],readnoise_cds[yy,xx])
        ct = c[np.ix_(idx,idx)]
        iw = linalg.solve_triangular(linalg.cholesky(ct,lower=True), np.eye(len(idx)), lower=True)
        hw = iw@h
        q, _ = linalg.qr(hw,mode="full")
        transform = q[:,h.shape[1]:].T@iw
        ainv = linalg.inv(hw.T@hw)
        projector = ainv@h.T@linalg.cho_solve(linalg.cho_factor(ct,lower=True), np.eye(len(idx)))
        info = dict(flat=int(flat), y=yy,x=xx,idx=idx,seg=seg,h=h,c=c,ct=ct,iw=iw,
                    transform=transform,ainv=ainv,projector=projector,
                    yt=g[idx,yy,xx], tm=transform@m[idx],tv=transform@v[idx])
        active.append(int(flat)); infos.append(info)
    active = np.asarray(active,int)
    if not len(active):
        raise ValueError("no estimable ROI histories")
    lrr = lap[active,:][:,active].tocsr()
    padding = physical_pilot.ravel().copy()
    padding[active] = 0
    bpad = np.asarray(lap[active,:]@np.square(padding)).ravel()
    fmax = max(1., 2*float(np.max(physical_pilot)))
    nrate = len(active)
    nbase = nrate if optical is None else 6
    if optical is not None:
        grid=np.asarray(optical["psf_grid"],float)
        spline=RectBivariateSpline(np.arange(grid.shape[1]),np.arange(grid.shape[0]),grid.T,kx=3,ky=3,s=0)
        os=float(optical["oversampling"])
        origin=np.asarray(optical.get("origin",[(grid.shape[1]-1)/2,(grid.shape[0]-1)/2]),float)
        flat_factor=np.asarray(optical["flat"],float).ravel()[active]/float(optical["photmjsr"])
        px=active%nx;py=active//nx
        ref=np.asarray(optical["initial_center_xy"],float)
        dx=px-ref[0];dy=py-ref[1]
        rx=float(np.max(np.abs(dx)));ry=float(np.max(np.abs(dy)))
        flux_scale=float(optical["initial_flux_SCI"])
        bg=np.asarray(optical["initial_background_SCI"],float)
        bgmin=max(1e-5,float(bg[0]-rx*abs(bg[1])-ry*abs(bg[2])))
        # Positive affine background parameterized by its minimum over a box;
        # no extra curvature or PSF freedom is introduced.
        gradmax=max(.03,5*float(np.max(np.abs(bg[1:]))))
        bgminmax=max(3.,5*abs(float(bg[0])))
        optical_lo=np.array([ref[0]-2,ref[1]-2,0.,0.,-gradmax,-gradmax])
        optical_hi=np.array([ref[0]+2,ref[1]+2,2.,bgminmax,gradmax,gradmax])
        optical_initial=np.array([*ref,1.,bgmin,bg[1],bg[2]])
        # B-spline coefficients bound its values inside support (positive basis
        # with partition of unity); this upper bound is fixed before fitting.
        pmax=max(0.,float(np.max(spline.get_coeffs())))
        fmax=max(fmax,float(np.max(flat_factor))*(2*flux_scale*pmax+bgminmax+2*rx*gradmax+2*ry*gradmax))
        def rate_map(z):
            cx,cy,fluxratio,bmin,bx,by=z[:6]
            gx=origin[0]+os*(px-cx);gy=origin[1]+os*(py-cy)
            inside=(gx>=0)&(gy>=0)&(gx<=grid.shape[1]-1)&(gy<=grid.shape[0]-1)
            p=np.zeros(nrate);pcx=np.zeros(nrate);pcy=np.zeros(nrate)
            p[inside]=spline.ev(gx[inside],gy[inside])
            pcx[inside]=-os*spline.ev(gx[inside],gy[inside],dx=1)
            pcy[inside]=-os*spline.ev(gx[inside],gy[inside],dy=1)
            ax=np.sqrt(bx*bx+1e-20);ay=np.sqrt(by*by+1e-20)
            b0=bmin+rx*ax+ry*ay
            flux=fluxratio*flux_scale
            rates=flat_factor*(flux*p+b0+bx*dx+by*dy)
            derivatives=flat_factor[:,None]*np.column_stack((flux*pcx,flux*pcy,flux_scale*p,
                                        np.ones(nrate),rx*bx/ax+dx,ry*by/ay+dy))
            return rates,derivatives
    else:
        def rate_map(z):
            return z[:nrate],None
    # Q_DN,p(t)>=t*f_p*(1-kappa_e*t*degree_p*g_p*f_p/2): this bound
    # protects every individual read, not merely the group average.
    kmax = 2./(4*float(np.max(full_gain))*fmax*float(np.max(read_times)))
    npar = nbase+int(redistribution)
    blocks = [len(i["tm"]) for i in infos]
    offsets = np.r_[0,np.cumsum(blocks)]
    yproj = np.concatenate([i["transform"]@i["yt"] for i in infos])
    def unpack(z):
        return rate_map(z)[0], (z[-1]*kmax if redistribution else 0.)
    def residual(z):
        f,k = unpack(z)
        b = lrr@np.square(f)+bpad
        return np.concatenate([i["tm"]*f[p]+.5*k*i["tv"]*b[p] for p,i in enumerate(infos)])-yproj
    def jacobian(z):
        f,k = unpack(z)
        b = lrr@np.square(f)+bpad
        rr,cc,dd=[],[],[]
        for p,i in enumerate(infos):
            lo,hi=offsets[p:p+2]
            for q,lval in zip(lrr.indices[lrr.indptr[p]:lrr.indptr[p+1]],lrr.data[lrr.indptr[p]:lrr.indptr[p+1]]):
                val=k*i["tv"]*lval*f[q]
                if q==p: val=val+i["tm"]
                rr.extend(range(lo,hi)); cc.extend([q]*(hi-lo)); dd.extend(val)
            if redistribution:
                rr.extend(range(lo,hi));cc.extend([nrate]*(hi-lo));dd.extend(.5*kmax*i["tv"]*b[p])
        rawjac=sparse.csr_matrix((dd,(rr,cc)),shape=(len(yproj),nrate+int(redistribution)))
        if optical is None: return rawjac
        mapped=rawjac[:,:nrate]@rate_map(z)[1]
        if redistribution: mapped=np.column_stack((mapped,rawjac[:,-1].toarray().ravel()))
        return mapped
    z0 = physical_pilot.ravel()[active].copy() if optical is None else optical_initial
    if redistribution: z0=np.r_[z0,.01]
    if optical is None:
        z0=np.maximum(z0,1e-10)
        lower=np.zeros(npar);upper=np.full(npar,fmax)
    else:
        lower=np.r_[optical_lo,0.] if redistribution else optical_lo
        upper=np.r_[optical_hi,1.] if redistribution else optical_hi
    if redistribution: upper[-1]=1
    fit=least_squares(residual,z0,jac=jacobian,bounds=(lower,upper),
                      x_scale="jac",ftol=1e-9,xtol=1e-9,gtol=1e-7,max_nfev=max_nfev,
                      tr_options={"atol":1e-12,"btol":1e-12,"maxiter":1500})
    f,k=unpack(fit.x)
    b=np.asarray(lrr@np.square(f)+bpad).ravel()
    jf=jacobian(fit.x)
    if sparse.issparse(jf): jf=jf.toarray()
    _,s,vt=linalg.svd(jf,full_matrices=False)
    tol=np.finfo(float).eps*max(jf.shape)*s[0]
    keep=s>tol
    cov=(vt[keep].T/np.square(s[keep]))@vt[keep]
    rank=int(np.sum(keep))
    held=[]; raw=[]; projected=[]; mean=[]; train_records=[]; drows=[]; vdiag=[]; independent_diag=[]; independent_drows=[]
    for p,i in enumerate(infos):
        idx=i["idx"]; yy=i["y"];xx=i["x"]
        source=m[idx]*f[p]+.5*k*v[idx]*b[p]
        intercept=i["projector"]@(i["yt"]-source)
        mu=i["h"]@intercept+source
        r=i["yt"]-mu
        raw.extend(r);mean.extend(mu);projected.extend(i["transform"]@r)
        train_records.extend([[i["flat"],int(j),int(sg),float(a),float(t),float(u),float(e)] for j,sg,a,t,u,e in zip(idx,i["seg"],i["yt"],mu,r,np.diag(i["ct"])**.5)])
        # Ninth-group holdout is prescribed for unchanged full nine-group wings.
        eligible=(ng==9 and train_stop==8 and not(int(pixeldq[yy,xx])&1)
                  and np.all((groupdq[:,yy,xx]&7)==0) and len(idx)==8)
        if not eligible: continue
        jh=8
        xt=np.zeros((len(idx),nrate));xh=np.zeros(nrate)
        for q,lval in zip(lrr.indices[lrr.indptr[p]:lrr.indptr[p+1]],lrr.data[lrr.indptr[p]:lrr.indptr[p+1]]):
            xt[:,q]=k*v[idx]*lval*f[q];xh[q]=k*v[jh]*lval*f[q]
        xt[:,p]+=m[idx];xh[p]+=m[jh]
        if optical is not None:
            rate_jac=rate_map(fit.x)[1]
            xt=xt@rate_jac;xh=xh@rate_jac
        if redistribution:
            xt=np.column_stack((xt,.5*kmax*v[idx]*b[p]));xh=np.r_[xh,.5*kmax*v[jh]*b[p]]
        cht=i["c"][jh,idx][None,:]
        rr=linalg.cho_solve(linalg.cho_factor(i["ct"],lower=True),cht.T).T
        mu_h=intercept[0]+m[jh]*f[p]+.5*k*v[jh]*b[p]
        prediction=float(mu_h+(rr@r)[0])
        dh=1-float((rr@np.ones(len(idx)))[0])
        slope_d=xh-(rr@xt)[0]-dh*(i["projector"]@xt)[0]
        noise=float(i["c"][jh,jh]-(rr@cht.T)[0,0])
        vdiag.append(noise+dh**2*i["ainv"][0,0]);drows.append(slope_d)
        independent_diag.append(float(i["c"][jh,jh]+i["ainv"][0,0]))
        independent_drows.append(xh-(i["projector"]@xt)[0])
        held.append([i["flat"],g[jh,yy,xx],prediction,g[jh,yy,xx]-prediction,mu_h,float((rr@r)[0])])
    held=np.asarray(held,float);drows=np.asarray(drows,float)
    predcov=np.diag(vdiag)+drows@cov@drows.T
    error=held[:,3]
    holdq=float(error@linalg.solve(predcov,error,assume_a="pos"))
    # Deliberately incorrect independence statistic retained to quantify its effect.
    independent_error=held[:,1]-held[:,4]
    independent_drows=np.asarray(independent_drows)
    independent_cov=np.diag(independent_diag)+independent_drows@cov@independent_drows.T
    independent_q=float(independent_error@linalg.solve(independent_cov,independent_error,assume_a="pos"))
    full_rates=physical_pilot.ravel().copy();full_rates[active]=f
    lapfull=np.asarray(lap@np.square(full_rates)).ravel()
    full_source=m[:,None]*full_rates[None,:]+.5*k*v[:,None]*lapfull[None,:]
    read_t=np.asarray(read_times).ravel()
    read_source=read_t[:,None]*full_rates[None,:]+.5*k*np.square(read_t[:,None])*lapfull[None,:]
    qtrain=float(np.dot(fit.fun,fit.fun))
    fullrank=rank+sum(i["h"].shape[1] for i in infos)
    dof=len(yproj)-rank
    ksigma=float(kmax*np.sqrt(cov[-1,-1])) if redistribution else 0.
    summary={
        "status":"CONDITIONAL_DIRECT_GROUP_DIAGNOSTIC_ONLY", "physical_calibration_pass":False,
        "mean_model":"a_p,segment + mean(t)_j*f_DN,p + kappa_e*mean(t^2)_j/2/g_p*L((g*f_DN)^2)_p",
        "rate_model":"free pixel rates" if optical is None else "flat/PHOTMJSR * (F*dated_STPSF(cx,cy)+b0+bx*dx+by*dy)",
        "rate_units":"DN/s", "kappa_units":"electron^-1", "model_redistribution":bool(redistribution),
        "optimizer_success":bool(fit.success),"optimizer_message":fit.message,"nfev":int(fit.nfev),"optimality":float(fit.optimality),
        "roi_bounds_local_y0_y1_x0_x1":list(roi),"roi_active_pixels":nrate,"training_observations":int(sum(len(i["idx"]) for i in infos)),
        "segment_intercepts":int(fullrank-rank),"profile_parameter_rank":rank,"profile_parameter_count":npar,
        "joint_parameter_rank":fullrank,"joint_parameter_count":npar+sum(i["h"].shape[1] for i in infos),
        "svd_tolerance":float(tol),"profile_condition_number":float(s[0]/s[-1]),"covariance_scale_from_residuals":False,
        "training_Q":qtrain,"training_dof":int(dof),"training_Q_per_dof":qtrain/dof,"training_conditional_working_p":float(chi2.sf(qtrain,dof)),
        "training_raw_residual_RMS_DN":float(np.sqrt(np.mean(np.square(raw)))),"training_max_abs_raw_residual_DN":float(np.max(np.abs(raw))),
        "kappa":float(k),"kappa_local_curvature_standard_error":ksigma,"kappa_upper_bound":float(kmax),
        "kappa_bound_active":bool(redistribution and (fit.x[-1]<1e-5 or fit.x[-1]>1-1e-5)),
        "parameter_coordinates":(["f_p DN/s in active_flat order"] if optical is None else ["cx local native pixel","cy local native pixel","F/initial_F","bmin SCI","bx SCI/pixel","by SCI/pixel"])+(["beta=kappa/kappa_upper_bound"] if redistribution else []),
        "held_group_independence":"internal temporal holdout only; already present in previously consumed native products; no new independent observation",
        "rate_upper_bound_DN_s":float(fmax),"rate_bound_active_count":int(np.sum((f<1e-7)|(f>fmax-1e-7))),
        "held_group":8,"held_pixels":len(held),"held_Q":holdq,"held_dof":len(held),"held_Q_per_dof":holdq/len(held),
        "held_conditional_working_p":float(chi2.sf(holdq,len(held))),"held_raw_RMS_DN":float(np.sqrt(np.mean(error**2))),
        "held_max_abs_raw_residual_DN":float(np.max(np.abs(error))),"held_covariance_includes_shared_reads":True,
        "deliberately_incorrect_independent_Q":independent_q,"deliberately_incorrect_independent_RMS_DN":float(np.sqrt(np.mean(independent_error**2))),
        "full_graph_sum_electron_charge_operator":float(full_gain@lapfull),
        "full_graph_electron_charge_operator_l1":float(np.sum(np.abs(full_gain*lapfull))),
        "relative_electron_charge_conservation_residual":float(abs(full_gain@lapfull)/max(np.sum(np.abs(full_gain*lapfull)),1e-300)),
        "full_graph_final_group_redistribution_electron":float(.5*k*v[-1]*(full_gain@lapfull)),
        "full_graph_final_group_redistribution_DN":float(.5*k*v[-1]*np.sum(lapfull)),
        "roi_final_group_boundary_exchange_DN":float(.5*k*v[-1]*np.sum(lapfull[active])),
        "roi_final_group_boundary_exchange_electron":float(.5*k*v[-1]*np.sum(full_gain[active]*lapfull[active])),
        "minimum_accumulated_source_mean_DN":float(np.min(full_source)),"positive_source_mean_verified":bool(np.min(full_source)>=-1e-9),
        "minimum_individual_read_source_DN":float(np.min(read_source)),"positive_individual_read_source_verified":bool(np.min(read_source)>=-1e-9),
        "padding_rates_fixed_training_pilot":True,"pilot":pilot_note,
        "gain_map_e_per_DN_min_max":[float(np.min(full_gain)),float(np.max(full_gain))],
        "conditional_limits":["frozen TRAIN-only Poisson plug-in covariance and fixed observed GROUPDQ masks",
            "CDS reference plus ideal independent-read noise; native SIRS spatial covariance omitted",
            "unit native post-linearity Jacobian approximation; exact nine-group Jacobians absent",
            "native correction of an averaged ADC group differs from averaging corrected individual reads; unresolved Jensen/history contribution",
            "late dark and calibration-reference covariance omitted; padded rates treated fixed",
            "spatial gain map fixed; gain uncertainty omitted; electron conservation is imposed before mixed-history slope selection",
            "stationary t^2 source surrogate does not identify charge-jump memory or a unique detector mechanism",
            "positive kappa boundary invalidates an unqualified symmetric confidence interval",
            "all working p values are conditional diagnostics; no physical response or lens validation"]}
    arrays=dict(active_flat=active,roi_rate_DN_s=f,pilot_rate_DN_s=pilot,full_rate_DN_s=full_rates.reshape(ny,nx),
                fitted_parameter_coordinates=fit.x,parameter_lower_bounds=lower,parameter_upper_bounds=upper,
                electron_L_f2_over_gain=lapfull.reshape(ny,nx),gain_e_per_DN=np.asarray(gain),parameter_covariance=cov,parameter_singular_values=s,
                projected_training_residual=np.asarray(projected),train_records=np.asarray(train_records),
                held_records=held,held_predictive_covariance=predcov,held_parameter_derivative=drows,
                held_conditional_base_variance=np.asarray(vdiag),group_mean_time_s=m,group_mean_squared_time_s2=v)
    if optical is not None:
        cx,cy,fr,bmin,bx,by=fit.x[:6]
        bg0=float(bmin+rx*np.sqrt(bx*bx+1e-20)+ry*np.sqrt(by*by+1e-20))
        summary["optical_fit"]={"centroid_local_xy_pixel":[float(cx),float(cy)],
            "flux_SCI":float(fr*flux_scale),"background_reference_local_xy_pixel":ref.tolist(),
            "background_plane_SCI":[bg0,float(bx),float(by)],"photmjsr":float(optical["photmjsr"]),
            "flux_parameter_scale_SCI":flux_scale,"minimum_background_SCI":float(bmin),
            "parameter_bound_active_indices":np.flatnonzero((fit.x-lower<1e-7)|(upper-fit.x<1e-7)).tolist(),
            "minimum_optical_rate_DN_s":float(np.min(f)),"maximum_optical_rate_DN_s":float(np.max(f)),
            "fixed_PSF_not_reestimated":True,"PSF_derivatives_propagated_through_electron_charge_operator":True,
            "mean_adequacy_at_working_alpha_0p001":bool(chi2.sf(qtrain,dof)>=.001),
            "physical_optical_response_validation_pass":False}
    return summary,arrays
