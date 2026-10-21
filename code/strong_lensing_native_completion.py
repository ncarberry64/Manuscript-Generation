"""Conditional full-history/IPC optical mean with exact branch derivatives.

The continuous native OLS source equations are optimized, and their float32
implementation is replayed independently. No statistical closure, full-well
identification, error scaling, physical lens objective or new PSF family enters.
"""
from __future__ import annotations
import time
import numpy as np
from scipy import sparse
from scipy.optimize import least_squares
import strong_lensing_joint_scene as js
import strong_lensing_response_completion as rc
from strong_lensing_native_calibration import CalibrationOrderedHistory,FixedHarmonicContext
from strong_lensing_native_ols import piecewise_ols_arrays


def right_harmonic(matrix,context):
    """Adjoint of the fixed-mask linear Dirichlet extension, retaining zeros."""
    a=sparse.csr_matrix(matrix)
    if not len(context.u):return a
    out=a.toarray()
    solve=context.factor.solve(out[:,context.u].T,trans='T').T
    out[:,context.k]-=solve@context.cross
    out[:,context.u]=0.
    return sparse.csr_matrix(out)


def native_rate_adjoint(history,state,selected,weights):
    """d slope/d effective-rate and d slope/d kappa by reverse read chain.

    This is algebraically the transpose of project_directions, avoiding nineteen
    dense context-by-optical-parameter read arrays. All nonlocal IPC and implicit
    first-group-intercept terms remain in the actual derivative.
    """
    selected=np.asarray(selected,int);w=np.asarray(weights,float)
    if w.shape!=(10,len(selected)) or not state['derivatives_available']:
        raise ValueError('Current branch nine-group plus ZEROFRAME weights required')
    hs=history.h[selected];ob=state['observation'];n=history.n;m=len(selected)
    beta=sparse.csr_matrix((m,n));read=[sparse.csr_matrix((m,n)) for _ in range(19)]
    for j in range(9):
        base=sparse.diags(w[j]*ob['native_linearity_derivative'][j,selected])@hs
        inv=ob['inverse_linearity_derivative'][j]
        beta+=base@sparse.diags(inv.mean(axis=0))
        for r in range(2):read[1+2*j+r]+=base@sparse.diags(inv[r]/(2*history.gain))
    base=sparse.diags(w[9]*state['native_zero_linearity_derivative'][selected]*history.nativezero_present[selected])@hs
    iz=state['inverse_zero_linearity_derivative'];beta+=base@sparse.diags(iz)
    read[0]+=base@sparse.diags(iz/history.gain)
    for r in range(2):read[1+r]+=beta@sparse.diags(state['intercept']['beta_signal_derivative'][r]/history.gain)
    physical=state['physical'];rate=sparse.csr_matrix((m,n));caps=sparse.csr_matrix((m,n));dk=np.zeros(m)
    for i,chain in zip(history.inverse,read):
        rate+=chain@physical['observed_rate_jacobian_s'][i]
        caps+=chain@physical['observed_cap_jacobian'][i]
        dk+=chain@physical['observed_kappa_jacobian_e2'][i]
    rate+=history.cap_multiplier*history.times[0].mean()*right_harmonic(caps,history.cap_context)
    rate=right_harmonic(rate,history.rate_context)
    effective=rate@sparse.diags(history.gain)@history.g@sparse.diags(1/history.gain)
    return effective,dk


class NativeCoupledScene(rc.CoupledScene):
    def __init__(self,target,fields,eta0,histories,wcss,origin,g,h,initial):
        super().__init__(target,fields,eta0)
        self.kappa_index=self.size;self.size+=1;self.scales=np.r_[self.scales,1e-7]
        self.contexts=[];self.mean_records=[];self.last_states=[]
        self.rejected_domain_trials=[]
        self.exposure_order=sorted(wcss);first=self.exposure_order[0]
        firstrows=[r for r in target.records if r['exposure']==first]
        # The sky anchor remains the original E1 coordinate, including held fits.
        from strong_lensing_pixel_geometry import tangent_plane_arcsec
        anchors=[tangent_plane_arcsec(np.asarray(wcss[first](*r['guess'])[:2]),origin) for r in firstrows]
        for i,record in enumerate(target.records+fields):
            field=i>=target.n;patch=4 if field else record['patch'];source=histories[record['exposure']]
            xx=source['x_pixel'][patch].ravel().astype(float);yy=source['y_pixel'][patch].ravel().astype(float)
            ctx={**record,'x':xx,'y':yy,'bounds':source['bounds_xyxy'][patch].tolist()}
            if not field:
                w=wcss[record['exposure']]
                def sky(x,y,w=w):
                    ra,dec=w(x,y)[:2];return tangent_plane_arcsec(np.stack([ra,dec],axis=-1),origin)
                js.build_scene_projection(ctx,sky,anchors[patch])
            hb=source['bounds_xyxy'][patch];vx=record['x'].astype(int)-hb[0];vy=record['y'].astype(int)-hb[1]
            ctx.update({'original':record,'field':field,'target_index':i if not field else None,'field_index':i-target.n if field else None,
                'patch':patch,'history':source,'selected':vy*65+vx,
                'conversion_e_per_s_per_SCI':(source['gain_e_per_DN'][patch]*source['flat'][patch]/float(source['photmjsr'])*(float(source['group_time_header_s'])/(5*float(source['frame_time_s'])))).ravel(),
                'native_SCI_factor':float(source['photmjsr'])/source['flat'][patch].ravel()[vy*65+vx]})
            self.contexts.append(ctx)
        for ctx in self.contexts:
            rate,_,_,_=self.optical(initial,ctx)
            cal=CalibrationOrderedHistory(ctx['history'],ctx['patch'],g,h,rate)
            ctx['calibration']=cal
            before=cal.evaluate(rate,initial[-1]);old=self.native_ols(ctx,before)[0]
            # The delivered compact template's outer support edge is not a
            # measured radiation boundary. Freeze an explicit finite trusted
            # source domain once; continue the exterior harmonically, never
            # reselect it during the optimizer or clip an interior rate.
            yy,xx=np.mgrid[:65,:65];radius=8 if ctx['field'] else 16
            trusted=(abs(xx-32)<=radius)&(abs(yy-32)<=radius)
            cal.rate_context=FixedHarmonicContext((65,65),trusted)
            continued=cal.rate_context.positive(cal.gain*(cal.g@(rate/cal.gain)))
            rawcaps=cal._raw_caps(continued)
            cal.cap_context=FixedHarmonicContext((65,65),np.isfinite(rawcaps)&(rawcaps>0))
            after=cal.evaluate(rate,initial[-1]);new=self.native_ols(ctx,after)[0]
            ctx['domain_prescription']={'trusted_preIPC_square_width_pixel':2*radius+1,
                'exterior':'Fixed positive harmonic Dirichlet continuation on full65 context; conditional boundary, not measured outer light.',
                'initial_native_rate_change_max_abs_DN_s':float(np.max(abs(new-old))),
                'initial_native_mean_change_max_abs_SCI':float(np.max(abs((new-old)*ctx['native_SCI_factor'])))}

    def bounds(self):
        lo,hi=super().bounds();return np.r_[lo,0.],np.r_[hi,1e-6]

    def diagnostics(self,x):
        r,j=self.evaluate(x,True);lo,hi=self.bounds();g=j.T@r
        # Bound proximity is dimensionless. A 1e-8 absolute kappa tolerance
        # would falsely classify an appreciable positive response as zero.
        tol=1e-8*self.scales
        atlo=np.isfinite(lo)&(abs(x-lo)<=tol);athi=np.isfinite(hi)&(abs(x-hi)<=tol)
        pg=g.copy();pg[atlo]=np.minimum(pg[atlo],0.);pg[athi]=np.maximum(pg[athi],0.);pg*=self.scales
        tcount=sum(len(rec['data']) for rec in self.target.records);fcount=sum(len(rec['data']) for rec in self.fields)
        return {'objective':float(r@r),'target_raw_chi2':float(r[:tcount]@r[:tcount]),'field_raw_chi2':float(r[tcount:tcount+fcount]@r[tcount:tcount+fcount]),
            'scaled_projected_KKT_gradient_max':float(abs(pg).max()),'scaled_projected_KKT_target_max':float(abs(pg[:self.n]).max()),
            'scaled_projected_KKT_eta_max':float(abs(pg[self.eta_start:self.field_start]).max()),
            'scaled_projected_KKT_field_nuisance_max':float(abs(pg[self.field_start:self.kappa_index]).max()),
            'scaled_projected_KKT_kappa':float(pg[-1]),'active_lower_indices':np.flatnonzero(atlo),'active_upper_indices':np.flatnonzero(athi),
            'function_evaluations':self.evaluations,'jacobian_evaluations':self.jacobians}

    @staticmethod
    def native_ols(ctx,state):
        h=ctx['history'];p=ctx['patch'];sel=ctx['selected']
        return piecewise_ols_arrays(state['groups_DN'].reshape(9,-1)[:,sel],h['groupdq'][p].reshape(9,-1)[:,sel],
            gain=h['gain_e_per_DN'][p].ravel()[sel],readnoise_cds_dn=h['readnoise_CDS_DN'][p].ravel()[sel],
            group_time=float(h['group_time_header_s']),frame_time=float(h['frame_time_s']),nframes=2,
            dark_dn_s=h['average_dark_current_DN_s'][p].ravel()[sel],zeroframe=state['zeroframe_DN'].ravel()[sel],
            pixeldq=h['pixeldq'][p].ravel()[sel],cast_input=False)

    def solve(self,initial,callback=None,max_nfev=400):
        base=self.evaluate(initial)[0];count=len(base);previous=[float(base@base)]
        started=time.perf_counter();lower,upper=self.bounds()
        def feasible_residual(x):
            try:return self.evaluate(x)[0]
            except ValueError as exc:
                # A finite barrier makes TRF reject and shorten an infeasible
                # trial; it does not clip photons or change the fixed mask.
                # Its Jacobian is never consumed at a rejected trial.
                message=str(exc)
                domain=('positive domain','negative','Excessive log-response','monotone','linearity','finite')
                if not any(word in message for word in domain):raise
                self.rejected_domain_trials.append({'error':message,'kappa':float(x[-1])})
                self.last_x=None;self.last_r=None;self.last_j=None
                return np.full(count,1e12)
        def progress(intermediate_result):
            d=self.diagnostics(intermediate_result.x);d['seconds']=time.perf_counter()-started
            d['relative_objective_change']=abs(d['objective']-previous[0])/max(1,previous[0]);previous[0]=d['objective']
            d['domain_rejected_trials']=len(self.rejected_domain_trials);self.history.append(d)
            if callback is not None:callback(intermediate_result.x,d)
        result=least_squares(feasible_residual,initial,jac=lambda x:self.evaluate(x,True)[1],
            bounds=(lower,upper),x_scale=self.scales,ftol=1e-14,xtol=1e-14,gtol=1e-8,
            max_nfev=max_nfev,callback=progress)
        return result,self.diagnostics(result.x)

    def optical(self,x,ctx,jacobian=False):
        rec=ctx['original'];resp=ctx['response']
        if ctx['field']:
            off=self.field_start+6*ctx['field_index'];q=x[off:off+6];sky=None
            plane_design=np.c_[np.ones(len(rec['x'])),rec['x']-rec['guess'][0],rec['y']-rec['guess'][1]];bgcols=np.arange(off+3,off+6)
        else:
            off=3*ctx['target_index'];q=x[off:off+3];sky=self.target.sky(x[:self.n])[9*ctx['patch']:9*ctx['patch']+9]
            bg=self.target.bg_start+3*self.target.eindex[rec['exposure']];bgcols=np.arange(bg,bg+3);plane_design=rec['detector_background']
        f=np.exp(q[2]);profile,dx,dy=resp.values(ctx['x']-q[0],ctx['y']-q[1],True);source=f*profile
        if sky is not None:source+=ctx['scene_matrix']@sky
        rate=source*ctx['conversion_e_per_s_per_SCI'];plane=plane_design@x[bgcols]
        if not jacobian:return rate,plane,None,None
        jl=np.zeros((len(rate),self.size));jl[:,off:off+3]=np.c_[-f*dx,-f*dy,f*profile]
        cal=f*resp.cal_values(ctx['x']-q[0],ctx['y']-q[1])
        if sky is not None:
            cal+=js.scene_calibration_columns(ctx,sky)
            if self.target.nsky:
                start=self.target.sky_start+9*ctx['patch'];jl[:,start:start+9]=ctx['scene_matrix']
        jl[:,self.eta_start:self.field_start]=cal;jl*=ctx['conversion_e_per_s_per_SCI'][:,None]
        jp=np.zeros((len(plane),self.size));jp[:,bgcols]=plane_design
        return rate,plane,jl,jp

    def evaluate(self,x,jacobian=False):
        if self.last_x is None or not np.array_equal(x,self.last_x):
            eta=x[self.eta_start:self.field_start]
            for ctx in self.contexts:
                ctx['response'].update(eta)
                if not ctx['field']:js.refresh_scene_projection(ctx)
            self.last_x=x.copy();self.last_r=None;self.last_j=None
        if self.last_r is not None and (not jacobian or self.last_j is not None):return self.last_r,self.last_j
        residual=[];matrices=[];records=[];states=[]
        for ctx in self.contexts:
            rate,plane,jl,jp=self.optical(x,ctx,jacobian);cal=ctx['calibration'];h=ctx['history'];p=ctx['patch'];sel=ctx['selected']
            state=cal.evaluate(rate,x[-1],derivatives=jacobian)
            val,weights,branches=self.native_ols(ctx,state)
            rec=ctx['original'];value=val*ctx['native_SCI_factor']+plane
            rr=(value-rec['data'])/rec['error'];residual.append(rr)
            if jacobian:
                chain,jk=native_rate_adjoint(cal,state,sel,weights)
                jj=(chain@jl)*ctx['native_SCI_factor'][:,None]+jp
                jj[:,-1]+=jk*ctx['native_SCI_factor'];matrices.append(jj/rec['error'][:,None])
            records.append({'name':rec['name'],'field':ctx['field'],'raw_chi2':float(rr@rr),'pixels':len(rr),
                'mean_SCI':value,'residual_SCI':value-rec['data'],'native_branch_nonunique_count':sum(not b.get('locally_unique_derivative',False) for b in branches),
                'native_SNR_boundary_count':sum(b.get('snr_boundary',False) for b in branches),'receipt':state['diagnostics'],'domain_prescription':ctx['domain_prescription']})
            states.append((state,plane,weights))
        rp,ap=self.target.prior(x[:self.n]);residual.extend([rp,x[self.eta_start:self.field_start]/self.prior_sd])
        if jacobian:
            matrices.append(np.c_[ap,np.zeros((len(rp),self.size-self.n))]);last=np.zeros((self.d,self.size));last[:,self.eta_start:self.field_start]=np.eye(self.d)/self.prior_sd;matrices.append(last)
        self.mean_records=records;self.last_states=states;self.last_r=np.concatenate(residual);self.evaluations+=1
        if jacobian:self.last_j=np.vstack(matrices);self.jacobians+=1
        return self.last_r,self.last_j
