"""Full nonlinear completion of the already-declared conditional optical mean.

No physical lens prediction, empirical error inflation or new PSF family enters
this optimizer. All four field-point nuisances are optimized jointly instead of
being held at a provisional profile application.
"""
from __future__ import annotations
import time
import numpy as np
from scipy.optimize import least_squares
import strong_lensing_joint_scene as js


def projected_gradient(parameters,gradient,lower,upper,scales,tolerance=1e-8):
    """Dimensionless half-objective KKT violation with correct bound signs."""
    lo=np.isfinite(lower)&np.isclose(parameters,lower,rtol=0,atol=tolerance)
    hi=np.isfinite(upper)&np.isclose(parameters,upper,rtol=0,atol=tolerance)
    answer=np.asarray(gradient).copy()
    answer[lo]=np.minimum(answer[lo],0.)
    answer[hi]=np.maximum(answer[hi],0.)
    return answer*np.asarray(scales),lo,hi


class CoupledScene:
    def __init__(self,target,fields,eta0,field_prior_sd=.2):
        self.target=target;self.fields=fields;self.d=len(eta0);self.n=target.size
        self.eta0=np.asarray(eta0);self.eta_start=self.n;self.field_start=self.n+self.d
        self.size=self.field_start+6*len(fields);self.prior_sd=field_prior_sd
        self.scales=np.r_[target.scales,np.full(self.d,.2),np.tile([1,1,1,10,10,10],len(fields))]
        self.last_x=None;self.last_r=None;self.last_j=None;self.evaluations=0;self.jacobians=0;self.history=[]

    def bounds(self):
        lo,hi=self.target.bounds();fl=[];fu=[]
        for r in self.fields:
            fl.extend([*np.asarray(r['guess'])-.7,-5.,-np.inf,-np.inf,-np.inf])
            fu.extend([*np.asarray(r['guess'])+.7,20.,np.inf,np.inf,np.inf])
        return np.r_[lo,np.full(self.d,-np.inf),fl],np.r_[hi,np.full(self.d,np.inf),fu]

    def field_application(self,parameters,record,calibration=False):
        cx,cy,lf=parameters[:3];f=np.exp(lf);x,y=record['x'],record['y'];origin=record['guess']
        p,dx,dy=record['response'].values(x-cx,y-cy,True)
        bg=np.c_[np.ones(len(x)),x-origin[0],y-origin[1]]
        mean=f*p+bg@parameters[3:];r=(mean-record['data'])/record['error']
        nuisance=np.c_[-f*dx,-f*dy,f*p,bg]/record['error'][:,None]
        cal=None if not calibration else f*record['response'].cal_values(x-cx,y-cy)/record['error'][:,None]
        return r,nuisance,cal

    def initial(self,p=None,eta=None):
        eta=self.eta0 if eta is None else np.asarray(eta)
        p=self.target.initial() if p is None else np.asarray(p)
        nuisance=[]
        for record in self.fields:
            record['response'].update(eta)
            first=np.r_[record['guess'],np.log(record['flux_guess']),np.median(record['data']),0,0]
            lo=np.r_[record['guess']-.7,-5.,[-np.inf]*3];hi=np.r_[record['guess']+.7,20.,[np.inf]*3]
            fit=least_squares(lambda q:self.field_application(q,record)[0],first,
                              jac=lambda q:self.field_application(q,record)[1],bounds=(lo,hi),
                              ftol=1e-11,xtol=1e-11,gtol=1e-8,max_nfev=100,x_scale='jac')
            nuisance.extend(fit.x)
        return np.r_[p,eta,nuisance]

    def evaluate(self,x,jacobian=False):
        if self.last_x is None or not np.array_equal(x,self.last_x):
            eta=x[self.eta_start:self.field_start]
            for rec in self.target.records:
                rec['response'].update(eta);js.refresh_scene_projection(rec)
            for rec in self.fields:rec['response'].update(eta)
            self.last_x=x.copy();self.last_j=None;self.last_r=None
        if self.last_r is not None and (not jacobian or self.last_j is not None):return self.last_r,self.last_j
        p=x[:self.n];eta=x[self.eta_start:self.field_start]
        rt,a,b,_=self.target.evaluate(p,with_cal=jacobian);rp,ap=self.target.prior(p)
        residuals=[rt];matrices=[]
        if jacobian:matrices.append(np.c_[a,b,np.zeros((len(rt),self.size-self.field_start))])
        for i,record in enumerate(self.fields):
            off=self.field_start+6*i;rr,jn,jc=self.field_application(x[off:off+6],record,jacobian)
            residuals.append(rr)
            if jacobian:
                jj=np.zeros((len(rr),self.size));jj[:,self.eta_start:self.field_start]=jc;jj[:,off:off+6]=jn;matrices.append(jj)
        residuals.extend([rp,eta/self.prior_sd])
        if jacobian:
            matrices.append(np.c_[ap,np.zeros((len(rp),self.size-self.n))])
            last=np.zeros((self.d,self.size));last[:,self.eta_start:self.field_start]=np.eye(self.d)/self.prior_sd;matrices.append(last)
        self.last_r=np.concatenate(residuals);self.evaluations+=1
        if jacobian:self.last_j=np.vstack(matrices);self.jacobians+=1
        return self.last_r,self.last_j

    def diagnostics(self,x):
        r,j=self.evaluate(x,True);lo,hi=self.bounds();g=j.T@r;pg,atlo,athi=projected_gradient(x,g,lo,hi,self.scales)
        tcount=sum(len(rec['data']) for rec in self.target.records);fcount=sum(len(rec['data']) for rec in self.fields)
        return {'objective':float(r@r),'target_raw_chi2':float(r[:tcount]@r[:tcount]),
                'field_raw_chi2':float(r[tcount:tcount+fcount]@r[tcount:tcount+fcount]),
                'scaled_projected_KKT_gradient_max':float(np.max(abs(pg))),
                'scaled_projected_KKT_target_max':float(np.max(abs(pg[:self.n]))),
                'scaled_projected_KKT_eta_max':float(np.max(abs(pg[self.eta_start:self.field_start]))),
                'scaled_projected_KKT_field_nuisance_max':float(np.max(abs(pg[self.field_start:]))),
                'active_lower_indices':np.flatnonzero(atlo),'active_upper_indices':np.flatnonzero(athi),
                'function_evaluations':self.evaluations,'jacobian_evaluations':self.jacobians}

    def solve(self,initial,callback=None,max_nfev=400):
        previous=[float(self.evaluate(initial)[0]@self.evaluate(initial)[0])];started=time.perf_counter();lower,upper=self.bounds()
        def progress(intermediate_result):
            d=self.diagnostics(intermediate_result.x);d['seconds']=time.perf_counter()-started
            d['relative_objective_change']=None if previous[0] is None else abs(d['objective']-previous[0])/max(1,previous[0]);previous[0]=d['objective']
            self.history.append(d)
            if callback is not None:callback(intermediate_result.x,d)
        result=least_squares(lambda x:self.evaluate(x)[0],initial,jac=lambda x:self.evaluate(x,True)[1],
             bounds=(lower,upper),x_scale=self.scales,ftol=1e-14,xtol=1e-14,gtol=1e-8,
             max_nfev=max_nfev,callback=progress,verbose=0)
        return result,self.diagnostics(result.x)
