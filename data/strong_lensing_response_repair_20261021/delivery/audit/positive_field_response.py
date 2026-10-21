"""Positive shared field ePSF calibration on valid early-group CAL rates.

Exploratory development, not independent validation. No physical lens score or
target geometry enters training. SAT alone is retained; DNU/nonfinite/ERR<=0 is
excluded. A common positive log-response is learned from 3/4 field exposures,
then checked on the fourth with only its point nuisance parameters fitted.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.linalg import null_space
from scipy.optimize import least_squares
from astropy.io import fits
from scipy.stats import chi2

HERE=Path(__file__).resolve().parent;REPO=HERE.parent/'repo';ART=REPO/'artifacts/strong_lensing_geometry_20261021'
sys.path.insert(0,str(HERE));from field_response_transfer import read_primary_fits

def get_fields():
    diag=json.loads((ART/'standard_psf_convention_diagnostic.json').read_text());out=[]
    for row in diag['records']:
        center=np.asarray(row['reproduced_global_centroid']);ix,iy=np.rint(center).astype(int)
        with fits.open(REPO/'data/strong_lensing_geometry_20261021/native'/row['exposure'],memmap=False) as h:
            d,e,q=[np.asarray(h[k].data[iy-8:iy+9,ix-8:ix+9]) for k in ['SCI','ERR','DQ']]
            yy,xx=np.mgrid[iy-8:iy+9,ix-8:ix+9];mask=np.isfinite(d)&np.isfinite(e)&(e>0)&((q.astype(np.uint64)&1)==0)
            g,_=read_primary_fits(ART/'standard_psf'/row['exposure'].replace('_cal.fits','_FIELD_standard_epsf.fits'))
            pedestal=max(0.,-float(g.min()))*1.1+1e-8
            out.append({'name':row['exposure'],'x':xx[mask]-center[0],'y':yy[mask]-center[1],
                'data':d[mask].astype(float),'error':e[mask].astype(float),'grid':g+pedestal,
                'pixel_dq':q[mask].astype(np.uint64),'pixel_x':xx[mask],'pixel_y':yy[mask],
                'signed_grid':g,'pedestal':pedestal,'origin':center,'sat_restored':int(np.sum(mask&((q&2)!=0))),
                'old_flux':row['flux_SCI_pixel_units'],'old_background':row['background_plane_SCI'],
                'interpolator':RectBivariateSpline(np.arange(101)/4-12.5,np.arange(101)/4-12.5,(g+pedestal).T)})
    return out

class Basis:
    def __init__(self,reference_grid,spacing=1.,width=.65):
        self.width=width;z=np.arange(-9.,9.0001,spacing);xx,yy=np.meshgrid(z,z);self.nodes=np.c_[xx.ravel(),yy.ravel()]
        self.nodes=self.nodes[np.linalg.norm(self.nodes,axis=1)<=10.5]
        g=np.arange(101)/4-12.5;gx,gy=np.meshgrid(g,g);self.gx,self.gy=gx,gy
        b,_,_=self.raw(gx.ravel(),gy.ravel());p=reference_grid.ravel()
        c=np.vstack([p@b,(p*gx.ravel())@b,(p*gy.ravel())@b])
        self.null=null_space(c);self.grid_basis=b@self.null
    def raw(self,x,y):
        dx=np.asarray(x)[:,None]-self.nodes[None,:,0];dy=np.asarray(y)[:,None]-self.nodes[None,:,1]
        b=np.exp(-(dx*dx+dy*dy)/(2*self.width**2))
        return b,-dx/self.width**2*b,-dy/self.width**2*b
    def evaluate(self,x,y,rec=None):return tuple(v@self.null for v in self.raw(x,y))
    def response_grid(self,rec,a):return rec['grid']*np.exp((self.grid_basis@a).reshape(101,101))

def fit(records,basis,prior_sd=.2,a_fixed=None,initial=None,max_nfev=70):
    free=a_fixed is None;nshape=basis.null.shape[1] if free else 0
    p=[];sizes=[9 if r.get('background_degree',1)==2 else 6 for r in records]
    offsets=nshape+np.r_[0,np.cumsum(sizes[:-1])]
    for r in records:
        bg=np.median(r['data']);flux=r['old_flux']
        p.extend([0.,0.,np.log(flux),bg-flux*r['pedestal'],0.,0.])
        if r.get('background_degree',1)==2:p.extend([0.,0.,0.])
    p=np.r_[np.zeros(nshape),p] if initial is None else initial.copy()
    cache={}
    def model_and_jac(params):
        if 'p' in cache and np.array_equal(params,cache['p']):return cache['v']
        a=params[:nshape] if free else a_fixed;rs=[];js=[]
        for i,r in enumerate(records):
            off=int(offsets[i]);cx,cy,lf,bg,bx,by=params[off:off+6]
            x,y=r['x']-cx,r['y']-cy
            b,dbx,dby=basis.evaluate(x,y,r);corr=np.exp(np.clip(b@a,-15,15))
            base=r['interpolator'].ev(x,y)
            px=r['interpolator'].ev(x,y,dx=1);py=r['interpolator'].ev(x,y,dy=1)
            f=np.exp(lf);point=f*base*corr;pred=point+bg+bx*r['x']+by*r['y']
            if sizes[i]==9:
                bxx,bxy,byy=params[off+6:off+9]
                pred=pred+bxx*r['x']**2+bxy*r['x']*r['y']+byy*r['y']**2
            resid=(pred-r['data'])/r['error'];jac=np.zeros((len(resid),len(params)))
            if free:jac[:,:nshape]=point[:,None]*b/r['error'][:,None]
            jac[:,off]=(-f*corr*(px+base*(dbx@a)))/r['error']
            jac[:,off+1]=(-f*corr*(py+base*(dby@a)))/r['error']
            jac[:,off+2]=point/r['error'];jac[:,off+3]=1/r['error']
            jac[:,off+4]=r['x']/r['error'];jac[:,off+5]=r['y']/r['error']
            if sizes[i]==9:
                jac[:,off+6]=r['x']**2/r['error'];jac[:,off+7]=r['x']*r['y']/r['error'];jac[:,off+8]=r['y']**2/r['error']
            rs.append(resid);js.append(jac)
        if free:
            rs.append(a/prior_sd);q=np.zeros((nshape,len(params)));q[:,:nshape]=np.eye(nshape)/prior_sd;js.append(q)
        result=(np.concatenate(rs),np.concatenate(js));cache['p']=params.copy();cache['v']=result;return result
    lo=np.r_[np.full(nshape,-np.inf),np.concatenate([[-.7,-.7,-5.]+[-np.inf]*(n-3) for n in sizes])]
    hi=np.r_[np.full(nshape,np.inf),np.concatenate([[.7,.7,20.]+[np.inf]*(n-3) for n in sizes])]
    solver_kwargs={}
    if getattr(basis,'iterative_solver',False):
        solver_kwargs={'tr_solver':'lsmr','tr_options':{'atol':1e-6,'btol':1e-6,'maxiter':250,'regularize':False}}
    res=least_squares(lambda p:model_and_jac(p)[0],p,jac=lambda p:model_and_jac(p)[1],bounds=(lo,hi),x_scale='jac',
                      max_nfev=max_nfev,ftol=3e-8,xtol=3e-8,gtol=1e-7,**solver_kwargs)
    a=res.x[:nshape] if free else a_fixed;rows=[];n=0
    for i,r in enumerate(records):
        ni=len(r['x']);v=res.fun[n:n+ni];j=res.jac[n:n+ni];off=int(offsets[i])
        jj=j[:,off:off+sizes[i]];jc,jn=jj[:,:2],jj[:,2:];prof=jc-jn@np.linalg.lstsq(jn,jc,rcond=1e-12)[0]
        cov=np.linalg.inv(prof.T@prof);q=float(v@v);dof=ni-sizes[i]
        # Moment convention is reported separately. Raw nuisance centroid is used
        # for detector-space fit; no physical lens coordinate enters the solve.
        grid=basis.response_grid(r,a);z=grid.sum()/16
        old=r['grid'];oldmean=np.array([(old*basis.gx).sum(),(old*basis.gy).sum()])/old.sum()
        newmean=np.array([(grid*basis.gx).sum(),(grid*basis.gy).sum()])/grid.sum()
        delta=newmean-oldmean
        center=res.x[off:off+2]+r['origin'];rows.append({'name':r['name'],'chi2':q,'dof_conditional_on_shape':dof,
            'p_conditional_on_shape':float(chi2.sf(q,dof)),'centroid_raw_pixel':center.tolist(),
            'centroid_gauge_pixel':(center+delta).tolist(),'moment_gauge_offset_pixel':delta.tolist(),
            'flux_normalization_grid_integral':float(z),'min_response_grid':float(grid.min()),
            'covariance_pixel2_conditional_on_shape':cov.tolist(),'restored_sat_pixels':r['sat_restored']})
        n+=ni
    return {'shape':a,'parameters':res.x,'fits':rows,'success':bool(res.success),'nfev':res.nfev,
            'message':res.message,'shape_penalty':float((a/prior_sd)@(a/prior_sd)) if free else None}

def serial(f):return {k:v for k,v in f.items() if k not in ['shape','parameters']}

def main():
    recs=get_fields();basis=Basis(recs[0]['grid']);report={'status':'DEVELOPMENT_FIELD_CALIBRATION_NO_PHYSICAL_VALIDATION',
        'mask':'Reject DNU/nonfinite/ERR<=0; SAT alone retained for valid calibrated early-group slope.',
        'positive_psf':'A constant pedestal makes the signed empirical template positive; baseline free flux/background absorbs this reparameterization. Shared smooth log-response keeps response positive.',
        'gauge':'Training log-shape removes first-order flux and two moments. Final reported centroids add the exact grid-moment difference and flux normalizes the positive grid; covariance here remains conditional on learned shape.',
        'shape_parameters':basis.null.shape[1],'shape_prior_sd':.2,'baseline':None,'folds':[]}
    zero=np.zeros(basis.null.shape[1]);baseline=fit(recs,basis,a_fixed=zero);report['baseline']=serial(baseline)
    print(json.dumps({'baseline':[r['chi2']/r['dof_conditional_on_shape'] for r in baseline['fits']]}),flush=True)
    for hold in range(4):
        train=fit([r for k,r in enumerate(recs) if k!=hold],basis);test=fit([recs[hold]],basis,a_fixed=train['shape'])
        row={'held_exposure':hold+1,'train':serial(train),'test':serial(test)};report['folds'].append(row)
        print(json.dumps({'held':hold+1,'train_success':train['success'],'train_nfev':train['nfev'],'test':test['fits']}),flush=True)
        (HERE/'positive_field_results.json').write_text(json.dumps(report,indent=2))
    allfit=fit(recs,basis,max_nfev=100);report['all_field_fit']=serial(allfit)
    np.savez(HERE/'positive_field_shape.npz',shape=allfit['shape'],nodes=basis.nodes,null=basis.null,width=basis.width)
    (HERE/'positive_field_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'all_field_fit':serial(allfit)}),flush=True)

if __name__=='__main__':main()
