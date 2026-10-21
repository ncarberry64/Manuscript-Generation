"""Consumed-data coherent sky and shared response uncertainty application.

The supplied positive field response is the expansion point. No lens prediction
enters this module. Native pixel-integrated ePSFs are sampled/convolved once.
All covariance/rank operations use one SVD convention on dimensionless columns.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares
from scipy.signal import fftconvolve
from scipy.stats import chi2

SVD_RCOND=1e-10


def svd_fit(matrix,rhs=None,scales=None,rcond=SVD_RCOND):
    j=np.asarray(matrix,float)
    scale=np.ones(j.shape[1]) if scales is None else np.asarray(scales,float)
    # Complete V is needed when the problem has fewer rows than columns;
    # otherwise a thin SVD silently omits exact null directions.
    u,s,v=np.linalg.svd(j*scale,full_matrices=j.shape[0]<j.shape[1])
    keep=s>rcond*(s[0] if len(s) else 1.)
    measured=v[:len(s)]
    inverse=(measured[keep].T/s[keep])@u[:,:len(s)][:,keep].T
    covariance=(measured[keep].T/s[keep]**2)@measured[keep]
    covariance=scale[:,None]*covariance*scale[None,:]
    result={'covariance':covariance,'rank':int(keep.sum()),'columns':j.shape[1],
            'singular_values':s,'null_vectors_dimensionless':np.vstack([measured[~keep],v[len(s):]]),
            'full_column_rank':bool(keep.sum()==j.shape[1]),
            'rcond':rcond,'column_scales':scale}
    if rhs is not None:result['solution']=scale*(inverse@np.asarray(rhs,float))
    return result


def profile_columns(nuisance,columns,residual=None,scales=None):
    n=np.asarray(nuisance,float);b=np.asarray(columns,float)
    scale=np.ones(n.shape[1]) if scales is None else np.asarray(scales,float)
    u,s,_=np.linalg.svd(n*scale,full_matrices=False)
    keep=s>SVD_RCOND*(s[0] if len(s) else 1.)
    q=u[:,keep]
    answer=b-q@(q.T@b)
    return (answer,None if residual is None else residual-q@(q.T@residual))


def propagated_cut_covariance(jpoint,jcal,calibration_covariance,prior_rows=None,scales=None):
    a=np.asarray(jpoint,float);b=np.asarray(jcal,float)
    prior=np.zeros((0,a.shape[1])) if prior_rows is None else np.asarray(prior_rows,float)
    fit=svd_fit(np.vstack([a,prior]),scales=scales)
    conditional=fit['covariance']
    sensitivity=-conditional@a.T@b
    marginal=conditional+sensitivity@calibration_covariance@sensitivity.T
    return conditional,marginal,sensitivity,fit


def conservative_charge_laplacian(charge):
    """Conserved neighbour exchange for electrons (or a common-gain unit).

    Raw DN with spatially varying gain must first be converted to electrons.
    This generic operator is not a fitted or transferred history response.
    """
    q=np.asarray(charge,float);result=np.zeros_like(q)
    dx=q[:,1:]-q[:,:-1];result[:,:-1]+=dx;result[:,1:]-=dx
    dy=q[1:]-q[:-1];result[:-1]+=dy;result[1:]-=dy
    return result


def history_rate(read_means,weights):
    """Explicit interface only: actual per-read means and slope weights required."""
    q,w=np.asarray(read_means,float),np.asarray(weights,float)
    if q.shape!=w.shape or q.ndim!=3 or not np.isfinite(q).all() or not np.isfinite(w).all():
        raise ValueError('Matching finite read,pixel-y,pixel-x arrays required')
    return np.sum(q*w,axis=0)


class PositiveResponse:
    def __init__(self,base_grid,nodes,null,width,eta,oversampling=4):
        self.over=oversampling;self.nodes=np.asarray(nodes);self.null=np.asarray(null)
        self.width=width;self.eta=np.asarray(eta)
        base=np.asarray(base_grid,float)
        self.pedestal=max(0.,-float(base.min()))*1.1+1e-8
        self.base=base+self.pedestal
        self.coord=(np.arange(base.shape[0])-(base.shape[0]-1)/2)/self.over
        yy,xx=np.meshgrid(self.coord,self.coord,indexing='ij');self.gx=xx;self.gy=yy
        self.base_spline=RectBivariateSpline(self.coord,self.coord,self.base.T)
        self.base_moment=np.array([(self.base*xx).sum(),(self.base*yy).sum()])/self.base.sum()
        raw=np.exp(-((xx.ravel()[:,None]-self.nodes[None,:,0])**2+
                     (yy.ravel()[:,None]-self.nodes[None,:,1])**2)/(2*width**2))
        self.basis=raw@self.null
        self.update(eta)

    def update(self,eta):
        self.eta=np.asarray(eta,float).copy();log=self.basis@self.eta
        if np.max(abs(log))>100:raise ValueError('Excessive log-response; reject coupled step')
        grid=self.base*np.exp(log.reshape(self.base.shape))
        self.normalization=float(grid.sum()/self.over**2)
        self.grid=grid/self.normalization
        self.node_coefficients=self.null@self.eta
        self.moment=np.array([(self.grid*self.gx).sum(),(self.grid*self.gy).sum()])/self.grid.sum()
        weighted=self.grid.ravel()[:,None]*self.basis
        dz=weighted.sum(axis=0)/self.grid.sum()
        self.dlog_normalization=dz
        self.derivative_grid=(weighted-self.grid.ravel()[:,None]*dz).T.reshape((-1,*grid.shape))
        self.moment_derivative=np.vstack([(self.gx.ravel()[:,None]*weighted).sum(axis=0),
                                         (self.gy.ravel()[:,None]*weighted).sum(axis=0)])/self.grid.sum()-self.moment[:,None]*dz

    def values(self,x,y,derivative=False):
        x,y=np.broadcast_arrays(np.asarray(x),np.asarray(y))
        inside=(abs(x)<=self.coord[-1])&(abs(y)<=self.coord[-1])
        dx=x.ravel()[:,None]-self.nodes[None,:,0]
        dy=y.ravel()[:,None]-self.nodes[None,:,1]
        raw=np.exp(-(dx*dx+dy*dy)/(2*self.width**2))
        correction=np.exp(raw@self.node_coefficients).reshape(x.shape)/self.normalization
        base=self.base_spline.ev(x,y)
        # Match the delivered family at the evaluation coordinates: spline of
        # the positive base TIMES exp(B(x,y) eta), not spline of their product.
        p=np.where(inside,base*correction,0.)
        if not derivative:return p
        logdx=((-dx/self.width**2*raw)@self.node_coefficients).reshape(x.shape)
        logdy=((-dy/self.width**2*raw)@self.node_coefficients).reshape(x.shape)
        return p,np.where(inside,correction*(self.base_spline.ev(x,y,dx=1)+base*logdx),0.),np.where(inside,correction*(self.base_spline.ev(x,y,dy=1)+base*logdy),0.)

    def cal_values(self,x,y):
        x,y=np.asarray(x),np.asarray(y)
        raw=np.exp(-((x.ravel()[:,None]-self.nodes[None,:,0])**2+
                     (y.ravel()[:,None]-self.nodes[None,:,1])**2)/(2*self.width**2))
        return self.values(x,y).ravel()[:,None]*(raw@self.null-self.dlog_normalization)


def build_scene_projection(record,sky_at_detector,anchor_arcsec,width_arcsec=.10,node_spacing_arcsec=.18):
    """Shared observed-sky Gaussian patch through exact GWCS and finite ePSF.

    The input callback returns exact tangent sky coordinates for detector points.
    SCI is the native surface-brightness product: constant sky remains constant
    under convolution. There is no extra PAM, pixel-area multiplication or second
    native-pixel integration. The ePSF already contains pixel integration.
    """
    x0,y0,x1,y1=record['bounds'];over=record['response'].over;pad=14
    fx=np.arange((x0-pad)*over,(x1-1+pad)*over+1)/over
    fy=np.arange((y0-pad)*over,(y1-1+pad)*over+1)/over
    yy,xx=np.meshgrid(fy,fx,indexing='ij');sky=np.asarray(sky_at_detector(xx,yy))
    centers=np.array([(a,b) for b in [-node_spacing_arcsec,0.,node_spacing_arcsec]
                      for a in [-node_spacing_arcsec,0.,node_spacing_arcsec]])+anchor_arcsec
    images=np.exp(-np.sum((sky[None]-centers[:,None,None,:])**2,axis=-1)/(2*width_arcsec**2))
    # Each Gaussian belongs to one explicitly finite observed-sky patch. There
    # are no data here fixing a full host/lens-light scene beyond these patches.
    images*=np.sum((sky-anchor_arcsec)**2,axis=-1)[None] <= .55**2
    record['sky_fine']=images;record['fine_shape']=xx.shape
    record['fine_ix']=np.rint((record['x']-fx[0])*over).astype(int)
    record['fine_iy']=np.rint((record['y']-fy[0])*over).astype(int)
    refresh_scene_projection(record)


def refresh_scene_projection(record):
    response=record['response'];fine=record['sky_fine']
    blurred=fftconvolve(fine,response.grid[None],mode='same',axes=(-2,-1))/response.over**2
    record['scene_matrix']=blurred[:,record['fine_iy'],record['fine_ix']].T


def scene_calibration_columns(record,coefficients):
    scene=np.tensordot(coefficients,record['sky_fine'],axes=(0,0))
    derivatives=record['response'].derivative_grid
    result=[]
    for start in range(0,len(derivatives),24):
        batch=derivatives[start:start+24]
        # scipy's mode='same' follows the FIRST input shape, including batch
        # axes. Explicitly broadcast the scene so no calibration modes vanish.
        scenes=np.broadcast_to(scene,(len(batch),*scene.shape))
        blurred=fftconvolve(scenes,batch,mode='same',axes=(-2,-1))/record['response'].over**2
        result.append(blurred[:,record['fine_iy'],record['fine_ix']].T)
    return np.column_stack(result)


class SceneProblem:
    def __init__(self,records,scene_fixed=None,scene_prior_sd=10.):
        self.records=records;self.n=len(records);self.exposures=sorted(set(r['exposure'] for r in records))
        self.eindex={e:k for k,e in enumerate(self.exposures)}
        self.scene_fixed=scene_fixed;self.scene_prior_sd=scene_prior_sd
        self.nsky=36 if scene_fixed is None else 0
        self.sky_start=3*self.n;self.bg_start=self.sky_start+self.nsky
        self.size=self.bg_start+3*len(self.exposures)
        self.scales=np.r_[np.tile([1.,1.,1.],self.n),np.full(self.nsky,10.),np.full(3*len(self.exposures),10.)]

    def initial(self):
        p=np.zeros(self.size)
        for i,r in enumerate(self.records):p[3*i:3*i+3]=[*r['guess'],np.log(r['flux_guess'])]
        if self.nsky:p[self.sky_start:self.bg_start]=.01
        for e,k in self.eindex.items():
            p[self.bg_start+3*k]=np.median(np.concatenate([r['data'] for r in self.records if r['exposure']==e]))
        return p

    def sky(self,p):
        return np.asarray(self.scene_fixed) if self.scene_fixed is not None else p[self.sky_start:self.bg_start]

    def bounds(self):
        lower=np.full(self.size,-np.inf);upper=np.full(self.size,np.inf)
        for i,r in enumerate(self.records):
            lower[3*i:3*i+3]=[*np.asarray(r['guess'])-.7,-5.]
            upper[3*i:3*i+3]=[*np.asarray(r['guess'])+.7,20.]
        if self.nsky:lower[self.sky_start:self.bg_start]=0.
        return lower,upper

    def evaluate(self,p,with_cal=False):
        rows=[];jrows=[];calrows=[];raw=[];sky=self.sky(p)
        for i,r in enumerate(self.records):
            cx,cy,lf=p[3*i:3*i+3];f=np.exp(lf);u,v=r['x']-cx,r['y']-cy
            profile,dx,dy=r['response'].values(u,v,True)
            coeff=sky[9*r['patch']:9*r['patch']+9]
            bg=self.bg_start+3*self.eindex[r['exposure']]
            mean=f*profile+r['scene_matrix']@coeff+r['detector_background']@p[bg:bg+3]
            residual=(mean-r['data'])/r['error'];j=np.zeros((len(mean),self.size))
            j[:,3*i:3*i+3]=np.c_[-f*dx,-f*dy,f*profile]/r['error'][:,None]
            if self.nsky:j[:,self.sky_start+9*r['patch']:self.sky_start+9*r['patch']+9]=r['scene_matrix']/r['error'][:,None]
            j[:,bg:bg+3]=r['detector_background']/r['error'][:,None]
            rows.append(residual);jrows.append(j);raw.append({'name':r['name'],'chi2':float(residual@residual),'pixels':len(mean),
                 'centroid_xy_pixel':[float(cx),float(cy)],'flux':float(f),'background':p[bg:bg+3].tolist(),
                 'response_grid_min':float(r['response'].grid.min()),'sampled_point_response_min':float(profile.min()),
                 'response_moment_pixel':r['response'].moment.tolist()})
            if with_cal:calrows.append((f*r['response'].cal_values(u,v)+scene_calibration_columns(r,coeff))/r['error'][:,None])
        return np.concatenate(rows),np.concatenate(jrows),None if not with_cal else np.concatenate(calrows),raw

    def prior(self,p):
        j=np.zeros((self.nsky,self.size))
        if self.nsky:j[:,self.sky_start:self.bg_start]=np.eye(self.nsky)/self.scene_prior_sd
        return p[self.sky_start:self.bg_start]/self.scene_prior_sd,j

    def fit(self,initial=None,max_nfev=80):
        def fun(p):return np.r_[self.evaluate(p)[0],self.prior(p)[0]]
        def jac(p):return np.vstack([self.evaluate(p)[1],self.prior(p)[1]])
        p=self.initial() if initial is None else np.asarray(initial).copy();lo,hi=self.bounds()
        return least_squares(fun,p,jac=jac,bounds=(lo,hi),x_scale=self.scales,
                             ftol=2e-8,xtol=2e-8,gtol=1e-6,max_nfev=max_nfev)


def point_covariance_block(problem,covariance):
    index=np.array([3*i+j for i in range(problem.n) for j in [0,1]])
    return covariance[np.ix_(index,index)]


def failure_summary(residual,jacobian):
    rank=svd_fit(jacobian)['rank'];q=float(residual@residual);dof=len(residual)-rank
    return {'chi2':q,'native_pixels':len(residual),'data_jacobian_rank':rank,
            'nominal_conditional_dof':dof,'conditional_p_value':float(chi2.sf(q,max(dof,1))),
            'absolute_adequacy_pass':bool(chi2.sf(q,max(dof,1))>=.001),
            'physical_validation_pass':False}
