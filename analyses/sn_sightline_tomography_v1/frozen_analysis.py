"""Frozen full-covariance SN validation and one scalar optical-path reconstruction.

Run: python run_validation.py. Inputs and analysis choices are bundled locally.
No axis, cutoff, fractions, background parameters, or transfer curves are fitted.
"""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import sys, json, hashlib, platform, re
import numpy as np
import pandas as pd
import scipy
import astropy
from scipy.linalg import cholesky, cho_solve, solve_triangular, svd, qr
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
from scipy.stats import chi2
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
IN=ROOT/'inputs'
sys.path.insert(0,str(IN))
import background
import bhsm_optical_background_reconstruction as optical
K=.018
AXIS=(211.48,-12.81)
PB=np.array([.03,.05,.10,.15,.30,.50,1.0])
DB=np.array([.10,.20,.40,.60,.80,1.20])
TB=np.array([.01,.03,.05,.10,.15,.30,.50,1.0])
EDGES=np.r_[0.,TB]
MODELS=['null','local','global','composite']
COLORS={'null':'#666666','local':'#187e68','global':'#c54b44','composite':'#7165ae'}
CHECKS={}
ARRAYS={}

def savejson(path,obj):
    def convert(x):
        if isinstance(x,np.ndarray): return x.tolist()
        if isinstance(x,np.generic): return x.item()
        raise TypeError(type(x).__name__)
    path.write_text(json.dumps(obj,indent=2,default=convert,allow_nan=False),encoding='utf-8')

def check(name,error,tolerance):
    CHECKS[name]=dict(error=float(error),tolerance=float(tolerance),passed=bool(error<tolerance))
    assert error<tolerance,(name,error,tolerance)

def vectors(ra,dec):
    a,d=np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])

def nuisance(z,edges):
    z=np.asarray(z)
    cols=[]
    for a,b in zip(edges[:-1],edges[1:]):
        take=(z>=a)&(z<b)
        if take.sum(): cols.extend([take.astype(float),take*(z-(a+b)/2)/(b-a)])
    return np.column_stack(cols)

class Gaussian:
    def __init__(self,C):
        self.C=np.asarray(C)
        self.L=cholesky(self.C,lower=True,check_finite=True)
    def white(self,a): return solve_triangular(self.L,a,lower=True,check_finite=False)
    def fit(self,X,y):
        Xw=self.white(X)
        scales=np.linalg.norm(Xw,axis=0)
        if np.any(scales==0): raise ValueError('Empty design column')
        u,s,vt=svd(Xw/scales,full_matrices=False,check_finite=False)
        if s[-1]<=s[0]*1e-12: raise ValueError('Rank deficient design')
        R=(vt.T/s)@u.T/scales[:,None]
        H=solve_triangular(self.L.T,R.T,lower=False,check_finite=False).T
        cov=R@R.T
        beta=H@y
        return dict(beta=beta,cov=cov,H=H,condition=float(s[0]/s[-1]),singular=s)
    def score(self,r,N):
        nw=self.white(N)
        q,rr=qr(nw,mode='economic',check_finite=False)
        rank=np.linalg.matrix_rank(rr)
        if rank!=N.shape[1]: raise ValueError('Nuisance rank deficiency')
        rw=self.white(r)
        rw=rw-q@(q.T@rw)
        score=float(rw@rw)
        df=len(r)-rank
        return dict(chi2=score,dof=df,p=float(chi2.sf(score,df)))

def contrast(q,C):
    g=Gaussian(C)
    w=g.white(q)
    score=float(w@w)
    return dict(chi2=score,dof=len(q),p=float(chi2.sf(score,len(q))))

def predictive_cov(C,f,var,cross):
    V=C+var*np.outer(f,f)-np.outer(cross,f)-np.outer(f,cross)
    return (V+V.T)/2

def read_data():
    p=pd.read_csv(IN/'pantheon_Pantheon+SH0ES.dat',sep=r'\s+',dtype={'CID':str})
    raw=np.loadtxt(IN/'pantheon_Pantheon+SH0ES_STAT+SYS.cov')
    n=int(raw[0]); assert n==len(p)==1701 and len(raw)==1+n*n
    C=raw[1:].reshape(n,n)
    check('Pantheon text covariance asymmetry within release rounding',np.max(abs(C-C.T)),1e-7)
    # The published decimal text differs by up to 3e-8 across the diagonal.
    # Symmetrize roundoff only; retain the exact downloaded file for provenance.
    C=(C+C.T)/2
    Gaussian(C)
    p['raw_row']=np.arange(n)
    h=pd.read_csv(IN/'des_DES-Dovekie_HD.csv',sep=r'\s+',comment='#',dtype={'CID':str})
    m=pd.read_csv(IN/'des_DES-Dovekie_Metadata.csv',sep=r'\s+',comment='#',dtype={'CID':str})
    assert h.CID.is_unique and m.CID.is_unique
    h['raw_row']=np.arange(len(h))
    h=h.merge(m[['CID','HOST_RA','HOST_DEC','HOST_ANGSEP']],on='CID',how='left',sort=False,validate='one_to_one')
    assert np.array_equal(h.raw_row,np.arange(len(h))) and h.HOST_RA.notna().all()
    positions=[]
    for f in sorted(IN.glob('*HEAD.FITS.gz')):
        data=fits.getdata(f,1)
        positions.extend(dict(CID=str(row['SNID']).strip(),RA=float(row['RA']),DEC=float(row['DEC']),position_file=f.name) for row in data)
    pos=pd.DataFrame(positions)
    # A SN header has one authoritative position; duplicate headers must agree.
    dup=pos.groupby('CID')[['RA','DEC']].agg(lambda v:np.ptp(v))
    assert (dup.to_numpy()<1e-7).all()
    pos=pos.drop_duplicates('CID')
    h=h.merge(pos,on='CID',how='left',sort=False,validate='one_to_one')
    assert np.array_equal(h.raw_row,np.arange(len(h)))
    assert h[['RA','DEC']].notna().all().all(),h.loc[h.RA.isna(),'CID'].tolist()
    assert ((h.RA>=0)&(h.RA<360)&(abs(h.DEC)<=90)).all()
    d=np.load(IN/'des_STAT+SYS.npz')
    nd=int(d[d.files[0]][0]); assert nd==len(h)==1820
    vals=d[d.files[1]]; assert len(vals)==nd*(nd+1)//2
    W=np.zeros((nd,nd)); W[np.triu_indices(nd)]=vals
    W[np.tril_indices(nd,-1)]=W.T[np.tril_indices(nd,-1)]
    L=cholesky(W,lower=True)
    CD=cho_solve((L,True),np.eye(nd),check_finite=False)
    CD=(CD+CD.T)/2
    check('DES full precision times covariance',np.max(abs(W@CD-np.eye(nd))),1e-8)
    # Demonstrate that marginal precision is the Schur complement, not W[keep,keep].
    a=np.arange(0,nd,40); b=np.setdiff1d(np.arange(nd),a)
    Schur=W[np.ix_(a,a)]-W[np.ix_(a,b)]@cho_solve((cholesky(W[np.ix_(b,b)],lower=True),True),W[np.ix_(b,a)])
    check('DES marginal Schur identity',np.max(abs(Schur@CD[np.ix_(a,a)]-np.eye(len(a)))),1e-8)
    # Resolve same-SN overlap across survey aliases using position as well as CID.
    def norm(s): return re.sub(r'^(sn|at)','',str(s).strip().lower())
    names=set(p.CID.map(norm))
    chord,nearest=cKDTree(vectors(p.RA,p.DEC)).query(vectors(h.RA,h.DEC))
    sep=np.rad2deg(2*np.arcsin(np.clip(chord/2,0,1)))*3600
    dz=abs(h.zHD.to_numpy()-p.zHD.to_numpy()[nearest])
    h['pantheon_overlap']=h.CID.map(norm).isin(names)|((sep<=5)&(dz<.01))
    h['nearest_pantheon_CID']=p.CID.to_numpy()[nearest]
    h['nearest_separation_arcsec']=sep
    h['nearest_delta_zHD']=dz
    h[['raw_row','CID','IDSURVEY','zHD','RA','DEC','position_file','pantheon_overlap','nearest_pantheon_CID','nearest_separation_arcsec','nearest_delta_zHD']].to_csv(ROOT/'DES_ROW_ORDER_AND_OVERLAP.csv',index=False)
    p[['raw_row','CID','IDSURVEY','zHD','RA','DEC']].to_csv(ROOT/'PANTHEON_ROW_ORDER.csv',index=False)
    return p,C,h,CD

def distances():
    sol=background.integrate(2.3986073)
    assert sol.success
    grid=np.linspace(0,2.1,42001)
    E=sol.sol(-np.log1p(grid))[2]
    chi=PchipInterpolator(grid,cumulative_trapezoid(1/E,grid,initial=0))
    def sk(x): return np.sin(np.sqrt(K)*x)/np.sqrt(K)
    def hk(x): return np.sqrt(K)/np.tan(np.sqrt(K)*x)
    den=float(hk(chi(.02))*chi(.02))
    table=pd.read_csv(IN/'R1_SI_preregistered_redshift_table.csv')
    az=table.A_mu_mag_dipole.to_numpy(); zz=table.z.to_numpy()
    anchor=np.interp(.02,zz,az)
    def transfer(z,model):
        z=np.asarray(z)
        local=hk(chi(z))*np.minimum(chi(z),chi(.03))/den
        glob=np.interp(z,zz,az)/anchor
        return dict(null=np.zeros_like(z),local=local,global_=glob,composite=(glob+5*local)/6)['global_' if model=='global' else model]
    def mu(z,zh): return 5*np.log10(299792.458/70*(1+zh)*sk(chi(z)))+25
    check('R1 present E matches supplied rounded V0',abs(sol.y[2,-1]-1),1e-7)
    check('Local transfer normalized at reference',abs(float(transfer(.02,'local'))-1),1e-12)
    check('Global transfer normalized at reference',abs(float(transfer(.02,'global'))-1),1e-12)
    return chi,hk,den,transfer,mu

def fixed_screen(name,frame,C,y,edges,A,varA,cross,transfer):
    z=frame.zHD.to_numpy(); cost=frame.cost.to_numpy()
    N=nuisance(z,edges); g=Gaussian(C)
    D=np.column_stack([((z>=a)&(z<b))*cost for a,b in zip(edges[:-1],edges[1:])])
    fit=g.fit(np.c_[N,D],y)
    P=fit['H'][N.shape[1]:]; obs=P@y; BC=P@C@P.T; cc=P@cross
    check(name+' compression covariance',np.max(abs(BC-fit['cov'][N.shape[1]:,N.shape[1]:])),1e-10)
    check(name+' compression nuisance removal',np.max(abs(P@N)),1e-9)
    check(name+' compression dipole identity',np.max(abs(P@D-np.eye(D.shape[1]))),1e-9)
    ARRAYS[name+'_compression']=P
    ARRAYS[name+'_bin_covariance']=BC
    ARRAYS[name+'_raw_rows']=frame.raw_row.to_numpy()
    ARRAYS[name+'_bin_training_crosscov']=cc
    out=dict(rows=len(frame),nuisance_columns=N.shape[1],bin_covariance=BC,bin_training_crosscov=cc,models={})
    records=[]
    rawnull=g.score(y,N)['chi2']; binnull=contrast(obs,BC)['chi2']
    for model in MODELS:
        f=transfer(z,model)*cost; pred=f*A; gp=P@f
        raw=g.score(y-pred,N)
        raw['delta_chi2_vs_null']=raw['chi2']-rawnull
        raw['log_likelihood_ratio_vs_null']=-raw['delta_chi2_vs_null']/2
        raw['likelihood_ratio_vs_null']=float(np.exp(np.clip(raw['log_likelihood_ratio_vs_null'],-700,700)))
        V=predictive_cov(C,f,varA,cross)
        predictive=Gaussian(V).score(y-pred,N)
        bp=gp*A; bv=predictive_cov(BC,gp,varA,cc)
        conditional=contrast(obs-bp,BC)
        conditional['delta_chi2_vs_null']=conditional['chi2']-binnull
        conditional['likelihood_ratio_vs_null']=float(np.exp(-conditional['delta_chi2_vs_null']/2))
        bs=contrast(obs-bp,bv)
        bs['p_bonferroni_three_alternatives']=min(1.,3*bs['p']) if model!='null' else bs['p']
        out['models'][model]=dict(raw_conditional=raw,raw_predictive=predictive,bin_conditional=conditional,bin_predictive=bs)
        ARRAYS[name+'_'+model+'_predictive_bin_covariance']=bv
        for j,(a,b) in enumerate(zip(edges[:-1],edges[1:])):
            take=(z>=a)&(z<b)
            records.append(dict(dataset=name,model=model,z_min=a,z_max=b,rows=int(take.sum()),z_mean=float(z[take].mean()),A_observed=obs[j],sigma_observed=np.sqrt(BC[j,j]),A_predicted_projected=bp[j],sigma_predictive_contrast=np.sqrt(bv[j,j]),bin_residual=obs[j]-bp[j],joint_bin_chi2=bs['chi2'],joint_bin_dof=bs['dof'],joint_bin_p=bs['p']))
    # Direct projected-precision check without constructing large inverse matrices.
    yn=g.white(y); nw=g.white(N); qrN=qr(nw,mode='economic')[0]
    check(name+' nuisance projection chi2 identity',abs(float(yn@yn-np.sum((qrN.T@yn)**2))-rawnull),1e-7)
    return out,records

def source_design(frame,chi,hk,den,edges):
    # Existing machinery, including its exact closed-S3 Green kernel, reused unchanged.
    tz,tl=optical.shell_operator(chi(frame.zHD),chi(edges),K)
    return frame.cost.to_numpy()[:,None]*hk(chi(frame.zHD))[:,None]*tz/den

def source_fit(frame,C,y,chi,hk,den,edges=EDGES):
    N=nuisance(frame.zHD,TB)
    T=source_design(frame,chi,hk,den,edges)
    g=Gaussian(C); nw=g.white(N); tw=g.white(T)
    q=qr(nw,mode='economic')[0]; tp=tw-q@(q.T@tw)
    spectrum=svd(tp,compute_uv=False)
    condition=float(spectrum[0]/spectrum[-1])
    if condition>=1e8:
        if len(edges)<=3: raise ValueError('Source is not identifiable even after deterministic coarsening')
        return source_fit(frame,C,y,chi,hk,den,np.delete(edges,-2))
    fit=g.fit(np.c_[N,T],y)
    H=fit['H'][N.shape[1]:]; S=fit['cov'][N.shape[1]:,N.shape[1]:]
    s=H@y
    check('source fit covariance '+str(len(frame)),np.max(abs(H@C@H.T-S)),1e-8)
    return dict(source=s,cov=S,H=H,T=T,N=N,spectrum=spectrum,condition=condition,edges=edges,fit_chi2=g.score(y-T@s,N),null_chi2=g.score(y,N))

def source_prediction(train,test,Ctest,ytest,cross,chi,hk,den):
    T=source_design(test,chi,hk,den,train['edges'])
    N=nuisance(test.zHD,TB)
    pred=T@train['source']
    CT=cross@train['H'].T
    V=Ctest+T@train['cov']@T.T-CT@T.T-T@CT.T
    V=(V+V.T)/2
    g=Gaussian(Ctest)
    baseline=g.score(ytest,N)
    plug=g.score(ytest-pred,N)
    plug['delta_chi2_vs_null']=plug['chi2']-baseline['chi2']
    predictive=Gaussian(V).score(ytest-pred,N)
    # Axis-sensitive bin-space transfer as well as all-row residual score.
    z=test.zHD.to_numpy(); columns=[]; used=[]
    for a,b in zip(TB[:-1],TB[1:]):
        take=(z>=a)&(z<b)
        if take.sum(): columns.append(take*test.cost.to_numpy()); used.append([a,b])
    D=np.column_stack(columns)
    bf=g.fit(np.c_[N,D],ytest)
    P=bf['H'][N.shape[1]:]
    bpred=P@pred; bobs=P@ytest
    bins=contrast(bobs-bpred,P@V@P.T)
    return dict(rows=len(test),raw_null=baseline,raw_plugin=plug,raw_predictive=predictive,bin_predictive=bins,bins=used,observed=bobs,predicted=bpred,observed_cov=P@Ctest@P.T,predictive_cov=P@V@P.T)

def tomography(p,Cp,yp,d,Cd,yd,chi,hk,den,transfer):
    # Only one observational source basis; CV fits are validation of this same estimator.
    pf=source_fit(p,Cp,yp,chi,hk,den)
    edges=pf['edges']; s=pf['source']; S=pf['cov']
    ARRAYS['source_covariance']=S; ARRAYS['source_estimator']=pf['H']; ARRAYS['source_raw_rows']=p.raw_row.to_numpy()
    ARRAYS['source_singular_values']=pf['spectrum']
    # Finite-volume equivalent source templates preserve integrated endpoint transfer.
    ce=chi(edges); I=np.zeros(len(edges))
    I[1:]=den*transfer(edges[1:],'global')/hk(ce[1:])
    glob=np.diff(I)/np.diff(ce)
    local=(edges[1:]<=.03+1e-12).astype(float)
    U=np.c_[glob,local]
    dec=Gaussian(S).fit(U,s)
    coeff=dec['beta']; persistent=glob*coeff[0]; localized=local*coeff[1]
    rem=s-persistent-localized
    rec=pd.DataFrame(dict(z_min=edges[:-1],z_max=edges[1:],chi_min=ce[:-1],chi_max=ce[1:],s_effective_mag=s,sigma=np.sqrt(np.diag(S)),global_equivalent_template=glob,local_template=local,coherent_component=persistent,localized_component=localized,residual_component=rem))
    rec.to_csv(ROOT/'BHSM_OPTICAL_BACKGROUND_RECONSTRUCTION_REAL.csv',index=False)
    pd.DataFrame(S,index=[f'{a:g}-{b:g}' for a,b in zip(edges[:-1],edges[1:])],columns=[f'{a:g}-{b:g}' for a,b in zip(edges[:-1],edges[1:])]).to_csv(ROOT/'BHSM_OPTICAL_BACKGROUND_SOURCE_COVARIANCE.csv')
    reconstruction=dict(rows=len(p),edges=edges,source=s,covariance=S,singular_values=pf['spectrum'],condition=pf['condition'],regularization_lambda=0.,fit_chi2=pf['fit_chi2'],null_chi2=pf['null_chi2'],coherent_local_coefficients=coeff,coherent_local_covariance=dec['cov'],decomposition_residual_chi2=Gaussian(S).score(s,U))
    # Correction: in-sample fit used eight source parameters in addition to nuisances.
    reconstruction['fit_chi2']['dof']-=len(s)
    reconstruction['fit_chi2']['p']=float(chi2.sf(reconstruction['fit_chi2']['chi2'],reconstruction['fit_chi2']['dof']))
    reconstruction['delta_chi2_fitted_source_vs_null']=pf['fit_chi2']['chi2']-pf['null_chi2']['chi2']
    reconstruction['nested_zero_source_p']=float(chi2.sf(-reconstruction['delta_chi2_fitted_source_vs_null'],len(s)))
    folds=[]
    sectors=np.floor((p.RA.to_numpy()%360)/90).astype(int)
    # Guard against splitting the same supernova across sky sectors.
    assert pd.DataFrame({'CID':p.CID,'sector':sectors}).groupby('CID').sector.nunique().max()==1
    for k in range(4):
        it=np.flatnonzero(sectors!=k); ih=np.flatnonzero(sectors==k)
        train=source_fit(p.iloc[it],Cp[np.ix_(it,it)],yp[it],chi,hk,den)
        pred=source_prediction(train,p.iloc[ih],Cp[np.ix_(ih,ih)],yp[ih],Cp[np.ix_(ih,it)],chi,hk,den)
        pred.update(ra_min=90*k,ra_max=90*(k+1),training_rows=len(it),training_condition=train['condition'],source=train['source'],source_covariance=train['cov'],edges=train['edges'])
        folds.append(pred)
    dt=source_prediction(pf,d,Cd,yd,np.zeros((len(d),len(p))),chi,hk,den)
    # Joint-source diagnostic; independent per-survey nuisance columns. No feedback into predictions.
    Td=source_design(d,chi,hk,den,edges); Nd=nuisance(d.zHD,TB)
    def projected_normal(C,T,N,y):
        g=Gaussian(C); tw=g.white(T); yw=g.white(y); nw=g.white(N); q=qr(nw,mode='economic')[0]
        tw=tw-q@(q.T@tw); yw=yw-q@(q.T@yw)
        return tw.T@tw,tw.T@yw
    Ip,bp=projected_normal(Cp,pf['T'],pf['N'],yp)
    Id,bd=projected_normal(Cd,Td,Nd,yd)
    Lj=cholesky(Ip+Id,lower=True); Sj=cho_solve((Lj,True),np.eye(len(s))); sj=cho_solve((Lj,True),bp+bd)
    delta=sj-s; delta_cov=S-Sj
    # Rank-aware statistic: DES may not identify individual foreground shells.
    ev,vec=np.linalg.eigh((delta_cov+delta_cov.T)/2); keep=ev>ev.max()*1e-9
    stat=float(np.sum((vec[:,keep].T@delta)**2/ev[keep])); df=int(keep.sum())
    reconstruction['joint_DES_source_diagnostic']=dict(source=sj,covariance=Sj,shift=delta,shift_covariance=delta_cov,shift_chi2=stat,dof=df,p=float(chi2.sf(stat,df)),assumption='Block-diagonal cross-survey covariance after same-SN removal; remaining common systematics unavailable.')
    reconstruction['angular_cross_validation']=folds
    reconstruction['DES_frozen_source_transfer']=dt
    return reconstruction,rec

def make_plots(bins,result,rec,transfer):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for ax,name,title in zip(axes,['Pantheon_heldout','DES_overlap_clean'],['Pantheon held-out','DES transfer (shared SNe removed)']):
        dat=bins[(bins.dataset==name)&(bins.model=='null')]
        ax.errorbar(dat.z_mean,dat.A_observed,yerr=dat.sigma_observed,fmt='o',color='black',label='Full-covariance bin GLS',capsize=3)
        for m in MODELS:
            mm=bins[(bins.dataset==name)&(bins.model==m)]
            ax.plot(mm.z_mean,mm.A_predicted_projected,'s-',markersize=3,color=COLORS[m],label=m)
        ax.set(xlabel='zHD',ylabel='Fixed-axis amplitude (mag)',title=title)
        ax.axhline(0,color='grey',lw=.5)
    axes[0].legend(fontsize=8); fig.savefig(ROOT/'FULL_COVARIANCE_AMPLITUDES_AND_DES_TRANSFER.png'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True)
    for j,name in enumerate(['Pantheon_heldout','DES_overlap_clean']):
        values=[result['fixed_tests'][name]['models'][m]['raw_conditional']['delta_chi2_vs_null'] for m in MODELS[1:]]
        ax.bar(np.arange(3)+(j-.5)*.35,values,width=.35,label=name.replace('_',' '))
    ax.axhline(0,color='black',lw=.8); ax.set_xticks(range(3),MODELS[1:]); ax.set(ylabel='Raw-row conditional chi-square minus null',title='Frozen predictions; negative favors model'); ax.legend(); fig.savefig(ROOT/'MODEL_VS_NULL_DELTA_CHI2.png'); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    x=np.arange(len(rec)); labels=[f'{a:g}-{b:g}' for a,b in zip(rec.z_min,rec.z_max)]
    axes[0].errorbar(x,rec.s_effective_mag,yerr=rec.sigma,fmt='o',color='black',capsize=3,label='GLS source ± 1 sigma')
    axes[0].plot(x,rec.coherent_component,label='Coherent equivalent')
    axes[0].plot(x,rec.localized_component,label='Localized')
    axes[0].plot(x,rec.residual_component,':',label='Residual')
    axes[0].set_xticks(x,labels,rotation=45); axes[0].set(ylabel='Effective source (mag)',xlabel='Redshift shell'); axes[0].legend(fontsize=8)
    cov=np.asarray(result['tomography']['covariance']); corr=cov/np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
    im=axes[1].imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r'); axes[1].set(title='Source correlation (full covariance in CSV)',xlabel='Shell index',ylabel='Shell index'); fig.colorbar(im,ax=axes[1]); fig.savefig(ROOT/'RECONSTRUCTED_SOURCE_AND_COVARIANCE.png'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True)
    t=result['tomography']; ax.semilogy(np.arange(1,len(t['singular_values'])+1),t['singular_values'],'o-',label=f"All Pantheon: condition {t['condition']:.1f}")
    ax.set(xlabel='Singular mode',ylabel='Whitened nuisance-projected singular value',title='Unregularized effective path-source inversion'); ax.legend(); fig.savefig(ROOT/'SOURCE_SINGULAR_VALUES.png'); plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4.2),constrained_layout=True)
    t=result['tomography']['DES_frozen_source_transfer']; x=np.mean(t['bins'],axis=1); sig=np.sqrt(np.diag(t['observed_cov'])); psig=np.sqrt(np.maximum(0,np.diag(t['predictive_cov'])-sig**2))
    ax.errorbar(x,t['observed'],yerr=sig,fmt='o',color='black',label='DES bin GLS')
    ax.errorbar(x,t['predicted'],yerr=psig,fmt='s-',capsize=3,label='Pantheon source prediction ± propagated uncertainty')
    ax.axhline(0,color='grey',lw=.7); ax.set(xlabel='zHD bin midpoint',ylabel='Projected fixed-axis amplitude (mag)',title='Real-source transfer: no DES source retuning'); ax.legend(fontsize=8); fig.savefig(ROOT/'TOMOGRAPHY_DES_TRANSFER.png'); plt.close(fig)

def main():
    p,Cp,d,Cd=read_data()
    chi,hk,den,transfer,mu=distances()
    axis=vectors(*[np.array([v]) for v in AXIS])[0]
    for frame in [p,d]: frame['cost']=vectors(frame.RA,frame.DEC)@axis
    yp=p.MU_SH0ES.to_numpy()-mu(p.zHD.to_numpy(),p.zHEL.to_numpy())
    yd=d.MU.to_numpy()-mu(d.zHD.to_numpy(),d.zHEL.to_numpy())
    it=np.flatnonzero((p.zHD>=.01)&(p.zHD<.03))
    ih=np.flatnonzero((p.zHD>=.03)&(p.zHD<1))
    X=np.c_[np.ones(len(it)),p.zHD.to_numpy()[it]-.02,p.cost.to_numpy()[it]]
    fit=Gaussian(Cp[np.ix_(it,it)]).fit(X,yp[it])
    A=float(fit['beta'][-1]); varA=float(fit['cov'][-1,-1]); w=fit['H'][-1]
    ARRAYS['training_amplitude_estimator']=w; ARRAYS['training_raw_rows']=it
    result=dict(schema='BHSM-official-full-covariance-validation-v1',analysis_plan_sha256=hashlib.sha256((ROOT/'ANALYSIS_PLAN.md').read_bytes()).hexdigest(),software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,pandas=pd.__version__,astropy=astropy.__version__,matplotlib=matplotlib.__version__),parameters=dict(Omega_m0=.31,Omega_r0=9e-5,Omega_k0=-.018,lam=1,q=3,V0=2.3986073,axis_RA_DEC_deg=AXIS,z_ref=.02,z_c=.03,global_fraction=1/6,local_fraction=5/6),input_summary=dict(pantheon_raw_rows=len(p),pantheon_unique_CIDs=p.CID.nunique(),des_raw_rows=len(d),missing_handoff_file='BHSM_real_data_prediction_screen.json',DES_overlaps_all=int(d.pantheon_overlap.sum())),training=dict(rows=len(it),unique_CIDs=p.iloc[it].CID.nunique(),A_mag=A,sigma_A_mag=np.sqrt(varA),monopole_and_slope=fit['beta'][:2],parameter_covariance=fit['cov']),fixed_tests={})
    bins=[]
    tests=[('Pantheon_heldout',p.iloc[ih],Cp[np.ix_(ih,ih)],yp[ih],PB,Cp[np.ix_(ih,it)]@w)]
    dbase=(d.IDSURVEY==10)&(d.zHD>=.1)&(d.zHD<1.2)
    for name,mask in [('DES_overlap_clean',dbase&~d.pantheon_overlap),('DES_all_release_survey10',dbase)]:
        ix=np.flatnonzero(mask)
        tests.append((name,d.iloc[ix],Cd[np.ix_(ix,ix)],yd[ix],DB,np.zeros(len(ix))))
    result['input_summary']['DES_highz_overlaps_removed']=int((dbase&d.pantheon_overlap).sum())
    for name,frame,C,y,edges,cross in tests:
        out,rows=fixed_screen(name,frame,C,y,edges,A,varA,cross,transfer)
        result['fixed_tests'][name]=out; bins.extend(rows)
    bins=pd.DataFrame(bins); bins.to_csv(ROOT/'BHSM_FULL_COVARIANCE_BIN_RESULTS.csv',index=False)
    # Seeded Monte Carlo verifies propagation of calibration uncertainty plus train/test correlation.
    P=ARRAYS['Pantheon_heldout_compression']; BC=ARRAYS['Pantheon_heldout_bin_covariance']; cc=ARRAYS['Pantheon_heldout_bin_training_crosscov']
    joint=np.block([[np.array([[varA]]),cc[None,:]],[cc[:,None],BC]])
    rng=np.random.default_rng(20260920); draw=rng.standard_normal((200000,len(cc)+1))@cholesky(joint,lower=True).T
    f=P@(p.cost.to_numpy()[ih]*transfer(p.zHD.to_numpy()[ih],'composite'))
    q=draw[:,1:]-draw[:,:1]*f
    exact=predictive_cov(BC,f,varA,cc); empirical=np.cov(q,rowvar=False)
    check('Seeded calibration-contrast covariance Monte Carlo relative Frobenius error',np.linalg.norm(empirical-exact)/np.linalg.norm(exact),.02)
    vals=np.sum(solve_triangular(cholesky(exact,lower=True),q.T,lower=True)**2,axis=0)
    check('Gaussian contrast nominal one-percent tail calibration',abs(float(np.mean(vals>chi2.isf(.01,len(cc))))-.01),.001)
    ip=np.flatnonzero((p.zHD>=.01)&(p.zHD<1))
    ides=np.flatnonzero(dbase&~d.pantheon_overlap&(d.zHD<1))
    t,rec=tomography(p.iloc[ip],Cp[np.ix_(ip,ip)],yp[ip],d.iloc[ides],Cd[np.ix_(ides,ides)],yd[ides],chi,hk,den,transfer)
    result['tomography']=t
    result['checks']=CHECKS
    result['limitations']=['Retrospective test of releases already used in the exploratory screen; not an independent prospective discovery.','No public Pantheon-DES cross-survey systematic covariance. DES results assume independence after explicit same-SN removal; no combined significance is claimed.','Published distances include survey-wide nuisance/bias calibration; held-out here means no BHSM amplitude/shape retuning, not independent retraining of the survey reductions.','One scalar distance observable cannot separate Z,F,G1,G2, prove time persistence, or establish the full n=2 angular pattern from a fixed cos(theta) template.','Eight fitted tomography shells are an observational reconstruction, not a BHSM prediction. Passing predictive goodness-of-fit with large source uncertainty is not evidence.']
    savejson(ROOT/'BHSM_FULL_COVARIANCE_VALIDATION.json',result)
    np.savez_compressed(ROOT/'GLS_OPERATORS_AND_COVARIANCES.npz',**ARRAYS)
    make_plots(bins,result,rec,transfer)
    print(json.dumps(dict(training=result['training'],tests={name:{m:entry['bin_predictive'] for m,entry in val['models'].items()} for name,val in result['fixed_tests'].items()},tomography=dict(condition=t['condition'],fit=t['fit_chi2'],source=t['source'].tolist(),source_sigma=np.sqrt(np.diag(t['covariance'])).tolist(),CV=[dict(rows=f['rows'],raw=f['raw_predictive'],bins=f['bin_predictive'],plugin_delta=f['raw_plugin']['delta_chi2_vs_null']) for f in t['angular_cross_validation']],DES=t['DES_frozen_source_transfer']['bin_predictive'],joint_shift=t['joint_DES_source_diagnostic']['p'])),default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item(),indent=2))

if __name__=='__main__': main()
