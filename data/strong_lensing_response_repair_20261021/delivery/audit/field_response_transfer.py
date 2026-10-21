"""Development-only shared ePSF residual calibration; no lens scores in fitting.

Archived consumed F115W fields only train the common correction. Four explicit
leave-one-exposure-out checks assess whether the response transfers between
dither phases. Frozen correction is then applied to already-consumed targets.
This is neither independent validation nor a claimed physical detector model.
"""
from pathlib import Path
import json, sys
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares
from scipy.linalg import null_space
from scipy.stats import chi2

HERE=Path(__file__).resolve().parent
REPO=HERE.parent/'repo'
ART=REPO/'artifacts/strong_lensing_geometry_20261021'

def read_primary_fits(path):
    raw=Path(path).read_bytes();hdr={};n=0
    while True:
        card=raw[n:n+80].decode('ascii');n+=80
        if card[:8].strip()=='END':break
        if card[8:10]=='= ':
            s=card[10:].split('/')[0].strip()
            try:hdr[card[:8].strip()]=float(s) if any(c in s for c in '.Ee') else int(s)
            except ValueError:hdr[card[:8].strip()]=s.strip("'")
    offset=((n+2879)//2880)*2880
    typ={-64:'>f8',-32:'>f4',32:'>i4',16:'>i2',8:'u1'}[hdr['BITPIX']]
    arr=np.frombuffer(raw,dtype=typ,count=hdr['NAXIS1']*hdr['NAXIS2'],offset=offset).reshape(hdr['NAXIS2'],hdr['NAXIS1']).astype(float)
    return arr,hdr

def spline(grid):
    return RectBivariateSpline(np.arange(grid.shape[1]),np.arange(grid.shape[0]),grid.T,kx=3,ky=3,s=0)

def eval_grid(interp,x,y,center):
    return interp.ev((x-center[0])*4+50,(y-center[1])*4+50)

def fit_arrays(rec,grid,degree=1):
    x,y,d,e=[np.asarray(rec[k],float) for k in ['x','y','data','error']]
    center=np.asarray(rec['center'],float);scale=float(np.median(e))
    bgcols=[np.ones(len(x)),x-center[0],y-center[1]]
    if degree==2:bgcols += [(x-center[0])**2,(x-center[0])*(y-center[1]),(y-center[1])**2]
    bg=np.column_stack(bgcols);interp=spline(grid)
    def predict(p):return np.exp(p[2])*eval_grid(interp,x,y,p[:2])+bg@p[3:]
    p0=np.r_[center,np.log(rec['flux']/scale),np.asarray(rec['background'],float)/scale]
    lo=np.r_[center-.7,-10.,np.full(bg.shape[1],-np.inf)]
    hi=np.r_[center+.7,30.,np.full(bg.shape[1],np.inf)]
    fit=least_squares(lambda p:(predict(p)-d/scale)/(e/scale),p0,bounds=(lo,hi),jac='3-point',x_scale='jac',max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-8)
    model=predict(fit.x)*scale
    jc,jn=fit.jac[:,:2],fit.jac[:,2:]
    jp=jc-jn@np.linalg.lstsq(jn,jc,rcond=1e-12)[0]
    cov=np.linalg.inv(jp.T@jp)
    return {'chi2':float(fit.fun@fit.fun),'dof':len(x)-len(fit.x),'center':fit.x[:2].tolist(),
            'centroid_covariance_pixel2':cov.tolist(),'flux':float(np.exp(fit.x[2])*scale),
            'background':(fit.x[3:]*scale).tolist(),'data_model':model,'success':bool(fit.success),'active':fit.active_mask.tolist()}

def fields():
    j=json.loads((ART/'standard_psf_convention_diagnostic.json').read_text())
    out=[]
    for row in j['records']:
        a=row['data_fit_arrays'];origin=np.asarray(row['field_stamp_bounds_xyxy'][:2]);cent=np.asarray(row['reproduced_global_centroid'])-origin
        g,_=read_primary_fits(ART/'standard_psf'/row['exposure'].replace('_cal.fits','_FIELD_standard_epsf.fits'))
        r={'name':row['exposure'],'x':a['pixel_x'],'y':a['pixel_y'],'data':a['pixel_data'],'error':a['pixel_error'],
           'center':cent.tolist(),'flux':row['flux_SCI_pixel_units'],'background':row['background_plane_SCI'],'grid':g}
        r['baseline']=fit_arrays(r,g);out.append(r)
    return out

class ShapeBasis:
    def __init__(self,spacing=1.,width=.65):
        vals=np.arange(-9.,9.00001,spacing);xx,yy=np.meshgrid(vals,vals)
        self.nodes=np.column_stack([xx.ravel(),yy.ravel()]);self.width=width
        radius=np.linalg.norm(self.nodes,axis=1)
        self.nodes=self.nodes[(radius>=2.)&(radius<=10.5)]
        g=np.arange(101)/4-12.5;gx,gy=np.meshgrid(g,g)
        self.gx,self.gy=gx,gy
        b=self.raw(gx.ravel(),gy.ravel())
        # Exact gauge on the 101x101 stored response: zero added flux and first moments.
        constraints=np.vstack([b.sum(axis=0),gx.ravel()@b,gy.ravel()@b])
        self.null=null_space(constraints)
        self.grid_basis=b@self.null
    def raw(self,x,y):
        r2=(np.asarray(x)[:,None]-self.nodes[None,:,0])**2+(np.asarray(y)[:,None]-self.nodes[None,:,1])**2
        return np.exp(-r2/(2*self.width**2))
    def matrix(self,x,y):return self.raw(x,y)@self.null
    def grid(self,coeff):return (self.grid_basis@coeff).reshape(101,101)

def learn(records,basis,reg,iterations=4):
    coeff=np.zeros(basis.null.shape[1]);results=[r['baseline'] for r in records]
    for iteration in range(iterations):
        lhs=[];rhs=[]
        for rec,fit in zip(records,results):
            x,y,d,e=[np.asarray(rec[k],float) for k in ['x','y','data','error']]
            c=np.array(fit['center']);f=fit['flux']
            xx=x-rec['center'][0];yy=y-rec['center'][1]
            back=np.column_stack([np.ones(len(x)),xx,yy])@np.array(fit['background'])
            rawpsf=eval_grid(spline(rec['grid']),x,y,c)
            lhs.append(f*basis.matrix(x-c[0],y-c[1])/e[:,None]);rhs.append((d-back-f*rawpsf)/e)
        a=np.concatenate(lhs);b=np.concatenate(rhs)
        coeff=np.linalg.solve(a.T@a+reg*np.eye(a.shape[1]),a.T@b)
        delta=basis.grid(coeff)
        results=[fit_arrays(r,r['grid']+delta) for r in records]
    return coeff,results

def clean_fit(f):return {k:v for k,v in f.items() if k!='data_model'}

def main():
    recs=fields();report={'scope':'Consumed-data development correction; field only training and selection; no lens target score in training.',
      'gauge':'Additive correction has exact zero grid integral and first moments; original PSF normalization/origin retained.',
      'baseline':[clean_fit(r['baseline']) for r in recs],'configurations':[],'not_independent_validation':True}
    # This finite comparison chooses a development correction on field data only.
    for spacing,width in [(1.5,.9),(1.,.65),(.75,.5)]:
        basis=ShapeBasis(spacing,width)
        for reg in [1e7,1e8,1e9]:
            row={'spacing':spacing,'width':width,'ridge':reg,'shape_parameters':basis.null.shape[1],'folds':[]}
            for hold in range(4):
                coef,train=learn([r for k,r in enumerate(recs) if k!=hold],basis,reg)
                corr=basis.grid(coef);test=fit_arrays(recs[hold],recs[hold]['grid']+corr)
                row['folds'].append({'held_exposure':hold+1,'train_chi2':sum(f['chi2'] for f in train),'test':clean_fit(test)})
            row['field_cv_chi2']=sum(f['test']['chi2'] for f in row['folds'])
            report['configurations'].append(row)
            print(json.dumps({k:v for k,v in row.items() if k!='folds'}),flush=True)
            (HERE/'field_response_results.json').write_text(json.dumps(report,indent=2))
    best=min(report['configurations'],key=lambda r:r['field_cv_chi2'])
    basis=ShapeBasis(best['spacing'],best['width']);coef,train=learn(recs,basis,best['ridge'],iterations=6);corr=basis.grid(coef)
    report['selected_field_configuration']={k:v for k,v in best.items() if k!='folds'}
    report['all_field_train_fits']=[clean_fit(r) for r in train]
    report['correction_checks']={'integral':float(corr.sum()/16),'first_x':float((corr*basis.gx).sum()/16),
                                'first_y':float((corr*basis.gy).sum()/16),'max_abs':float(np.max(abs(corr))),
                                'psf_negative_pixels':[int(np.sum(r['grid']+corr<0)) for r in recs]}
    np.save(HERE/'field_response_correction.npy',corr)
    report['target_transfer']=[]
    for path in sorted((ART/'jwst_target_diagnostic/exposures').glob('*.json')):
        ex=json.loads(path.read_text());entries=[]
        for p in ex['points']:
            grid,_=read_primary_fits(REPO/p['field_psf_path'])
            rec={'x':p['pixel_x'],'y':p['pixel_y'],'data':p['pixel_data'],'error':p['pixel_error'],
                 'center':p['centroid_xy_pixel'],'flux':p['flux_in_SCI_pixel_units'],'background':p['background_polynomial_SCI']}
            # Background basis origin differs in original fit; rerun baseline consistently first.
            base=fit_arrays(rec,grid,degree=2);rec.update({k:base[k] for k in ['center','flux','background']})
            fit=fit_arrays(rec,grid+corr,degree=2)
            entries.append({'label':p['label'],'original_chi2':p['chi2'],'baseline_refit':clean_fit(base),'corrected':clean_fit(fit),
                            'centroid_shift_pixel':(np.array(fit['center'])-p['centroid_xy_pixel']).tolist(),
                            'original_tangent':p['tangent_east_north_arcsec'],'GWCS_jacobian':p['GWCS_jacobian_arcsec_per_pixel']})
        report['target_transfer'].append({'exposure':path.name,'points':entries})
    (HERE/'field_response_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'selected':report['selected_field_configuration'],'gauge':report['correction_checks'],
                      'target_chi2':[[p['corrected']['chi2']/p['corrected']['dof'] for p in e['points']] for e in report['target_transfer']]}),flush=True)

if __name__=='__main__':main()
