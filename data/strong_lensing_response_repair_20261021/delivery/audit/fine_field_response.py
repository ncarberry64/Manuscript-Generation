"""One prespecified finer-resolution positive field-response experiment.

RBF spacing0.5pixel,width0.35pixel replaces spacing1,width0.65. Same priorSD0.2,
positivity,pedestal convention,native errors and4fold field checks. All target
pixels and physical predictions are absent. No resolution sweep is performed.
"""
from pathlib import Path
import json,sys
import numpy as np
from scipy.linalg import svdvals
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
from positive_field_response import get_fields,fit,serial
from spatial_field_response import residual_summary

class FineBasis:
    iterative_solver=True
    def __init__(self,reference_grid):
        self.width=.35;z=np.arange(-9.,9.0001,.5);xx,yy=np.meshgrid(z,z);nodes=np.c_[xx.ravel(),yy.ravel()]
        self.nodes=nodes[np.linalg.norm(nodes,axis=1)<=10.5];self.null=np.zeros((1,len(self.nodes)))
        z=np.arange(101)/4-12.5;self.gx,self.gy=np.meshgrid(z,z)
        b,_,_=self.raw(self.gx.ravel(),self.gy.ravel());p=reference_grid.ravel()
        c=np.vstack([p@b,(p*self.gx.ravel())@b,(p*self.gy.ravel())@b])
        self.gauge=np.linalg.qr(c.T,mode='reduced')[0]
        self.grid_basis=b-(b@self.gauge)@self.gauge.T
    def raw(self,x,y):
        dx=np.asarray(x)[:,None]-self.nodes[None,:,0];dy=np.asarray(y)[:,None]-self.nodes[None,:,1]
        b=np.exp(-(dx*dx+dy*dy)/(2*self.width**2))
        return b,-dx/self.width**2*b,-dy/self.width**2*b
    def evaluate(self,x,y,rec=None):
        return tuple(b-(b@self.gauge)@self.gauge.T for b in self.raw(x,y))
    def response_grid(self,rec,a):return rec['grid']*np.exp((self.grid_basis@a).reshape(101,101))

def data_information(recs,basis,result):
    a=result['shape'];nshape=len(a);rows=[]
    for i,r in enumerate(recs):
        cx,cy,lf,bg,bx,by=result['parameters'][nshape+6*i:nshape+6*(i+1)]
        x,y=r['x']-cx,r['y']-cy;b,dbx,dby=basis.evaluate(x,y,r)
        corr=np.exp(b@a);p=r['interpolator'].ev(x,y);px=r['interpolator'].ev(x,y,dx=1);py=r['interpolator'].ev(x,y,dy=1)
        f=np.exp(lf);point=f*p*corr
        jshape=point[:,None]*b/r['error'][:,None]
        jn=np.column_stack([-f*corr*(px+p*(dbx@a)),-f*corr*(py+p*(dby@a)),point,np.ones(len(x)),r['x'],r['y']])/r['error'][:,None]
        prof=jshape-jn@np.linalg.lstsq(jn,jshape,rcond=1e-12)[0]
        rows.append(prof)
    s=svdvals(np.concatenate(rows));rank=int(np.sum(s>s[0]*1e-10))
    return {'native_data_rows':sum(len(r['x']) for r in recs),'shape_parameters_before3_gauge':len(a),
      'profiled_shape_numerical_rank_relative1e-10':rank,'unidentified_directions_including3_gauge':len(a)-rank,
      'linearized_effective_shape_df_with_prior':float(np.sum(s*s/(s*s+25))),
      'singular_values_max_min':[float(s[0]),float(s[-1])],
      'scope':'Rank/profile is local linearization. Proper ridge identifies coefficients by prior, not extra calibration data.'}

def main():
    recs=get_fields();basis=FineBasis(recs[0]['grid']);report={'status':'EXPLORATORY_FIELD_ONLY_FINER_RESPONSE',
      'spacing_pixel':.5,'width_pixel':.35,'shape_prior_sd':.2,'coefficient_count':basis.null.shape[1],
      'native_errors_unchanged':True,'target_pixels_used':False,'lens_predictions_used':False,'folds':[]}
    for hold in range(4):
        trainrecs=[r for i,r in enumerate(recs) if i!=hold]
        tr=fit(trainrecs,basis,max_nfev=70);te=fit([recs[hold]],basis,a_fixed=tr['shape'],max_nfev=70)
        info=data_information(trainrecs,basis,tr)
        row={'held':hold+1,'train':serial(tr),'test':serial(te),'training_identifiability':info};report['folds'].append(row)
        (HERE/'fine_field_results.json').write_text(json.dumps(report,indent=2))
        print(json.dumps({'held':hold+1,'train_success':tr['success'],'nfev':tr['nfev'],'held_chi2':te['fits'][0]['chi2'],
          'held_chi2_per_dof':te['fits'][0]['chi2']/283,'rank':info['profiled_shape_numerical_rank_relative1e-10'],
          'effective_df':info['linearized_effective_shape_df_with_prior']}),flush=True)
    allfit=fit(recs,basis,max_nfev=100);report['all_field_fit']=serial(allfit);report['residuals']=residual_summary(recs,basis,allfit)
    report['all_training_identifiability']=data_information(recs,basis,allfit)
    report['held_exposure_chi2_total']=sum(r['test']['fits'][0]['chi2'] for r in report['folds'])
    np.savez(HERE/'fine_field_shape.npz',shape=allfit['shape'],parameters=allfit['parameters'],nodes=basis.nodes,gauge=basis.gauge,width=basis.width)
    (HERE/'fine_field_results.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'held_total':report['held_exposure_chi2_total'],'all_fit':serial(allfit),'rank':report['all_training_identifiability']}),flush=True)

if __name__=='__main__':main()
