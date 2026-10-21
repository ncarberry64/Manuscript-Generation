"""Replay saved new joint-scene means and retain native target residual maps.

No optimizer, field nuisance refit, new response or unused observation is used.
Existing result files are read-only and are checked byte-for-byte before/after.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_joint_scene as run
js=run.js
OUT=run.OUT/'raw_residuals'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def export(label,problem,p,eta,full,expected):
    run.set_eta(problem.records,eta)
    residual,j,_,rows=problem.evaluate(p)
    arrays={k:[] for k in ['name','bounds','x','y','SCI','ERR','DQ','mask','mean','raw_residual','standardized_residual']}
    total=0.;stats=[];start=0
    for k,rec in enumerate(problem.records):
        x0,y0,x1,y1=rec['bounds'];yy,xx=np.mgrid[y0:y1,x0:x1]
        sci,err,dq,mask=full[rec['name']]
        ix=rec['x'].astype(int)-x0;iy=rec['y'].astype(int)-y0
        recreated=np.zeros(mask.shape,bool);recreated[iy,ix]=True
        if not np.array_equal(mask,recreated):raise AssertionError('Native mask differs from original saved-mean pixel set')
        if not np.array_equal(sci[mask],rec['data']) or not np.array_equal(err[mask],rec['error']):raise AssertionError('Native SCI/ERR row order mismatch')
        cx,cy,lf=p[3*k:3*k+3];sky=problem.sky(p)[9*rec['patch']:9*rec['patch']+9]
        bg=problem.bg_start+3*problem.eindex[rec['exposure']]
        value=np.exp(lf)*rec['response'].values(rec['x']-cx,rec['y']-cy)+rec['scene_matrix']@sky+rec['detector_background']@p[bg:bg+3]
        raw=value-rec['data'];standardized=raw/rec['error'];n=len(raw)
        if not np.array_equal(standardized,residual[start:start+n]):raise AssertionError('Original residual application disagrees')
        start+=n;mean=np.full(mask.shape,np.nan);rr=np.full(mask.shape,np.nan);zz=np.full(mask.shape,np.nan)
        mean[mask]=value;rr[mask]=raw;zz[mask]=standardized
        for key,val in [('name',rec['name']),('bounds',rec['bounds']),('x',xx),('y',yy),('SCI',sci),('ERR',err),('DQ',dq),('mask',mask),('mean',mean),('raw_residual',rr),('standardized_residual',zz)]:arrays[key].append(val)
        q=float(standardized@standardized);total+=q
        stats.append({'name':rec['name'],'pixels':n,'Q':q,'standardized_residual_max_abs':float(abs(standardized).max()),'raw_residual_mean_SCI':float(raw.mean()),'raw_residual_RMS_SCI':float(np.sqrt(np.mean(raw**2)))})
    if not np.isclose(total,expected,rtol=3e-14,atol=1e-9):raise AssertionError(f'{label} Q replay mismatch: {total} vs {expected}')
    path=OUT/f'{label}.npz';np.savez_compressed(path,**{k:np.asarray(v) for k,v in arrays.items()})
    return {'file':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':sha(path),'Q_replayed':total,'Q_original':expected,'difference':total-expected,'records':stats}


def plot_final(receipt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    saved=np.load(OUT/'all_data_final.npz');fig,axes=plt.subplots(4,4,figsize=(10,10),constrained_layout=True)
    for i,ax in enumerate(axes.ravel()):
        data=np.ma.masked_where(~saved['mask'][i],saved['standardized_residual'][i])
        im=ax.imshow(data,origin='lower',cmap='RdBu_r',vmin=-10,vmax=10)
        q=receipt['records'][i]['Q'];mx=receipt['records'][i]['standardized_residual_max_abs']
        name=str(saved['name'][i]);ax.set_title(f'E{i//4+1} {name[-1]}: Q={q:.0f}, max|r|={mx:.1f}',fontsize=9)
        ax.set_xticks([]);ax.set_yticks([])
    fig.colorbar(im,ax=axes,label='(mean - SCI) / original ERR; display clipped at +/-10',shrink=.8)
    fig.suptitle('Saved one-step joint-scene mean: consumed F115W development\nAll residuals fail native-pixel adequacy; no error rescaling or validation',fontsize=12)
    path=OUT/'all_data_final_standardized_residuals.png';fig.savefig(path,dpi=170);plt.close(fig)
    return {'file':str(path.relative_to(ROOT)).replace('\\','/'),'sha256':sha(path),'display_only_clip':[-10,10],'raw_arrays_unclipped':True}


def main():
    OUT.mkdir(exist_ok=True)
    files=[run.OUT/f'{name}{suffix}' for name in ['all_data','held_exposure_training'] for suffix in ['.json','_linearized_arrays.npz']]
    files.extend(run.OUT/p for p in ['held_exposure_application.json','execution_summary.json'])
    before={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files}
    records,fields,eta0,wcss,origin=run.load_inputs()
    # Recover full17x17 stamps, including rejected pixels for explicit mask
    # retention. Model/residual values at excluded pixels remain declared NaN.
    full={}
    for exposure in sorted({r['exposure'] for r in records}):
        path=ROOT/'data/strong_lensing_geometry_20261021/native'/exposure
        with fits.open(path,memmap=False) as h:
            for rec in records:
                if rec['exposure']!=exposure:continue
                x0,y0,x1,y1=rec['bounds'];s,e,d=[np.array(h[k].data[y0:y1,x0:x1]) for k in ['SCI','ERR','DQ']]
                mask,_=run.jwst_native_rate_mask(s,e,d,primary_header=h[0].header);full[rec['name']]=(s,e,d,mask)
    applications={};allproblem=js.SceneProblem(records)
    for label,problem in [('all_data',allproblem),('held_exposure_training',js.SceneProblem([r for r in records if r['exposure']!=allproblem.exposures[-1]]))]:
        saved=np.load(run.OUT/f'{label}_linearized_arrays.npz');report=json.loads((run.OUT/f'{label}.json').read_text())
        for suffix,key in [('baseline','initial'),('final','updated')]:
            expected=report['baseline_known_ERR_fit']['chi2'] if key=='initial' else report['coupled_GN']['original_nonlinear_known_ERR_fit']['chi2']
            applications[f'{label}_{suffix}']=export(f'{label}_{suffix}',problem,saved[f'parameters_{key}'],saved[f'eta_{key}'],full,expected)
        if label=='held_exposure_training':training=problem;trainp=saved['parameters_updated'].copy();traineta=saved['eta_updated'].copy()
    held=[r for r in records if r['exposure']==allproblem.exposures[-1]];problem=js.SceneProblem(held,scene_fixed=training.sky(trainp))
    report=json.loads((run.OUT/'held_exposure_application.json').read_text());p=problem.initial()
    rowmap={r['name']:r for r in report['records']}
    for k,rec in enumerate(problem.records):
        row=rowmap[rec['name']];p[3*k:3*k+3]=[*row['centroid_xy_pixel'],np.log(row['flux'])]
    p[problem.bg_start:]=report['records'][0]['background']
    applications['held_exposure']=export('held_exposure',problem,p,traineta,full,report['known_ERR_fit']['chi2'])
    plot=plot_final(applications['all_data_final'])
    after={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files}
    if before!=after:raise AssertionError('Original saved result bytes changed during read-only mean replay')
    receipt={'scope':'Read-only saved-mean replay on already-consumed native target stamps. No optimizer, field nuisance refit, changed mean, response update, ERR scaling, physical lens fit or unused data application.',
             'array_conventions':'Arrays retain full17x17 native x/y, SCI, ERR, DQ and mask. Model/raw/standardized residuals are defined on exactly the original valid rows; excluded pixels are NaN in these three arrays. Raw SCI/ERR/DQ at excluded pixels are retained unchanged.',
             'applications':applications,'plot':plot,'original_results_sha256_before':before,'original_results_sha256_after':after,'original_results_byte_identical':before==after,
             'source_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [Path(__file__).resolve(),ROOT/'code/strong_lensing_joint_scene.py',ROOT/'scripts/run_he0435_joint_scene.py']},
             'input_provenance':'../execution_summary.json records exact native receipts, response NPZ,20localPSFs and4archivedGWCS hashes. No new observation is opened.',
             'field_map_scope':'Only target means are exported. Field nuisance flux/background values were not retained in the original fit receipt; reconstructing those would require an additional nuisance fit, which this no-refit exporter does not perform.',
             'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10}
    run.dump(OUT/'execution_receipt.json',receipt)
    for name,a in applications.items():print(name,a['Q_replayed'],a['difference'])


if __name__=='__main__':main()
