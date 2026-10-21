"""No-fit native raw residual partition by radius, jumps and group support."""
from pathlib import Path
import argparse,hashlib,re,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_he0435_native_completion as run


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--label',choices=['all_data','training','held_target'],default='held_target');args=parser.parse_args()
    dest=run.OUT/args.label;path=dest/'native_replay_maps.npz';maps=np.load(path);rows=[]
    histories={}
    for key in maps.files:
        if not key.endswith('_raw_residual'):continue
        name=key.removesuffix('_raw_residual');stem=re.match(r'(jw\d+_\d+_\d+_nrcb3)',name).group(1)
        if stem not in histories:
            with np.load(run.BASE/'histories'/f'{stem}_histories.npz') as f:histories[stem]={k:f[k] for k in f.files}
        h=histories[stem];label=name[-1] if name[-1] in 'ABCD' else 'FIELD';p=list(h['labels']).index(label)
        x=maps[name+'_x'].astype(int);y=maps[name+'_y'].astype(int);bx,by,_,_=h['bounds_xyxy'][p];ix=x-bx;iy=y-by
        z=maps[name+'_standardized_residual'];raw=maps[name+'_raw_residual'];dq=h['groupdq'][p,:,iy,ix]
        # Advanced indexing returns selected-pixel,group ordering here.
        if dq.shape!=(len(z),9):dq=dq.T
        jump=np.any((dq&4)!=0,axis=1);full=np.all((dq&3)==0,axis=1)
        radius=np.hypot(x-h['prior_guess_native_xy'][p,0],y-h['prior_guess_native_xy'][p,1])
        def summary(mask):
            return {'pixels':int(mask.sum()),'Q_original_native_ERR':float(z[mask]@z[mask]),'raw_RMS_SCI':float(np.sqrt(np.mean(raw[mask]**2))) if mask.any() else None}
        grid=np.full((17,17),np.nan);grid[y-y.min(),x-x.min()]=z
        corr={}
        for axis,a,b in [('x',grid[:,:-1],grid[:,1:]),('y',grid[:-1],grid[1:])]:
            use=np.isfinite(a)&np.isfinite(b);corr[axis]=float(np.corrcoef(a[use],b[use])[0,1])
        rows.append({'name':name,'label':label,'all':summary(np.ones(len(z),bool)),
            'native_jump_history':summary(jump),'no_native_jump_history':summary(~jump),
            'all9groups_DNU_SAT_eligible':summary(full),'fewer_eligible_groups':summary(~full),
            'radial':[{'radius_pixel':[lo,hi],**summary((radius>=lo)&(radius<hi))} for lo,hi in [(0,1.5),(1.5,3),(3,5),(5,13)]],
            'adjacent_standardized_residual_correlation':corr})
    run.prior.dump(dest/'spatial_residual_diagnosis.json',{'scope':'No-fit partition of unchanged official-native saved residuals. Jump flags identify potentially unmodeled deposited charge/segment-state response; they are not discarded to improve fit.',
        'rows':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'saved_mean_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'physical_validation_pass':False,'physical_admission':False,'physical_guard':2e-10})
    print([(r['name'],r['all']['Q_original_native_ERR'],r['native_jump_history']['Q_original_native_ERR']) for r in rows],flush=True)


if __name__=='__main__':main()
