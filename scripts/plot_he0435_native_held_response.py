"""Show original-ERR held residuals before/after the converged native response."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artifacts/strong_lensing_response_completion_20261021'


def main():
    a=OUT/'static_diagnostics/held_target_strict.npz'
    b=OUT/'native_coupled/held_target/native_replay_maps.npz'
    old=np.load(a);new=np.load(b)
    fig,axes=plt.subplots(2,4,figsize=(14.8,7.5),layout='constrained')
    checks=[]
    for p,label in enumerate('ABCD'):
        prefix=f'jw01198007001_02101_00004_nrcb3_cal_{label}'
        grid=np.full((17,17),np.nan);x0,y0,_,_=old['bounds'][p]
        x=new[prefix+'_x'].astype(int)-x0;y=new[prefix+'_y'].astype(int)-y0
        grid[y,x]=new[prefix+'_standardized_residual']
        previous=np.where(old['mask'][p],old['standardized_residual'][p],np.nan)
        for row,(value,name) in enumerate([(previous,'Static'),(grid,'Native coupled')]):
            q=float(np.nansum(value**2))
            im=axes[row,p].imshow(np.ma.masked_invalid(value),origin='lower',cmap='RdBu_r',vmin=-10,vmax=10,interpolation='nearest')
            axes[row,p].set_title(f'{name} {label}: Q={q:.1f}',fontsize=12)
            axes[row,p].set_xticks([]);axes[row,p].set_yticks([])
        checks.append({'label':label,'static_Q':float(np.nansum(previous**2)),
                       'native_original_ERR_Q':float(np.nansum(grid**2)),
                       'source_array_Q':float(np.sum(new[prefix+'_standardized_residual']**2))})
    fig.suptitle('Held E4 target: the converged native response retains the broad mean failure\nE1–E3 target training; all four field stars calibrate; original ERR unchanged',fontsize=14)
    fig.colorbar(im,ax=axes,label='(mean − SCI) / original native ERR; colors clipped at ±10',shrink=.86)
    path=ROOT/'manuscript/figures/he0435_converged_held_residuals_20261021.png'
    fig.savefig(path,dpi=180);plt.close(fig)
    dest=OUT/'paper_review';dest.mkdir(exist_ok=True)
    (dest/'held_figure_receipt.json').write_text(json.dumps({'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [a,b]},'checks':checks,'unclipped_statistics':True,'error_rescaling':False},indent=2)+'\n',encoding='utf-8',newline='\n')


if __name__=='__main__':main()
