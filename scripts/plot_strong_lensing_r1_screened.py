"""Generate the scientific figure from the frozen local prediction."""
from pathlib import Path
import json
import sys
import os
import tempfile
os.environ.setdefault('MPLCONFIGDIR',str(Path(tempfile.gettempdir())/'bhsm_lensing_r1_mplconfig_20261020'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from strong_lensing_r1_screened import R1Coefficients,R1ScreenedLens
from strong_lensing_local_family import TolmanVIILens

def main(output_dir=None, formats=('pdf', 'png', 'svg')):
    result=json.loads((ROOT/'artifacts/strong_lensing_r1_screened_20261020/results.json').read_text())
    base=TolmanVIILens(mass_Msun=result['body']['mass_Msun'],radius_Mpc=result['body']['radius_Mpc'])
    lens=R1ScreenedLens(base,R1Coefficients(**result['coefficients']))
    radius=np.geomspace(.1,1000.,180)
    force=np.array([lens.extra_force(r/1000) for r in radius])*1e6
    projected=np.array([lens.extra_alpha(r/1000)/base.application([r/1000,0])[1][0] for r in radius])*1e6
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                        'axes.spines.right':False,'axes.labelcolor':'#243447','text.color':'#243447',
                        'axes.edgecolor':'#8192a3','xtick.color':'#425467','ytick.color':'#425467',
                        'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(8.4,3.35),gridspec_kw={'width_ratios':[1.3,1]})
    ax=axes[0]
    ax.loglog(radius,force,color='#205a9c',lw=2,label='Local force')
    ax.loglog(radius,projected,color='#008a83',lw=2,label='Projected deflection')
    ax.axvline(3,color='#9aa4b0',ls=':',lw=1)
    ax.axvline(lens.rV_Mpc*1000,color='#9aa4b0',ls='--',lw=1)
    ax.text(3.5,25,'Body radius',fontsize=7.5,rotation=90,va='bottom')
    ax.text(lens.rV_Mpc*1000*1.13,2.05,r'$r_V \simeq 212$ kpc',fontsize=7.5,rotation=90,va='bottom')
    ax.set(xlabel='Radius or impact parameter (kpc)',ylabel='Additional response (ppm)',ylim=(2,2400),xlim=(.1,1000))
    ax.set_title('(a) Radial response',loc='left',fontsize=10,fontweight='bold')
    ax.grid(axis='y',which='major',color='#e7ecf0',lw=.65)
    ax.legend(loc='upper left',frameon=False,fontsize=8)
    ax=axes[1]
    values=result['changes']['relative_delays_seconds']
    bars=ax.bar([0,1],values,width=.52,color=['#205a9c','#008a83'])
    ax.set_xticks([0,1],['Central minus\npositive outer','Negative outer minus\npositive outer'])
    ax.tick_params(axis='x',length=0,labelsize=8)
    ax.set(ylabel='Additional relative delay (s)',ylim=(0,31))
    ax.set_title('(b) Delay changes',loc='left',fontsize=10,fontweight='bold')
    ax.grid(axis='y',color='#e7ecf0',lw=.65)
    ax.set_axisbelow(True)
    for b,v in zip(bars,values):
        ax.text(b.get_x()+b.get_width()/2,v+.8,f'{v:.2f} s',ha='center',fontsize=10,fontweight='bold')
    fig.subplots_adjust(left=.09,right=.985,bottom=.25,top=.89,wspace=.38)
    fig.text(.09,.035,r'Frozen R1 cubic model: $M=10^{11}M_\odot$, $R=3$ kpc, $z_l=0.3$, $z_s=1$; matched leading-order comparison.',fontsize=7.7)
    out = output_dir or ROOT/'manuscript/figures'
    out.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        path=out/f'r1_screening_prediction.{ext}'
        fig.savefig(path,dpi=240,facecolor='white',metadata={'Creator':'Reproducible matplotlib figure'} if ext=='pdf' else None)
        print(path)
    plt.close(fig)
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, help='Write replay figures separately from the supplied publication assets.')
    parser.add_argument('--formats', nargs='+', choices=('pdf', 'png', 'svg'), default=('pdf', 'png', 'svg'), help='Select replay formats without changing the figure data.')
    args = parser.parse_args()
    main(args.output_dir, args.formats)
