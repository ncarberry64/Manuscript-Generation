"""Render publication figures from frozen tables; no fit or scientific rerun."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'manuscript/figures'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                     'pdf.fonttype': 42, 'savefig.dpi': 220})

def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(OUT / f'{name}.{ext}', bbox_inches='tight')
    plt.close(fig)

def main():
    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.0, 7.0))
    ax.set(xlim=(0, 10), ylim=(0, 10)); ax.axis('off')
    boxes = [
      (0.3, 8.3, 9.4, 1.2, 'CONDITIONAL: closed $S^3$ parent and action-owned sector\nOPEN: background normalization and local matrix elements', '--'),
      (0.3, 6.6, 9.4, 1.1, 'DERIVED on a regular common domain\nConstraints / gauge quotient → physical Hessian → Schur reduction', '-'),
      (0.3, 4.9, 9.4, 1.1, 'CONDITIONAL: select one normalized physical $n=2$ profile\nDERIVED: $X_2=(q_2,\\Pi_2)$; exact dipole $\\sin(2\\psi)\\,\\hat n\\cdot\\hat p$', '--'),
      (0.3, 3.1, 9.4, 1.2, 'DERIVED for fixed response: one state → induced metric\nGrowth / Weyl-lensing / fixed-redshift distance / Jacobi–Sachs optics', '-'),
      (0.3, 1.5, 4.45, 1.0, 'CONDITIONAL prospective gate A\nOne shared two-component state', '--'),
      (5.25, 1.5, 4.45, 1.0, 'CONDITIONAL prospective gate B\nOne coherent global amplitude', '--'),
      (0.3, 0.0, 9.4, 0.9, 'OPEN: state direction, seam export, source statistics and joint covariance\nDERIVED: isolated covariance rank ≤ m for m coherent amplitudes', ':'),
    ]
    for x,y,w,h,label,style in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.07',
                                   fc='white', ec='black', lw=1.2, linestyle=style))
        ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=9.3,linespacing=1.5)
    for a,b in [(8.25,7.8),(6.55,6.1),(4.85,4.4)]:
        ax.annotate('',(5,b),(5,a),arrowprops={'arrowstyle':'->'})
    for x in (2.5,7.5):
        ax.annotate('',(x,2.6),(5,3.0),arrowprops={'arrowstyle':'->'})
    save(fig,'constraint_prediction_chain')
    with (ROOT/'analyses/sn_sightline_tomography_v1/interpretation/LIKELIHOODS.csv').open() as f:
        rows={r['test']:r for r in csv.DictReader(f)}
    keys=['P_original_heldout_release_only','P_original_training_to_heldout',
          'P_training_to_heldout_fixed_external_PV_sensitivity',
          'DES_original_heldout_release_only','DES_Pantheon_transfer']
    labels=['Pantheon+ released distances [PRIMARY]', 'Pantheon+ source/environment [CONDITIONAL]',
            'Pantheon+ fixed external PV [SENSITIVITY ONLY]',
            'DES released distances [PRIMARY]', 'DES nuisance transfer [CONDITIONAL]']
    fig,ax=plt.subplots(figsize=(8.8,4.5))
    for i,k in enumerate(keys):
        v=float(rows[k]['delta_chi2']); y=4-i
        ax.plot([0,v],[y,y],color='0.5',lw=1.1)
        ax.scatter(v,y,s=72,marker='D' if i==2 else 'o',
                   facecolors='white' if i==2 else 'black',edgecolors='black',zorder=3)
        ax.annotate(f'{v:+.3f}',(v,y),xytext=(0,11),textcoords='offset points',ha='center',fontsize=9)
    ax.axvline(0,color='black',lw=.9); ax.set_yticks(range(5),labels[::-1],fontsize=9)
    ax.set(xlim=(-6,4),ylim=(-.6,4.6),xlabel=r'$\Delta\chi^2=\chi^2_{\rm frozen}-\chi^2_{\rm null}$')
    ax.set_title('Negative: favors frozen template     |     Positive: favors null',fontsize=10,pad=16)
    for sp in ('top','right','left'): ax.spines[sp].set_visible(False)
    ax.tick_params(axis='y',length=0); ax.grid(axis='x',alpha=.15)
    fig.text(.5,.005,'External PV: coefficient fixed at 1; complete joint uncertainty unavailable. No error bars inferred.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.055,1,1)); save(fig,'frozen_sn_comparison')

if __name__=='__main__': main()
