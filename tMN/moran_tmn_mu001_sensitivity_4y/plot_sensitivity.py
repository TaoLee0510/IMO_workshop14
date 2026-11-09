import json,csv,math,os
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_joint_mpl')
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
ROOT=Path(__file__).resolve().parent
CONFIG=json.loads((ROOT/'config.json').read_text());COHORT=json.loads((ROOT/'data/cohort.json').read_text())
MUS=CONFIG['sensitivity']['mu'];DS=CONFIG['sensitivity']['D_init']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
def load():
    return [json.loads(p.read_text()) for p in (ROOT/'results/fits').glob('*.json') if not p.stem.endswith('_evaluations')]
def save(fig,name):
    path=ROOT/'figures'/name;path.parent.mkdir(exist_ok=True,parents=True)
    fig.savefig(path.with_suffix('.pdf'),bbox_inches='tight');fig.savefig(path.with_suffix('.png'),dpi=180,bbox_inches='tight');plt.close(fig)
def matrix(ax,values,title,log=False,fmt='.2f',bounds=None):
    masked=np.ma.masked_invalid(values);kwargs={'cmap':'viridis'}
    if log:
        positive=values[np.isfinite(values)&(values>0)]
        kwargs['norm']=LogNorm(vmin=float(positive.min()),vmax=max(float(positive.max()),float(positive.min())*1.01))
    elif bounds:kwargs.update(vmin=bounds[0],vmax=bounds[1])
    cmap=plt.colormaps['viridis'].copy();cmap.set_bad('#dedede');kwargs['cmap']=cmap
    im=ax.imshow(masked,aspect='auto',**kwargs)
    ax.set_xticks(np.arange(len(MUS)),['0' if x==0 else f'{x:g}' for x in MUS],rotation=45,ha='right')
    ax.set_yticks(np.arange(len(DS)),[str(x) for x in DS]);ax.set(xlabel='C-to-D probability per C division (mu)',ylabel='Measured-baseline D cells',title=title)
    for i in range(len(DS)):
        for j in range(len(MUS)):
            if not np.isfinite(values[i,j]):ax.text(j,i,'n.a.',ha='center',va='center',fontsize=6,color='#444444');continue
            tone=float(im.norm(values[i,j]));color='white' if tone<.55 else '#151515'
            ax.text(j,i,format(values[i,j],fmt),ha='center',va='center',fontsize=10.5,color=color)
    return im
def main():
    fits=load();good=[r for r in fits if r['status']=='fitted'];assert len(fits)==79 and len(good)==78
    grid={(r['parameters']['mu'],r['parameters']['D_init']):r for r in good if r['parameters']['b']==CONFIG['reference']['b'] and r['parameters']['s_D']==CONFIG['reference']['s_D']}
    best=min(r['joint_nll'] for r in grid.values());nll=np.full((8,9),np.nan);alpha=nll.copy();theta=nll.copy()
    for i,d in enumerate(DS):
        for j,mu in enumerate(MUS):
            r=grid.get((mu,d))
            if r:nll[i,j]=r['joint_nll']-best;alpha[i,j]=r['alpha'];theta[i,j]=r['theta_days']
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained');im=matrix(ax,nll,'Joint likelihood after re-fitting alpha and theta',fmt='.1f');fig.colorbar(im,ax=ax,label='Delta NLL from best grid scenario');save(fig,'mu_D_joint_likelihood')
    fig,axes=plt.subplots(1,2,figsize=(12,6.5),layout='constrained')
    for ax,vals,title,log in [(axes[0],alpha,'Fitted alpha',False),(axes[1],theta,'Fitted theta (days)',True)]:
        im=matrix(ax,vals,title,log=log,fmt='.0f' if log else '.1f');fig.colorbar(im,ax=ax,shrink=.8)
    save(fig,'mu_D_fitted_parameters')
    risks=list(csv.DictReader((ROOT/'results/patient_risks.csv').open()));bykey={(r['scenario'],r['patient_id']):r for r in risks}
    for patient in COHORT:
        fig,axes=plt.subplots(2,1,figsize=(10,11),layout='constrained')
        for ax,year in zip(axes,[1,4]):
            vals=np.full((8,9),np.nan)
            for i,d in enumerate(DS):
                for j,mu in enumerate(MUS):
                    fit=grid.get((mu,d))
                    if fit:vals[i,j]=float(bykey[(fit['scenario'],patient['patient_id'])][f'year{year}_probability'])
            im=matrix(ax,vals,f'Year {year} tMN probability',log=True,fmt='.1e');fig.colorbar(im,ax=ax,label='Probability (log color)',shrink=.8)
        fig.suptitle(patient['patient_id'].replace('_',' ')+'   '+patient['gene']+'   Own baseline state',fontsize=12)
        save(fig,'patients/'+patient['patient_id']+'_mu_D_risk')
    reference=next(r for r in fits if r['scenario']=='reference')
    for key,label in [('b','Turnover b (per day)'),('s_D','Extra D fitness s_D')]:
        rows=sorted([r for r in good if (r['scenario']=='reference' or r['scenario'].startswith(key+'_'))],key=lambda r:r['parameters'][key])
        x=[r['parameters'][key] for r in rows];fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
        for column,color,title in [('joint_nll','#173f5e','Joint'),('VAF_nll','#1765b0','VAF'),('clinical_given_VAF_nll','#c82d36','Clinical given VAF')]:
            axes[0].plot(x,[r[column]-reference[column] for r in rows],marker='o',ms=4,label=title,color=color)
        axes[0].axhline(0,color='gray',lw=.8);axes[0].set(xlabel=label,ylabel='NLL change from reference');axes[0].legend(fontsize=8)
        for pid in ['Patient_103','Patient_197','Patient_235']:
            axes[1].plot(x,[float(bykey[(r['scenario'],pid)]['year1_probability']) for r in rows],marker='o',ms=4,label=pid.replace('_',' '))
        axes[1].set(yscale='log',xlabel=label,ylabel='One-year tMN probability');axes[1].legend(fontsize=8)
        if key=='b':
            for ax in axes:ax.set_xscale('log')
        save(fig,key+'_sensitivity')
    print('Sensitivity figures complete.')
if __name__=='__main__':main()
