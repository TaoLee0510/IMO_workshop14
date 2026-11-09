"""Native simulation batching, exact-trajectory output and transparent plots."""
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path
from dataclasses import asdict
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_joint_mpl')
ROOT=Path(__file__).resolve().parent


def simulate_job(task):
    label,values,alpha,theta,d_fraction,settings=task
    n=values['N'];c0=round(settings['initial_C_fraction']*n);d0=round(d_fraction*n)
    folder=ROOT/'results/trajectories';folder.mkdir(parents=True,exist_ok=True)
    name=f'{label}_d{d0}'
    output=folder/f'{name}.npz';meta=folder/f'{name}.json'
    simulator_sha=hashlib.sha256((ROOT/'simulator.cpp').read_bytes()).hexdigest()
    signature=json.dumps([label,values,alpha,theta,d0,settings['master_seed'],simulator_sha],sort_keys=True)
    job_seed=int.from_bytes(hashlib.sha256(signature.encode()).digest()[:8],'little')
    if output.exists() and meta.exists():
        saved=json.loads(meta.read_text())
        if saved['job_seed']==str(job_seed) and saved['replicates']==settings['replicates']:return saved
    prefix=folder/name
    command=[str(ROOT/'bin/moran_simulator'),str(prefix),str(n),str(values['b']),str(values['s_D']),str(values['mu']),
             str(values['area_days']),str(alpha),str(theta),str(c0),str(d0),str(settings['horizon_days']),
             str(settings['replicates']),str(job_seed),str(values['threshold_fraction'])]
    completed=subprocess.run(command,capture_output=True,text=True,check=True)
    timing=json.loads(completed.stdout)
    raw=Path(str(prefix)+'.bin')
    values_array=np.fromfile(raw,dtype=np.int32).reshape(settings['replicates'],settings['horizon_days']+1,2)
    assert np.all(values_array>=0) and np.all(values_array.sum(axis=2)<=n)
    assert np.all(values_array[:,0,0]==c0) and np.all(values_array[:,0,1]==d0)
    np.savez_compressed(output,counts=values_array,days=np.arange(settings['horizon_days']+1,dtype=np.int32))
    raw.unlink() # owned temporary raw output, replaced by lossless compressed counts
    events=list(csv.DictReader(Path(str(prefix)+'.csv').open()))
    assert len(events)==settings['replicates'] and len({r['seed'] for r in events})==settings['replicates']
    first=[float(r['first_hit_day']) if r['first_hit_day'] else math.inf for r in events]
    result={'name':name,'label':label,'parameters':values,'alpha':alpha,'theta_days':theta,
            'D_initial_fraction':d_fraction,'C_initial_fraction':settings['initial_C_fraction'],'replicates':settings['replicates'],
            'job_seed':str(job_seed),'ever_D_fraction':sum(int(r['ever_D']) for r in events)/len(events),
            'extant_D_at_four_years_fraction':float(np.mean(values_array[:,-1,1]>0)),
            'mean_C_four_year_fraction':float(values_array[:,-1,0].mean()/n),
            'mean_D_four_year_fraction':float(values_array[:,-1,1].mean()/n),
            **{f'MC_risk_year_{year}':sum(t<=365*year for t in first)/len(first) for year in range(1,5)},
            **timing,'simulator_source_sha256':simulator_sha,'simulation_initial_counts':{'W':n-c0-d0,'C':c0,'D':d0},
            'count_columns':['C','D'],'count_dtype':'int32','algorithm':'exact time-inhomogeneous Moran process; state-specific future-gamma thinning bound; null replacements skipped'}
    meta.write_text(json.dumps(result,indent=2))
    print(f"SIM {name}: {len(events)} trajectories; {timing['seconds']:.1f}s",flush=True)
    return result


def plot_grid(config,metadata):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=ROOT/'results/figures';folder.mkdir(exist_ok=True)
    settings=config['simulation'];n=config['reference']['N']
    for d0 in settings['initial_D_fractions']:
        for typ,index,color in [('single_mutant',0,'#1765b0'),('double_hit',1,'#c82d36')]:
            fig,axes=plt.subplots(3,3,figsize=(11,8),sharex=True,sharey=True,layout='constrained')
            for ia,alpha in enumerate(settings['alpha']):
                for it,theta in enumerate(settings['theta_days']):
                    ax=axes[ia,it]
                    record=next(r for r in metadata if r['label']==f'grid_a{alpha:g}_t{theta:g}' and r['D_initial_fraction']==d0)
                    data=np.load(ROOT/'results/trajectories'/f"{record['name']}.npz")
                    years=data['days']/365;counts=data['counts'][:,:,index]/n*100
                    ax.plot(years,counts.T,color=color,alpha=settings['trajectory_opacity'],lw=.45,rasterized=True)
                    ax.plot(years,counts.mean(axis=0),color=color,lw=2.1,label='Mean of all 200')
                    ax.set_title(fr'$\alpha={alpha:g},\ \theta={theta:g}$ days',fontsize=10)
                    ax.set_xlim(0,4);ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.15)
                    if ia==2:ax.set_xlabel('Years after CAR infusion')
                    if it==0:ax.set_ylabel('Cell fraction (%)')
            title='Single mutant C' if index==0 else 'Double hit D'
            fig.suptitle(f'{title}: 200 exact Moran trajectories per panel | initial D={d0:.0%}\nFixed gamma advantage area = 1 day; C(0)/N = 2%',fontsize=12)
            from matplotlib.lines import Line2D
            fig.legend(handles=[Line2D([0],[0],color=color,alpha=.2,lw=.6,label='Individual trajectory'),
                                Line2D([0],[0],color=color,lw=2.1,label='Mean of all 200 (including zeros)')],
                       loc='outside lower center',ncol=2,frameon=False,fontsize=9)
            name=f'{typ}_initial_D_{d0:g}'
            fig.savefig(folder/f'{name}.pdf');fig.savefig(folder/f'{name}.png',dpi=180);plt.close(fig)
    # Early-time means in one compact comparison, without obscuring rare zero-D trajectories.
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ir,d0 in enumerate(settings['initial_D_fractions']):
        for j,r in enumerate(metadata):
            if not r['label'].startswith('grid_') or r['D_initial_fraction']!=d0:continue
            data=np.load(ROOT/'results/trajectories'/f"{r['name']}.npz")
            for col in range(2):
                avg=data['counts'][:,:,col].mean(axis=0)/n*100
                axes[ir,col].plot(data['days'][:181],avg[:181],label=f"a={r['alpha']:g},t={r['theta_days']:g}",lw=1,
                                  color=plt.cm.Blues(.4+.5*(j%9)/8) if col==0 else plt.cm.Reds(.4+.5*(j%9)/8))
        for col in range(2):axes[ir,col].set(xlabel='Days after CAR infusion',ylabel='Mean cell fraction (%)',title=f"{'C' if col==0 else 'D'}; initial D={d0:.0%}")
    axes[0,1].legend(fontsize=7,ncol=2)
    fig.savefig(folder/'early_means.pdf');fig.savefig(folder/'early_means.png',dpi=180);plt.close(fig)
