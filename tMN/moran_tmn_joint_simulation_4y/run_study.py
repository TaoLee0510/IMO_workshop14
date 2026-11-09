"""Resumable full study: validation, joint fits, exact simulations, risks, reports."""
import argparse
import csv
import itertools
import json
import math
import os
import platform
import subprocess
import sys
import time
from dataclasses import replace,asdict
from datetime import datetime
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from zoneinfo import ZoneInfo
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_joint_mpl')
import numpy as np
import scipy
from model import Parameters,risk_curve,cohort_loglik
from source_data import prepare
from fitting import fit_scenario
from simulations import simulate_job,plot_grid
from validate import validate

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'


def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def scenarios(config):
    ref=Parameters(**config['reference'])
    answer=[('reference',asdict(ref))]
    for key,values in config['sensitivity'].items():
        for value in values:
            answer.append((f'{key}_{value:g}',asdict(replace(ref,**{key:value}))))
    answer.append(('old_rates_fixed_area',asdict(replace(ref,**config['old_rates_comparator']))))
    answer.append(('vaf_only_comparator',asdict(ref)))
    return answer


def run_fits(config,cohort,workers=3,reference_only=False):
    definitions=scenarios(config)
    if reference_only:definitions=definitions[:1]
    refpath=RESULTS/'fits/reference.json'
    reference=fit_scenario((definitions[0][0],definitions[0][1],cohort,{}))
    tasks=[]
    for name,values in definitions[1:]:
        fit_cohort=[dict(r,outcome={'kind':'unknown'}) for r in cohort] if name=='vaf_only_comparator' else cohort
        tasks.append((name,values,fit_cohort,{'starts':[(reference['alpha'],reference['theta_days']),(2.,2.)],'maxiter':70}))
    fits=[reference]
    if tasks:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for f in as_completed([pool.submit(fit_scenario,task) for task in tasks]):fits.append(f.result())
    fits.sort(key=lambda x:next(i for i,(name,v) in enumerate(definitions) if name==x['scenario']))
    simple=[{k:v for k,v in fit.items() if not isinstance(v,(dict,list))}|fit['parameters'] for fit in fits]
    write_csv(RESULTS/'joint_fits.csv',simple)
    details=[{'scenario':fit['scenario'],**{k:v for k,v in row.items() if k!='VAF_predictions'}} for fit in fits for row in fit.get('patient_details',[])]
    write_csv(RESULTS/'patient_likelihoods.csv',details)
    predictions=[{'scenario':fit['scenario'],'patient_id':row['patient_id'],'gene':row['gene'],**prediction}
                 for fit in fits for row in fit.get('patient_details',[]) for prediction in row.get('VAF_predictions',[])]
    write_csv(RESULTS/'VAF_predictions.csv',predictions)
    return fits


def run_simulations(config,fits,workers=3,grid_only=False):
    tasks=[];ref=config['reference'];settings=config['simulation']
    for a,theta,d in itertools.product(settings['alpha'],settings['theta_days'],settings['initial_D_fractions']):
        tasks.append((f'grid_a{a:g}_t{theta:g}',ref,a,theta,d,settings))
    if not grid_only:
        for fit in fits:
            if fit['status'] not in ('fitted','numerical_refinement_required'):continue
            for d in settings['initial_D_fractions']:
                # Clinical D_init is a baseline-count assumption, distinct from prescribed simulation D(0).
                tasks.append((f"joint_{fit['scenario']}",fit['parameters'],fit['alpha'],fit['theta_days'],d,settings))
    metadata=[]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for f in as_completed([pool.submit(simulate_job,task) for task in tasks]):metadata.append(f.result())
    metadata.sort(key=lambda x:x['name'])
    write_csv(RESULTS/'simulation_summary.csv',[{k:v for k,v in r.items() if not isinstance(v,dict)}|r['parameters'] for r in metadata])
    (RESULTS/'simulation_metadata.json').write_text(json.dumps(metadata,indent=2))
    plot_grid(config,metadata)
    plot_joint_simulations(config,metadata)
    return metadata


def plot_joint_simulations(config,metadata):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=RESULTS/'figures';folder.mkdir(exist_ok=True)
    reference=[r for r in metadata if r['label']=='joint_reference']
    if not reference:return
    fig,axes=plt.subplots(2,2,figsize=(10,7),sharex=True,layout='constrained')
    for ir,d in enumerate(config['simulation']['initial_D_fractions']):
        r=next(x for x in reference if x['D_initial_fraction']==d)
        data=np.load(RESULTS/'trajectories'/f"{r['name']}.npz")
        years=data['days']/365
        for col,color in [(0,'#1765b0'),(1,'#c82d36')]:
            fraction=data['counts'][:,:,col]/r['parameters']['N']*100
            axes[ir,col].plot(years,fraction.T,color=color,alpha=.035,lw=.45,rasterized=True)
            axes[ir,col].plot(years,fraction.mean(axis=0),color=color,lw=2)
            axes[ir,col].set(title=f"{'Single mutant C' if col==0 else 'Double hit D'}; initial D={d:.0%}",ylabel='Cell fraction (%)',xlabel='Years after CAR infusion',ylim=(0,None))
    fit=reference[0]
    fig.suptitle(f"Joint-fit parameters: alpha={fit['alpha']:.3g}, theta={fit['theta_days']:.3g} days; area fixed at 1 day")
    fig.savefig(folder/'joint_optimum_trajectories.pdf');fig.savefig(folder/'joint_optimum_trajectories.png',dpi=180);plt.close(fig)


def risk_task(task):
    fit,cohort=task;name=fit['scenario'];path=RESULTS/'risks'/f'{name}.json'
    path.parent.mkdir(exist_ok=True)
    meta_path=path.with_suffix('.meta.json')
    signature={'alpha':fit.get('alpha'),'theta_days':fit.get('theta_days'),'parameters':fit['parameters']}
    if path.exists() and meta_path.exists() and json.loads(meta_path.read_text())==signature:return json.loads(path.read_text())
    if 'alpha' not in fit:return []
    p=Parameters(**fit['parameters']);rows=[]
    for row in cohort:
        curve=risk_curve(fit['alpha'],fit['theta_days'],p,2*row['baseline_vaf'],row['baseline_day'])
        rows.append({'scenario':name,'population':'patient','patient_id':row['patient_id'],'baseline_founder_fraction':2*row['baseline_vaf'],
                     **{f'year_{year}_probability':curve[float(365*year)] for year in range(1,5)}})
    # Standardized simulations hold single mutant fraction at 0.02.
    for d in (0.,.01):
        pp=replace(p,D_init=round(d*p.N))
        curve=risk_curve(fit['alpha'],fit['theta_days'],pp,.02+d)
        rows.append({'scenario':name,'population':'standardized','initial_C_fraction':.02,'initial_D_fraction':d,
                     **{f'year_{year}_probability':curve[float(365*year)] for year in range(1,5)}})
    path.write_text(json.dumps(rows,indent=2,allow_nan=False))
    meta_path.write_text(json.dumps(signature,indent=2))
    print(f'RISKS {name}: {len(rows)} rows',flush=True)
    return rows


def profile_task(task):
    kind,value,reference,cohort=task
    p=Parameters(**reference['parameters']);a=reference['alpha'];theta=reference['theta_days']
    from scipy.optimize import minimize_scalar
    if kind=='alpha_profile':
        def objective(z):return -cohort_loglik(cohort,value,math.exp(z),p,reference['cap'],1)[0]
        fitted=minimize_scalar(objective,bounds=(math.log(.25),math.log(1460)),method='bounded',options={'maxiter':28,'xatol':.01})
        candidates=[(float(fitted.fun),float(fitted.x))]+[(objective(z),z) for z in (math.log(.25),math.log(1460),math.log(theta))]
        val,z=min(candidates);alpha,th=value,math.exp(z)
    else:
        def objective(z):return -cohort_loglik(cohort,1+math.exp(z),value,p,reference['cap'],1)[0]
        fitted=minimize_scalar(objective,bounds=(math.log(.01),math.log(49)),method='bounded',options={'maxiter':28,'xatol':.01})
        candidates=[(float(fitted.fun),float(fitted.x))]+[(objective(z),z) for z in (math.log(.01),math.log(49),math.log(a-1))]
        val,z=min(candidates);alpha,th=1+math.exp(z),value
    row={'profile':kind,'fixed_value':value,'alpha':alpha,'theta_days':th,'joint_nll':float(val),
         'delta_nll_from_joint_optimum':float(val)-reference['joint_nll'],'profile_optimizer_success':bool(fitted.success)}
    print(f'PROFILE {kind} {value:.3g}: deltaNLL={val-reference["joint_nll"]:.4g}',flush=True)
    return row


def likelihood_profiles(config,cohort,reference):
    target=RESULTS/'likelihood_profiles.csv'
    meta=RESULTS/'profile_reference.json'
    signature={'alpha':reference['alpha'],'theta_days':reference['theta_days'],'joint_nll':reference['joint_nll']}
    if target.exists() and meta.exists() and json.loads(meta.read_text())==signature:return
    a=reference['alpha'];theta=reference['theta_days']
    points=[('alpha_profile',v) for v in sorted(set([1.01,2.,10.,50.,a]))]
    points +=[('theta_profile',v) for v in sorted(set([.25,2.,21.,1460.,theta]))]
    rows=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        for future in as_completed([pool.submit(profile_task,(kind,value,reference,cohort)) for kind,value in points]):rows.append(future.result())
    rows.sort(key=lambda r:(r['profile'],r['fixed_value']))
    write_csv(target,rows)
    meta.write_text(json.dumps(signature,indent=2))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['validate','fit-reference','fit','simulate-grid','simulate','risks','profile','report','all'],default='all')
    parser.add_argument('--workers',type=int,default=3)
    args=parser.parse_args()
    config=json.loads((ROOT/'config.json').read_text());RESULTS.mkdir(exist_ok=True)
    (ROOT/'bin').mkdir(exist_ok=True)
    if not (ROOT/'bin/moran_simulator').exists() or (ROOT/'simulator.cpp').stat().st_mtime>(ROOT/'bin/moran_simulator').stat().st_mtime:
        subprocess.run(['clang++','-O3','-std=c++17',str(ROOT/'simulator.cpp'),'-o',str(ROOT/'bin/moran_simulator')],check=True)
    cohort_path=ROOT/'data/cohort.json'
    cohort=json.loads(cohort_path.read_text()) if cohort_path.exists() and cohort_path.stat().st_mtime>(ROOT/'config.json').stat().st_mtime else prepare()
    began=datetime.now(ZoneInfo('America/New_York')).isoformat()
    if args.stage in ('all','validate'):validate(config,cohort)
    if args.stage in ('all','fit','fit-reference'):
        fits=run_fits(config,cohort,args.workers,args.stage=='fit-reference')
    else:
        fits=[json.loads(p.read_text()) for p in sorted((RESULTS/'fits').glob('*.json')) if not p.stem.endswith('_evaluations')]
    if args.stage in ('all','simulate','simulate-grid'):run_simulations(config,fits,args.workers,args.stage=='simulate-grid')
    if args.stage in ('all','risks'):
        allrows=[]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            for f in as_completed([pool.submit(risk_task,(fit,cohort)) for fit in fits if 'alpha' in fit]):allrows.extend(f.result())
        write_csv(RESULTS/'four_year_risks.csv',allrows)
        reference=next(f for f in fits if f['scenario']=='reference')
        likelihood_profiles(config,cohort,reference)
    if args.stage=='profile':likelihood_profiles(config,cohort,next(f for f in fits if f['scenario']=='reference'))
    if args.stage in ('all','report'):
        from reporting import make_report
        make_report(config,cohort,fits)
    status={'last_stage':args.stage,'started_at':began,'finished_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
            'fixed_area_days':config['reference']['area_days'],'gamma_normalization':'area',
            'joint_fit_exists':(RESULTS/'fits/reference.json').exists(),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
            'platform':platform.platform(),'clinical_policy':'confirmed source-record joint inference; unknown outcomes not negative'}
    (RESULTS/f'stage_{args.stage}_status.json').write_text(json.dumps(status,indent=2))
    (RESULTS/'last_run_status.json').write_text(json.dumps(status,indent=2))


if __name__=='__main__':main()
