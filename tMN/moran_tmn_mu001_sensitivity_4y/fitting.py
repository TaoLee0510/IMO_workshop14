"""Joint alpha/theta optimization and conditional likelihood diagnostics."""
import json
import os
import hashlib
import math
import time
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import minimize
from model import Parameters,cohort_loglik,peak_extra
if os.environ.get('TMN_ISOLATE_NUMERICS')=='1':
    from isolated_numerics import cohort_loglik

ROOT=Path(__file__).resolve().parent
def save_json(path,value):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False));temp.replace(path)


def fit_scenario(task):
    name,values,cohort,settings=task
    p=Parameters(**values)
    fit_bounds=settings.get('fit_bounds',json.loads((ROOT/'config.json').read_text())['fit_bounds'])
    amin,amax=map(float,fit_bounds['alpha']);tmin,tmax=map(float,fit_bounds['theta_days'])
    if not (1<amin<amax and 0<tmin<tmax):raise ValueError('Invalid gamma shape/scale bounds')
    bounds=[(math.log(amin-1),math.log(amax-1)),(math.log(tmin),math.log(tmax))]
    signature={'parameters':values,'fit_bounds':fit_bounds,
               'cohort_sha256':hashlib.sha256(json.dumps(cohort,sort_keys=True).encode()).hexdigest(),
               'model_sha256':hashlib.sha256((ROOT/'model.py').read_bytes()).hexdigest()}
    target=ROOT/'results'/settings.get('output_dir','fits')/f'{name}.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    saved=json.loads(target.read_text()) if target.exists() else None
    if saved and saved.get('optimizer_version',0)>=4 and saved.get('fit_signature')==signature and saved.get('status') in ('fitted','structurally_incompatible') and (p.sigma_vaf>=.15 or saved.get('precision_version',0)>=2) and not settings.get('force',False):return saved
    start=time.time()
    if p.mu==0 and p.D_init==0 and any(r['outcome']['kind']=='exact' for r in cohort):
        out={'scenario':name,'status':'structurally_incompatible','reason':'mu=0 and D_init=0 cannot produce observed tMN events','parameters':values,'optimizer_version':4,'fit_bounds':fit_bounds,'fit_signature':signature}
        save_json(target,out);return out
    cap=settings.get('cap',512);refinement=settings.get('refinement',1)
    cache={};errors=[];evaluations=[]
    def objective(z):
        z=np.clip(np.asarray(z,dtype=float),[v[0] for v in bounds],[v[1] for v in bounds])
        key=tuple(np.round(z,9))
        if key not in cache:
            alpha=1+math.exp(float(z[0]));theta=math.exp(float(z[1]))
            try:
                ll,_=cohort_loglik(cohort,alpha,theta,p,cap,refinement)
                value=-ll if math.isfinite(ll) else 1e100
            except (ArithmeticError,ValueError,FloatingPointError) as exc:
                errors.append({'alpha':alpha,'theta_days':theta,'error':str(exc)});value=1e100
                print(f"EVALUATION_ERROR {name}: alpha={alpha:.12g} theta={theta:.12g} {exc}",flush=True)
            cache[key]=value;evaluations.append({'alpha':alpha,'theta':theta,'nll':value})
            if len(evaluations)%20==0:
                print(f"EVAL {name}: n={len(evaluations)} bestNLL={min(cache.values()):.5g} elapsed={time.time()-start:.1f}s",flush=True)
        return cache[key]
    # Multiple starting times are important under a fixed-area pulse.
    anchors=settings.get('starts',[(2.,2.),(4.,7.),(4.,30.)])+[(amin,tmin),(2.,tmin),(8.,tmin),(amax,tmin),(8.,21.),(amax,tmax)]
    if saved and 'alpha' in saved:anchors.append((saved['alpha'],saved['theta_days']))
    evaluation_path=target.parent/f'{name}_evaluations.json'
    if evaluation_path.exists():
        try:
            previous_evaluations=json.loads(evaluation_path.read_text())
            finite=[e for e in previous_evaluations if math.isfinite(e['nll']) and e['nll']<1e90]
            if finite:
                previous_best=min(finite,key=lambda e:e['nll'])
                anchors.append((previous_best['alpha'],previous_best['theta']))
                # Previous values only supply a starting point; every likelihood
                # is evaluated again under the current inputs and source.
        except (ValueError,KeyError,TypeError):pass
    anchors=[(min(max(float(a),amin),amax),min(max(float(t),tmin),tmax)) for a,t in anchors]
    ranked=sorted(set(anchors),key=lambda pair:objective(np.log([pair[0]-1,pair[1]])))
    starts=ranked[:settings.get('nstarts',2)]
    solutions=[]
    for alpha,theta in starts:
        z=np.log([alpha-1,theta])
        simplex=np.tile(z,(3,1))
        for j in range(2):simplex[j+1,j]+=.35 if z[j]+.35<=bounds[j][1] else -.35
        result=minimize(objective,z,method='Nelder-Mead',bounds=bounds,
                        options={'maxiter':settings.get('maxiter',90),'xatol':.003,'fatol':.002,'adaptive':False,'initial_simplex':simplex})
        solutions.append(result)
    best=min(solutions,key=lambda x:x.fun)
    if best.fun>=1e90:raise ArithmeticError(f'All joint evaluations failed in {name}: {errors[-3:]}')
    alpha=1+math.exp(float(best.x[0]));theta=math.exp(float(best.x[1]))
    save_json(target.parent/f'{name}_evaluations.json',evaluations)
    print(f"OPTIMIZED {name}: alpha={alpha:.7g} theta={theta:.7g} NLL={best.fun:.7g}; precision checks starting",flush=True)
    ll,details=cohort_loglik(cohort,alpha,theta,p,cap,refinement,True)
    # Numerical checks operate on log likelihood, not only total probability mass.
    cap_ll,_=cohort_loglik(cohort,alpha,theta,p,min(p.N,cap*2),refinement)
    refined_ll,_=cohort_loglik(cohort,alpha,theta,p,min(p.N,cap*2),refinement*5)
    cap_delta=abs(cap_ll-ll);tolerance_delta=abs(refined_ll-cap_ll)
    if cap_delta>.005 or tolerance_delta>.005:
        if cap<p.N and settings.get('expansion_count',0)<8:
            expanded=dict(settings,cap=min(p.N,cap*2),refinement=refinement*2,expansion_count=settings.get('expansion_count',0)+1,force=True,starts=[(alpha,theta)]+starts)
            return fit_scenario((name,values,cohort,expanded))
    full_N_delta=None
    if name=='reference':
        full_ll,_=cohort_loglik(cohort,alpha,theta,p,p.N,refinement)
        full_N_delta=abs(full_ll-ll)
    out={'scenario':name,'optimizer_version':4,'fit_bounds':fit_bounds,'fit_signature':signature,'precision_version':2 if p.sigma_vaf<.15 else 1,'status':'fitted' if cap_delta<.005 and tolerance_delta<.005 and (full_N_delta is None or full_N_delta<.005) else 'numerical_refinement_required',
         'parameters':values,'alpha':alpha,'theta_days':theta,'peak_day':(alpha-1)*theta,
         'peak_additional_fitness':peak_extra(alpha,theta,p),'area_days':p.area_days,
         'joint_nll':-ll,'VAF_nll':-sum(r['VAF_loglik'] for r in details),
         'clinical_given_VAF_nll':-sum(r['clinical_given_VAF_loglik'] for r in details),
         'clinical_events_used':sum(r['outcome']['kind']=='exact' for r in cohort),'censor_records_used':sum(r['outcome']['kind']=='censor' for r in cohort),
         'unknown_outcomes':sum(r['outcome']['kind']=='unknown' for r in cohort),'patients':len(cohort),'VAF_observations':sum(len(r['observations']) for r in cohort),
         'optimizer_success':bool(best.success),'start_nlls':[float(r.fun) for r in solutions],
         'boundary_optimum':any(min(abs(best.x[j]-bounds[j][0]),abs(best.x[j]-bounds[j][1]))<.02 for j in range(2)),
         'cap':cap,'cap_refinement_loglik_difference':cap_delta,'tolerance_refinement_loglik_difference':tolerance_delta,
         'full_N_loglik_difference':full_N_delta,
         'patient_details':details,'evaluation_count':len(evaluations),'numerical_errors':errors[:10],
         'seconds':time.time()-start,'inference':'shared-state hybrid joint likelihood; source-based exploratory clinical inference'}
    save_json(target,out)
    save_json(target.parent/f'{name}_evaluations.json',evaluations)
    print(f"FIT {name}: alpha={alpha:.5g} theta={theta:.5g} jointNLL={-ll:.5g} cap_delta={cap_delta:.2g} seconds={out['seconds']:.1f}",flush=True)
    return out
