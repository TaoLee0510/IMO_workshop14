"""Joint alpha/theta optimization and conditional likelihood diagnostics."""
import json
import math
import time
from pathlib import Path
from dataclasses import asdict
import numpy as np
from scipy.optimize import minimize
from model import Parameters,cohort_loglik,peak_extra

ROOT=Path(__file__).resolve().parent


def fit_scenario(task):
    name,values,cohort,settings=task
    p=Parameters(**values)
    target=ROOT/'results'/settings.get('output_dir','fits')/f'{name}.json'
    target.parent.mkdir(parents=True,exist_ok=True)
    saved=json.loads(target.read_text()) if target.exists() else None
    if saved and saved.get('optimizer_version',0)>=2 and saved.get('status') in ('fitted','structurally_incompatible') and (p.sigma_vaf>=.15 or saved.get('precision_version',0)>=2) and not settings.get('force',False):return saved
    start=time.time()
    if p.mu==0 and p.D_init==0 and any(r['outcome']['kind']=='exact' for r in cohort):
        out={'scenario':name,'status':'structurally_incompatible','reason':'mu=0 and D_init=0 cannot produce observed tMN events','parameters':values,'optimizer_version':2}
        target.write_text(json.dumps(out,indent=2));return out
    cap=settings.get('cap',512);refinement=settings.get('refinement',1)
    cache={};errors=[];evaluations=[]
    def objective(z):
        key=tuple(np.round(z,9))
        if key not in cache:
            alpha=1+math.exp(float(z[0]));theta=math.exp(float(z[1]))
            try:
                ll,_=cohort_loglik(cohort,alpha,theta,p,cap,refinement)
                value=-ll if math.isfinite(ll) else 1e100
            except (ArithmeticError,ValueError,FloatingPointError) as exc:
                errors.append(str(exc));value=1e100
            cache[key]=value;evaluations.append({'alpha':alpha,'theta':theta,'nll':value})
            if len(evaluations)%20==0:
                print(f"EVAL {name}: n={len(evaluations)} bestNLL={min(cache.values()):.5g} elapsed={time.time()-start:.1f}s",flush=True)
        return cache[key]
    bounds=[(math.log(.01),math.log(49)),(math.log(.25),math.log(1460))]
    # Multiple starting times are important under a fixed-area pulse.
    anchors=settings.get('starts',[(2.,2.),(4.,7.),(4.,30.)])+[(1.01,.25),(2.,.25),(8.,.25),(50.,.25),(8.,21.),(50.,1460.)]
    if saved and 'alpha' in saved:anchors.append((saved['alpha'],saved['theta_days']))
    ranked=sorted(set(anchors),key=lambda pair:objective(np.log([pair[0]-1,pair[1]])))
    starts=ranked[:2]
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
    ll,details=cohort_loglik(cohort,alpha,theta,p,cap,refinement,True)
    # Numerical checks operate on log likelihood, not only total probability mass.
    cap_ll,_=cohort_loglik(cohort,alpha,theta,p,min(p.N,cap*2),refinement)
    refined_ll,_=cohort_loglik(cohort,alpha,theta,p,min(p.N,cap*2),refinement*5)
    cap_delta=abs(cap_ll-ll);tolerance_delta=abs(refined_ll-cap_ll)
    if cap_delta>.03 or tolerance_delta>.03:
        if cap<4096 and not settings.get('expanded',False):
            expanded=dict(settings,cap=cap*2,refinement=refinement*2,expanded=True,force=True,starts=[(alpha,theta)]+starts)
            return fit_scenario((name,values,cohort,expanded))
    full_N_delta=None
    if name=='reference':
        full_ll,_=cohort_loglik(cohort,alpha,theta,p,p.N,refinement)
        full_N_delta=abs(full_ll-ll)
    out={'scenario':name,'optimizer_version':2,'precision_version':2 if p.sigma_vaf<.15 else 1,'status':'fitted' if cap_delta<.03 and tolerance_delta<.03 and (full_N_delta is None or full_N_delta<.03) else 'numerical_refinement_required',
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
    target.write_text(json.dumps(out,indent=2,allow_nan=False))
    (target.parent/f'{name}_evaluations.json').write_text(json.dumps(evaluations))
    print(f"FIT {name}: alpha={alpha:.5g} theta={theta:.5g} jointNLL={-ll:.5g} cap_delta={cap_delta:.2g} seconds={out['seconds']:.1f}",flush=True)
    return out
