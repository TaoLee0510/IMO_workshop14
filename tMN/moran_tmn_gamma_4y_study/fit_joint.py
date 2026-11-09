"""Readiness check and cohort-wide alpha/theta fit of the forward joint likelihood.

Default: validate clinical input without running an optimizer. Use --run only
with reconciled records. No conversion of absent event flags into negatives.
The initial BE likelihood must undergo further convergence checks before use.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from gamma_model import Parameters, joint_loglik
from run_initial import load_vaf, json_write

ROOT=Path(__file__).resolve().parent


def yes(value):
    return str(value).lower() in ('true','1','yes')


def check_data(path):
    pids,baselines,requested,observed,source=load_vaf()
    records=list(csv.DictReader(path.open()))
    if len({r['patient_id'] for r in records})!=len(records):
        raise ValueError('Clinical outcomes must have exactly one row per patient.')
    byid={r['patient_id']:r for r in records}
    issues=[]; outcomes={}
    for pid in pids:
        row=byid.get(pid,{})
        reasons=[]
        for key in ('index_infusion_verified','ascertainment_verified','clinical_likelihood_ready'):
            if not yes(row.get(key,'')):reasons.append(key)
        event=row.get('verified_event_day','')
        censor=row.get('verified_censor_day','')
        if event and censor:
            reasons.append('event_and_censor_both_set_choose_first_outcome')
        elif event:
            if not yes(row.get('first_event_verified','')):reasons.append('first_event_verified')
            day=float(event)
            outcomes[pid]={'kind':'exact','day':day} if day<=1460 else {'kind':'censor','day':1460}
        elif censor:
            outcomes[pid]={'kind':'censor','day':min(1460,float(censor))}
        else:
            reasons.append('verified_event_or_event_free_censor_day')
        if row.get('verified_death_day','') and float(row['verified_death_day'])<=1460:
            reasons.append('death_within_horizon_requires_competing_event_extension')
        if pid in outcomes:
            end=outcomes[pid]['day']
            obs_days=[float(r['day']) for r in source if r['patient_id']==pid]
            if end<=0 or any(day>end or (day==end and outcomes[pid]['kind']=='exact') for day in obs_days):
                reasons.append('observation_after_clinical_outcome_requires_rebuilding_analysis_cohort')
        if reasons:issues.append({'patient_id':pid,'missing_or_unresolved':reasons})
    result={'clinical_joint_fit_ready':not issues and bool(pids),'patients':len(pids),
            'patients_requiring_reconciliation':len(issues),'issues':issues,
            'note':'No clinical optimizer was run by this readiness check.'}
    return result,(pids,baselines,requested,observed,source,outcomes)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--clinical',type=Path,default=ROOT/'data/clinical_reconciliation.csv')
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--max-step',type=float,default=.25)
    args=parser.parse_args()
    results=ROOT/'results';results.mkdir(exist_ok=True)
    report,inputs=check_data(args.clinical)
    json_write(results/'clinical_readiness.json',report)
    print(json.dumps(report,indent=2))
    if not args.run:return
    if not report['clinical_joint_fit_ready']:
        raise SystemExit('Clinical joint fitting deferred: source reconciliation is incomplete; see results/clinical_readiness.json.')
    config=json.loads((ROOT/'config.json').read_text())
    p=Parameters(**config['reference'])
    pids,baselines,requested,observed,source,outcomes=inputs
    obs_by_patient={pid:[] for pid in pids}
    for r,v in zip(source,observed):obs_by_patient[r['patient_id']].append((float(r['day']),float(v)))
    def objective(z):
        alpha,theta=1+float(np.exp(z[0])),float(np.exp(z[1]))
        total=0.
        for pid,(day,vaf) in zip(pids,baselines):
            ll=joint_loglik(alpha=alpha,theta=theta,p=p,baseline_day=day,baseline_vaf=vaf,
                observations=obs_by_patient[pid],outcome=outcomes[pid],max_step=args.max_step)
            if not math.isfinite(ll):return math.inf
            total-=ll
        return total
    fitted=[]
    for alpha,theta in [(2,10),(3,50),(6,15)]:
        fit=minimize(objective,np.log([alpha-1,theta]),method='Nelder-Mead',
            bounds=[(math.log(.01),math.log(49)),(math.log(.25),math.log(1460))],
            options={'maxiter':200,'xatol':1e-4,'fatol':1e-4})
        fitted.append(fit)
    best=min(fitted,key=lambda x:x.fun)
    # Do not call an optimization result validated until timestep refinement agrees.
    out={'analysis_type':'HYBRID_FORWARD_JOINT_FIT_PENDING_PRODUCTION_VALIDATION',
         'alpha':1+float(np.exp(best.x[0])),'theta_days':float(np.exp(best.x[1])),
         'negative_log_likelihood':float(best.fun),'optimizer_success':bool(best.success),
         'patients':len(pids),'max_step_days':args.max_step,
         'competing_death_included':False,'likelihood_convergence_and_uncertainty_pending':True}
    json_write(results/'joint_fit_unvalidated.json',out)
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
