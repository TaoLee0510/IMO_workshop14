"""Start the study: validate the engine, fit VAF diagnostics, report readiness.

Initial VAF-only results and standardized risk demonstrations are explicitly
unadjusted model diagnostics. This script never manufactures clinical labels.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
import os
import platform
import sys
from dataclasses import replace
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_gamma_mpl')
import numpy as np
import scipy
from scipy.optimize import least_squares
from scipy.special import logit
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from gamma_model import Parameters, fitness, background, first_passage, first_passage_adaptive, joint_loglik, observation_log_weights, meanfield_predictions
from prepare_data import write_csv

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'
CONFIG=json.loads((ROOT/'config.json').read_text())


def json_write(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False))


def numerical_validation():
    p=Parameters(**CONFIG['reference'])
    assert fitness(-3,3.,10.,p)==(1.,1.)
    fm,fd=fitness(20.,3.,10.,p)
    assert abs(fm-1.1)<1e-14 and abs(fd-fm-.021)<1e-14
    assert fitness(0.,3.,10.,p)[0]==1
    zero=first_passage(alpha=3,theta=10,p=replace(p,mu=0),baseline_vaf=.01)
    assert all(v==0 for v in zero.values())
    assert replace(p,N=10000).K==101 and replace(p,N=1000000).K==10001
    impossible=joint_loglik(alpha=3,theta=10,p=replace(p,N=200,mu=0),baseline_day=0,baseline_vaf=.1,
        observations=[],outcome={'kind':'exact','day':50},max_step=.5)
    assert impossible == -math.inf
    # Independent exact matrix calculation for the one-dimensional filter:
    # constant fM=1, mu=0 makes q constant, with a stochastic D and extra fD.
    small=replace(p,N=100,D_init=1,mu=0,amplitude=0,threshold_fraction=.05,r_D=.01,sigma_vaf=.2)
    q=(.2-small.D_init/small.N)/(1-small.D_init/small.N)
    d=np.arange(small.K,dtype=float)
    fd=1+small.r_D/small.b
    pd=fd*d/(small.N+(fd-1)*d)
    birth=small.b*(small.N-d)*pd
    death=small.b*d*(1-pd)
    generator=diags([birth[:-1],-(birth+death),death[1:]],[-1,0,1],format='csr')
    state=np.zeros(small.K); state[1]=1
    state=expm_multiply(generator*5,state)
    state*=np.exp(observation_log_weights(.101,q,d,small))
    state=expm_multiply(generator*15,state)
    exact_log=math.log(birth[-1]*state[-1])
    comparisons=[]
    for step in (.5,.25,.125):
        got=joint_loglik(alpha=3,theta=10,p=small,baseline_day=0,baseline_vaf=.1,
            observations=[(5,.101)],outcome={'kind':'exact','day':20},max_step=step)
        comparisons.append(dict(step_days=step,log_likelihood=got,absolute_log_error=abs(got-exact_log)))
    assert comparisons[-1]['absolute_log_error']<.02
    assert comparisons[-1]['absolute_log_error']<comparisons[0]['absolute_log_error']
    # Unknown outcomes must integrate over all D, including those above K.
    ll_unknown=joint_loglik(alpha=3,theta=10,p=small,baseline_day=0,baseline_vaf=.1,
        observations=[(20,.11)],outcome={'kind':'unknown'},max_step=.125)
    ll_censor=joint_loglik(alpha=3,theta=10,p=small,baseline_day=0,baseline_vaf=.1,
        observations=[(20,.11)],outcome={'kind':'censor','day':20},max_step=.125)
    assert ll_unknown >= ll_censor-1e-9
    # Full two-dimensional small-N Moran vs hybrid under a gamma pulse.
    # This comparison assesses the M-background approximation, not just BE error.
    n=80; m0=16; pp=replace(p,N=n,mu=.002,r_D=.003,threshold_fraction=.05)
    states=[(m,d) for d in range(pp.K) for m in range(n-d+1)]
    lookup={x:i for i,x in enumerate(states)}
    from scipy.integrate import solve_ivp
    def exact_rhs(t,probs):
        fm,fd=fitness(t,3,5,pp)
        out=np.zeros_like(probs)
        for j,(m,d) in enumerate(states):
            w=n-m-d; total=w+fm*m+fd*d
            pw=w/total; pm=(1-pp.mu)*fm*m/total; pd=(fd*d+pp.mu*fm*m)/total
            for dest,rate in [((m+1,d),pp.b*w*pm),((m-1,d),pp.b*m*pw),((m,d+1),pp.b*w*pd),
                              ((m,d-1),pp.b*d*pw),((m+1,d-1),pp.b*d*pm),((m-1,d+1),pp.b*m*pd)]:
                if rate>0:
                    out[j]-=rate*probs[j]
                    if dest in lookup: out[lookup[dest]]+=rate*probs[j]
        return out
    initial=np.zeros(len(states)); initial[lookup[(m0,0)]]=1
    sol=solve_ivp(exact_rhs,(0,30),initial,t_eval=[10,20,30],rtol=1e-7,atol=1e-11,max_step=.5)
    assert sol.success
    exact_risks=1-sol.y.sum(axis=0)
    hybrid=first_passage(alpha=3,theta=5,p=pp,baseline_vaf=.1,targets=[10,20,30],max_step=.125)
    hybrid_difference=max(abs(exact_risks[j]-hybrid[t]) for j,t in enumerate([10,20,30]))
    assert hybrid_difference<.01
    # Structural recovery of the mean-field initializer under known noiseless data.
    bases=[(0,.01),(0,.05)]
    requested=[(i,t) for i in range(2) for t in [3,7,15,30,60,120,365]]
    true_alpha,true_theta=3.2,18.
    truth=meanfield_predictions(true_alpha,true_theta,p,bases,requested)
    def residual(z):
        pred=meanfield_predictions(1+np.exp(z[0]),np.exp(z[1]),p,bases,requested)
        return logit(2*pred)-logit(2*truth)
    fitted=least_squares(residual,np.log([1.5,30]),bounds=(np.log([.01,.25]),np.log([49,1460])),xtol=1e-9,ftol=1e-9,gtol=1e-9)
    recovered=[1+float(np.exp(fitted.x[0])),float(np.exp(fitted.x[1]))]
    assert max(abs(recovered[0]-true_alpha),abs(recovered[1]-true_theta))<.02
    # Production-sized first-passage timestep refinement at a standardized baseline.
    risk_refinement=[]
    for step in (2.,1.,.5):
        risks=first_passage(alpha=3,theta=20,p=p,baseline_vaf=.01,max_step=step)
        values=list(risks.values())
        assert all(0<=v<=1 for v in values) and all(a<=b for a,b in zip(values,values[1:]))
        risk_refinement.append(dict(step_days=step,risks=risks))
    error=max(abs(risk_refinement[-1]['risks'][t]-risk_refinement[-2]['risks'][t]) for t in CONFIG['risk_days'])
    assert error<.001
    adaptive_risks=first_passage_adaptive(alpha=3,theta=20,p=p,baseline_vaf=.01)
    refined_adaptive_risks=first_passage_adaptive(alpha=3,theta=20,p=p,baseline_vaf=.01,refinement=10)
    adaptive_error=max(abs(adaptive_risks[t]-refined_adaptive_risks[t]) for t in CONFIG['risk_days'])
    assert adaptive_error<1e-7
    result=dict(status='passed',scope='Initial validation, not a comprehensive production inference certification',
        pre_CAR_neutrality=True, gamma_peak=True, zero_second_hit_invariant=True,
        impossible_event_zero_likelihood=True, unknown_outcome_not_censor=True,
        exact_1D_joint_log_likelihood=exact_log,joint_filter_refinements=comparisons,
        small_N_full_Moran_risks=dict(zip([10,20,30],map(float,exact_risks))),
        small_N_hybrid_risks=hybrid,small_N_max_absolute_difference=float(hybrid_difference),
        noiseless_meanfield_parameter_recovery=dict(true_alpha=true_alpha,true_theta=true_theta,
            recovered_alpha=recovered[0],recovered_theta=recovered[1],not_a_joint_parameter_recovery_test=True),
        full_N_risk_refinement=risk_refinement,full_N_last_refinement_max_absolute_difference=float(error),
        adaptive_full_N_risks=adaptive_risks,adaptive_tolerance_refinement_max_absolute_difference=adaptive_error,
        note='BE has appreciable relative error at tiny early probabilities; demonstration risks use adaptive LSODA.')
    json_write(RESULTS/'validation.json',result)
    print('Validation passed.',flush=True)


def load_vaf():
    rows=list(csv.DictReader((ROOT/'data/analysis_longitudinal.csv').open()))
    pids=sorted({r['patient_id'] for r in rows},key=lambda x:int(x.split('_')[1]))
    baselines=[]; requested=[]; observed=[]; source=[]
    for i,pid in enumerate(pids):
        rr=[r for r in rows if r['patient_id']==pid]
        baseline=next(r for r in rr if r['role']=='initial_condition')
        baselines.append((float(baseline['day']),float(baseline['vaf_percent'])/100))
        for r in rr:
            if r['role']=='likelihood_observation':
                requested.append((i,float(r['day'])))
                observed.append(float(r['vaf_percent'])/100)
                source.append(r)
    return pids,baselines,requested,np.array(observed),source


def fit_vaf_diagnostic(run_sensitivity=False):
    pids,baselines,requested,observed,source=load_vaf()
    reference=Parameters(**CONFIG['reference'])
    scenarios=[('reference',reference)]
    if run_sensitivity:
        for key,values in {'amplitude':[.03,.3],'b':[1/30,1.], 'mu':[0.,1e-4], 'r_D':[0.,.01]}.items():
            for value in values:
                scenarios.append((f'{key}_{value:g}',replace(reference,**{key:value})))
    summary=[]
    best_ref=None
    for name,p in scenarios:
        def residual(z):
            pred=meanfield_predictions(1+np.exp(z[0]),np.exp(z[1]),p,baselines,requested)
            # Numerical protection for saturated mean-field fractions, not an event-probability floor.
            return (logit(np.clip(2*pred,1e-14,1-1e-14))-logit(2*observed))/p.sigma_vaf
        candidates=[]
        for alpha,theta in [(2,10),(3,50),(6,15)]:
            fitted=least_squares(residual,np.log([alpha-1,theta]),
                bounds=(np.log([.01,.25]),np.log([49,1460])),max_nfev=150,
                diff_step=1e-4,ftol=2e-7,xtol=2e-7,gtol=2e-7)
            candidates.append(fitted)
        best=min(candidates,key=lambda x:np.dot(x.fun,x.fun))
        alpha,theta=1+float(np.exp(best.x[0])),float(np.exp(best.x[1]))
        singular=np.linalg.svd(best.jac,compute_uv=False)
        condition=float(singular[0]/singular[-1]) if singular[-1]>0 else None
        boundary=bool(np.any(np.abs(best.x-np.log([.01,.25]))<.02) or np.any(np.abs(best.x-np.log([49,1460]))<.02))
        record=dict(scenario=name,analysis_type='PRELIMINARY_VAF_ONLY_MEANFIELD_NOT_CLINICAL_JOINT_FIT',
            patients=len(pids),observations=len(observed),alpha=alpha,theta_days=theta,peak_day=(alpha-1)*theta,
            amplitude=p.amplitude,b=p.b,mu=p.mu,r_D=p.r_D,D_init=p.D_init,
            residual_sum_squares=float(np.dot(best.fun,best.fun)),
            negative_log_likelihood_z=float(.5*np.dot(best.fun,best.fun)+len(observed)*math.log(p.sigma_vaf*math.sqrt(2*math.pi))),
            maximum_abs_standardized_residual=float(np.abs(best.fun).max()),jacobian_condition=condition,
            optimizer_success=bool(best.success),boundary_optimum=boundary,
            start_cost_min=float(min(np.dot(x.fun,x.fun) for x in candidates)),
            start_cost_max=float(max(np.dot(x.fun,x.fun) for x in candidates)),
            clinical_outcomes_used=0,interpretation='diagnostic only; missing outcomes, detection selection and approximation remain unresolved')
        summary.append(record)
        print(json.dumps(record),flush=True)
        if name=='reference':
            best_ref=(alpha,theta)
            pred=meanfield_predictions(alpha,theta,p,baselines,requested)
            predictions=[dict(patient_id=r['patient_id'],gene=r['gene'],variant_id=r['variant_id'],sample_type=r['sample_type'],
                day=r['day'],observed_vaf_percent=100*observed[j],predicted_vaf_percent=100*pred[j],
                standardized_logit_residual=float(best.fun[j])) for j,r in enumerate(source)]
            write_csv(RESULTS/'preliminary_vaf_predictions.csv',predictions)
    write_csv(RESULTS/'preliminary_vaf_fits.csv',summary)
    alpha,theta=best_ref
    risk_rows=[]
    # Standardized illustration only; not per-patient clinical incidence.
    for vaf in (.005,.01,.05):
        risks=first_passage_adaptive(alpha=alpha,theta=theta,p=reference,baseline_vaf=vaf)
        refined=first_passage_adaptive(alpha=alpha,theta=theta,p=reference,baseline_vaf=vaf,refinement=10)
        risk_rows.append(dict(analysis_type='STANDARDIZED_MODEL_DEMONSTRATION_NOT_CLINICALLY_CALIBRATED',
            baseline_vaf_percent=vaf*100,alpha=alpha,theta_days=theta,
            **{f'year_{year}_probability':risks[day] for year,day in enumerate(CONFIG['risk_days'],1)},
            endpoint='aggregate D first passage; competing death omitted',integrator='adaptive_banded_LSODA',
            tolerance_refinement_max_absolute_difference=max(abs(risks[t]-refined[t]) for t in CONFIG['risk_days'])))
    write_csv(RESULTS/'standardized_risk_demonstration.csv',risk_rows)
    # Profile grid is a diagnostic of shape/scale dependence, not a confidence interval.
    profile=[]
    for a in [1.02,1.2,1.5,2.,3.,5.,10.,25.,49.]:
        for th in [.25,1.,3.,10.,30.,100.,300.,1000.,1460.]:
            pred=meanfield_predictions(a,th,reference,baselines,requested)
            rr=(logit(np.clip(2*pred,1e-14,1-1e-14))-logit(2*observed))/reference.sigma_vaf
            profile.append(dict(alpha=a,theta_days=th,negative_log_likelihood_z=float(.5*rr@rr+len(rr)*math.log(reference.sigma_vaf*math.sqrt(2*math.pi)))))
    write_csv(RESULTS/'preliminary_likelihood_grid.csv',profile)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    days=np.linspace(0,365,1000)
    fm,_=fitness(days,alpha,theta,reference)
    axes[0].plot(days,fm-1,color='#235c82');axes[0].set(xlabel='Days after CAR infusion',ylabel='Founder fitness advantage',title='Preliminary VAF-only fitted pulse')
    axes[1].scatter(observed*100,predictions and np.array([r['predicted_vaf_percent'] for r in predictions]),color='#235c82',s=20)
    lim=max(observed.max()*100,max(r['predicted_vaf_percent'] for r in predictions))
    axes[1].plot([0,lim],[0,lim],color='#aaa',linestyle='--')
    axes[1].set(xlabel='Observed VAF (%)',ylabel='Predicted VAF (%)',title=f'Blood-cell index markers; {len(pids)} patients')
    fig.suptitle('Diagnostic only: no clinical outcomes used',fontsize=11)
    fig.savefig(RESULTS/'preliminary_diagnostics.pdf');fig.savefig(RESULTS/'preliminary_diagnostics.png',dpi=160);plt.close(fig)
    return summary


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--mode',choices=['validate','vaf','all'],default='all')
    parser.add_argument('--initial-sensitivity',action='store_true')
    args=parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    start=datetime.now(ZoneInfo('America/New_York')).isoformat()
    if args.mode in ('validate','all'):numerical_validation()
    fits=fit_vaf_diagnostic(args.initial_sensitivity) if args.mode in ('vaf','all') else []
    audit=json.loads((ROOT/'data/audit.json').read_text())
    status=dict(started_at=start,finished_at=datetime.now(ZoneInfo('America/New_York')).isoformat(),
        phase='initial_implementation_and_diagnostics',protocol_pdf='main.pdf',
        validation_completed=(RESULTS/'validation.json').exists(),preliminary_VAF_scenarios=len(fits),
        clinical_joint_fit_completed=False,clinical_joint_fit_ready=audit['full_clinical_joint_fit_ready'],
        clinical_missing_fields=audit['clinical_missing_fields'],
        completed=['English protocol and compiled PDF','source audit','longitudinal index-marker dataset',
                   'hybrid shared-state forward-likelihood engine','initial numerical validation',
                   'VAF-only diagnostic fits if requested','standardized model risk demonstrations if requested'],
        pending=['clinical time-origin/outcome/ascertainment reconciliation','detection and compartment review',
                 'production joint parameter-recovery and likelihood convergence checks',
                 'clinical joint fit','patient uncertainty intervals','complete refitted clinical sensitivity analysis',
                 'clinical 1-4-year probabilities and final XLSX'],
        python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__)
    json_write(RESULTS/'run_status.json',status)


if __name__=='__main__':main()
