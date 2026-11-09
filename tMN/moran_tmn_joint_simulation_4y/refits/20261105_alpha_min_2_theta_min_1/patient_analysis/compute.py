"""Checkpointed patient-specific simulations, shared-state diagnostics and risks."""
import argparse,subprocess,time,math,hashlib,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from common import *
from scipy.special import logit,expit,ndtr,ndtri_exp,log_ndtr,logsumexp
from scipy.optimize import brentq,minimize_scalar
from model import Backgrounds,ForwardFilter,tilt_for_horizon,cohort_loglik

RISK_DAYS=sorted(set(range(60,1461,60))|{30,127,365,589,730,1095,1460})

def simulate(row):
    out=ROOT/'results/trajectories'/row['patient_id'];meta=out.with_suffix('.json')
    out.parent.mkdir(exist_ok=True,parents=True)
    signature=sha(ROOT/'patient_simulator.cpp')+MANIFEST['fit_sha256']+json.dumps(row,sort_keys=True)
    seed=int.from_bytes(hashlib.sha256((signature+str(MANIFEST['master_seed'])).encode()).digest()[:8],'little')
    if meta.exists() and json.loads(meta.read_text()).get('signature')==signature:return json.loads(meta.read_text())
    init=initial(row);reps=MANIFEST['replicates_per_patient'];horizon=MANIFEST['horizon_days']
    command=[str(ROOT/'bin/patient_simulator'),str(out),str(P.N),str(P.b),str(P.s_D),str(P.mu),str(P.area_days),str(A),str(THETA),
             str(init['initial_C']),str(P.D_init),str(horizon),str(reps),str(seed),str(P.threshold_fraction),str(row['baseline_day'])]
    result=subprocess.run(command,capture_output=True,text=True,check=True)
    raw=out.with_suffix('.bin');cd=np.fromfile(raw,dtype=np.int32).reshape(reps,horizon+1,2)
    counts=np.concatenate(((P.N-cd.sum(axis=2))[:,:,None],cd),axis=2).astype(np.int32)
    assert np.all(counts>=0) and np.all(counts.sum(axis=2)==P.N)
    np.savez_compressed(out.with_suffix('.npz'),counts=counts,days=np.arange(horizon+1,dtype=np.int32),columns=np.array(['W','C','D']))
    raw.unlink()
    events=list(csv.DictReader(out.with_suffix('.csv').open()));assert len(events)==reps
    assert all(int(e['car_C'])==counts[j,0,1] and int(e['car_D'])==counts[j,0,2] for j,e in enumerate(events))
    hits=[float(e['first_hit_day']) if e['first_hit_day'] else math.inf for e in events]
    rec={**init,'signature':signature,'job_seed':str(seed),'replicates':reps,**json.loads(result.stdout),
         'mean_C_year4':float(counts[:,-1,1].mean()),'mean_D_year4':float(counts[:,-1,2].mean()),
         'fraction_ever_D':sum(int(e['ever_D']) for e in events)/reps,
         **{f'MC_probability_year{y}':sum(t<=365*y for t in hits)/reps for y in range(1,5)}}
    write_json(meta,rec);print(f"SIM {row['patient_id']}: {reps} paths, {rec['seconds']:.1f}s",flush=True);return rec

def predictive(f,time,vaf=None):
    d,w,marked=weights(f,time);means=logit(marked)
    def cdf(z):return float(w@ndtr((z-means)/P.sigma_vaf))
    lo=float(means.min()-10*P.sigma_vaf);hi=float(means.max()+10*P.sigma_vaf)
    result={'day':float(time),'mean_VAF':float(w@marked/2)}
    for label,prob in [('lower_VAF',.025),('upper_VAF',.975)]:
        result[label]=float(expit(brentq(lambda z:cdf(z)-prob,lo,hi))/2)
    if vaf is not None:
        z=(logit(2*vaf)-means)/P.sigma_vaf
        logw=np.full_like(w,-np.inf);positive=w>0;logw[positive]=np.log(w[positive])
        logcdf=float(logsumexp(logw+log_ndtr(z)));logsf=float(logsumexp(logw+log_ndtr(-z)))
        residual=float(ndtri_exp(logcdf)) if logcdf<math.log(.5) else -float(ndtri_exp(logsf))
        pit=cdf(logit(2*vaf))
        result.update(observed_VAF=vaf,PIT=pit,log_predictive_CDF=logcdf,log_predictive_survival=logsf,predictive_residual=residual,
                      approximate_z_residual=float((logit(2*vaf)-logit(2*result['mean_VAF']))/P.sigma_vaf))
    return result

def vaf_task(row):
    path=ROOT/'results/patients'/f"{row['patient_id']}_VAF.json"
    if path.exists():
        saved=json.loads(path.read_text())
        if saved.get('diagnostic_version')==2:return saved
    known=row['outcome']['kind']!='unknown';end=row['outcome'].get('day',1460.)
    bg=Backgrounds([row],A,THETA,P,1)
    eta=tilt_for_horizon(A,THETA,P,end,row['baseline_day']) if known else 0.
    f=ForwardFilter(P,A,THETA,bg.for_patient(0),killed=known,cap=FIT['cap'],eta=eta,refinement=1)
    if row['baseline_day']<0:
        f.propagate(row['baseline_day'],0.)
        if known:f.logscale-=f.physical_logmass();f.state[-1]=0.
    obs={x['day']:x['vaf'] for x in row['observations']}
    times=sorted(set([0.,float(end)]+list(np.arange(1,min(end,60)+1,2))+list(np.arange(90,end,60))+list(obs)))
    curves=[];predictions=[];start=0.
    for t in times:
        f.propagate(start,t);point=predictive(f,t,obs.get(t));curves.append(point)
        if t in obs:
            predictions.append(point.copy());f.observation(obs[t],t)
        start=t
    ll=f.event_logdensity(end) if row['outcome']['kind']=='exact' else f.physical_logmass()
    saved=next(x for x in FIT['patient_details'] if x['patient_id']==row['patient_id'])
    delta=abs(ll-saved['joint_loglik'])
    assert delta<.01,(row['patient_id'],delta)
    for pt,old in zip(predictions,saved['VAF_predictions']):
        assert abs(pt['mean_VAF']-old['predicted_vaf'])<1e-5
    rec={'patient_id':row['patient_id'],'diagnostic_version':2,'curves':curves,'predictions':predictions,'recomputed_joint_loglik':ll,'saved_loglik_difference':delta}
    write_json(path,rec);print(f"VAF {row['patient_id']}: {len(predictions)} observations, logL difference {delta:.2g}",flush=True);return rec

def risk_at(row,t,refinement=2):
    bg=Backgrounds([row],A,THETA,P,refinement)
    eta=tilt_for_horizon(A,THETA,P,t,row['baseline_day'])
    f=ForwardFilter(P,A,THETA,bg.for_patient(0),killed=True,cap=P.K-1,eta=eta,refinement=refinement)
    if row['baseline_day']<0:
        f.propagate(row['baseline_day'],0.);f.logscale-=f.physical_logmass();f.state[-1]=0.
    f.propagate(0.,t);mass=f.physical_logmass(include_absorbed=True)
    assert abs(mass)<2e-5,(row['patient_id'],t,mass)
    lp=f.log_cdf()-mass
    return {'patient_id':row['patient_id'],'day':t,'probability':math.exp(lp) if lp>-745 else None,
            'log10_probability':lp/math.log(10),'numeric_underflow':lp<=-745}

def risk_task(row):
    path=ROOT/'results/patients'/f"{row['patient_id']}_risks.json"
    if path.exists():return json.loads(path.read_text())
    points=[]
    for t in RISK_DAYS:
        points.append(risk_at(row,t))
    assert all(a['log10_probability']<=b['log10_probability']+1e-5 for a,b in zip(points,points[1:]))
    rec={'patient_id':row['patient_id'],'conditioning':'own baseline and event-free CAR day 0; no follow-up VAF or eventual outcome used','points':points}
    write_json(path,rec);print(f"RISK {row['patient_id']}: P365={next(r['probability'] for r in points if r['day']==365):.4g}",flush=True);return rec

def likelihood_task(task):
    kind,x,y=task;tag=f'{kind}_{x:.9g}_{y:.9g}'.replace('.','p')
    path=ROOT/'results/likelihood'/f'{tag}.json'
    if path.exists():return json.loads(path.read_text())
    if kind=='surface':
        a,t=x,y;nll=-cohort_loglik(COHORT,a,t,P,FIT['cap'],1)[0]
        rec={'kind':kind,'alpha':a,'theta_days':t,'joint_nll':nll,'delta_nll':nll-FIT['joint_nll']}
    else:
        if kind=='alpha_profile':
            fun=lambda z:-cohort_loglik(COHORT,x,math.exp(z),P,FIT['cap'],1)[0]
            bounds=(0.,math.log(1460));edges=[0.,math.log(1460)]
        else:
            fun=lambda z:-cohort_loglik(COHORT,1+math.exp(z),x,P,FIT['cap'],1)[0]
            bounds=(0.,math.log(49));edges=[0.,math.log(7),math.log(A-1),math.log(49)]
        sol=minimize_scalar(fun,bounds=bounds,method='bounded',options={'maxiter':24,'xatol':.01})
        nll,z=min([(float(sol.fun),float(sol.x))]+[(fun(z),z) for z in edges])
        a,t=(x,math.exp(z)) if kind=='alpha_profile' else (1+math.exp(z),x)
        rec={'kind':kind,'fixed_value':x,'alpha':a,'theta_days':t,'joint_nll':nll,'delta_nll':nll-FIT['joint_nll'],
             'optimizer_success':bool(sol.success),'method':'bounded nuisance search plus explicit candidate anchors'}
    write_json(path,rec);print(f"LIK {kind}: a={a:.4g}, t={t:.4g}, deltaNLL={rec['delta_nll']:.4g}",flush=True);return rec

def pool_run(function,tasks,workers):
    with ProcessPoolExecutor(max_workers=workers) as pool:return [f.result() for f in as_completed([pool.submit(function,t) for t in tasks])]

def consolidate():
    patients=ROOT/'results/patients'
    if all((patients/f"{r['patient_id']}_VAF.json").exists() for r in COHORT):
        obs=[{'patient_id':r['patient_id'],'gene':r['gene'],**p} for r in COHORT for p in json.loads((patients/f"{r['patient_id']}_VAF.json").read_text())['predictions']]
        write_csv(ROOT/'results/VAF_diagnostics.csv',obs)
    if all((patients/f"{r['patient_id']}_risks.json").exists() for r in COHORT):
        rows=[];risk_rows=[]
        for r in COHORT:
            points=json.loads((patients/f"{r['patient_id']}_risks.json").read_text())['points'];risk_rows+=points
            rows.append({**initial(r),**{f'year{y}_probability':next(p['probability'] for p in points if p['day']==365*y) for y in range(1,5)}})
        write_csv(ROOT/'results/patient_risks.csv',rows);write_csv(ROOT/'results/risk_curves.csv',risk_rows)
    records=[json.loads(p.read_text()) for p in (ROOT/'results/trajectories').glob('*.json')]
    if len(records)==14:
        write_csv(ROOT/'results/simulation_summary.csv',[{k:v for k,v in r.items() if k!='signature'} for r in records])
        daily=[];yearly=[]
        for row in COHORT:
            data=np.load(ROOT/'results/trajectories'/f"{row['patient_id']}.npz");counts=data['counts'];days=data['days']
            for day in days:
                vals=counts[:,day,:];lo,hi=np.quantile(vals,[.025,.975],axis=0)
                rec={'patient_id':row['patient_id'],'day':int(day),**{f'mean_{v}_cells':float(vals[:,j].mean()) for j,v in enumerate(['W','C','D'])},
                     'lower_C_cells':float(lo[1]),'upper_C_cells':float(hi[1]),'lower_D_cells':float(lo[2]),'upper_D_cells':float(hi[2])}
                daily.append(rec)
                if day in [0,365,730,1095,1460]:yearly.append(rec)
        write_csv(ROOT/'results/trajectory_daily_summary.csv',daily);write_csv(ROOT/'results/trajectory_yearly_summary.csv',yearly)
    lik=[json.loads(p.read_text()) for p in (ROOT/'results/likelihood').glob('*.json')]
    if lik:
        write_csv(ROOT/'results/likelihood_surface.csv',[r for r in lik if r['kind']=='surface'])
        write_csv(ROOT/'results/likelihood_profiles.csv',[r for r in lik if r['kind']!='surface'])
    write_csv(ROOT/'results/initial_states.csv',[initial(r) for r in COHORT])

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['simulate','vaf','risks','likelihood','consolidate','all'],default='all');parser.add_argument('--workers',type=int,default=3);args=parser.parse_args()
    verify_inputs()
    for d in ['patients','trajectories','likelihood']:(ROOT/'results'/d).mkdir(exist_ok=True,parents=True)
    if args.stage in ['simulate','all']:
        subprocess.run(['clang++','-O3','-std=c++17',str(ROOT/'patient_simulator.cpp'),'-o',str(ROOT/'bin/patient_simulator')],check=True)
        pool_run(simulate,COHORT,args.workers)
    if args.stage in ['vaf','all']:pool_run(vaf_task,COHORT,args.workers)
    if args.stage in ['risks','all']:pool_run(risk_task,COHORT,args.workers)
    if args.stage in ['likelihood','all']:
        aa=[2.,4.,8.,12.,A,20.,30.,40.,50.];tt=[1.,1.5,2.,4.,7.,14.,30.,90.,365.,1460.]
        tasks=[('surface',a,t) for a in aa for t in tt]+[('alpha_profile',a,0.) for a in aa if a!=40]+[('theta_profile',t,0.) for t in [1.,1.5,2.,7.,30.,365.,1460.]]
        pool_run(likelihood_task,tasks,args.workers)
    consolidate();print(f'STAGE {args.stage} complete',flush=True)
if __name__=='__main__':main()
