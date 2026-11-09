"""79 fixed-parameter scenarios, each re-optimizing the shared Gamma shape/scale."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,csv,hashlib,json,math,time,shutil,resource,fcntl
from functools import wraps
from pathlib import Path
from dataclasses import replace,asdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply
from model import Parameters,cohort_loglik,ForwardFilter,rates,risk_curve
from risk_backend import risk_curve
from fitting import fit_scenario
import memory_control
ROOT=Path(__file__).resolve().parent
CONFIG=json.loads((ROOT/'config.json').read_text());COHORT=json.loads((ROOT/'data/cohort.json').read_text())
P=Parameters(**CONFIG['reference'])
def resource_checkpoint(stage,scenario=None):
    usage=shutil.disk_usage(ROOT)
    if usage.free<512*1024**2:
        raise OSError('Less than 512 MiB free disk space; checkpoint preserved, execution stopped before the next task.')
    checkpoint={'updated_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
        'pid':os.getpid(),'stage':stage,'scenario':scenario,'mode':'isolated','workers':1,
        'pool_jobs':int(os.environ.get('TMN_POOL_JOBS','1')),'numerical_threads':1,
        'peak_process_memory_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        'free_disk_GiB':usage.free/1024**3}
    write_json(ROOT/'results/workers'/f'{os.getpid()}.json',checkpoint)
    write_json(ROOT/'results/runtime.json',checkpoint)
def write_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_suffix(path.suffix+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)
def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r));path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+f'.{os.getpid()}.tmp')
    with temp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    temp.replace(path)
def collection_lock(function):
    @wraps(function)
    def locked(*args,**kwargs):
        with (ROOT/'results/collection.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            return function(*args,**kwargs)
    return locked
def definitions():
    rows=[]
    for mu in CONFIG['sensitivity']['mu']:
        for d in CONFIG['sensitivity']['D_init']:
            name='reference' if mu==.01 and d==0 else f'mu_{mu:g}_D_{d}'
            rows.append((name,asdict(replace(P,mu=mu,D_init=d)),{'family':'mu_D_grid','mu':mu,'D_init':d}))
    for key in ['b','s_D']:
        for value in CONFIG['sensitivity'][key]:
            if value==getattr(P,key):continue
            rows.append((f'{key}_{value:.9g}',asdict(replace(P,**{key:value})),{'family':key+'_OAT','mu':P.mu,'D_init':0}))
    assert len(rows)==79 and len(set(name for name,p,m in rows))==79
    return rows
def benchmark():
    result=[];d=np.arange(P.K,dtype=float);old=json.loads((ROOT/'data/previous_reference.json').read_text())
    for q,mu,d0,end in [(.00544,.01,0,127.),(.98,.1,0,5.),(.04,0.,500,127.)]:
        p=replace(P,mu=mu,D_init=d0,area_days=0);qfun=lambda t:q
        birth,death=rates(d,q,0.,2.,1.,p)
        matrix=diags([birth[:-1],-(birth+death),death[1:]],[-1,0,1],shape=(p.K,p.K),format='csc')
        initial=np.zeros(p.K);initial[d0]=1;exact=expm_multiply(matrix*end,initial)
        density=exact[-1]*birth[-1];assert density>0
        f=ForwardFilter(p,2.,1.,qfun,killed=True,cap=p.K-1,eta=1.,refinement=5)
        f.propagate(0.,end);error=abs(math.log(density)-f.event_logdensity(end))
        assert error<.002,(q,mu,d0,error)
        result.append({'q':q,'mu':mu,'D_init':d0,'day':end,'independent_density_log_error':error})
        print('GENERATOR_CHECK',result[-1],flush=True)
    times=[]
    for cap,refine in [(512,1),(1024,1),(1024,5),(P.N,1)]:
        start=time.time();ll,details=cohort_loglik(COHORT,old['alpha'],old['theta_days'],P,cap,refine,True)
        item={'minimum_cap':cap,'refinement':refine,'joint_nll':-ll,'seconds':time.time()-start,'patient_caps':{r['patient_id']:r['cap'] for r in details}}
        times.append(item);print('BENCHMARK',json.dumps(item),flush=True)
    assert max(abs(r['joint_nll']-times[0]['joint_nll']) for r in times)<.005
    write_json(ROOT/'results/validation/benchmark.json',{'status':'passed','generator_checks':result,'likelihood_convergence':times})
@collection_lock
def collect():
    definitions_by_name={n:(v,m) for n,v,m in definitions()}
    fits=[json.loads(p.read_text()) for p in (ROOT/'results/fits').glob('*.json') if not p.stem.endswith('_evaluations')]
    rows=[];patients=[];pred=[]
    for f in fits:
        values,meta=definitions_by_name[f['scenario']]
        rows.append({k:v for k,v in f.items() if not isinstance(v,(dict,list))}|f['parameters']|meta)
        for r in f.get('patient_details',[]):
            patients.append({'scenario':f['scenario'],**{k:v for k,v in r.items() if k!='VAF_predictions'}})
            pred += [{'scenario':f['scenario'],'patient_id':r['patient_id'],**v} for v in r.get('VAF_predictions',[])]
    write_csv(ROOT/'results/joint_fits.csv',sorted(rows,key=lambda r:r['scenario']!='reference'))
    write_csv(ROOT/'results/patient_likelihoods.csv',patients);write_csv(ROOT/'results/VAF_predictions.csv',pred)
    write_json(ROOT/'results/progress.json',{'updated_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),'expected_scenarios':79,'saved_scenarios':len(fits),'fitted':sum(f['status']=='fitted' for f in fits),'incompatible':sum(f['status']=='structurally_incompatible' for f in fits)})
    return fits
def run_reference():
    old=json.loads((ROOT/'data/previous_reference.json').read_text())
    return fit_scenario(('reference',asdict(P),COHORT,{'starts':[(old['alpha'],old['theta_days']),(2.,1.),(8.,1.),(16.,1.)],'maxiter':120,'nstarts':3}))
def run_fits(workers):
    ref=run_reference();tasks=[]
    for name,p,meta in definitions():
        if name=='reference':continue
        tasks.append((name,p,COHORT,{'starts':[(ref['alpha'],ref['theta_days']),(2.,1.),(8.,1.)],'maxiter':100,'nstarts':2}))
    if workers==1:
        for task in tasks:
            resource_checkpoint('fit',task[0]);out=fit_scenario(task);collect()
            print(f"COMPLETE {out['scenario']}: {out['status']}",flush=True)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for f in as_completed([pool.submit(fit_scenario,t) for t in tasks]):
                out=f.result();collect();print(f"COMPLETE {out['scenario']}: {out['status']}",flush=True)
    collect()
def risk_task(fit):
    path=ROOT/'results/risks'/f"{fit['scenario']}.json"
    signature={'parameters':fit['parameters'],'alpha':fit['alpha'],'theta_days':fit['theta_days'],'model_sha256':hashlib.sha256((ROOT/'model.py').read_bytes()).hexdigest(),'risk_source_sha256':hashlib.sha256((ROOT/'risk_backend.py').read_bytes()).hexdigest()}
    if path.exists() and json.loads(path.read_text()).get('signature')==signature:return json.loads(path.read_text())['rows']
    p=Parameters(**fit['parameters']);rows=[]
    for row in COHORT:
        if os.environ.get('TMN_ISOLATE_NUMERICS')=='1':
            from isolated_numerics import patient_risk
            curve=patient_risk(fit['alpha'],fit['theta_days'],p,2*row['baseline_vaf'],row['baseline_day'])
        else:curve=risk_curve(fit['alpha'],fit['theta_days'],p,2*row['baseline_vaf'],row['baseline_day'])
        rows.append({'scenario':fit['scenario'],'patient_id':row['patient_id'],'gene':row['gene'],'baseline_VAF':row['baseline_vaf'],
                     'baseline_day':row['baseline_day'],'initial_D':p.D_init,'mu':p.mu,
                     **{f'year{y}_probability':curve[float(365*y)] for y in range(1,5)}})
    write_json(path,{'signature':signature,'rows':rows});print(f"RISK {fit['scenario']} complete",flush=True);return rows
def run_risks(workers):
    fits=[f for f in collect() if f['status']=='fitted'];rows=[]
    if workers==1:
        for fit in fits:
            resource_checkpoint('risks',fit['scenario']);rows+=risk_task(fit)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for f in as_completed([pool.submit(risk_task,t) for t in fits]):rows+=f.result()
    write_csv(ROOT/'results/patient_risks.csv',rows)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['benchmark','reference','fit','risks','collect','all'],default='all');parser.add_argument('--workers',type=int,default=1);parser.add_argument('--scenario');args=parser.parse_args()
    if args.scenario:
        if args.workers!=1 or args.stage not in ['fit','risks']:parser.error('--scenario requires a serial fit or risks stage')
        case=next((r for r in definitions() if r[0]==args.scenario),None)
        if case is None:parser.error('Unknown scenario')
        resource_checkpoint(args.stage,args.scenario)
        if args.stage=='fit':
            ref=run_reference()
            out=ref if args.scenario=='reference' else fit_scenario((case[0],case[1],COHORT,{'starts':[(ref['alpha'],ref['theta_days']),(2.,1.),(8.,1.)],'maxiter':100,'nstarts':2}))
            collect();print(f"COMPLETE {out['scenario']}: {out['status']}",flush=True)
        else:
            fit=json.loads((ROOT/'results/fits'/f'{args.scenario}.json').read_text())
            if fit['status']=='fitted':risk_task(fit)
        return
    if args.stage in ['benchmark','all']:benchmark()
    if args.stage=='reference':run_reference();collect()
    if args.stage in ['fit','all']:run_fits(args.workers)
    if args.stage in ['risks','all']:run_risks(args.workers)
    if args.stage=='collect':collect()
if __name__=='__main__':main()
