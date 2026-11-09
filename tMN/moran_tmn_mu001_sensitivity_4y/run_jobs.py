"""Bounded scenario jobs, with one isolated numerical child per job.

The controller uses only the standard library. Every likelihood evaluation or
patient forecast exits before the same job starts another numerical calculation.
Scenario stages run in order: fitting first, then risks. Use --jobs 3 for the
user-authorized three-job configuration. Shared CSV collection is file-locked.
"""
import argparse
from collections import deque
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from zoneinfo import ZoneInfo
from run_serial import names

ROOT=Path(__file__).resolve().parent
def now():return datetime.now(ZoneInfo('America/New_York')).isoformat()
def write_state(**values):
    target=ROOT/'results/serial_controller.json';temp=target.with_suffix(f'.{os.getpid()}.json.tmp')
    temp.write_text(json.dumps({'updated_at':now(),'controller_pid':os.getpid(),'controller_source':'run_jobs.py',**values},indent=2)+'\n');temp.replace(target)
def terminate(jobs):
    for job in jobs:
        try:os.killpg(job['process'].pid,signal.SIGTERM)
        except ProcessLookupError:pass
    for job in jobs:
        try:job['process'].wait(timeout=10)
        except subprocess.TimeoutExpired:
            try:os.killpg(job['process'].pid,signal.SIGKILL)
            except ProcessLookupError:pass
            job['process'].wait()
        job['log'].close()
def run_pool(stage,scenarios,limit,environment):
    pending=deque(scenarios);active={};completed=0;last_update=0.
    (ROOT/'logs').mkdir(exist_ok=True)
    try:
        while pending or active:
            for pid,job in list(active.items()):
                code=job['process'].poll()
                if code is None:continue
                job['log'].close();del active[pid]
                if code:raise RuntimeError(f"{stage} {job['scenario']} failed ({code}); see {job['log_path']}")
                completed+=1
                print(f"{now()} DONE {stage} {job['scenario']} ({completed}/{len(scenarios)})",flush=True)
            while pending and len(active)<limit:
                if shutil.disk_usage(ROOT).free<512*1024**2:
                    raise OSError('Less than 512 MiB free; checkpoints preserved')
                scenario=pending.popleft()
                command=[sys.executable,'-u','run_sensitivity.py','--stage',stage,'--scenario',scenario,'--workers','1']
                if stage=='risk_consolidation':command=[sys.executable,'-u','run_sensitivity.py','--stage','risks','--workers','1']
                path=ROOT/'logs'/f'{stage}_{scenario}.log';log=path.open('a')
                log.write(f'\n{now()} JOB START\n');log.flush()
                child=subprocess.Popen(command,cwd=ROOT,env=environment,start_new_session=True,stdout=log,stderr=subprocess.STDOUT)
                active[child.pid]={'process':child,'scenario':scenario,'log':log,'log_path':str(path)}
                print(f'{now()} START {stage} {scenario} PID={child.pid}',flush=True)
            if time.monotonic()-last_update>=2 or not active:
                write_state(status='running',stage=stage,workers=limit,numerical_threads=1,
                    pending_scenarios=len(pending),completed_in_stage=completed,total_in_stage=len(scenarios),
                    active_jobs=[{'pid':pid,'scenario':job['scenario'],'log':job['log_path']} for pid,job in active.items()])
                last_update=time.monotonic()
            if active:time.sleep(.2)
    finally:
        terminate(list(active.values()))
def main():
    def stop(*_):raise KeyboardInterrupt('Controller stopped by SIGTERM')
    signal.signal(signal.SIGTERM,stop)
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['fit','risks','numerical'],default='numerical')
    parser.add_argument('--jobs',type=int,default=1)
    parser.add_argument('--scenarios',nargs='+')
    args=parser.parse_args()
    if not 1<=args.jobs<=3:parser.error('Use 1 to 3 jobs')
    environment=dict(os.environ,TMN_ISOLATE_NUMERICS='1',TMN_POOL_JOBS=str(args.jobs))
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:
        environment[key]='1'
    config=json.loads((ROOT/'config.json').read_text());scenarios=names(config)
    if args.scenarios:
        if not set(args.scenarios)<=set(scenarios):parser.error('Unknown scenario')
        scenarios=list(dict.fromkeys(args.scenarios))
    with (ROOT/'results/serial_controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            for stage in (['fit','risks'] if args.stage=='numerical' else [args.stage]):
                run_pool(stage,scenarios,args.jobs,environment)
            if args.stage in ['risks','numerical'] and not args.scenarios:
                run_pool('risk_consolidation',['all'],1,environment)
            write_state(status='complete',stage=args.stage,workers=args.jobs,numerical_threads=1,active_jobs=[])
        except BaseException as exc:
            write_state(status='stopped',stage=args.stage,workers=args.jobs,reason=str(exc),active_jobs=[])
            raise
if __name__=='__main__':main()
