"""Resume numerical work with one fresh numerical process at a time.

Only the standard library is imported by this lightweight controller. Each
scenario runs to its atomic checkpoint, then exits and returns its memory to
the operating system before the next scenario starts. No parallel pool is used.
"""
import argparse
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

ROOT=Path(__file__).resolve().parent
def now():return datetime.now(ZoneInfo('America/New_York')).isoformat()
def write_state(**values):
    target=ROOT/'results/serial_controller.json';temp=target.with_suffix('.json.tmp')
    temp.write_text(json.dumps({'updated_at':now(),'controller_pid':os.getpid(),**values},indent=2)+'\n');temp.replace(target)
def names(config):
    result=['reference']
    for mu in config['sensitivity']['mu']:
        for d in config['sensitivity']['D_init']:
            if not (mu==.01 and d==0):result.append(f'mu_{mu:g}_D_{d}')
    for key in ['b','s_D']:
        for value in config['sensitivity'][key]:
            if value!=config['reference'][key]:result.append(f'{key}_{value:.9g}')
    assert len(result)==79 and len(set(result))==79
    return result
def run(command,stage,scenario,environment):
    if shutil.disk_usage(ROOT).free<512*1024**2:raise OSError('Less than 512 MiB free; checkpoints retained')
    print(f'{now()} START {stage} {scenario}',flush=True)
    child=subprocess.Popen(command,cwd=ROOT,env=environment,start_new_session=True)
    write_state(status='running',stage=stage,scenario=scenario,child_pid=child.pid,workers=1,numerical_threads=1)
    try:
        code=child.wait()
    except BaseException:
        try:os.killpg(child.pid,signal.SIGTERM)
        except ProcessLookupError:pass
        try:child.wait(timeout=30)
        except subprocess.TimeoutExpired:
            try:os.killpg(child.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            child.wait()
        raise
    if code:raise RuntimeError(f'{stage} {scenario} failed with return code {code}')
    print(f'{now()} DONE {stage} {scenario}',flush=True)
def main():
    def stop(*_):raise KeyboardInterrupt('Controller stopped by SIGTERM')
    signal.signal(signal.SIGTERM,stop)
    parser=argparse.ArgumentParser()
    parser.add_argument('--stage',choices=['fit','risks','numerical'],default='numerical')
    args=parser.parse_args()
    environment=dict(os.environ)
    environment['TMN_ISOLATE_NUMERICS']='1'
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:
        environment[key]='1'
    config=json.loads((ROOT/'config.json').read_text())
    with (ROOT/'results/serial_controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            for stage in (['fit','risks'] if args.stage=='numerical' else [args.stage]):
                for name in names(config):
                    # The worker verifies signatures before reusing checkpoints.
                    run([sys.executable,'-u','run_sensitivity.py','--stage',stage,'--scenario',name,'--workers','1'],stage,name,environment)
            if args.stage in ['risks','numerical']:
                run([sys.executable,'-u','run_sensitivity.py','--stage','risks','--workers','1'],'risk_consolidation',None,environment)
            write_state(status='complete',stage=args.stage,workers=1,numerical_threads=1,child_pid=None)
        except BaseException as exc:
            write_state(status='stopped',stage=args.stage,reason=str(exc),child_pid=None)
            raise
if __name__=='__main__':main()
