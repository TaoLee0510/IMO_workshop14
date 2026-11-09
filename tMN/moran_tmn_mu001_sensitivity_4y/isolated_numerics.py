"""Execute one likelihood or one patient's risk forecast, then release memory.

The child imports the same owning model and receives the same arguments. Only
process lifetime changes; no equations, tolerances, random draws or fitted
parameter definitions are altered. Children never launch numerical children.
"""
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import resource
import time
from datetime import datetime
from zoneinfo import ZoneInfo

ROOT=Path(__file__).resolve().parent
PREFIX='TMN_ISOLATED_RESULT '
def invoke(payload):
    environment=dict(os.environ)
    environment['TMN_ISOLATE_NUMERICS']='0'
    result=subprocess.run([sys.executable,str(Path(__file__).resolve())],
        input=json.dumps(payload),capture_output=True,text=True,env=environment,cwd=ROOT)
    records=[line[len(PREFIX):] for line in result.stdout.splitlines() if line.startswith(PREFIX)]
    if not records:
        raise ArithmeticError(f'Isolated calculation failed ({result.returncode}): {result.stderr[-2000:]} {result.stdout[-1000:]}')
    answer=json.loads(records[-1])
    if not answer['ok']:raise ArithmeticError(answer['error'])
    return answer['result']

def cohort_loglik(cohort,alpha,theta,p,cap=512,refinement=1,decompose=False):
    from dataclasses import asdict
    ll,details=invoke({'operation':'likelihood','cohort':cohort,'alpha':alpha,'theta':theta,
        'parameters':asdict(p),'cap':cap,'refinement':refinement,'decompose':decompose})
    return (-math.inf if ll is None else ll),details

def patient_risk(alpha,theta,p,founder,baseline):
    from dataclasses import asdict
    values=invoke({'operation':'risk','alpha':alpha,'theta':theta,'parameters':asdict(p),
                   'founder':founder,'baseline':baseline})
    return {float(k):v for k,v in values.items()}

def main():
    start=time.time()
    payload=json.loads(sys.stdin.read())
    try:
        import model
        import memory_control
        p=model.Parameters(**payload['parameters'])
        if payload['operation']=='likelihood':
            ll,details=model.cohort_loglik(payload['cohort'],payload['alpha'],payload['theta'],p,
                payload['cap'],payload['refinement'],payload['decompose'])
            result=[ll if math.isfinite(ll) else None,details]
        elif payload['operation']=='risk':
            from risk_backend import risk_curve
            result=risk_curve(payload['alpha'],payload['theta'],p,payload['founder'],payload['baseline'])
        else:raise ValueError('Unknown isolated operation')
        answer={'ok':True,'result':result}
    except Exception as exc:
        answer={'ok':False,'error':f'{type(exc).__name__}: {exc}'}
    record={'time':datetime.now(ZoneInfo('America/New_York')).isoformat(),'pid':os.getpid(),
        'operation':payload['operation'],'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,
        'seconds':time.time()-start,'ok':answer['ok']}
    with (ROOT/'results/resource_observations.jsonl').open('a') as log:log.write(json.dumps(record)+'\n')
    print(PREFIX+json.dumps(answer,allow_nan=False),flush=True)
if __name__=='__main__':main()
