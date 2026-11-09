import csv,json,hashlib,math,os
from pathlib import Path
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_patient_mpl')
import numpy as np
from model import Parameters
import memory_control
ROOT=Path(__file__).resolve().parent
FIT=json.loads((ROOT/'inputs/fit.json').read_text())
COHORT=json.loads((ROOT/'inputs/cohort.json').read_text())
MANIFEST=json.loads((ROOT/'analysis_manifest.json').read_text())
P=Parameters(**FIT['parameters']);A=FIT['alpha'];THETA=FIT['theta_days']
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temp.replace(path)
def write_csv(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def initial(row):
    c=round(2*P.N*row['baseline_vaf'])-P.D_init
    return {'patient_id':row['patient_id'],'gene':row['gene'],'sample_type':row['sample_type'],
            'baseline_day':row['baseline_day'],'baseline_VAF':row['baseline_vaf'],
            'baseline_founder_fraction':2*row['baseline_vaf'],'initial_W':P.N-c-P.D_init,'initial_C':c,'initial_D':P.D_init,
            'rounding_C_cells':c+P.D_init-2*P.N*row['baseline_vaf'],'outcome':row['outcome']['kind'],'outcome_day':row['outcome'].get('day'),
            'followup_VAF_count':len(row['observations'])}
def verify_inputs():
    assert sha(ROOT/'inputs/fit.json')==MANIFEST['fit_sha256']
    assert sha(ROOT/'inputs/cohort.json')==MANIFEST['cohort_sha256']
    assert sha(ROOT/'model.py')==MANIFEST['model_sha256']
    assert len(COHORT)==14 and P.area_days==1 and FIT['fit_bounds']=={'alpha':[2.,50.],'theta_days':[1.,1460.]}
def weights(f,time):
    positive=f.state[:len(f.d)]>0
    d=f.d[positive];logs=np.log(f.state[:len(f.d)][positive])-f.eta*d
    from scipy.special import logsumexp
    w=np.exp(logs-logsumexp(logs));q=f.qfun(time)
    return d,w,np.clip(q+(1-q)*d/P.N,1e-15,1-1e-15)
