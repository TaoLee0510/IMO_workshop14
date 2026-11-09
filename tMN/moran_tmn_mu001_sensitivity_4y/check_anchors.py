import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from model import Parameters,Backgrounds,patient_loglik
ROOT=Path(__file__).resolve().parent
COHORT=json.loads((ROOT/'data/cohort.json').read_text())
P=Parameters(**json.loads((ROOT/'config.json').read_text())['reference'])
def check(pair):
    a,t=pair;bg=Backgrounds(COHORT,a,t,P);out={'alpha':a,'theta_days':t,'patients':[]}
    for i,row in enumerate(COHORT):
        try:
            start=time.time();ll,d=patient_loglik(row,i,bg,a,t,P,512,1)
            out['patients'].append({'patient':row['patient_id'],'ll':ll,'seconds':time.time()-start})
        except Exception as e:out['error_patient']=row['patient_id'];out['error']=str(e);break
    print(json.dumps({k:v for k,v in out.items() if k!='patients'}),flush=True);return out
def main():
    pairs=[(2.,1.),(8.,1.),(50.,1.),(8.,21.),(50.,1460.)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        out=[f.result() for f in as_completed([pool.submit(check,p) for p in pairs])]
    (ROOT/'results/validation/anchor_checks.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':main()
