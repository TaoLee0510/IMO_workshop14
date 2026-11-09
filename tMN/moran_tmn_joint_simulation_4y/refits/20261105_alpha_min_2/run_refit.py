"""Reference shared-state joint refit, preserving the preceding study outputs."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import csv
import hashlib
import json
import math
import platform
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import scipy
from fitting import fit_scenario

ROOT=Path(__file__).resolve().parent


def write_rows(path,rows):
    keys=list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def main():
    config=json.loads((ROOT/'config.json').read_text())
    cohort=json.loads((ROOT/'data/cohort.json').read_text())
    old=json.loads((ROOT/'data/previous_reference.json').read_text())
    assert config['reference']==old['parameters']
    assert config['fit_bounds']=={'alpha':[2.0,50.0],'theta_days':[0.25,1460.0]}
    started=datetime.now(ZoneInfo('America/New_York')).isoformat()
    print('START: alpha=[2,50], theta=[0.25,1460] days; reference fixed parameters unchanged; full shared-state joint likelihood.',flush=True)
    fit=fit_scenario(('reference',config['reference'],cohort,{'maxiter':150}))
    evaluations=json.loads((ROOT/'results/fits/reference_evaluations.json').read_text())
    assert all(2<=x['alpha']<=50+1e-12 and .25-1e-12<=x['theta']<=1460+1e-9 for x in evaluations)
    assert 2<=fit['alpha']<=50 and .25<=fit['theta_days']<=1460
    assert fit['status']=='fitted' and fit['optimizer_success']
    assert fit['clinical_events_used']==2 and fit['censor_records_used']==1 and fit['unknown_outcomes']==11
    assert math.isclose(fit['joint_nll'],fit['VAF_nll']+fit['clinical_given_VAF_nll'],abs_tol=1e-7)
    assert not fit['numerical_errors'],fit['numerical_errors']
    write_rows(ROOT/'results/joint_fits.csv',[{k:v for k,v in fit.items() if not isinstance(v,(dict,list))}|fit['parameters']])
    write_rows(ROOT/'results/patient_likelihoods.csv',[{k:v for k,v in r.items() if k!='VAF_predictions'} for r in fit['patient_details']])
    write_rows(ROOT/'results/VAF_predictions.csv',[{'patient_id':r['patient_id'],'gene':r['gene'],**p} for r in fit['patient_details'] for p in r['VAF_predictions']])
    comparison=[]
    for name,f,lower in [('previous_alpha_min_1.01',old,1.01),('new_alpha_min_2',fit,2.0)]:
        comparison.append({'fit':name,'alpha_lower_bound':lower,**{k:f[k] for k in ['alpha','theta_days','peak_day','peak_additional_fitness','joint_nll','VAF_nll','clinical_given_VAF_nll']}})
    write_rows(ROOT/'results/comparison.csv',comparison)
    summary={'started_at':started,'finished_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
             'fit_bounds':config['fit_bounds'],'status':fit['status'],'alpha':fit['alpha'],'theta_days':fit['theta_days'],
             'joint_nll':fit['joint_nll'],'delta_nll_vs_previous':fit['joint_nll']-old['joint_nll'],
             'all_evaluated_points_within_bounds':True,'evaluated_points':len(evaluations),'fixed_parameters_unchanged':True,
             'clinical_events_used':fit['clinical_events_used'],'censor_records_used':fit['censor_records_used'],
             'cap_refinement_loglik_difference':fit['cap_refinement_loglik_difference'],
             'tolerance_refinement_loglik_difference':fit['tolerance_refinement_loglik_difference'],
             'full_N_loglik_difference':fit['full_N_loglik_difference'],
             'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),
             'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir() if p.suffix in ['.py','.json']}}
    (ROOT/'results/refit_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':main()
