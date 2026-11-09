"""Final provenance, numerical, trajectory and artifact reconciliation."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from model import Parameters,cohort_loglik

ROOT=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--artifacts',action='store_true');args=parser.parse_args()
    config=json.loads((ROOT/'config.json').read_text());cohort=json.loads((ROOT/'data/cohort.json').read_text())
    fits=[json.loads(p.read_text()) for p in sorted((ROOT/'results/fits').glob('*.json')) if not p.stem.endswith('_evaluations')]
    assert len(fits)==21
    compatible=[r for r in fits if r['status']=='fitted'];assert len(compatible)==20
    assert all(r.get('optimizer_version')==2 for r in fits)
    assert all(r['optimizer_success'] for r in compatible)
    assert next(r for r in fits if r['scenario']=='sigma_vaf_0.1')['precision_version']==2
    consistency=[]
    for fit in compatible:
        data=[dict(r,outcome={'kind':'unknown'}) for r in cohort] if fit['scenario']=='vaf_only_comparator' else cohort
        ll,_=cohort_loglik(data,fit['alpha'],fit['theta_days'],Parameters(**fit['parameters']),fit['cap'],1)
        delta=abs(-ll-fit['joint_nll'])
        assert delta<.03,(fit['scenario'],delta)
        assert abs(fit['joint_nll']-fit['VAF_nll']-fit['clinical_given_VAF_nll'])<1e-7
        consistency.append({'scenario':fit['scenario'],'saved_vs_current_NLL_difference':delta})
    metadata=json.loads((ROOT/'results/simulation_metadata.json').read_text());assert len(metadata)==58
    assert sum(r['replicates'] for r in metadata)==11600
    assert sum(r['replicates'] for r in metadata if r['label'].startswith('grid_'))==3600
    cppsha=hashlib.sha256((ROOT/'simulator.cpp').read_bytes()).hexdigest();seeds=[]
    for rec in metadata:
        assert rec['replicates']==200 and rec['simulator_source_sha256']==cppsha
        data=np.load(ROOT/'results/trajectories'/f"{rec['name']}.npz")
        counts=data['counts'];assert counts.shape==(200,1461,2)
        n=rec['parameters']['N'];d0=round(rec['D_initial_fraction']*n)
        assert np.all(counts>=0) and np.all(counts.sum(axis=2)<=n)
        assert np.all(counts[:,0,0]==2000) and np.all(counts[:,0,1]==d0)
        assert abs(counts[:,-1,0].mean()/n-rec['mean_C_four_year_fraction'])<1e-12
        assert abs(counts[:,-1,1].mean()/n-rec['mean_D_four_year_fraction'])<1e-12
        rows=list(csv.DictReader((ROOT/'results/trajectories'/f"{rec['name']}.csv").open()))
        assert len(rows)==200;seeds.extend(r['seed'] for r in rows)
        K=math.floor(rec['parameters']['threshold_fraction']*n)+1
        if d0>=K:assert all(float(r['first_hit_day'])==0 for r in rows)
    assert len(set(seeds))==len(seeds)
    risks=list(csv.DictReader((ROOT/'results/four_year_risks.csv').open()));assert len(risks)==320
    for row in risks:
        values=[float(row[f'year_{year}_probability']) for year in range(1,5)]
        assert all(0<=v<=1 for v in values) and all(a<=b+1e-10 for a,b in zip(values,values[1:])),row
    profile_ref=json.loads((ROOT/'results/profile_reference.json').read_text())
    reference=next(r for r in fits if r['scenario']=='reference')
    assert profile_ref['alpha']==reference['alpha'] and profile_ref['theta_days']==reference['theta_days']
    profile=list(csv.DictReader((ROOT/'results/likelihood_profiles.csv').open()))
    assert min(float(r['delta_nll_from_joint_optimum']) for r in profile)>-.03
    validation=json.loads((ROOT/'results/validation.json').read_text());assert validation['status']=='passed'
    artifact_info={}
    if args.artifacts:
        from pypdf import PdfReader
        from openpyxl import load_workbook
        pdf=ROOT/'main.pdf';reader=PdfReader(pdf);assert len(reader.pages)>=10
        log=(ROOT/'main.log').read_text();assert 'Overfull' not in log
        xlsx=ROOT/'outputs/20261105_fixed_area/Fixed_area_joint_Moran_study.xlsx'
        book=load_workbook(xlsx,read_only=True,data_only=True)
        def populated_rows(sheet):
            return [row for row in book[sheet].values if any(value is not None for value in row)]
        assert len(populated_rows('Joint fits'))==22 and len(populated_rows('Four-year risks'))==321
        assert len(populated_rows('Simulations'))==59 and len(populated_rows('Patients'))==15
        mu_headers=list(book['Joint fits'].values)[0];mu_index=mu_headers.index('mu')
        assert any(r[mu_index]==1e-7 for r in list(book['Joint fits'].values)[1:])
        def records(sheet):
            values=iter(book[sheet].values);headers=next(values)
            return [dict(zip(headers,row)) for row in values]
        saved_fits={r['scenario']:r for r in records('Joint fits')}
        for r in fits:
            row=saved_fits[r['scenario']]
            assert row['status']==r['status']
            for key in ['alpha','theta_days','joint_nll','VAF_nll','clinical_given_VAF_nll']:
                expected=r.get(key)
                assert row[key] is None if expected is None else math.isclose(row[key],expected,rel_tol=1e-12,abs_tol=0)
        def risk_key(r):return (r['scenario'],r['population'],r.get('patient_id') or '',float(r.get('initial_D_fraction') or 0))
        source_risks={risk_key(r):r for r in risks}
        for row in records('Four-year risks'):
            source=source_risks[risk_key(row)]
            for year in range(1,5):
                key=f'year_{year}_probability';assert math.isclose(row[key],float(source[key]),rel_tol=1e-12,abs_tol=0)
        patients={r['patient_id']:r for r in records('Patients')}
        assert patients['Patient_25']['outcome_day']==1460 and patients['Patient_25']['listed_first_event_day']==1561
        assert patients['Patient_103']['outcome']=='exact' and patients['Patient_197']['outcome']=='exact'
        assert sum(r['outcome']=='unknown' for r in patients.values())==11
        for sheet in book:
            for row in sheet.values:
                assert not any(isinstance(v,str) and v in ['#REF!','#VALUE!','#DIV/0!','#NAME?','#NUM!'] for v in row)
        artifact_info={'pdf_pages':len(reader.pages),'workbook_sheets':book.sheetnames,'xlsx_bytes':xlsx.stat().st_size,
                       'xlsx_numeric_values_reconciled_to_saved_fits_and_risks':True}
        book.close()
    report={'status':'passed','checked_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
            'fit_scenarios':21,'compatible_joint_or_comparator_fits':20,'structurally_incompatible_scenarios':1,
            'primary_simulations':3600,'all_reported_simulations':11600,'unique_replicate_seeds':len(seeds),
            'risk_rows':320,'risk_units':'probability on [0,1]','saved_fit_reconciliation':consistency,
            'source_code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir() if p.suffix in ['.py','.cpp','.json','.mjs','.tex']},
            **artifact_info}
    target=ROOT/'results'/('delivery_audit.json' if args.artifacts else 'computation_audit.json')
    target.write_text(json.dumps(report,indent=2))
    print(json.dumps({k:v for k,v in report.items() if k not in ['source_code_sha256','saved_fit_reconciliation']},indent=2))


if __name__=='__main__':main()
