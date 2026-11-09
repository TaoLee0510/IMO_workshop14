"""Prepare the confirmed 14-patient source-record joint cohort without invented negatives."""
import csv
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def prepare():
    config=json.loads((ROOT/'config.json').read_text())
    bundled=ROOT/'data/source_inputs'
    previous=bundled if (bundled/'analysis_longitudinal.csv').exists() else ROOT.parent/'moran_tmn_gamma_4y_study/data'
    workbooks=bundled if (bundled/'41588_2026_2795_MOESM7_ESM.xlsx').exists() else ROOT.parent
    rows=list(csv.DictReader((previous/'analysis_longitudinal.csv').open()))
    inventory={r['patient_id']:r for r in csv.DictReader((previous/'index_marker_inventory.csv').open())}
    cohort=[]
    ids=sorted(set(r['patient_id'] for r in rows)|{'Patient_103','Patient_25'},key=lambda x:int(x.split('_')[1]))
    for pid in ids:
        rec=inventory[pid]
        observations=[{'day':float(r['day']),'vaf':float(r['vaf_percent'])/100,
                       'source_excel_row':int(r['source_excel_row'])}
                      for r in rows if r['patient_id']==pid and r['role']=='likelihood_observation']
        outcome=config['clinical'].get(pid,{'kind':'unknown'}).copy()
        cohort.append({'patient_id':pid,'gene':rec['gene'],'variant_id':rec['variant_id'],
                       'sample_type':rec['sample_type'],'baseline_day':float(rec['baseline_day']),
                       'baseline_vaf':float(rec['baseline_vaf_percent'])/100,
                       'observations':observations,'outcome':outcome,
                       'clinical_evidence':'source_record_assumption_confirmed' if outcome['kind']!='unknown' else 'unknown_not_negative'})
    assert len(cohort)==14 and sum(len(r['observations']) for r in cohort)==42
    assert sum(r['outcome']['kind']=='exact' for r in cohort)==2
    assert sum(r['outcome']['kind']=='censor' for r in cohort)==1
    (ROOT/'data').mkdir(exist_ok=True)
    (ROOT/'data/cohort.json').write_text(json.dumps(cohort,indent=2))
    fields=['patient_id','gene','variant_id','sample_type','baseline_day','baseline_vaf','followup_observations','outcome','outcome_day','clinical_evidence']
    with (ROOT/'data/cohort.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in cohort:
            w.writerow({**{k:r[k] for k in fields if k in r},'followup_observations':len(r['observations']),
                        'outcome':r['outcome']['kind'],'outcome_day':r['outcome'].get('day','')})
    inputs=[workbooks/'41588_2026_2795_MOESM7_ESM.xlsx',workbooks/'41588_2026_2795_MOESM9_ESM.xlsx',
            previous/'analysis_longitudinal.csv',previous/'index_marker_inventory.csv']
    audit={'patients':14,'VAF_followup_observations':42,'exact_events':2,'four_year_censor_records':1,'unknown_outcomes':11,
           'clinical_policy':'Use the user-confirmed source times, not blanket manual verification flags.',
           'source_assumptions':['Listed tMN time is first diagnosis from the relevant index CAR infusion.','Patient_25 first diagnosis at day 1561 supports survival through day 1460.'],
           'excluded_clinical_cases':{'Patient_127':'no eligible pre-CAR blood-cell initial state','Patient_18':'unresolved CAR time origin and no eligible initial state',
                                      'Patient_69':'unresolved multiple CAR origins','Patient_275':'diagnosis time missing; postdiagnosis VAF cannot be identified'},
           'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}}
    (ROOT/'data/audit.json').write_text(json.dumps(audit,indent=2))
    return cohort


if __name__=='__main__':
    print(json.dumps(prepare(),indent=2))
