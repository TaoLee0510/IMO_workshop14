"""Prepare typed, source-preserving tables for the Artifact Tool workbook."""
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'workbook_support'

def table(name, file, first=()):
    path=ROOT/file
    records=list(csv.DictReader(path.open()))
    columns=list(records[0])
    columns=[k for k in first if k in columns]+[k for k in columns if k not in first]
    kinds={}
    for key in columns:
        values=[r[key] for r in records if r[key]!='']
        kind='text'
        if key not in ('first_vaf_source_value','last_vaf_source_value','tmn_source_excel_rows'):
            if values and all(x in ('True','False') for x in values): kind='boolean'
            elif values:
                try:
                    if all(math.isfinite(float(x)) for x in values): kind='number'
                except ValueError: pass
        kinds[key]=kind
    def typed(value,key):
        if value=='': return None
        if kinds[key]=='number': return float(value)
        if kinds[key]=='boolean': return value=='True'
        return value
    return {'name':name,'source':file,'columns':columns,'kinds':kinds,'rows':[[typed(r[k],k) for k in columns] for r in records]}

def main():
    specs=[
      ('Summary','results/fit_gene_summary.csv',('gene','n_fitted_patient_gene_groups','median_s2','median_cluster_bootstrap_95_low','median_cluster_bootstrap_95_high','n_zero_boundary','n_mismatch_flag','exploratory_small_n')),
      ('Patient Gene Fits','results/fit_reference.csv',('patient_id','gene','s2','eligible_variants','total_variants','status','boundary','mismatch_flag','tmn_source_flag')),
      ('Gene Contrasts','results/fit_gene_contrasts.csv',()),
      ('Sensitivity Summary','results/fit_sensitivity_summary.csv',('scenario','gene','n_fitted_patient_gene_groups','median_s2','n_zero_boundary','n_mismatch_flag','median_cluster_bootstrap_95_low','median_cluster_bootstrap_95_high')),
      ('Sensitivity Fits','results/fit_all_scenarios.csv',('scenario','patient_id','gene','s2','eligible_variants','status','boundary','mismatch_flag','tmn_source_flag')),
      ('Patient Risk Scenarios','results/risk_patient_gene_scenarios.csv',('patient_id','gene','mu','r_per_day','D0','day','P_first_passage','status','s2','founder_count')),
      ('Standardized Risk','results/risk_standardized_scenarios.csv',('gene','mu','r_per_day','D0','day','P_first_passage','status','s2','founder_count')),
      ('Variant Predictions','results/fit_variant_predictions.csv',('patient_id','gene','variant_id','first_day','last_day','first_vaf_percent','last_vaf_percent','predicted_last_vaf_percent','logit_residual')),
      ('Gene Eligibility','data/gene_eligibility.csv',()),
      ('Endpoints','data/all_endpoints.csv',()),
    ]
    tables=[table(*spec) for spec in specs]
    outcomes=json.loads((ROOT/'data/outcome_evidence.json').read_text())
    columns=['patient_id','genes','tmn_source_flag','label_interpretation','event_free_followup_available','clinical_record_present','clinical_source_excel_row','clinical_cohort','time_cell_therapy_to_tMN_days']
    rows=[]
    for o in outcomes:
        for c in o['clinical_records'] or [None]:
            rows.append([o['patient_id'],';'.join(o['genes']),o['tmn_source_flag'],o['label_interpretation'],o['event_free_followup_available'],bool(c),c['source_excel_row'] if c else None,c['cohort'] if c else None,c['time_cell_therapy_to_tMN'] if c else None])
    tables.append({'name':'Outcome Evidence','source':'data/outcome_evidence.json','columns':columns,'rows':rows,'kinds':{k:('number' if k in ('tmn_source_flag','clinical_source_excel_row','time_cell_therapy_to_tMN_days') else 'boolean' if k in ('event_free_followup_available','clinical_record_present') else 'text') for k in columns}})
    config=json.loads((ROOT/'approved_scenarios.json').read_text())
    settings=[['N',config['N'],'cells','Fixed progenitor pool'],['first passage threshold',config['threshold'],'double-hit cells','Strictly above 1% of N'],['s2 lower bound',0,'relative fitness','All patient/gene estimates are nonnegative']]
    units={'b':'per day','tau':'days','Tmax':'days after CAR T','rise_doubling_days':'days','decay_half_life_days':'days','s1':'relative fitness'}
    for key,value in config['reference'].items(): settings.append([key,value,units[key],'Approved reference assumption'])
    for key,values in config['one_at_a_time_growth_sensitivity'].items():
        for value in values:settings.append([key,value,units[key],'One-at-a-time refit sensitivity'])
    for key,values in config['risk_scenarios'].items():
        for value in values:settings.append([key,value,'per division' if key.startswith('mu') else 'per day' if 'growth' in key else 'cells','Exploratory risk scenario'])
    tables.append({'name':'Model Settings','source':'approved_scenarios.json','columns':['parameter','value','units','scenario_role'],'rows':settings,'kinds':{'parameter':'text','value':'number','units':'text','scenario_role':'text'}})
    payload={'tables':tables,'audit':json.loads((ROOT/'data/audit.json').read_text()),'fit_validation':json.loads((ROOT/'results/fit_validation.json').read_text()),'config':config}
    (OUT/'workbook_payload.json').write_text(json.dumps(payload,allow_nan=False,separators=(',',':')))
    print(json.dumps({t['name']:len(t['rows']) for t in tables}))

if __name__=='__main__': main()
