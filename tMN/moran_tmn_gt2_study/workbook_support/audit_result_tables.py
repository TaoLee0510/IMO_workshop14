"""Independent table integrity checks, without fitting or simulation."""
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def read(file): return list(csv.DictReader((ROOT/file).open()))

def main():
    endpoints=read('data/all_endpoints.csv')
    selected=[r for r in endpoints if r['gene_retained']=='True']
    eligible=[r for r in selected if r['analysis_included']=='True']
    reference=read('results/fit_reference.csv')
    fitted=[r for r in reference if r['s2']!='']
    assert len(reference)==108 and len(fitted)==83 and len(eligible)==110
    assert {(r['patient_id'],r['gene']) for r in selected}=={(r['patient_id'],r['gene']) for r in reference}
    assert {(r['patient_id'],r['gene']) for r in eligible}=={(r['patient_id'],r['gene']) for r in fitted}
    for fit in reference:
        ee=[r for r in eligible if r['patient_id']==fit['patient_id'] and r['gene']==fit['gene']]
        rr=[r for r in selected if r['patient_id']==fit['patient_id'] and r['gene']==fit['gene']]
        assert len(ee)==int(fit['eligible_variants'])
        assert len(rr)==int(fit['total_variants'])
        assert int(fit['tmn_source_flag'])==max(int(r['tmn_source_flag']) for r in rr)
        if fit['s2']!='': assert float(fit['s2'])>=0
        else: assert fit['status']=='not_estimable_from_positive_endpoints'
    predictions=read('results/fit_variant_predictions.csv')
    assert len(predictions)==110
    assert {(r['patient_id'],r['gene'],r['variant_id']) for r in predictions}=={(r['patient_id'],r['gene'],r['variant_id']) for r in eligible}
    for r in predictions:
        assert 0<float(r['predicted_last_vaf_percent'])<50
        pred=2*float(r['predicted_last_vaf_percent'])/100
        observed=2*float(r['last_vaf_percent'])/100
        assert abs((math.log(pred/(1-pred))-math.log(observed/(1-observed)))-float(r['logit_residual']))<1e-10
    summaries=read('results/fit_gene_summary.csv')
    medians={}
    for s in summaries:
        rr=[r for r in fitted if r['gene']==s['gene']]
        medians[s['gene']]=statistics.median(float(r['s2']) for r in rr)
        assert medians[s['gene']]==float(s['median_s2'])
        assert len(rr)==int(s['n_fitted_patient_gene_groups'])
        assert sum(r['boundary']=='True' for r in rr)==int(s['n_zero_boundary'])
        assert sum(r['mismatch_flag']=='True' for r in rr)==int(s['n_mismatch_flag'])
        assert int(s['bootstrap_requested_draws'])==int(s['bootstrap_effective_draws'])+int(s['bootstrap_gene_absent_draws'])
    contrasts=read('results/fit_gene_contrasts.csv')
    assert len(contrasts)==15
    for c in contrasts:
        left,right=c['contrast'].split(' minus ')
        assert float(c['median_difference_s2'])==medians[left]-medians[right]
        assert float(c['bonferroni_familywise_95_low'])<=float(c['cluster_bootstrap_95_low'])<=float(c['cluster_bootstrap_95_high'])<=float(c['bonferroni_familywise_95_high'])
        if int(c['paired_patients'])<2: assert c['paired_bootstrap_95_low']==c['paired_bootstrap_95_high']==''
    allfits=read('results/fit_all_scenarios.csv')
    assert len(allfits)==1616
    assert len(read('results/fit_sensitivity_summary.csv'))==90
    result={'fit_tables_passed':True,'reference_groups':len(reference),'fitted_groups':len(fitted),'fitted_unique_patients':len({r['patient_id'] for r in fitted}),'eligible_variants':len(eligible),'contrasts':len(contrasts),'all_scenario_rows':len(allfits),'risk_tables_checked':False}
    riskfile=ROOT/'results/risk_patient_gene_scenarios.csv'
    if riskfile.exists():
        risk=read('results/risk_patient_gene_scenarios.csv')
        groups=defaultdict(list)
        for r in risk:
            if r['status']!='computed':continue
            p=float(r['P_first_passage']);assert 0<=p<=1
            if float(r['mu'])==0 and float(r['D0'])==0:assert p==0
            key=tuple(r[k] for k in ('patient_id','gene','mu','r_per_day','D0'))
            groups[key].append((float(r['day']),p))
        for series in groups.values():
            ordered=sorted(series)
            assert all(a[1]<=b[1]+1e-12 for a,b in zip(ordered,ordered[1:]))
        assert {(r['patient_id'],r['gene']) for r in risk}=={(r['patient_id'],r['gene']) for r in fitted}
        standard=read('results/risk_standardized_scenarios.csv');assert len(standard)==1440
        result.update(risk_tables_checked=True,patient_risk_rows=len(risk),standardized_risk_rows=len(standard),patient_risk_scenarios=len(groups))
    (ROOT/'workbook_support/result_reconciliation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result))

if __name__=='__main__':main()
