"""Reconcile completed scenarios, source signatures, risk tables and artifacts."""
import csv,hashlib,json,math
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from openpyxl import load_workbook
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parent
def read(path):return json.loads((ROOT/path).read_text())
def rows(path):return list(csv.DictReader((ROOT/path).open()))
def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def close(a,b):return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=0)
def sheet_rows(book,name):
    vals=list(book[name].values);return [dict(zip(vals[0],v)) for v in vals[1:] if any(x is not None for x in v)]
def main():
    cohort=read('data/cohort.json');config=read('config.json')
    fits=[json.loads(p.read_text()) for p in (ROOT/'results/fits').glob('*.json') if not p.stem.endswith('_evaluations')]
    good={r['scenario']:r for r in fits if r['status']=='fitted'}
    assert len(fits)==79 and len(good)==78
    incompatible=[f for f in fits if f['status']=='structurally_incompatible']
    assert len(incompatible)==1 and incompatible[0]['parameters']['mu']==0 and incompatible[0]['parameters']['D_init']==0
    assert sum(r['outcome']['kind']=='exact' for r in cohort)==2
    assert sum(r['outcome']['kind']=='censor' for r in cohort)==1
    assert sum(r['outcome']['kind']=='unknown' for r in cohort)==11
    cohort_hash=hashlib.sha256(json.dumps(cohort,sort_keys=True).encode()).hexdigest()
    for fit in good.values():
        assert fit['fit_signature']['model_sha256']==sha('model.py')
        assert fit['fit_signature']['cohort_sha256']==cohort_hash
        assert fit['fit_signature']['fit_bounds']==config['fit_bounds']
        assert 2<=fit['alpha']<=50+1e-12 and 1<=fit['theta_days']<=1460+1e-10
        assert fit['parameters']['area_days']==1 and fit['parameters']['N']==100000
        assert fit['cap_refinement_loglik_difference']<.005
        assert fit['tolerance_refinement_loglik_difference']<.005
        assert not fit['numerical_errors']
        assert abs(fit['joint_nll']-fit['VAF_nll']-fit['clinical_given_VAF_nll'])<1e-8
        assert all(2*100000*r['baseline_vaf']>=fit['parameters']['D_init'] for r in cohort)
    assert good['reference']['full_N_loglik_difference']<.005
    risks=rows('results/patient_risks.csv');assert len(risks)==1092
    lookup={(r['scenario'],r['patient_id']):r for r in risks};assert len(lookup)==1092
    for scenario,fit in good.items():
        saved=read(f'results/risks/{scenario}.json');s=saved['signature']
        assert s['parameters']==fit['parameters'] and s['alpha']==fit['alpha'] and s['theta_days']==fit['theta_days']
        assert s['model_sha256']==sha('model.py') and s['risk_source_sha256']==sha('risk_backend.py')
        assert len(saved['rows'])==14
        for patient in cohort:
            row=lookup[scenario,patient['patient_id']]
            probabilities=[float(row[f'year{y}_probability']) for y in range(1,5)]
            assert all(0<=p<=1 for p in probabilities)
            assert all(a<=b+1e-10 for a,b in zip(probabilities,probabilities[1:]))
    child=read('patient_analysis/results/validation.json');assert child['status']=='passed'
    assert child['simulations']==2800 and child['unique_seeds']==2800 and child['risk_precision_checks']==17
    assert abs(read('patient_analysis/results/gamma_check.json')['integral_0_infinity']-1)<1e-9
    comparison=read('results/validation/Moran_comparison.json');assert comparison['patient_year_comparisons']==56
    assert comparison['hybrid_probabilities_in_MC_95_intervals']==56
    assert len(list((ROOT/'figures').rglob('*.pdf')))==18
    assert len(list((ROOT/'patient_analysis/figures').rglob('*.pdf')))==64
    assert len(list((ROOT/'figures').rglob('*.png')))==18
    assert len(list((ROOT/'patient_analysis/figures').rglob('*.png')))==64
    workbook=ROOT/'outputs/20261105_mu001_sensitivity/Mu001_joint_sensitivity_results.xlsx'
    book=load_workbook(workbook,read_only=True,data_only=True)
    assert len(book.sheetnames)==8
    assert book['Parameters']['B12'].value==1 and book['Parameters']['B13'].value==1001 and book['Parameters']['B14'].value==0
    for cell,value in {'B5':good['reference']['alpha'],'B6':1,'B7':.01,'B8':0,'B9':good['reference']['joint_nll'],
       'B10':good['reference']['VAF_nll'],'B11':good['reference']['clinical_given_VAF_nll'],'B12':79,'B13':78,'B14':1,'B15':42,'B17':22/42,'B18':78}.items():
        actual=book['Summary'][cell].value
        assert actual is not None and math.isclose(actual,value,rel_tol=1e-12,abs_tol=1e-10),(cell,actual,value)
    saved_risks=sheet_rows(book,'Patient risks');assert len(saved_risks)==1092
    for row in saved_risks:
        expected=lookup[row['scenario'],row['patient_id']]
        for year in range(1,5):assert close(row[f'year{year}_probability'],expected[f'year{year}_probability'])
    jf=sheet_rows(book,'Joint fits');assert len(jf)==79
    for row in jf:
        if row['status']=='fitted':
            expected=good[row['scenario']]
            for key in ['alpha','theta_days','joint_nll','VAF_nll','clinical_given_VAF_nll']:assert close(row[key],expected[key])
            assert row['boundary_optimum']=='At bound'
        else:assert row['alpha'] is None and row['theta_days'] is None
    for name,count in [('Reference VAF',42),('Reference sims',14),('Reference yearly',70),('Initial states',112)]:
        assert len(sheet_rows(book,name))==count
    errors={'#REF!','#VALUE!','#DIV/0!','#NAME?','#N/A','#NUM!','#NULL!','#SPILL!','#CALC!'}
    for sh in book:
        for row in sh.values:assert not any(isinstance(v,str) and v in errors for v in row)
    book.close()
    formulas=load_workbook(workbook,read_only=True,data_only=False)
    assert formulas['Summary']['B5'].value=="='Joint fits'!E2"
    assert formulas['Summary']['B18'].value=='=COUNTIF(\'Joint fits\'!M2:M80,"At bound")'
    formulas.close()
    pdf=ROOT/'output/pdf/Mu001_joint_sensitivity_report.pdf';reader=PdfReader(pdf)
    assert 30<=len(reader.pages)<=45
    assert sum(int(page.get('/Rotate',0))%180==90 for page in reader.pages)==2
    text=' '.join(page.extract_text() or '' for page in reader.pages)
    assert '455.808076' in text and '52.4' in text
    log=(pdf.with_suffix('.log')).read_text()
    assert 'Overfull' not in log and 'Undefined control sequence' not in log
    assert len(list((ROOT/'qa/pdf_pages').glob('page-*.png')))==len(reader.pages)
    report={'status':'passed','validated_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
        'fitted_scenarios':78,'structurally_incompatible_scenarios':1,'patient_risk_records':1092,'yearly_probability_values':4368,
        'reference_trajectories':2800,'figures_PDF_and_PNG_pairs':82,'workbook_sheets':8,'PDF_pages':len(reader.pages),
        'fit_source_signatures_verified':True,'risk_source_signatures_verified':True,'workbook_probabilities_and_NLLs_reconciled':True,
        'fixed_Gamma_integral':1,'reference_model_and_full_Moran_comparison':comparison,
        'outputs_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [pdf,workbook]}}
    (ROOT/'results/validation/delivery.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
