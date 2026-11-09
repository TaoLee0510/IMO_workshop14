import argparse,json,csv,math,subprocess
from common import *
from compute import risk_at,pool_run
def precision(row):
    saved=json.loads((ROOT/'results/patients'/f"{row['patient_id']}_risks.json").read_text())['points'];checks=[]
    for day in ([365,1460] if row['patient_id'] in ['Patient_103','Patient_235','Patient_318'] else [365]):
        current=next(p for p in saved if p['day']==day);refined=risk_at(row,day,10)
        relative=abs(math.expm1((refined['log10_probability']-current['log10_probability'])*math.log(10)))
        assert relative<1e-3,(row['patient_id'],day,relative)
        checks.append({'patient_id':row['patient_id'],'day':day,'refinement_2_vs_10_relative_difference':relative})
    return checks
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--artifacts',action='store_true');args=parser.parse_args();verify_inputs()
    source=Path(MANIFEST['fit_source']);assert sha(source)==MANIFEST['fit_sha256']
    seeds=[];simulations=0
    for row in COHORT:
        pid=row['patient_id'];data=np.load(ROOT/'results/trajectories'/f'{pid}.npz');c=data['counts']
        assert c.shape==(200,1461,3) and np.all(c>=0) and np.all(c.sum(axis=2)==P.N)
        events=list(csv.DictReader((ROOT/'results/trajectories'/f'{pid}.csv').open()));assert len(events)==200
        seeds += [r['seed'] for r in events];simulations+=len(events)
        assert all(int(r['car_C'])==c[j,0,1] and int(r['car_D'])==c[j,0,2] for j,r in enumerate(events))
        if row['baseline_day']==0:assert np.all(c[:,0,1]==initial(row)['initial_C']) and np.all(c[:,0,2]==P.D_init)
        vaf=json.loads((ROOT/'results/patients'/f'{pid}_VAF.json').read_text());assert vaf['saved_loglik_difference']<.01
        assert len(vaf['predictions'])==len(row['observations'])
        for p,o in zip(vaf['predictions'],row['observations']):assert p['day']==o['day'] and p['observed_VAF']==o['vaf'] and 0<p['lower_VAF']<p['upper_VAF']<.5
    assert len(set(seeds))==2800
    target=ROOT/'results/risk_precision.json'
    if not target.exists():write_json(target,[r for checks in pool_run(precision,COHORT,3) for r in checks])
    prefix=ROOT/'qa/native_neutral'
    command=[str(ROOT/'bin/patient_simulator'),str(prefix),'80',str(1/28),'0','0','1',str(A),'1','20','0','3','200','123456','0.5','-10']
    subprocess.run(command,check=True,capture_output=True,text=True)
    native_rows=list(csv.DictReader(prefix.with_suffix('.csv').open()));car=np.array([int(r['car_C']) for r in native_rows])
    assert all(int(r['car_D'])==0 and int(r['final_D'])==0 for r in native_rows)
    assert abs(car.mean()-20)<=5*car.std(ddof=1)/math.sqrt(200)+.05
    report={'status':'passed','patients':14,'followup_VAF_observations':42,'exact_events':2,'censor_records':1,'unknown_outcomes':11,
       'simulations':simulations,'unique_seeds':len(set(seeds)),'pre_CAR_neutral_mean_C':float(car.mean()),'pre_CAR_expected_C':20,
       'fit_and_cohort_hashes_verified':True,'full_population_conservation':True,'risk_precision_checks':len(json.loads(target.read_text())),
       'fixed_Gamma_PDF_integral':json.loads((ROOT/'results/gamma_check.json').read_text())['integral_0_infinity']}
    if args.artifacts:
        from pypdf import PdfReader
        from openpyxl import load_workbook
        pdf=ROOT/'output/pdf/Patient_specific_Moran_analysis.pdf';reader=PdfReader(pdf);assert len(reader.pages)>=20
        assert 'Overfull' not in (pdf.parent/'Patient_specific_Moran_analysis.log').read_text()
        book=load_workbook(ROOT/'outputs/20261105_patient_analysis/Patient_specific_Moran_results.xlsx',read_only=True,data_only=True);assert len(book.sheetnames)==8
        for cell,key in [('D16','VAF_nll'),('E16','clinical_given_VAF_nll'),('F16','joint_nll')]:
            assert math.isclose(book['Likelihoods'][cell].value,FIT[key],rel_tol=1e-12,abs_tol=1e-9),(cell,book['Likelihoods'][cell].value)
        data=[r for r in book['Patient risks'].values if any(v is not None for v in r)]
        header=data[1];rows=[dict(zip(header,r)) for r in data[2:16]]
        source_rows=list(csv.DictReader((ROOT/'results/patient_risks.csv').open()))
        for row,expected in zip(rows,source_rows):
            assert row['patient_id']==expected['patient_id']
            for y in range(1,5):assert math.isclose(row[f'year{y}_probability'],float(expected[f'year{y}_probability']),rel_tol=1e-12)
        for sh in book:
            for row in sh.values:assert not any(isinstance(v,str) and v in ['#REF!','#VALUE!','#DIV/0!','#NAME?','#NUM!'] for v in row)
        report.update(pdf_pages=len(reader.pages),workbook_sheets=book.sheetnames,workbook_risks_reconciled=True);book.close()
    report['source_code_sha256']={p.name:sha(p) for p in ROOT.iterdir() if p.suffix in ['.py','.cpp','.mjs','.tex','.json']}
    write_json(ROOT/'results'/('delivery_validation.json' if args.artifacts else 'validation.json'),report)
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
