"""Archive the reviewed outputs, scientific inputs, code and raw trajectories."""
import hashlib,json,zipfile
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
ROOT=Path(__file__).resolve().parent
def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while chunk:=f.read(1024*1024):h.update(chunk)
    return h.hexdigest()
def main():
    validation=json.loads((ROOT/'results/validation/delivery.json').read_text());assert validation['status']=='passed'
    visual=json.loads((ROOT/'results/validation/visual_review.json').read_text());assert visual['status']=='passed'
    selected=set()
    for pattern in ['*.py','*.mjs','*.tex','README.md','requirements.txt','config.json','provenance.json',
        'data/*','results/*.csv','results/fits/*.json','results/risks/*.json','results/comparison/*.json','results/validation/*.json',
        'figures/**/*.pdf','figures/**/*.png','patient_analysis/*.py','patient_analysis/*.cpp','patient_analysis/*.json',
        'patient_analysis/README.md','patient_analysis/requirements.txt','patient_analysis/inputs/*.json',
        'patient_analysis/figures/**/*.pdf','patient_analysis/figures/**/*.png','patient_analysis/results/*.csv',
        'patient_analysis/results/*.json','patient_analysis/results/patients/*.json','patient_analysis/results/trajectories/*.npz',
        'patient_analysis/results/trajectories/*.csv','patient_analysis/results/trajectories/*.json']:
        selected.update(p for p in ROOT.glob(pattern) if p.is_file())
    selected={p for p in selected if p.name not in ['native_kernel.py','native_derivative.json']}
    pdf=ROOT/'output/pdf/Mu001_joint_sensitivity_report.pdf'
    workbook=ROOT/'outputs/20261105_mu001_sensitivity/Mu001_joint_sensitivity_results.xlsx'
    selected.update([pdf,workbook])
    files=sorted(selected);prefix=ROOT.name+'/'
    manifest={'created_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
        'study':'Four-year joint Moran CAR-T/tMN fit; reference mu=0.01',
        'validation':validation,'files':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':digest(p)} for p in files}}
    output=ROOT/'outputs/20261105_mu001_sensitivity/TMN_mu001_complete_study.zip';temp=output.with_suffix('.zip.tmp')
    with zipfile.ZipFile(temp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as archive:
        for p in files:archive.write(p,prefix+str(p.relative_to(ROOT)))
        archive.writestr(prefix+'DELIVERY_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
    with zipfile.ZipFile(temp) as archive:
        assert archive.testzip() is None
        for name,record in manifest['files'].items():
            assert hashlib.sha256(archive.read(prefix+name)).hexdigest()==record['sha256'],name
    temp.replace(output)
    status={'status':'complete','completed_at':manifest['created_at'],'fitted_scenarios':78,'patient_risk_records':1092,
        'PDF_pages':validation['PDF_pages'],'workbook_sheets':8,'figure_pairs':82,'raw_Moran_paths':2800,
        'archive_files':len(files)+1,'archive_MiB':output.stat().st_size/1024**2,
        'outputs':{str(p.relative_to(ROOT)):{'bytes':p.stat().st_size,'sha256':digest(p)} for p in [pdf,workbook,output]}}
    (ROOT/'results/delivery_status.json').write_text(json.dumps(status,indent=2)+'\n');print(json.dumps(status,indent=2))
if __name__=='__main__':main()
