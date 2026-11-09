"""Compare the inference approximation with 200 independent full Moran paths.

One patient is evaluated per fresh process. This validation does not refit the
model or alter the frozen trajectories. Monte Carlo intervals describe simulation
sampling error at fixed parameters, not uncertainty in the fitted parameters.
"""
import argparse,csv,hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def records(path):return list(csv.DictReader(path.open()))
def patient(pid):
    import numpy as np
    from scipy.stats import beta
    from model import Parameters,Backgrounds,ForwardFilter
    import memory_control
    fit=json.loads((ROOT/'results/fits/reference.json').read_text());p=Parameters(**fit['parameters'])
    row=next(r for r in json.loads((ROOT/'data/cohort.json').read_text()) if r['patient_id']==pid)
    signature=hashlib.sha256(Path(__file__).read_bytes()+(ROOT/'model.py').read_bytes()+(ROOT/'results/fits/reference.json').read_bytes()).hexdigest()
    target=ROOT/'results/comparison'/f'{pid}.json';target.parent.mkdir(exist_ok=True)
    if target.exists() and json.loads(target.read_text()).get('signature')==signature:return
    a,t=fit['alpha'],fit['theta_days'];bg=Backgrounds([row],a,t,p,2);qfun=bg.for_patient(0)
    filt=ForwardFilter(p,a,t,qfun,killed=False,cap=p.N,window=True,refinement=2)
    if row['baseline_day']<0:
        pre=ForwardFilter(p,a,t,qfun,killed=True,cap=p.K-1,refinement=2)
        pre.propagate(row['baseline_day'],0.);pre.condition_event_free()
        filt.state[:p.K]=pre.state[:p.K]*math.exp(pre.logscale)
    counts=np.load(ROOT/'patient_analysis/results/trajectories'/f'{pid}.npz')['counts']
    hits=records(ROOT/'patient_analysis/results/trajectories'/f'{pid}.csv')
    risks=next(r for r in records(ROOT/'results/patient_risks.csv') if r['scenario']=='reference' and r['patient_id']==pid)
    output=[];start=0.
    for year in range(1,5):
        day=year*365;filt.propagate(start,float(day));start=float(day)
        mass_error=abs(filt.physical_logmass());assert mass_error<2e-5,(pid,day,mass_error)
        weights=filt.state/filt.state.sum();mean_d=float(weights@filt.d);mean_c=qfun(day)*(p.N-mean_d)
        vals=counts[:,day,:];mc=vals.mean(axis=0);se=vals.std(axis=0,ddof=1)/math.sqrt(200)
        h=sum(bool(e['first_hit_day']) and float(e['first_hit_day'])<=day for e in hits)
        lower=float(beta.ppf(.025,h,201-h)) if h else 0.
        upper=float(beta.ppf(.975,h+1,200-h)) if h<200 else 1.
        predicted=float(risks[f'year{year}_probability'])
        output.append({'patient_id':pid,'day':day,'hybrid_mean_C_cells':mean_c,'hybrid_mean_D_cells':mean_d,
            'MC_mean_C_cells':float(mc[1]),'MC_mean_D_cells':float(mc[2]),'MC_SE_C_cells':float(se[1]),'MC_SE_D_cells':float(se[2]),
            'C_mean_error_percentage_points':float((mean_c-mc[1])/p.N*100),
            'D_mean_error_percentage_points':float((mean_d-mc[2])/p.N*100),
            'hybrid_tMN_probability':predicted,'MC_tMN_probability':h/200,'MC_events':int(h),
            'MC_probability_lower_95':lower,'MC_probability_upper_95':upper,
            'hybrid_probability_in_MC_interval':lower-1e-10<=predicted<=upper+1e-10,
            'forward_mass_log_error':mass_error,'window_edge_probability':filt.window_edge_probability})
    target.write_text(json.dumps({'signature':signature,'rows':output},indent=2)+'\n')
    print('COMPARE',pid,flush=True)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--patient');args=parser.parse_args()
    if args.patient:patient(args.patient);return
    cohort=json.loads((ROOT/'data/cohort.json').read_text())
    for row in cohort:subprocess.run([sys.executable,str(Path(__file__).resolve()),'--patient',row['patient_id']],check=True)
    rows=[r for p in cohort for r in json.loads((ROOT/'results/comparison'/f"{p['patient_id']}.json").read_text())['rows']]
    with (ROOT/'results/Moran_comparison.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    summary={'status':'complete','patient_year_comparisons':len(rows),'replicates_per_patient':200,
        'max_C_mean_error_percentage_points':max(abs(r['C_mean_error_percentage_points']) for r in rows),
        'max_D_mean_error_percentage_points':max(abs(r['D_mean_error_percentage_points']) for r in rows),
        'max_absolute_probability_difference':max(abs(r['hybrid_tMN_probability']-r['MC_tMN_probability']) for r in rows),
        'hybrid_probabilities_in_MC_95_intervals':sum(r['hybrid_probability_in_MC_interval'] for r in rows),
        'interval_definition':'Pointwise exact binomial 95% interval for 200 full-process simulations; not fitted-parameter uncertainty',
        'max_forward_mass_log_error':max(r['forward_mass_log_error'] for r in rows)}
    (ROOT/'results/validation/Moran_comparison.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
