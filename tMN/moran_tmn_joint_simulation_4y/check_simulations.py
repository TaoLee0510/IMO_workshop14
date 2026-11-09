"""Compare exact Monte Carlo endpoint moments with the inference approximation."""
import csv
import json
import math
from pathlib import Path
from dataclasses import replace
import numpy as np
from scipy.stats import beta
from model import Parameters,Backgrounds,ForwardFilter,risk_curve

ROOT=Path(__file__).resolve().parent


def check():
    config=json.loads((ROOT/'config.json').read_text());reference=Parameters(**config['reference'])
    rows=[]
    for alpha,theta in [(2.,2.),(4.,7.),(8.,21.)]:
        for initial_D in (0,1000):
            label=f'grid_a{alpha:g}_t{theta:g}_d{initial_D}'
            data=np.load(ROOT/'results/trajectories'/f'{label}.npz')['counts']
            p=replace(reference,D_init=initial_D)
            founder=.02+initial_D/p.N
            row={'baseline_vaf':founder/2,'baseline_day':0.,'observations':[],'outcome':{'kind':'unknown'}}
            background=Backgrounds([row],alpha,theta,p)
            f=ForwardFilter(p,alpha,theta,background.for_patient(0),killed=False,cap=10000,eta=0.,refinement=5)
            f.propagate(0.,1460.)
            physical=f.state[:len(f.d)]/f.state[:len(f.d)].sum()
            mean_D=float(np.sum(f.d*physical)/p.N)
            variance_D=float(np.sum((f.d/p.N-mean_D)**2*physical))
            mean_C=background.for_patient(0)(1460)*(1-mean_D)
            for index,typ,hybrid in [(0,'C',mean_C),(1,'D',mean_D)]:
                final=data[:,-1,index]/p.N;mean=float(final.mean());se=float(final.std(ddof=1)/math.sqrt(len(final)))
                rows.append({'alpha':alpha,'theta_days':theta,'initial_D_fraction':initial_D/p.N,'type':typ,
                             'Monte_Carlo_mean_fraction':mean,'Monte_Carlo_mean_SE':se,'hybrid_mean_fraction':hybrid,
                             'difference_in_MC_SE':(mean-hybrid)/se if se>0 else None,
                             'hybrid_predicted_MC_mean_SE':math.sqrt(variance_D/len(final)) if typ=='D' else None,
                             'difference_in_hybrid_predicted_MC_SE':(mean-hybrid)/math.sqrt(variance_D/len(final)) if typ=='D' and variance_D>0 else None,
                             'note':'Rare D tails make 200-replicate mean SE unstable when few nonzero paths exist.' if typ=='D' and initial_D==0 else 'Monte Carlo error and approximation error are distinct.'})
            events=list(csv.DictReader((ROOT/'results/trajectories'/f'{label}.csv').open()))
            k=sum(bool(r['first_hit_day']) and float(r['first_hit_day'])<=1460 for r in events);n=len(events)
            lo=float(beta.ppf(.025,k,n-k+1)) if k else 0.
            hi=float(beta.ppf(.975,k+1,n-k)) if k<n else 1.
            theory=risk_curve(alpha,theta,p,founder,targets=[1460.])[1460.]
            rows.append({'alpha':alpha,'theta_days':theta,'initial_D_fraction':initial_D/p.N,'type':'four_year_first_passage',
                         'Monte_Carlo_risk':k/n,'MC_binomial_lower_95':lo,'MC_binomial_upper_95':hi,'hybrid_risk':theory,
                         'note':'A binomial interval measures Monte Carlo uncertainty, not clinical risk uncertainty.'})
            print(f'CROSSCHECK alpha={alpha:g}, theta={theta:g}, D0={initial_D}',flush=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (ROOT/'results/simulation_crosscheck.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    (ROOT/'results/simulation_crosscheck.json').write_text(json.dumps(rows,indent=2,allow_nan=False))


if __name__=='__main__':check()
