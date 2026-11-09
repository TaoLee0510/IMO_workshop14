"""Derived interpretation of existing diagnostic results; no parameter refitting."""
import csv
import json
import math
import os
from dataclasses import replace
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/tmn_gamma_mpl')
import numpy as np
from scipy.stats import gamma
from scipy.optimize import brentq
from scipy.special import logit
from gamma_model import Parameters, meanfield_predictions
from run_initial import load_vaf

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/interpretation'
OUT.mkdir(exist_ok=True)
config=json.loads((ROOT/'config.json').read_text())
p=Parameters(**config['reference'])
fits=list(csv.DictReader((ROOT/'results/preliminary_vaf_fits.csv').open()))
ref=fits[0]; alpha=float(ref['alpha']); theta=float(ref['theta_days'])
pids,bases,requested,obs,source=load_vaf()
predictions=list(csv.DictReader((ROOT/'results/preliminary_vaf_predictions.csv').open()))
mode=(alpha-1)*theta
area=p.amplitude/gamma.pdf(mode,alpha,scale=theta)
half=lambda t:(alpha-1)*math.log(t/mode)-(t-mode)/theta+math.log(2)
left=brentq(half,mode/1000,mode);right=brentq(half,mode,mode+20*theta)
null=meanfield_predictions(alpha,theta,replace(p,amplitude=0),bases,requested)
residual=(logit(2*null)-logit(2*obs))/p.sigma_vaf
rss=float(ref['residual_sum_squares']);nullrss=float(residual@residual)
contributions=[]
for pid in pids:
    rr=[r for r in predictions if r['patient_id']==pid]
    sse=sum(float(r['standardized_logit_residual'])**2 for r in rr)
    contributions.append(dict(patient_id=pid,gene=rr[0]['gene'],observations=len(rr),rss=sse,rss_share_percent=100*sse/rss))
contributions.sort(key=lambda r:-r['rss'])
summary=dict(analysis_scope='Derived interpretation of existing preliminary VAF-only mean-field results; no clinical joint fit or new parameter fit',
    alpha=alpha,theta_days=theta,peak_day=mode,peak_fitness=1+p.amplitude,
    half_maximum_interval_days=[left,right],pulse_fwhm_days=right-left,
    central_95pct_integrated_advantage_interval_days=gamma.ppf([.025,.975],alpha,scale=theta).tolist(),
    integrated_advantage_days=area,two_type_no_second_hit_odds_ratio_upper_bound=math.exp(p.b*area),
    RMS_standardized_residual=math.sqrt(rss/len(obs)),reference_rss=rss,
    null_comparison='Set founder CAR amplitude to zero, keeping all other fixed parameters unchanged; D extra fitness is retained',
    no_founder_CAR_pulse_rss=nullrss,rss_reduction_percent=100*(nullrss-rss)/nullrss,
    null_comparison_is_not_a_formal_hypothesis_test=True,patient_rss_contributions=contributions)
(OUT/'interpretation_calculations.json').write_text(json.dumps(summary,indent=2))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
selected=['Patient_44','Patient_235','Patient_180']
for ax,pid in zip(axes,selected):
    i=pids.index(pid); day0,vaf0=bases[i]
    rr=[r for r in predictions if r['patient_id']==pid]
    times=np.array([float(r['day']) for r in rr]);values=np.array([float(r['observed_vaf_percent']) for r in rr])
    grid=np.linspace(0,max(times),400)
    curve=meanfield_predictions(alpha,theta,p,[bases[i]],[(0,t) for t in grid])
    ax.plot(grid,100*curve,color='#bc5c23',lw=2,label='Shared-pulse prediction')
    ax.plot([day0,0],[100*vaf0,100*vaf0],color='#bc5c23',lw=2)
    ax.scatter(times,values,color='#235c82',s=30,zorder=3,label='Observed follow-up')
    ax.scatter([day0],[100*vaf0],marker='D',facecolor='white',edgecolor='#235c82',s=38,zorder=3,label='Conditioned baseline')
    ax.set(title=f"{pid.replace('_',' ')} | {rr[0]['gene']}",xlabel='Days relative to CAR infusion',ylabel='VAF (%)')
    ax.set_ylim(bottom=0 if pid!='Patient_235' else 38)
    ax.grid(axis='y',alpha=.2)
handles,labels=axes[0].get_legend_handles_labels()
fig.legend(handles,labels,loc='outside lower center',ncol=3,frameon=False)
fig.suptitle('Preliminary VAF-only mean-field fit: expansion and decline are poorly captured',fontsize=12)
fig.savefig(OUT/'patient_trajectory_diagnostics.png',dpi=180)
fig.savefig(OUT/'patient_trajectory_diagnostics.pdf')
plt.close(fig)
print(json.dumps(summary,indent=2))
