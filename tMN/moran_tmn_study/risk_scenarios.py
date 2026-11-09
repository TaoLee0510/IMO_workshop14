"""Aggregate double-hit first-passage scenarios; no clinical-outcome calibration.

The D birth/death process retains exact fixed-N Moran replacement rates conditional
on a deterministic mutant proportion q(t) among the remaining N-D cells. This
approximates stochastic fluctuations of M, not the aggregate threshold event.
The killed forward equation on D=0,...,1000 is solved by implicit Euler.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import hashlib
from pathlib import Path
import time
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/tmn_risk_matplotlib')
import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
N = 100_000
B = 1 / 7
MU_GRID = [0, 1e-7, 1e-6, 1e-5, 1e-4]
R_GRID = [0, .001, .003, .01]
D0_GRID = [0, 1, 10, 100]
GENES = ['DNMT3A', 'TP53', 'PPM1D']


def fitness(t, s2):
    ld = .1 * math.exp(-(t + 5) / 7) if t >= -5 else 0.
    car = 0.
    if t >= 0:
        car = s2 * math.exp((math.log(2) / 2) * (t - 10)) if t < 10 else s2 * math.exp(-(math.log(2) / 30) * (t - 10))
    return 1. + ld + car


def background_derivative(q, fm, mu, b):
    # Composition among non-D cells, ignoring D in the mean fitness denominator.
    return b * q * (1 - q) * ((1 - mu) * fm - 1) / (1 + q * (fm - 1))


def next_background(q, t, step, s2, mu, b, constant_fm=None):
    f0 = constant_fm if constant_fm is not None else fitness(t + 1e-10, s2)
    fh = constant_fm if constant_fm is not None else fitness(t + step / 2, s2)
    f1 = constant_fm if constant_fm is not None else fitness(t + step - 1e-10, s2)
    k1 = background_derivative(q, f0, mu, b)
    k2 = background_derivative(q + step * k1 / 2, fh, mu, b)
    k3 = background_derivative(q + step * k2 / 2, fh, mu, b)
    k4 = background_derivative(q + step * k3, f1, mu, b)
    new = q + step * (k1 + 2*k2 + 2*k3 + k4) / 6
    if not -1e-10 <= new <= 1 + 1e-10:
        raise ValueError('Background RK4 left [0,1]; reduce timestep.')
    return min(1., max(0., new))


def rates(d, q, fm, sd, mu, n=N, b=B):
    remaining = n - d
    # q represents M/(N-D), so baseline D is removed from founder burden.
    m = remaining * q
    fd = fm + sd
    denominator = remaining + (fm - 1)*m + fd*d
    pd = (fd*d + mu*fm*m) / denominator
    return b * remaining * pd, b * d * (1-pd)


def first_passage(*, s2, founder_count, d0, mu, r, t0, targets,
                  n=N, b=B, refinement=1., constant_fm=None, uniform_step=None):
    """Return cumulative first-passage probabilities conditional on baseline.

    founder_count is M0+D0. Noninteger total clone size is retained in the
    deterministic M background; D0 and the absorbing threshold are integers.
    """
    k = math.floor(.01*n) + 1
    if not 0 <= d0 <= founder_count <= n:
        raise ValueError('Incompatible founder burden and D0.')
    targets = sorted(set(float(x) for x in targets if x >= t0))
    if not targets:
        return {}, {'max_mass_error': 0.}
    if d0 >= k:
        return {x: 1. for x in targets}, {'max_mass_error': 0.}
    if mu == 0 and d0 == 0:
        return {x: 0. for x in targets}, {'max_mass_error': 0.}
    d = np.arange(k, dtype=float)
    p = np.zeros(k); p[d0] = 1.
    q = (founder_count - d0)/(n-d0)
    sd = r/b
    band = np.zeros((3, k))
    t = float(t0)
    absorbed = 0.
    max_error = 0.
    result = {}
    # Split at treatment discontinuities and the peak; output times are exact.
    boundaries = sorted(set(targets + [v for v in (-5., 0., 10., 180., 365.) if t0 < v <= targets[-1]]))
    target_set = set(targets)
    if t in target_set:
        result[t] = 0.
    for end in boundaries:
        while t < end - 1e-9:
            nominal = uniform_step if uniform_step is not None else (5. if t < -5 or t >= 365 else 2. if t >= 180 else 1.)
            step = min(nominal/refinement, end-t)
            q_next = next_background(q,t,step,s2,mu,b,constant_fm)
            fm = constant_fm if constant_fm is not None else fitness(t+step/2,s2)
            birth, death = rates(d,(q+q_next)/2,fm,sd,mu,n,b)
            band[1] = 1 + step*(birth+death)
            band[0, 1:] = -step*death[1:]
            band[2, :-1] = -step*birth[:-1]
            p = solve_banded((1,1),band,p,overwrite_b=False,check_finite=False)
            absorbed += step*birth[-1]*p[-1]
            t += step
            q = q_next
        max_error = max(max_error, abs(float(np.sum(p))+absorbed-1))
        if p.min() < -1e-10 or not -1e-10 <= absorbed <= 1+1e-9:
            raise ArithmeticError('Probability conservation/positivity failure.')
        if end in target_set:
            result[end] = min(1.,max(0.,float(absorbed)))
    return result, {'max_mass_error': max_error, 'final_background_nonD_mutant_fraction': q}


def exact_small_n(n, m0, d0, fm, sd, mu, b, times):
    """Exact full two-dimensional Moran forward equation at small N, constant fitness."""
    k = math.floor(.01*n)+1
    states = [(m,d) for d in range(k) for m in range(n-d+1)]
    lookup = {state:i for i,state in enumerate(states)}
    ii=[]; jj=[]; vv=[]
    for j,(m,d) in enumerate(states):
        w=n-m-d; f=w+fm*m+(fm+sd)*d
        pw=w/f; pm=(1-mu)*fm*m/f; pd=((fm+sd)*d+mu*fm*m)/f
        events=[((m+1,d),b*w*pm),((m-1,d),b*m*pw),((m,d+1),b*w*pd),((m,d-1),b*d*pw),((m+1,d-1),b*d*pm),((m-1,d+1),b*m*pd)]
        total=0.
        for dest,rate in events:
            if rate <= 0: continue
            total += rate
            if dest in lookup:
                ii.append(lookup[dest]); jj.append(j); vv.append(rate)
        ii.append(j); jj.append(j); vv.append(-total)
    generator=coo_matrix((vv,(ii,jj)),shape=(len(states),len(states))).tocsr()
    p=np.zeros(len(states)); p[lookup[(m0,d0)]]=1.
    answer={}; previous=0.
    for t in times:
        p=expm_multiply(generator*(t-previous),p)
        answer[t]=1-float(p.sum()); previous=t
    return answer


def validate():
    records=[]
    # Exact invariants and an independent full Moran comparison.
    assert first_passage(s2=1,founder_count=1000,d0=0,mu=0,r=.01,t0=-10,targets=[365])[0][365] == 0
    n=200; times=[25.,50.,100.]
    for m0,d0,mu,r,fm in [(80,0,.001,.003,1.05),(79,1,.001,.003,1.05),(80,0,.01,0,1.),(19,1,.001,.01,1.1)]:
        exact=exact_small_n(n,m0,d0,fm,r/B,mu,B,times)
        approx,diag=first_passage(s2=0,founder_count=m0+d0,d0=d0,mu=mu,r=r,t0=0,targets=times,n=n,constant_fm=fm,uniform_step=.025)
        records.append({'kind':'exact_small_N_comparison','N':n,'M0':m0,'D0':d0,'mu':mu,'r':r,'fM':fm,'exact':exact,'hybrid':approx,'max_absolute_difference':max(abs(exact[t]-approx[t]) for t in times),'diagnostics':diag})
    for s2,d0,mu,r,founder in [(.3,0,1e-6,.01,1000),(1.,100,1e-5,.003,1000),(0.,10,1e-4,.01,1000),(2.,1,1e-4,.01,200),(.5,1,1e-6,.01,20)]:
        series=[]
        for refinement in [1.,2.,4.]:
            risk,diag=first_passage(s2=s2,founder_count=founder,d0=d0,mu=mu,r=r,t0=-6,targets=[365.,730.,1825.],refinement=refinement)
            series.append({'refinement':refinement,'risk':risk,'diagnostics':diag})
        records.append({'kind':'time_step_convergence','s2':s2,'D0':d0,'mu':mu,'r':r,'founder_count':founder,'series':series,'max_difference_coarse_vs_fine':max(abs(series[0]['risk'][t]-series[2]['risk'][t]) for t in [365.,730.,1825.]),'max_difference_half_vs_fine':max(abs(series[1]['risk'][t]-series[2]['risk'][t]) for t in [365.,730.,1825.])})
    ROOT.joinpath('results').mkdir(exist_ok=True)
    ROOT.joinpath('results/risk_validation.json').write_text(json.dumps(records,indent=2))
    print(json.dumps({'validation':[{'kind':x['kind'],'max_abs_error':x.get('max_absolute_difference',x.get('max_difference_coarse_vs_fine'))} for x in records]}),flush=True)
    return records


def read_fit_records(path):
    payload=json.loads(Path(path).read_text())
    # The endpoint fitter owns the schema; accept explicit common container names.
    records=payload if isinstance(payload,list) else payload.get('patient_gene_fits',payload.get('fits',payload.get('results',[])))
    if not records:
        raise ValueError(f'No patient-gene records in {path}; top-level keys={list(payload)}')
    return payload,records


def s2_value(row):
    for key in ('s2','s2_hat','s2_estimate'):
        if row.get(key) is not None:
            return float(row[key])
    raise KeyError(f'No s2 field: {row.keys()}')


def one_unit(task):
    rows=[]; curves=[]
    targets=task['targets']
    for mu,r,d0 in itertools.product(MU_GRID,R_GRID,D0_GRID):
        base={k:task[k] for k in ('scope','patient_id','gene','variant_id','s2','first_day','last_day','founder_count')}
        base.update({k:task.get(k,'') for k in ('fit_status','growth_fit_boundary','growth_fit_mismatch','quality_flags')})
        base.update(mu=mu,r_per_day=r,s_D=r/B,D0=d0)
        if d0 > task['founder_count']:
            rows.append({**base,'status':'incompatible_D0_exceeds_founder','day':None,'P_first_passage':None})
            continue
        risk,diag=first_passage(s2=task['s2'],founder_count=task['founder_count'],d0=d0,mu=mu,r=r,t0=task['first_day'],targets=targets,refinement=task.get('refinement',1.))
        for day in targets:
            if day not in risk: continue
            row={**base,'status':'computed','day':day,'P_first_passage':risk[day],'max_mass_error':diag['max_mass_error']}
            if day in task['report_targets']:
                rows.append(row)
            if task['scope']=='standardized':
                curves.append(row)
    return rows,curves


def write_csv(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader(); writer.writerows(rows)


def tasks_from_fit(fit_path,mode,refinement):
    payload,records=read_fit_records(fit_path)
    usable=[]
    for row in records:
        try:s2=s2_value(row)
        except (KeyError,ValueError,TypeError):continue
        if math.isfinite(s2):usable.append((row,s2))
    tasks=[]
    if mode in ('standardized','all'):
        for gene in GENES:
            vals=[s for row,s in usable if row.get('gene')==gene]
            if not vals:continue
            targets=sorted(set([0.,365.,730.,1825.]+list(np.arange(30.,1825.,30.))))
            tasks.append(dict(scope='standardized',patient_id='',gene=gene,variant_id='',s2=float(np.median(vals)),first_day=0.,last_day=1825.,founder_count=1000.,targets=targets,report_targets=[365.,730.,1825.],refinement=refinement))
    if mode in ('patients','all'):
        for row,s2 in usable:
            pid=row.get('patient_id',row.get('Patient_ID'))
            gene=row.get('gene')
            candidates=[v for v in row.get('variant_initial_states',[]) if 0<float(v['first_vaf_percent'])<50 and 0<float(v['last_vaf_percent'])<50 and float(v['last_day'])>float(v['first_day'])]
            if not candidates:continue
            index=max(candidates,key=lambda v:(float(v['first_vaf_percent']),v['variant_id']))
            first,last=float(index['first_day']),float(index['last_day'])
            targets=sorted(set(t for t in [365.,730.,1825.,last] if t>=first))
            tasks.append(dict(scope='patient_gene_index_variant',patient_id=pid,gene=gene,variant_id=index['variant_id'],s2=s2,first_day=first,last_day=last,founder_count=2*float(index['first_vaf_percent'])/100*N,targets=targets,report_targets=targets,refinement=refinement,fit_status=row.get('status',''),growth_fit_boundary=row.get('boundary',False),growth_fit_mismatch=row.get('mismatch_flag',False),quality_flags=index.get('quality_flags','')))
    return tasks


def make_figure(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.titlesize':10})
    fig,axes=plt.subplots(3,4,figsize=(12,8),constrained_layout=True)
    for i,gene in enumerate(GENES):
        for j,r in enumerate(R_GRID):
            ax=axes[i,j]; matrix=np.full((4,5),np.nan)
            for row in rows:
                if row['gene']==gene and row.get('day')==1825 and row['r_per_day']==r and row['status']=='computed':
                    matrix[D0_GRID.index(row['D0']),MU_GRID.index(row['mu'])]=100*row['P_first_passage']
            im=ax.imshow(matrix,origin='lower',aspect='auto',cmap='Blues',norm=Normalize(0,100))
            ax.set_xticks(range(5),['0',r'$10^{-7}$',r'$10^{-6}$',r'$10^{-5}$',r'$10^{-4}$'])
            ax.set_yticks(range(4),[str(x) for x in D0_GRID])
            ax.set_title(f'{gene} | r = {r:g}/day')
            if i==2:ax.set_xlabel('Second-hit probability per division')
            if j==0:ax.set_ylabel('Baseline double-hit cells')
            for yy in range(4):
                for xx in range(5):
                    value=matrix[yy,xx]
                    if np.isfinite(value):
                        label='<0.1' if 0<value<.1 else f'{value:.1f}'
                        ax.text(xx,yy,label,ha='center',va='center',fontsize=8,color='white' if value>55 else '#202020')
    fig.suptitle('Five-year aggregate double-hit first-passage scenarios\nCommon founder burden: 1% of cells at day 0; gene median fitted s₂',fontsize=12)
    fig.colorbar(im,ax=axes,shrink=.7,label='Model first-passage probability (%)')
    fig.get_layout_engine().set(rect=(0,.045,1,.955))
    fig.text(.02,.014,'Approximate scenarios; reference s₂ held fixed, not jointly refitted to VAF endpoints; no clinical calibration.',fontsize=9)
    ROOT.joinpath('figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures/risk_standardized_5year.png',dpi=180)
    fig.savefig(ROOT/'figures/risk_standardized_5year.pdf')
    plt.close(fig)


def make_curves(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    from matplotlib.lines import Line2D
    colors={'DNMT3A':'#3568A6','TP53':'#D77B29','PPM1D':'#7A5BA7'}
    styles={0:':',.001:'--',.003:'-.',.01:'-'}
    fig,axes=plt.subplots(1,3,figsize=(12,4.8),sharex=True,sharey=True,layout='constrained')
    for ax,gene in zip(axes,GENES):
        selected=[r for r in rows if r['gene']==gene and r['status']=='computed']
        days=sorted(set(r['day'] for r in selected))
        lower=[min(r['P_first_passage'] for r in selected if r['day']==d) for d in days]
        upper=[max(r['P_first_passage'] for r in selected if r['day']==d) for d in days]
        ax.fill_between(np.array(days)/365,lower,upper,color=colors[gene],alpha=.12)
        for advantage in R_GRID:
            curve=sorted([r for r in selected if r['mu']==1e-6 and r['D0']==10 and r['r_per_day']==advantage],key=lambda r:r['day'])
            ax.plot([r['day']/365 for r in curve],[r['P_first_passage'] for r in curve],color=colors[gene],linestyle=styles[advantage],linewidth=1.8)
        ax.set_title(gene,color=colors[gene],fontweight='bold')
        ax.set(xlim=(0,5),ylim=(0,1),xlabel='Years since CAR T (day 0)')
        ax.set_xticks(range(6));ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#dddddd',linewidth=.6)
    axes[0].set_ylabel('Aggregate first-passage probability')
    fig.suptitle('Double-hit threshold crossing over time\nLines: μ = 10⁻⁶ per division, D₀ = 10; shading: all 80 approved scenarios (not a confidence interval)',fontsize=11)
    handles=[Line2D([0],[0],color='#333333',linestyle=styles[r],label=f'r = {r:g}/day') for r in R_GRID]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.075),ncol=4,frameon=False)
    fig.get_layout_engine().set(rect=(0,.18,1,.82))
    fig.text(.025,.025,'Common founder burden: 1% of cells; gene median reference s₂. Approximate, retrospective scenarios; no clinical calibration.',fontsize=8.5)
    fig.savefig(ROOT/'figures/risk_standardized_curves.png',dpi=180)
    fig.savefig(ROOT/'figures/risk_standardized_curves.pdf')
    plt.close(fig)


def render_existing():
    for filename,renderer in [('risk_standardized_scenarios.csv',make_figure),('risk_standardized_curves.csv',make_curves)]:
        with ROOT.joinpath('results',filename).open() as handle: rows=list(csv.DictReader(handle))
        for row in rows:
            for key in ('day','r_per_day','mu','D0','P_first_passage'):
                if row.get(key) not in (None,''):row[key]=float(row[key])
        renderer(rows)


def coupled_endpoint(key):
    """Diagnostic only: full mean-field propagation without refitting s2."""
    s2,founder,d0,mu,r,t0,last=key
    state=np.array([(founder-d0)/N,d0/N],dtype=float)
    def rhs(t,xy):
        x,y=xy; fm=fitness(t,s2); fd=fm+r/B
        mean=1-x-y+fm*x+fd*y
        return [B*((1-mu)*fm*x/mean-x),B*((fd*y+mu*fm*x)/mean-y)]
    start=t0
    for end in sorted(set([last]+[v for v in (-5.,0.,10.) if t0<v<last])):
        if end<=start:continue
        sol=solve_ivp(rhs,(start,end),state,rtol=1e-8,atol=1e-11,max_step=15.)
        if not sol.success:raise ArithmeticError(sol.message)
        state=sol.y[:,-1]; start=end
    return 50.*float(state.sum()),float(state[1])


def endpoint_check(workers):
    path=ROOT/'results/risk_patient_gene_scenarios.csv'
    with path.open() as handle:rows=list(csv.DictReader(handle))
    with ROOT.joinpath('data/endpoints.csv').open() as handle:
        observed={(r['patient_id'],r['gene'],r['variant_id']):float(r['last_vaf_percent']) for r in csv.DictReader(handle)}
    keys=[]; rowkeys={}
    for i,row in enumerate(rows):
        if row['status']!='computed':continue
        key=tuple(float(row[col]) for col in ('s2','founder_count','D0','mu','r_per_day','first_day','last_day'))
        rowkeys[i]=key
        keys.append(key)
    unique=list(dict.fromkeys(keys))
    print(f'Checking coupled endpoint for {len(unique)} unique parameter scenarios.',flush=True)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        predictions=dict(zip(unique,executor.map(coupled_endpoint,unique,chunksize=20)))
    residuals=[]
    for i,key in rowkeys.items():
        row=rows[i]
        prediction,last_d=predictions[key]
        truth=observed[(row['patient_id'],row['gene'],row['variant_id'])]
        delta=prediction-truth
        row.update(coupled_endpoint_predicted_vaf_percent=prediction,observed_last_vaf_percent=truth,coupled_endpoint_residual_percentage_points=delta,coupled_endpoint_D_fraction=last_d,scenario_endpoint_refitted=False)
        residuals.append(abs(delta))
    write_csv(path,rows)
    summary={'unique_parameter_scenarios':len(unique),'diagnostic':'Full coupled three-state mean-field endpoint under each risk scenario; reference s2 was NOT re-fitted.','absolute_residual_percentage_points_median':float(np.median(residuals)),'absolute_residual_percentage_points_max':float(np.max(residuals)),'interpretation':'Large residuals identify scenarios incompatible with the observed endpoint at the fixed reference s2. No observation-error threshold, likelihood, exclusion, or re-fitting has been imposed.'}
    ROOT.joinpath('results/risk_endpoint_diagnostics.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary),flush=True)


def audit_outputs():
    with ROOT.joinpath('results/risk_patient_gene_scenarios.csv').open() as handle:rows=list(csv.DictReader(handle))
    groups={}
    for row in rows:
        if row['status']!='computed':continue
        probability=float(row['P_first_passage'])
        assert 0<=probability<=1
        if float(row['mu'])==0 and int(row['D0'])==0:assert probability==0
        key=tuple(row[k] for k in ('patient_id','gene','mu','r_per_day','D0'))
        groups.setdefault(key,[]).append((float(row['day']),probability))
    for values in groups.values():
        values=sorted(values)
        assert all(a[1]<=b[1]+1e-12 for a,b in zip(values,values[1:]))
    probability,diag=first_passage(s2=.1,founder_count=1000,d0=100,mu=1e-6,r=.01,t0=-133,targets=[-5.,0.,365.])
    preld_vaf,preld_d=coupled_endpoint((.1,1000,100,1e-6,.01,-133.,-5.))
    audit={'computed_patient_horizon_rows':len(rows),'unique_patient_scenarios':len(groups),'probability_bounds_and_time_monotonicity_pass':True,'mu0_D0zero_gives_exactly_zero_risk':True,'max_mass_error':max(float(row['max_mass_error']) for row in rows if row['status']=='computed'),'minimum_index_founder_cells':min(float(row['founder_count']) for row in rows),'source_sha256':{name:hashlib.sha256(ROOT.joinpath(name).read_bytes()).hexdigest() for name in ['data/endpoints.csv','results/fit_reference.json','risk_scenarios.py']},'pre_LD_propagation_diagnostic':{'note':'Illustrative approved-grid D0=100,mu=1e-6,r=.01 with 1%founder cells at day-133; this is a code diagnostic, not a fitted patient.','initial_vaf_percent':.5,'meanfield_vaf_percent_at_day_minus5':preld_vaf,'meanfield_D_fraction_at_day_minus5':preld_d,'conditional_first_passage_probability':probability,'probability_diagnostics':diag}}
    ROOT.joinpath('results/risk_audit.json').write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit),flush=True)


def run(mode,fit_path,workers,refinement):
    tasks=tasks_from_fit(fit_path,mode,refinement)
    print(f'Computing {len(tasks)} units, 80 scenarios each; workers={workers}',flush=True)
    start=time.monotonic(); rows=[]; curves=[]
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for index,(new_rows,new_curves) in enumerate(executor.map(one_unit,tasks),start=1):
            rows.extend(new_rows); curves.extend(new_curves)
            print(f'Completed {index}/{len(tasks)} units ({time.monotonic()-start:.1f}s)',flush=True)
    ROOT.joinpath('results').mkdir(exist_ok=True)
    standardized=[r for r in rows if r['scope']=='standardized']
    patients=[r for r in rows if r['scope']!='standardized']
    if standardized:
        write_csv(ROOT/'results/risk_standardized_scenarios.csv',standardized)
        write_csv(ROOT/'results/risk_standardized_curves.csv',curves)
        make_figure(standardized)
        make_curves(curves)
    if patients:write_csv(ROOT/'results/risk_patient_gene_scenarios.csv',patients)
    summary={'model':'fixed-N aggregate D first passage with deterministic non-D composition','threshold_cells':1001,'N':N,'b_per_day':B,'mu_grid':MU_GRID,'r_per_day_grid':R_GRID,'D0_grid':D0_GRID,'units':tasks,'elapsed_seconds':time.monotonic()-start,'computed_rows':sum(r['status']=='computed' for r in rows),'incompatible_scenarios':sum(r['status']!='computed' for r in rows),'limitations':['Approximate stochastic M background; conditional D first passage solved by a killed birth-death master equation.','Risk is conditional on no prior threshold crossing at first observation and the specified baseline D0.','Fitted s2 uses the later VAF, so these are retrospective conditional scenarios, not baseline-only prospective forecasts.','Reference fitted s2 is held fixed across mutation/fitness/D0 scenarios; scenarios are not all re-fitted to endpoint data.','Index variant selected by largest eligible baseline VAF; unphased variants are not summed.','Clinical tMN probabilities are not calibrated; absence from tMN_VAF is not verified event-free follow-up.','Scenario variation is not a confidence interval and does not include fitted-s2 or kinetic uncertainty.','Standardized comparisons use a common 1% founder cell fraction, not each gene cohort’s observed baseline distribution.']}
    ROOT.joinpath(f'results/risk_{mode}_metadata.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps({'output_rows':len(rows),'elapsed_seconds':summary['elapsed_seconds']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['validate','standardized','patients','all','endpoint-check','figures','audit'],default='validate')
    parser.add_argument('--fit',type=Path,default=ROOT/'results/fit_reference.json')
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--refinement',type=float,default=1.)
    args=parser.parse_args()
    if args.mode=='validate':validate()
    elif args.mode=='endpoint-check':endpoint_check(args.workers)
    elif args.mode=='figures':render_existing()
    elif args.mode=='audit':audit_outputs()
    else:run(args.mode,args.fit,args.workers,args.refinement)
