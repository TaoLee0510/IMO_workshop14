"""English result tables, TeX results appendix and scientific diagnostic plots."""
import csv
import json
import math
import os
import subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parent


def tex(text):
    return str(text).replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')


def number(value):
    if value==0:return '$0$'
    if abs(value)<.001 or abs(value)>10000:
        mantissa,exponent=f'{value:.3e}'.split('e')
        return f'${float(mantissa):g}\\times10^{{{int(exponent)}}}$'
    return f'${value:.5g}$'


def probability_number(value):
    if .999<value<1:
        return '$1-'+number(1-value).strip('$')+'$'
    return number(value)


def make_report(config,cohort,fits):
    results=ROOT/'results';out=ROOT/'outputs/20261105_fixed_area';out.mkdir(parents=True,exist_ok=True)
    ref=next(r for r in fits if r['scenario']=='reference')
    risks=list(csv.DictReader((results/'four_year_risks.csv').open()))
    simulation=json.loads((results/'simulation_metadata.json').read_text())
    validation=json.loads((results/'validation.json').read_text())
    from scipy.special import gammainc
    fitrows=[]
    for r in sorted(fits,key=lambda x:(x['scenario']!='reference',x['scenario'])):
        fitrows.append({'scenario':r['scenario'],'status':r['status'],
                       **{k:r.get(k) for k in ['alpha','theta_days','peak_day','peak_additional_fitness','joint_nll','VAF_nll','clinical_given_VAF_nll','boundary_optimum','cap_refinement_loglik_difference','tolerance_refinement_loglik_difference']},
                       **r['parameters'],
                       'area_within_four_years_days':r['parameters']['area_days']*float(gammainc(r['alpha'],1460/r['theta_days'])) if 'alpha' in r else None})
    def numeric_row(row):
        answer={}
        for k,v in row.items():
            if isinstance(v,str):
                try:v=float(v) if v else None
                except ValueError:pass
            answer[k]=v
        return answer
    riskrows=[numeric_row(r) for r in risks]
    riskrows.sort(key=lambda r:(r['scenario']!='reference',r['scenario'],r['population']!='standardized',str(r.get('patient_id') or ''),r.get('initial_D_fraction') or 0))
    simrows=[{k:r[k] for k in ['name','label','alpha','theta_days','C_initial_fraction','D_initial_fraction','replicates','ever_D_fraction','extant_D_at_four_years_fraction','mean_C_four_year_fraction','mean_D_four_year_fraction','MC_risk_year_1','MC_risk_year_2','MC_risk_year_3','MC_risk_year_4','job_seed']} for r in simulation]
    patientrows=[{k:r[k] for k in ['patient_id','gene','variant_id','sample_type','baseline_day','baseline_vaf']}|
                 {'followup_VAF_count':len(r['observations']),'outcome':r['outcome']['kind'],'outcome_day':r['outcome'].get('day'),
                  'listed_first_event_day':r['outcome'].get('first_event_record_day',r['outcome'].get('day') if r['outcome']['kind']=='exact' else None),
                  'clinical_source_sheet':r['outcome'].get('source_sheet'),'clinical_source_excel_row':r['outcome'].get('source_excel_row'),
                  'source_evidence':r['clinical_evidence']} for r in cohort]
    parameters=[{'parameter':k,'reference':v,'units':{'b':'per day','s_D':'dimensionless','mu':'per C division','area_days':'days','sigma_vaf':'logit SD'}.get(k,'cells' if k in ['N','D_init'] else 'fraction'),
                 'definition':{'area_days':'Integral of f_C(t)-1 over post-CAR time; peak is not fixed.','s_D':'Fixed extra D fitness above C after CAR; zero before CAR.','b':'Moran replacement scale; total event rate N*b.'}.get(k,k),
                 'source':'Model rate reference: https://pubmed.ncbi.nlm.nih.gov/33636128/' if k=='b' else 'User-confirmed modeling assumption'} for k,v in config['reference'].items()]
    likelihood_rows=[{'scenario':f['scenario'],**{k:v for k,v in r.items() if k!='VAF_predictions'}} for f in fits for r in f.get('patient_details',[])]
    likelihood_rows.sort(key=lambda r:(r['scenario']!='reference',r['scenario'],int(r['patient_id'].split('_')[1])))
    profiles=list(csv.DictReader((results/'likelihood_profiles.csv').open())) if (results/'likelihood_profiles.csv').exists() else []
    crosschecks=json.loads((results/'simulation_crosscheck.json').read_text()) if (results/'simulation_crosscheck.json').exists() else []
    payload={'reference':ref,'fit_rows':fitrows,'risk_rows':riskrows,'simulation_rows':simrows,'patient_rows':patientrows,'parameter_rows':parameters,
             'likelihood_rows':likelihood_rows,'profile_rows':[numeric_row(r) for r in profiles],'crosscheck_rows':crosschecks,
             'notes':['Fixed area, not fixed peak. Main additional-fitness area is 1 day.','2 exact tMN events and 1 source-supported four-year censor record enter the joint likelihood.',
                      'Reference alpha and theta reach lower search bounds; the precise early peak is not resolved by the VAF observations.',
                      '11 clinical outcomes remain unknown, not negative.','Inference is a deterministic-C-background/stochastic-D approximation checked against exact Moran simulation.',
                      'Risks are model threshold probabilities; clinical ascertainment, marker causality and competing death are not calibrated.',
                      'Area is fixed on 0 to infinity; pulses beyond four years can contribute less area during observed follow-up.']}
    (results/'workbook_payload.json').write_text(json.dumps(payload,indent=2,allow_nan=False))
    # Figure: joint parameter fits and observed-data residuals.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from model import Parameters,fitness
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    end=min(365.,max(7.,ref['peak_day']+10*ref['theta_days']))
    days=np.unique(np.r_[0.,np.geomspace(1e-6,end,2000),np.linspace(0,end,1000)]);p=Parameters(**ref['parameters'])
    axes[0].plot(days,[fitness(float(t),ref['alpha'],ref['theta_days'],p)[0]-1 for t in days],color='#1765b0')
    axes[0].set(xlabel='Days after CAR infusion',ylabel='Additional C fitness',title='Joint fitted pulse; full area fixed at 1 day')
    predictions=[dict(patient_id=r['patient_id'],**v) for r in ref['patient_details'] for v in r.get('VAF_predictions',[])]
    axes[1].scatter([100*r['observed_vaf'] for r in predictions],[100*r['predicted_vaf'] for r in predictions],s=20,color='#1765b0')
    axes[1].plot([0,50],[0,50],color='#aaa',ls='--');axes[1].set(xlabel='Observed VAF (%)',ylabel='One-step filter prediction (%)',title='Predictions before current VAF observation')
    fig.savefig(results/'figures/joint_fit_diagnostics.pdf');fig.savefig(results/'figures/joint_fit_diagnostics.png',dpi=180);plt.close(fig)
    text=r'''\clearpage\section{Executed results}
\textbf{Analysis scope.} The source-record joint cohort contains 14 patients, 42 follow-up VAF observations, 2 exact tMN times (days 127 and 589), 1 source-supported day-1460 censor record and 11 unknown clinical outcomes. The latter contribute VAF likelihood only. These are exploratory, source-cohort conditional estimates, not calibrated population incidence.

\subsection{Reference joint fit}
'''
    text+=f"The reference fit gives $\\alpha={ref['alpha']:.5g}$ and $\\theta={ref['theta_days']:.5g}$ days; the peak occurs at day ${ref['peak_day']:.5g}$ with additional fitness ${ref['peak_additional_fitness']:.5g}$. The integrated additional advantage is fixed at 1 day. "
    text+=f"The joint negative log likelihood is ${ref['joint_nll']:.6g}$, decomposed into VAF ${ref['VAF_nll']:.6g}$ and clinical-given-VAF ${ref['clinical_given_VAF_nll']:.6g}$. "
    text+=f"The cap-refinement difference is ${ref['cap_refinement_loglik_difference']:.3g}$ log units; the tolerance-refinement difference is ${ref['tolerance_refinement_loglik_difference']:.3g}$. "
    text+=f"The optimum {'reaches' if ref['boundary_optimum'] else 'does not reach'} a numerical search boundary. Optimizer status: {tex(ref['optimizer_success'])}.\n"
    if ref['alpha']<1.011 and ref['theta_days']<.251:
        text+=r'Both parameters reach their lower search limits. The resulting peak precedes the earliest post-CAR VAF observation (day 1); its precise sub-day timing is not empirically resolved and should not be interpreted as a measured biological peak.'+'\n'
    text+=r'\begin{center}\includegraphics[width=\textwidth]{results/figures/joint_fit_diagnostics.pdf}\end{center}'+'\n'
    text+=r'''\subsection{Clinical likelihood contributions}
\begin{center}\begin{tabular}{lrrr}\toprule
Patient & Outcome & VAF NLL & Clinical-given-VAF NLL\\\midrule
'''
    for row in ref['patient_details']:
        if row['outcome']!='unknown':text+=f"{tex(row['patient_id'])} & {row['outcome']} & {number(-row['VAF_loglik'])} & {number(-row['clinical_given_VAF_loglik'])} \\\\\n"
    text+=r'\bottomrule\end{tabular}\end{center}'+'\n'
    text+=r'''\begin{samepage}\subsection{Standardized four-year probabilities}
Both scenarios start with $C/N=0.02$; $D/N$ is stated separately. Numbers are probabilities, not percentages. The $D/N=0.01$ scenario starts one cell below the strict 1\% threshold and therefore describes existing double-hit dynamics.
\begin{center}\begin{tabular}{lrrrr}\toprule
Initial $D/N$ & Year 1 & Year 2 & Year 3 & Year 4\\\midrule
'''
    for row in riskrows:
        if row['scenario']=='reference' and row['population']=='standardized':
            text+=f"{row['initial_D_fraction']:g} & "+' & '.join(probability_number(row[f'year_{y}_probability']) for y in range(1,5))+r' \\'+'\n'
    text+=r'\bottomrule\end{tabular}\end{center}\end{samepage}'+'\n'
    text+=r'''\subsection{Fixed-parameter sensitivity}
Every compatible scenario re-estimates the same shared alpha/theta pair using the full joint likelihood. The no-mutation/no-initial-D scenario is structurally incompatible with observed tMN and is recorded as such, without a probability floor. Varying $b$ holds $s_D$ fixed. Full numeric results, including likelihood components and numerical checks, are in the workbook and CSV files.
Changing the observation SD changes the observation model and its likelihood normalization; do not rank biological assumptions solely by cross-SD NLL. The fixed gamma area is over $[0,\infty)$; the area actually falling inside four years is recorded separately.
\begin{longtable}{lrrrr}\toprule
Scenario & $\alpha$ & $\theta$ (days) & Joint NLL & Boundary\\\midrule\endhead
'''
    for row in fitrows:
        if row['alpha'] is None:text+=f"{tex(row['scenario'])} & \\multicolumn{{4}}{{l}}{{Structurally incompatible}} \\\\\n"
        else:text+=f"{tex(row['scenario'])} & {row['alpha']:.4g} & {row['theta_days']:.4g} & {row['joint_nll']:.5g} & {row['boundary_optimum']} \\\\\n"
    text+=r'\bottomrule\end{longtable}'+'\n'
    text+=r'''\clearpage\subsection{Uncertainty and validation}
Profile likelihood results optimize the other parameter at each fixed value. A flat or boundary-limited profile means that a finite two-sided interval is not established by these data. Fixed-parameter scenario variation is not a confidence interval. Initial validation includes fixed-area integrals, an independently exponentiated tiny-event generator, a small-N full Moran comparison and neutral Monte Carlo moment checks. Numerical agreement does not establish the clinical meaning of the threshold or representative sampling.
'''
    if (results/'likelihood_profiles.csv').exists():
        profiles=list(csv.DictReader((results/'likelihood_profiles.csv').open()))
        text+=r'\begin{center}\begin{tabular}{lrr}\toprule Profile & Fixed value & $\Delta$ NLL\\\midrule'+'\n'
        for row in profiles:text+=f"{tex(row['profile'])} & {float(row['fixed_value']):.4g} & {float(row['delta_nll_from_joint_optimum']):.4g} \\\\\n"
        text+=r'\bottomrule\end{tabular}\end{center}'+'\n'
    if (results/'simulation_crosscheck.json').exists():
        text+=r'''\paragraph{Monte Carlo cross-check.} Exact Moran endpoint means and first-passage proportions are compared with the inference approximation at three grid points and both initial-D settings. Rare double-hit means with $D(0)=0$ have heavy sampling tails; an empirical SE from 200 paths can substantially understate uncertainty when large clones are not sampled. The cross-check therefore also supplies the D-mean Monte Carlo SE predicted by the forward distribution. This is distinct from clinical parameter uncertainty. See \texttt{simulation\_crosscheck.csv}.
'''
    text+=r'\clearpage\section{Exact stochastic Moran trajectories}'+'\n'
    for d in (.01,0.):
        for typ,title in [('single_mutant','Single mutant C: blue'),('double_hit','Double hit D: red')]:
            text+=f"\\subsection*{{{title}; initial $D/N={d:g}$}}\n"
            text+=f"\\includegraphics[width=\\textwidth]{{results/figures/{typ}_initial_D_{d:g}.pdf}}\n"
            text+=r'Each panel contains 200 complete random trajectories and their arithmetic mean. Area, not peak height, is held fixed. Trajectories continue after threshold crossing.\clearpage'+'\n'
    text+=r'\includegraphics[width=\textwidth]{results/figures/joint_optimum_trajectories.pdf}'+'\n'
    (ROOT/'executed_results.tex').write_text(text)
    # Artifact Tool authoring uses bundled runtime and libraries, not repository dependencies.
    dependencies=Path('/Users/4482173/.cache/codex-runtimes/codex-primary-runtime/dependencies')
    node=str(dependencies/'node/bin/node');modules=ROOT/'node_modules'
    if not modules.exists():modules.symlink_to(dependencies/'node/node_modules',target_is_directory=True)
    marker=results/'artifact_authoring_started.txt'
    if not marker.exists():
        mark='/Users/4482173/.codex/plugins/cache/openai-primary-runtime/spreadsheets/26.1030.10813/skills/spreadsheets/container_tools/mark_artifact_operation_started.mjs'
        subprocess.run([node,mark,'--operation-kind','create','--expected-output-count','1','--output-format','xlsx'],check=True)
        marker.write_text('create xlsx marker emitted once\n')
    subprocess.run([node,str(ROOT/'build_workbook.mjs')],check=True,cwd=ROOT)
    for _ in range(2):subprocess.run(['pdflatex','-interaction=batchmode','-halt-on-error','main.tex'],check=True,cwd=ROOT)
    print('REPORT AND WORKBOOK COMPLETE',flush=True)
