"""English result narrative and typed workbook payload after all scenarios."""
import json,csv,math
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def rows(name):
    values=list(csv.DictReader((ROOT/name).open()))
    for row in values:
        for key,v in list(row.items()):
            if v=='':row[key]=None
            elif v in ['True','False']:row[key]=v=='True'
            elif key not in ['scenario','status','patient_id','gene','sample_type','outcome','family','variant_id','source','reason','inference']:
                try:row[key]=float(v)
                except ValueError:pass
    return values
def esc(s):return str(s).replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
def sci(v):
    if v is None:return 'n.a.'
    if v==0:return '0'
    power=math.floor(math.log10(abs(v)))
    return '$'+f'{v/10**power:.3f}'+r'\times10^{'+str(power)+'}$'
def main():
    config=json.loads((ROOT/'config.json').read_text());cohort=json.loads((ROOT/'data/cohort.json').read_text())
    fits=rows('results/joint_fits.csv');risks=rows('results/patient_risks.csv')
    assert len(fits)==79 and sum(r['status']=='fitted' for r in fits)==78
    ref=json.loads((ROOT/'results/fits/reference.json').read_text());old=json.loads((ROOT/'data/previous_reference.json').read_text())
    best=min((r for r in fits if r['status']=='fitted'),key=lambda r:r['joint_nll'])
    diagnostics=json.loads((ROOT/'patient_analysis/results/diagnostic_summary.json').read_text())
    comparison=json.loads((ROOT/'results/validation/Moran_comparison.json').read_text())
    parts=[r'\clearpage\section{Executed joint reference fit}',
       f"The reference fixes mu=0.01 and initial D=0. Its fitted shape is {ref['alpha']:.8g}, scale {ref['theta_days']:.8g} days, and derived mode {ref['peak_day']:.8g} days. "
       f"Joint NLL is {ref['joint_nll']:.6f}, decomposed into VAF NLL {ref['VAF_nll']:.6f} and clinical-given-VAF NLL {ref['clinical_given_VAF_nll']:.6f}. "
       f"The preceding mu=1e-6 reference had NLL {old['joint_nll']:.6f}. Parameters were re-estimated after changing mu.",
       f"The lowest NLL across compatible sensitivity scenarios is {best['joint_nll']:.6f}, in "+esc(best['scenario'])+
       ". This is a conditional scenario comparison; the requested reference remains mu=0.01 and initial D=0.",
       f"All 78 conditional optimizations ended at at least one search boundary. At the reference, the 95\\% one-step observation interval covers {round(42*diagnostics['predictive_interval_coverage_fraction'])}/42 follow-up VAFs ({100*diagnostics['predictive_interval_coverage_fraction']:.1f}\\%). VAF RMSE is {diagnostics['VAF_RMSE_percentage_points']:.3f} percentage points. These diagnostics show substantial model/data mismatch; the fitted timescale and clinical probabilities should be interpreted conditionally, rather than as validated physiological estimates.",
       r'\begin{landscape}\section{Mu and initial-D sensitivity}\begin{center}\includegraphics[width=.98\linewidth,height=.78\textheight,keepaspectratio]{figures/mu_D_joint_likelihood.pdf}\end{center}',
       "Every compatible cell re-optimizes alpha and theta. The gray zero-mutation/zero-D corner is structurally incompatible with observed positive tMN. No probability floor is used.",
       r'\clearpage\begin{center}\includegraphics[width=\linewidth]{figures/mu_D_fitted_parameters.pdf}\end{center}',
       "Parameter maps show conditional optima within the common bounds. A value at a search boundary is not an identified physiological timescale. At alpha=50 and theta=1460 days, the Gamma mode is 71540 days. The fraction of the normalized Gamma area occurring within the four-year observation window therefore varies across scenarios, while its full integral remains 1.",
       r'\end{landscape}\clearpage\section{Turnover and double-hit fitness sensitivity}',
       r'\begin{center}\includegraphics[width=\linewidth]{figures/b_sensitivity.pdf}\\[10pt]\includegraphics[width=\linewidth]{figures/s_D_sensitivity.pdf}\end{center}',
       r'\clearpage\section{New-reference patient diagnostics}',
       r'\begin{center}\includegraphics[width=.98\linewidth]{patient_analysis/figures/fit_diagnostics.pdf}\end{center}',
       r'\clearpage\section{New-reference Gamma pulse and one-year risk}',
       r'\begin{center}\includegraphics[width=\linewidth]{patient_analysis/figures/gamma_and_fitness.pdf}\\[10pt]\includegraphics[width=.88\linewidth,height=.65\textheight,keepaspectratio]{patient_analysis/figures/one_year_tMN_heatmap.pdf}\end{center}',
       r'\clearpage\subsection{Reference patient risk table}',
       r'\begin{center}\small\begin{tabular}{lrrrr}\toprule Patient & Year 1 & Year 2 & Year 3 & Year 4\\\midrule']
    reference_risks=[r for r in risks if r['scenario']=='reference']
    order={r['patient_id']:i for i,r in enumerate(cohort)};reference_risks.sort(key=lambda r:order[r['patient_id']])
    for r in reference_risks:parts.append(esc(r['patient_id'])+' & '+' & '.join(sci(r[f'year{y}_probability']) for y in range(1,5))+r'\\')
    parts += [r'\bottomrule\end{tabular}\end{center}',
       "Risk forecasts condition on each patient's measured baseline and event-free CAR day 0. They do not use the subsequent individual VAF or eventual outcome to set the forecast. The fitted common parameters were learned in this cohort, so these are in-sample conditional predictions.",
       r'\subsection{Full-process simulation comparison}',
       f"Across 56 patient/year comparisons at the new reference, the maximum difference between the inference approximation and 200 full Moran paths is {comparison['max_C_mean_error_percentage_points']:.4f} percentage points for mean C and {comparison['max_D_mean_error_percentage_points']:.4f} percentage points for mean D. The largest absolute tMN probability difference is {comparison['max_absolute_probability_difference']:.4f}. All 56 model probabilities fall within pointwise exact 95\\% binomial intervals for the simulated event proportions. These intervals measure Monte Carlo error at fixed parameters, not fitted-parameter uncertainty. With 200 paths, empirical event fractions have increments of 0.005 and cannot verify extremely small tail probabilities. Agreement with the simulator verifies the approximation in this reference setting; it does not establish agreement with clinical incidence.",
       r'\subsection{Patient figure guide}',
       "Open diamonds are baseline VAF; black points are follow-up measurements. The blue VAF band is a 95\\% observation prediction interval. Blue C and red D trajectory panels contain all 200 full Moran paths with opacity 0.035 and their arithmetic mean. Dotted lines show known diagnosis or censoring times. Process paths continue after diagnosis without modeling postdiagnosis treatment. Unknown outcomes remain unknown."]
    for r in cohort:parts += [r'\clearpage\begin{center}\includegraphics[width=.97\linewidth,height=.9\textheight,keepaspectratio]{patient_analysis/figures/patients/'+r['patient_id']+r'_dashboard.pdf}\end{center}']
    parts += [r'\clearpage\section{Patient-specific mu--D risk responses}',
       "Each patient's sensitivity panel holds the measured baseline founder burden fixed and subtracts D from that total when initializing C. Colors are probabilities on logarithmic scales. Year-one and year-four maps are supplied for every patient."]
    for i,r in enumerate(cohort):parts += [(r'\clearpage' if i else '')+r'\begin{center}\includegraphics[width=.98\linewidth,height=.85\textheight,keepaspectratio]{figures/patients/'+r['patient_id']+r'_mu_D_risk.pdf}\end{center}']
    (ROOT/'executed_results.tex').write_text('\n\n'.join(parts)+'\n')
    p=ref['parameters'];params=[{'parameter':k,'value':v,'definition':definition} for k,v,definition in [
        ('alpha',ref['alpha'],'Shared fitted Gamma shape'),('theta_days',ref['theta_days'],'Shared fitted Gamma scale, days'),
        ('N',p['N'],'Competing population size, cells'),('b',p['b'],'Moran turnover per day'),('s_D',p['s_D'],'Extra post-CAR D fitness'),
        ('mu',p['mu'],'Effective C-to-D probability per C division'),('D_init',p['D_init'],'Measured-baseline D cells'),
        ('Gamma_PDF_integral',1.,'Integral of normalized Gamma PDF on [0,infinity)'),('threshold_fraction',p['threshold_fraction'],'Strict aggregate D fraction threshold'),
        ('sigma_vaf',p['sigma_vaf'],'Observation SD on logit(2 VAF) scale'),
        ('fitness_multiplier',p['area_days'],'Fixed numeric multiplier on normalized Gamma density'),
        ('K_cells',math.floor(p['threshold_fraction']*p['N'])+1,'First-passage barrier: floor(N times threshold fraction) plus 1'),
        ('diagnosis_delay_days',0,'No delay between threshold crossing and modeled tMN endpoint')]]
    states=[]
    for d in config['sensitivity']['D_init']:
        for r in cohort:
            total=2*p['N']*r['baseline_vaf'];states.append({'patient_id':r['patient_id'],'gene':r['gene'],'baseline_day':r['baseline_day'],'baseline_VAF':r['baseline_vaf'],
              'D_init':d,'initial_C_continuous':total-d,'initial_C_integer':round(total)-d,'initial_W_integer':p['N']-round(total),'D_fraction_of_founder':d/total,'sample_type':r['sample_type']})
    compare_by_key={(r['patient_id'],r['day']):r for r in rows('results/Moran_comparison.csv')}
    yearly=[r|compare_by_key.get((r['patient_id'],r['day']),{}) for r in rows('patient_analysis/results/trajectory_yearly_summary.csv')]
    payload={'fits':fits,'risks':sorted(risks,key=lambda r:(r['scenario']!='reference',r['scenario'],order[r['patient_id']])),
      'reference':ref,'reference_risks':reference_risks,'parameters':params,'initial_states':states,
      'diagnostics':diagnostics,'comparison':comparison,
      'observations':rows('patient_analysis/results/VAF_diagnostics.csv'),
      'simulations':rows('patient_analysis/results/simulation_summary.csv'),
      'yearly':yearly}
    (ROOT/'results/workbook_payload.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    print('Executed TeX pages and workbook payload created.')
if __name__=='__main__':main()
