"""Build executed-result TeX pages and a typed, traceable workbook payload."""
from common import *
from plot_figures import read_csv
def esc(s):
    return str(s).replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
def sci(v):
    if v==0:return '0'
    power=math.floor(math.log10(abs(v)))
    return '$'+f'{v/10**power:.3f}'+r'\times10^{'+str(power)+'}$'
def main():
    risks=read_csv('patient_risks.csv');obs=read_csv('VAF_diagnostics.csv');sim=read_csv('simulation_summary.csv');yearly=read_csv('trajectory_yearly_summary.csv')
    diag=json.loads((ROOT/'results/diagnostic_summary.json').read_text());lik=json.loads((ROOT/'results/likelihood_diagnostic_summary.json').read_text())
    sections=[r'\clearpage\section{Executed diagnostics}',
      f"The frozen joint NLL is {FIT['joint_nll']:.6f}, comprising VAF NLL {FIT['VAF_nll']:.6f} and clinical-given-VAF NLL {FIT['clinical_given_VAF_nll']:.6f}. "
      f"The observed one-step VAF RMSE is {diag['VAF_RMSE_percentage_points']:.3f} percentage points. "
      f"The empirical coverage of the 95\\% observation prediction intervals is {100*diag['predictive_interval_coverage_fraction']:.1f}\\% across 42 observations. "
      "These are training-cohort diagnostics at fixed fitted parameters, not held-out validation.",
      r'\begin{center}\includegraphics[width=.98\linewidth]{figures/fit_diagnostics.pdf}\end{center}',
      "Quantile residuals use stable log-CDF/log-survival inversion of the predictive mixture, preserving extreme normal scores without probability clipping. Raw PIT values, stable log tails and approximate standardized residuals are retained.",
      r'\clearpage\subsection{Residuals by patient}\begin{center}\includegraphics[width=\linewidth]{figures/residuals_by_patient.pdf}\end{center}',
      r'\clearpage\subsection{Likelihood support}',
      r'\begin{center}\includegraphics[width=.91\linewidth]{figures/likelihood_surface.pdf}\\[8pt]\includegraphics[width=.98\linewidth]{figures/likelihood_profiles.pdf}\end{center}',
      "The profile searches optimize the other parameter within the confirmed bounds and evaluate explicit boundary/anchor candidates. A boundary-limited or flat curve does not establish a finite confidence interval."]
    if lik['better_than_frozen_by_more_than_0.01']:sections.append("A diagnostic search found a lower NLL than the frozen fit. This report retains the frozen fit, and the alternative point is recorded in the diagnostic data.")
    sections += [r'\clearpage\section{Gamma pulse and fitness functions}',r'\begin{center}\includegraphics[width=\linewidth]{figures/gamma_and_fitness.pdf}\end{center}',
      f"The normalized density has shape {A:.6f}, scale {THETA:g} day and mode {(A-1)*THETA:.6f} days. Its integral on zero to infinity is 1. Peak height is derived.",
      r'\clearpage\section{Patient one-year tMN probability}',
      r'\begin{center}\includegraphics[width=\linewidth,height=.77\textheight,keepaspectratio]{figures/one_year_tMN_heatmap.pdf}\end{center}',
      "Rows are ordered by baseline founder fraction. The heatmap uses a logarithmic color scale and annotates actual probabilities. Each forecast conditions on the patient's baseline state and event-free CAR day 0; subsequent VAF and eventual clinical outcome are not used to set that patient's risk.",
      r'\clearpage\subsection{Year 1--4 probability table}',
      r'\begin{center}\small\begin{tabular}{lrrrrr}\toprule Patient & Baseline VAF (\%) & Year 1 & Year 2 & Year 3 & Year 4\\\midrule']
    for r in risks:sections.append(esc(r['patient_id'])+f" & {100*r['baseline_VAF']:.3f} & "+' & '.join(sci(r[f'year{y}_probability']) for y in range(1,5))+r'\\')
    sections += [r'\bottomrule\end{tabular}\end{center}',
       "Probabilities refer to first crossing of the aggregate double-hit threshold, not a clinically calibrated incidence model. With 200 simulations, an empirical event fraction has increments of 0.005. Rare risks therefore come from the shared-state forward model, with exact simulations used for checks.",
       r'\clearpage\section{Patient overview figures}']
    for kind,caption in [('VAF','VAF observation and sequential prediction'),('C','Single-mutant C: 200 trajectories and their arithmetic mean'),('D','Double-hit D: 200 trajectories and their arithmetic mean'),('tMN','Baseline-initialized cumulative tMN probability on a log10 scale')]:
        sections += [r'\begin{center}\includegraphics[width=.99\linewidth]{figures/overview_'+kind+r'.pdf}\end{center}',esc(caption),r'\clearpage']
    sections += [r'\section{Individual patient dashboards}',
      "Black points are follow-up VAF observations; open diamonds are baseline values. Blue VAF bands are 95\\% observation prediction intervals. "
      "The VAF curve is a sequential predictive mean and ends at diagnosis for known events. C and D panels show 200 independent full Moran paths (opacity 0.035) and their arithmetic mean, including zero/extinct paths. "
      "Trajectories continue as an idealized process after an observed diagnosis; postdiagnosis treatment is not modeled. "
      "Red dotted lines mark observed tMN and gray dotted lines mark censoring. Unknown outcomes are not negatives. The 14 individual VAF, trajectory and tMN figures are also supplied separately."]
    for row in COHORT:sections += [r'\clearpage\begin{center}\includegraphics[width=.97\linewidth,height=.90\textheight,keepaspectratio]{figures/patients/'+row['patient_id']+r'_dashboard.pdf}\end{center}']
    (ROOT/'results_pages.tex').write_text('\n\n'.join(sections)+'\n')
    params=[{'parameter':k,'value':v,'definition':definition} for k,v,definition in [
      ('alpha',A,'Shared fitted Gamma shape'),('theta_days',THETA,'Shared fitted Gamma scale, days'),('N',P.N,'Competing population, cells'),
      ('b',P.b,'Moran turnover scale per day'),('s_D',P.s_D,'Extra post-CAR D fitness'),('mu',P.mu,'C to D mutation probability per C division'),
      ('Gamma_PDF_integral',1.,'Integral of normalized Gamma density over [0,infinity)'),('D_init',P.D_init,'D count at each measured baseline'),
      ('threshold_fraction',P.threshold_fraction,'Strict aggregate D fraction threshold'),('sigma_vaf',P.sigma_vaf,'SD on logit(2 VAF) observation scale')]]
    initial_rows=[{**initial(r),'variant_id':r['variant_id'],'source':'MOESM9 longitudinal extraction; MOESM7 clinical source; exact rows in cohort input'} for r in COHORT]
    likelihoods=[{'patient_id':r['patient_id'],'gene':r['gene'],'outcome':r['outcome'],'VAF_nll':-r['VAF_loglik'],'clinical_given_VAF_nll':-r['clinical_given_VAF_loglik'],'joint_nll':-r['joint_loglik']} for r in FIT['patient_details']]
    payload={'risks':risks,'observations':obs,'simulations':sim,'yearly_trajectories':yearly,'risk_curves':read_csv('risk_curves.csv'),'initial_states':initial_rows,'parameters':params,'likelihoods':likelihoods,'fit':FIT,
      'notes':['One founder marker per patient; baseline founder fraction = 2 x VAF under a diploid heterozygous assumption.',
               'Risk forecasts condition on baseline state and event-free CAR day 0; follow-up VAF and eventual outcome are not used to set individual risks.',
               'Parameters are frozen at the latest joint fit. Prediction intervals do not include fitted-parameter uncertainty.',
               'Eleven unknown clinical outcomes are not coded as negatives. Patient 25 is censored at day 1460 based on its listed first event at day 1561.']}
    write_json(ROOT/'results/workbook_payload.json',payload);print('TeX result pages and workbook payload written.')
if __name__=='__main__':main()
