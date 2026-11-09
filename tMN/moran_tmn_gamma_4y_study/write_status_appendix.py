"""Generate an English execution-status appendix from inspected output tables."""
import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
audit=json.loads((ROOT/'data/audit.json').read_text())
validation=json.loads((ROOT/'results/validation.json').read_text())
fits=list(csv.DictReader((ROOT/'results/preliminary_vaf_fits.csv').open()))
reference=next(r for r in fits if r['scenario']=='reference')
risks=list(csv.DictReader((ROOT/'results/standardized_risk_demonstration.csv').open()))
text=r'''\clearpage
\appendix
\section{Initial execution record: 5 November 2026}
\textbf{Status.} Protocol authoring, data extraction, initial engine validation and preliminary VAF diagnostics have been executed. Clinical joint inference is pending verified outcome and follow-up data. No missing clinical outcome was coded as a negative.

\textbf{Source and cohort audit.} Scanned 135,806 mutation rows and matched 1,696 longitudinal records to the selected six-gene endpoint inventory: 973 cfDNA, 353 PBL, 320 PBMC and 50 tumor observations. All matched endpoint allele pairs use REF-to-TUMOR orientation. The audit flags 98 duplicate-call rows, 287 values outside the positive diploid mapping, and 320 sample-type naming conflicts. These flag categories can overlap. There are no sample-day metadata conflicts among these matched observations.

The initial patient-level blood-cell diagnostic includes 12 patients with 54 observations: 12 baseline initial conditions and 42 likelihood observations. Exclusions among 67 source patients are 21 without a unique eligible pre-CAR blood-cell baseline, 31 whose baseline-selected marker has no eligible follow-up, 2 with unresolved multiple CAR time origins, and 1 source-positive patient without a diagnosis time. These stricter rules materially narrow the earlier endpoint cohort. The 320 raw PBMC records are labeled PBL in sample metadata; the mutation-table label is retained and compartments are not pooled within a patient.

\textbf{Validation.} Initial invariants passed. An independent constant-rate matrix exponential agrees with the refined shared-state joint log likelihood to 0.0068 log units at a 0.125-day backward-Euler step. One small-$N$ full Moran comparison differs from the hybrid risk by approximately $1.11\times10^{-4}$ in absolute probability. Noiseless synthetic mean-field data recover known shape/scale parameters. This is not a clinical joint parameter-recovery experiment. Adaptive full-$N$ risk tolerance checks also passed. These are initial numerical checks, not broad certification of the hybrid approximation.

\textbf{Preliminary fit.} One reference and eight individually refitted VAF-only scenarios have completed. The reference has
'''
text+=f"$\\alpha={float(reference['alpha']):.3f}$, $\\theta={float(reference['theta_days']):.4f}$ days, and peak day ${float(reference['peak_day']):.2f}$. "
text+=f"The standardized residual sum of squares is ${float(reference['residual_sum_squares']):.2f}$ across {reference['observations']} observations, with maximum absolute residual ${float(reference['maximum_abs_standardized_residual']):.2f}$ measurement SDs. "
text+=r'''The shared nonnegative pulse leaves substantial lack of fit; some sensitivity fits also reach the numerical search boundary. Optimizer convergence does not establish biological adequacy or parameter certainty. No clinical outcomes were used and no confidence intervals have been established.

\textbf{Risk calculation demonstration only.} The following are model threshold probabilities under the preliminary VAF fit and fixed reference parameters, starting at day 0 with the indicated VAF. They are \textbf{not clinically calibrated tMN probabilities}, omit competing death, and are not patient-specific estimates. Entries are probabilities on a 0--1 scale, not percentages.

\begin{center}
\begin{tabular}{@{}lrrrr@{}}
\toprule
Baseline VAF & Year 1 & Year 2 & Year 3 & Year 4\\\midrule
'''
def tex_number(value):
    if 0 < abs(value) < 0.001:
        mantissa,exponent=f'{value:.2e}'.split('e')
        return f'${float(mantissa):g}\\times10^{{{int(exponent)}}}$'
    return f'${value:.3g}$'


for row in risks:
    values=' & '.join(tex_number(float(row[f'year_{year}_probability'])) for year in range(1,5))
    text+=f"{float(row['baseline_vaf_percent']):g}\\% & {values} \\\\\n"
text+=r'''\bottomrule
\end{tabular}
\end{center}
Standalone risks use adaptive banded LSODA with explicit absorbed flux and a tenfold tolerance comparison. The maximum absolute difference across these demonstrations is below $3\times10^{-10}$. This numerical agreement does not resolve the clinical data and model assumptions.

\textbf{Remaining clinical requirement.} Verify index infusion, first diagnosis or event-free censor date, death information and ascertainment. The readiness check currently finds no analysis patient with a complete verified clinical record. The shared-state likelihood implementation is available, but the clinical optimizer, full refitted clinical sensitivity study, uncertainty intervals and final clinical 1--4-year estimates have not been run. See the CSV reconciliation table and machine-readable readiness/status reports.
'''
(ROOT/'initial_status.tex').write_text(text)
print('Wrote initial_status.tex')
