ENDPOINT-ONLY MORAN tMN STUDY
DNMT3A, TP53, PPM1D

Start with main.tex or main.pdf for the two-page mathematical specification.
approved_scenarios.json records the settings explicitly approved by the user.

INPUTS
The two original MOESM Excel files live in the parent directory. Only
MOESM9/Gene_Level_Expansion_Data supplies VAF observations. MOESM7/tMN_VAF
supplies the requested source-list flag; MOESM7/Auto_CAR_tMN supplies clinical
timing discrepancy flags. No longitudinal VAF observations are used.
data/endpoints.csv preserves all 130 variant trajectories and their source rows.
data/audit.json records source hashes, counts, and quality flags.

MAIN RESULTS
results/fit_reference.csv: all 85 patient/gene groups, including 18 with no
estimable positive diploid-compatible endpoint pair. There are 67 conditional
fits across 52 distinct patients, using 93 variant endpoint pairs.
results/fit_variant_predictions.csv: observed/predicted endpoints and residuals.
results/fit_gene_summary.csv and fit_gene_contrasts.csv: gene medians and
patient-cluster bootstrap intervals. All three reference gene contrasts include
zero. These comparisons apply to the estimable subset, not the complete cohort.
results/fit_all_scenarios.csv: reference, 13 one-at-a-time parameter changes,
and a separate exclusion sensitivity for endpoints after a listed diagnosis.
figures/fit_all_patient_gene_trajectories.pdf: all 85 patient/gene panels.
Only endpoint dots are data; lines are model interpolations, not observations.

FITTING INTERPRETATION
The primary stage fixes mutation and baseline double-hit burden at zero and
fits a nonnegative CAR selection amplitude per patient/gene. It conditions on
positive observed baseline VAF and uses a pseudo-Huber loss in marked-cell
log-odds. It is not a calibrated sequencing likelihood, a hierarchical posterior,
or a joint fit of transformation parameters. No read counts or zero-VAF
pseudocounts are invented. Zero endpoints and non-diploid observations remain
flagged and excluded from this conditional numerical objective.
Nineteen fitted groups hit s2=0; 22 have model-mismatch flags. Such a boundary
estimate does not demonstrate no CAR effect or validate the model for a
declining trajectory. Unphased variants share a patient/gene parameter but are
not assumed to be independent clones or summed into a single observed burden.
Gene bootstrap intervals account for patient clustering, not uncertainty in
fixed kinetics, baseline measurement, clonal structure, or sample composition.
The initial proposed hierarchical measurement model was not estimated because
observation/detection parameters are unavailable; no shrinkage or individual
posterior confidence intervals are claimed.

RISK INTERPRETATION
results/risk_standardized_scenarios.csv and risk_standardized_curves.csv use
gene-median s2, a common 1% founder-marked cell fraction at day 0, and all 80
approved mutation/fitness/initial-D scenarios per gene.
results/risk_patient_gene_scenarios.csv uses the largest-baseline eligible
variant for each fitted patient/gene as a prespecified index marker. It retains
baseline compatibility and model-mismatch flags and coupled mean-field
endpoint predictions/residuals under each scenario.
The event is FIRST crossing of aggregate D >= 1001, not an individual lineage
above threshold at the final day. Double hits may exist initially. The hybrid
solver approximates the stochastic single-hit background while retaining a
fixed-N absorbing double-hit birth/death process. See risk_validation.json for
timestep and exact small-N comparisons; these are not full-N accuracy guarantees.
The deterministic single-hit background omits founder extinction and variability,
which can matter for small founder populations over long follow-up.
An additional full-N exact Moran simulation with 1,500 replicates tests one
rare-founder scenario (M0=125, D0=0, s2=0.1, mu=1e-4, b*s_D=0.01/day).
At five years its first-passage estimate is 0.2753 (Monte Carlo 95% interval
0.2533-0.2985), compared with 0.2829 from the approximation. The approximation
also lies within the simulation intervals at one and two years. This is one
case-specific check, not an accuracy guarantee for the full scenario grid.
See results/risk_full_N_validation.json; its intervals quantify simulation
error, not patient-risk uncertainty.

Risk is conditional on being below threshold at the first observed day and on
the specified baseline D0. Patient curves are retrospective scenarios using
s2 fitted with the last VAF, not prospective baseline-only predictions. They
are not all jointly refitted to the endpoint: large coupled endpoint residuals
show scenarios that do not reproduce the observed final VAF at that fixed s2.
No observation-error cutoff was invented to accept/reject scenarios. Do not
interpret all scenario-grid probabilities as equally supported by the data.
Do not sum variant/gene risks: clonal overlap and malignant lineage are unknown.
No clinical event likelihood is fitted: absence from tMN_VAF is not verified
negative status, event-free censoring dates are unavailable, clinical infusion
origins are unresolved, and the 1% operational threshold is not clinically
calibrated. Scenario envelopes are not confidence intervals.

REPRODUCTION
Use Python 3.12+ and the pinned packages in requirements.txt. Run from this folder:

  python prepare_endpoints.py
  python fit_endpoints.py --workers 2
  python risk_scenarios.py --mode validate
  python risk_scenarios.py --mode standardized --workers 3
  python risk_scenarios.py --mode patients --workers 3
  python risk_scenarios.py --mode endpoint-check --workers 3
  python risk_scenarios.py --mode audit
  pdflatex -interaction=nonstopmode -halt-on-error main.tex

Optional full-N stochastic validation requires a C++17 compiler:

  python results/risk_full_N_validate.py

The ZIP includes source code and results but excludes compiled executables,
temporary rendering files, and the original Excel inputs. Keep the two original
MOESM files in the parent directory to rerun extraction; the included extracted
endpoint table is sufficient to rerun fitting and scenario calculations.

The original run used the bundled Python runtime with additional scientific
packages installed in /private/tmp/tmn_scientific_runtime. That temporary path
is not required when the requirements are installed in your own environment.

SOURCE CONTEXT
Study: https://www.nature.com/articles/s41588-026-02795-5
CAR peak timing: https://doi.org/10.3390/biology13050350
CAR multiphasic kinetics: https://pubmed.ncbi.nlm.nih.gov/39240498/
Only the approximate CAR peak window has direct exposure-based support here;
the other fixed quantities are approved modeling scenarios, not measured
progenitor-selection parameters. CAR exposure is not itself progenitor fitness.
