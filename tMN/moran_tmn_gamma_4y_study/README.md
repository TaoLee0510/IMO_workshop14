# Four-year gamma-fitness tMN study

The English protocol is `main.tex`; its compiled version is `main.pdf`.
`config.json` records the proposed reference assumptions adopted when the user
instructed implementation to begin on 5 November 2026. Values are assumptions,
not measured biological constants. Existing studies and source workbooks are
retained separately.

The primary specification removes all pre-CAR mutant fitness advantages and
uses a fixed-peak gamma pulse after infusion. Only shape alpha and scale theta
are free parameters. There is one shared pair across the cohort. Extra D
fitness begins at day 0. The modeled threshold is aggregate D exceeding 1%.

## Execution

From this directory, with Python 3.12 and `requirements.txt` installed:

```sh
python prepare_data.py
python run_initial.py --mode all --initial-sensitivity
python fit_joint.py
python write_status_appendix.py
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The initial run used the bundled Python runtime and the pre-existing scientific
packages at `/private/tmp/tmn_scientific_runtime` via `PYTHONPATH`; versions are
recorded in `results/run_status.json`. The paths are an execution convenience,
not a portable environment specification.

## What has started

- Source hashing and row-level longitudinal joins, with allele orientation,
  duplicates, compartment, timing and clinical-outcome audits.
- Baseline-only selection of one index marker per patient. PBL/PBMC are kept
  separate; cfDNA and tumor are not used in the primary blood-cell mapping.
- A hybrid forward filter in which VAF weights and tMN outcomes act on the
  same stochastic D distribution, with deterministic M/(N-D) background.
- Unknown outcomes use the unrestricted physical D process; they are not
  interpreted as right-censored or negative.
- Initial validations against an independent matrix-exponential joint
  likelihood, a small-N full Moran process, known synthetic mean-field
  parameters, and production-size risk tolerance refinements.
- A preliminary mean-field VAF-only diagnostic fit and eight individually
  refitted sensitivity scenarios. These are NOT clinical joint estimates.
- Standardized 1-4-year threshold-risk demonstrations, NOT clinical incidence.

The primary VAF diagnostic excludes unresolved multiple CAR time origins and
source-positive cases lacking a diagnosis time. Some PBMC records are called
PBL in the metadata; the raw mutation-table label is retained, and the mismatch
is recorded. No sample-day metadata conflicts are accepted. Missing filtered
mutation rows do not establish zero VAF. Positive-detection selection, baseline
measurement conditioning, unphased marker choice and fixed fitness assumptions
remain limitations.

## Clinical data readiness

`data/clinical_reconciliation.csv` contains one row per patient and no invented
negatives. Its empty verified fields require source adjudication. Review index
infusion, first tMN diagnosis, last confirmed event-free contact, death and
cohort ascertainment. A VAF sample date is not a censor date. Diagnoses beyond
four years must use the verified first-event/time-origin information.
Re-extraction writes `clinical_evidence_inventory.csv` and preserves an existing
reconciliation file. Source time-origin/exclusion rules must also be reviewed
when rebuilding a cohort after clinical adjudication.

`python fit_joint.py` checks these fields and writes
`results/clinical_readiness.json`. It does not run an optimizer.
After reconciliation, `python fit_joint.py --run` provides the initial joint
fitting interface. It refuses unresolved data and records results as pending
production validation. Death within four years requires an additional
competing-event model. The full joint uncertainty and sensitivity analyses
remain future study stages; initial diagnostic tables must not be relabeled
as completed joint analyses.

The numerical joint filter initially uses backward Euler and needs timestep
refinement. Standalone risk demonstrations use adaptive banded LSODA and
explicit absorption flux, with tolerance comparisons. Only one small-N
full-process comparison is an initial check, not broad validation of the
deterministic M approximation. No clinical joint parameter-recovery claim is
made by the noiseless mean-field initialization test.

## Main outputs

- `data/audit.json`: source and analysis eligibility audit.
- `data/longitudinal_matched.csv`: all matched longitudinal source observations.
- `data/index_marker_inventory.csv`: inclusion/exclusion by patient.
- `data/analysis_longitudinal.csv`: selected baseline and VAF observations.
- `data/clinical_reconciliation.csv`: clinical fields requiring verification.
- `results/validation.json`: initial validation results and their limits.
- `results/preliminary_vaf_fits.csv`: diagnostic fit and initial sensitivities.
- `results/preliminary_vaf_predictions.csv`: individual diagnostic residuals.
- `results/preliminary_likelihood_grid.csv`: shape/scale diagnostic grid.
- `results/preliminary_diagnostics.pdf`: standalone diagnostic figure.
- `results/standardized_risk_demonstration.csv`: explicitly uncalibrated risks.
- `results/clinical_readiness.json`: missing fields by analyzed patient.
- `results/run_status.json`: completed and pending study stages, versions, times.

Final clinical 1-4-year risks, patient uncertainty intervals and the complete
refitted joint sensitivity analysis are pending verified clinical follow-up.
