# Fixed-area gamma Moran study

This is the user-confirmed four-year joint-fit and simulation study. The final
correction is **fixed area, never fixed peak**. With time in days,
`f_C(t)=1+area_days*GammaPDF(t; alpha, theta)` after CAR infusion, and
`f_D=f_C+s_D`. Before day 0 all relative fitnesses equal 1. The reference
advantage area is 1 day. `b=1/28` per day, `s_D=0.02`, `mu=1e-6` per C division.

The prescribed simulation states at day 0 are C=2,000 and D=1,000 or 0,
with W making up the fixed total 100,000. C excludes D. Full stochastic
simulations use all three types; inference explicitly uses a deterministic
C/(N-D) background and stochastic D approximation.

## Reproducible execution

After the code and protocol have been written, compile the protocol, then run:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
python run_study.py --stage all --workers 3
```

Resumable stages: `validate`, `fit-reference`, `fit`, `simulate-grid`,
`simulate`, `risks`, `profile`, and `report`. Completed fits and simulation jobs are
checkpointed. Rebuild the native binary automatically if its source changes.
Do not reuse checkpoints after changing the biological configuration; archive
or isolate the previous run first. The local execution uses bundled Python
with `PYTHONPATH=/private/tmp/tmn_scientific_runtime`; requirements list versions.
The compiler is clang++ with C++17 and O3. The report/workbook stage uses the
bundled Node and Artifact Tool through a local dependency symlink.

## Sources and clinical analysis

The supplied MOESM7/MOESM9 files and audited endpoint-matched longitudinal data
are read without modifying the originals. `source_data.py` records hashes and
builds the user-confirmed 14-person cohort. Patient 103 day 127 and Patient 197
day 589 enter as events. Patient 25 is censored at 1460 under the confirmed
assumption that its listed day-1561 event is the first diagnosis. Eleven
unknown outcomes contribute VAF only and are not treated as negatives.
Baseline-only patients with usable clinical outcomes are included.
Copies of the four hashed source inputs are included in `data/source_inputs`,
so rebuilding the cohort does not require the earlier study directory.

The likelihood applies VAF weights and clinical first-passage information on
the same stochastic D distribution. It is not a product of independent VAF
and event marginal likelihoods. Unknown outcomes use a reflecting numerical
D cap with cap refinement, rather than a clinical absorbing threshold. The
reference also checks the full unrestricted D=0..N chain. Exponential tilting
is a similarity transformation of the actual finite-N generator, used only
to resolve tiny event densities; no event probability floor is used.

All compatible sensitivity scenarios refit alpha/theta. Holding s_D fixed
while varying b changes derived b*s_D. The mu=0/D_init=0 scenario is
structurally incompatible with positive events. A separate VAF-only fit is
retained solely as a comparator to the joint fit.

## Simulations and figures

The nine alpha/theta combinations (2,4,8) x (2,7,21 days), each with two initial
D fractions and 200 replicates, produce 3,600 primary trajectories. Fitted
reference/sensitivity parameter sets get 200 replicates per initial D fraction.
The exact time-inhomogeneous Moran algorithm skips type-preserving events
analytically and uses a proven state-specific thinning bound. Daily integer
counts, replicate seeds, first-passage times and metadata are retained.
Trajectories continue after crossing. All-zero/extinct replicates are included
in means. Four main figures use blue for C and red for D; 200 faint curves and
a bold arithmetic mean appear in each panel.

Risk outputs use forward first-passage calculations rather than claiming rare
probability precision from 200 replicates. The D=1% scenario starts one cell
below the strict >1% threshold and describes existing-D behavior. It is not a
healthy baseline incidence cohort. Source selection, missing outcomes,
founder causality, diagnosis mapping and competing death remain uncalibrated.
Profile likelihoods quantify shape/scale support conditional on fixed
parameters; sensitivity ranges are not confidence intervals.

## Deliverables

- `main.tex` and `main.pdf`: English protocol and executed results.
- `results/joint_fits.csv`, `patient_likelihoods.csv`, `VAF_predictions.csv`.
- `results/four_year_risks.csv`, `likelihood_profiles.csv`.
- `results/simulation_summary.csv`, `simulation_metadata.json`.
- `results/trajectories/*.npz` and matching replicate CSV/metadata JSON.
- `results/figures/*.pdf` and PNG, with four main 3x3 trajectory grids.
- `outputs/20261105_fixed_area/Fixed_area_joint_Moran_study.xlsx`.
- `results/validation.json`, input hashes, numeric checks, logs and seeds.

Previous peak-normalized work remains separate and is superseded for the
current requested study.
