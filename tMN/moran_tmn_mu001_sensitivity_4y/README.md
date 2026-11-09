# Four-year joint CAR-T/tMN refit: mu = 0.01

The completed delivery contains 78 compatible joint fits, one structurally
incompatible scenario, 1092 patient/scenario risk records (four yearly values
each), 2800 reference Moran paths, and 82 PDF/PNG figure pairs. Main outputs are
`output/pdf/Mu001_joint_sensitivity_report.pdf` and
`outputs/20261105_mu001_sensitivity/Mu001_joint_sensitivity_results.xlsx`.
`results/validation/delivery.json` records the source, numerical and export checks.
The report states the model's limited VAF coverage and boundary optima explicitly.

The reference and sensitivity design were confirmed before execution. Scientific
outputs are in English. The preceding study is preserved in its original folder.

## Low-memory execution (resumed 5 November 2026)

The user increased the authorized concurrency to **three jobs** on 5 November
2026. The completed numerical run used `run_jobs.py --stage numerical --jobs 3`. It uses
only the standard library and limits the pool to three independent scenarios.
Each job has **one numerical thread**, with all five thread settings fixed to 1.
`isolated_numerics.py` runs each likelihood evaluation or patient's risk
calculation in a fresh process, which exits and returns its memory to the operating
system. Within a job the optimizer waits before starting the next evaluation.
Each scenario optimizer also exits after completion. `memory_control.py` collects
completed solver cycles and releases unused macOS allocator pages. These controls
do not change the model equations, likelihood, bounds or numerical tolerances.

Shared CSV/progress collection is protected by a file lock. JSON and CSV exports
use atomic replacement with process-specific temporary names. Each scenario has
its own log; worker resource checkpoints are stored under `results/workers/`.
`run_serial.py` remains available for the older one-job mode, and both controllers
share a lock to prevent simultaneous controller launches.

At the completed reference, the isolated likelihood reproduced joint NLL
455.8080763255657 exactly and its numerical process peaked at 177 MiB. This is a
measured reference test, not a guaranteed peak for every sensitivity case.

Completed scenario files are reused only after checking fixed parameters, fit
bounds, cohort hash, model hash and numerical validation status. Interrupted work
is restarted; previously evaluated likelihood values are never trusted as a
numerical cache for a new fit. Atomic checkpoint writes protect completed files.
Execution stops before a new task if free disk space is below 512 MiB.

The numerical command finishes all fits, then calculates scenario risks with the
same three-job limit. It does not start plotting or reporting while fitting or
risk calculation is active.

Progress and process identifiers:

- `results/serial_controller.json`: active controller and numerical worker.
- `results/progress.json`: saved, fitted and structurally incompatible counts.
- `results/runtime.json`: most recent numerical task and resource checkpoint.
- `parallel_resume.log`: starts and completions in the current three-job pool.
- `logs/fit_*.log` and `logs/risks_*.log`: detailed per-scenario output.
- `serial_resume.log`: historical output from the preceding one-job run.
- `results/resource_observations.jsonl`: peak memory and runtime of every completed
  isolated likelihood/patient-risk calculation.

To pause, send SIGTERM to the **controller PID recorded in this study's
`serial_controller.json`**. Its handler terminates only the active process group
it started for this study, including the numerical child, and preserves saved
scenario files. Do not use a system-wide Python kill.

## Scientific definition

The normalized Gamma PDF integrates to 1 over [0, infinity) for every alpha and
theta. The peak is derived. The legacy `area_days=1` field is a numerical fitness
multiplier; it does not mean a one-day Gamma density integral.

Shared fitted parameters are alpha in [2, 50] and theta in [1, 1460] days. Fixed
reference parameters are N = 100000 cells, mu = 0.01 per C division, D_init = 0
cells at measured baseline, b = 1/28 per day, s_D = 0.02, strict D/N > 0.01
(1001-cell barrier), zero diagnosis delay and observation SD = 0.25 on
logit(2 VAF). Before CAR all fitnesses equal 1; mutation remains active. After
CAR, f_W = 1, f_C = 1 + Gamma PDF and f_D = f_C + s_D.

The hybrid inference retains deterministic q=C/(N-D) and a stochastic finite-N
D chain. The VAF and clinical terms use the same latent state. Conditional
clinical log likelihood is a diagnostic decomposition, not a second independent
fit. The 14-patient cohort has 42 follow-up VAFs, 2 exact events, 1 censor record
and 11 unknown clinical outcomes. Baseline-only patients do not contribute
fabricated follow-up observations. The aggregate D threshold is a model endpoint,
not an independently calibrated clinical diagnosis model.

`config.json` records the approved sensitivity levels. There are 79 scenarios:
72 mu-by-initial-D cells and 7 additional turnover/fitness cases. Every compatible
case refits alpha and theta. The mu=0, D_init=0 corner is structurally incompatible
with observed events, so 78 cases can be fitted. Founder burden C+D=2N*baseline
VAF is held fixed when D_init changes.

## Files and downstream stages

- `results/fits/reference.json`: new reference parameter estimates and validation.
- `results/fits/*.json`: scenario fits and separate objective evaluation records.
- `results/joint_fits.csv`: collected scenario estimates.
- `results/patient_risks.csv`: all scenario patient year-1 to year-4 probabilities
  after risk consolidation.
- `patient_analysis/`: frozen new-reference fit, 2800 complete Moran simulations,
  patient VAF diagnostics, risk curves and figure pipeline.
- `results/validation/`: independent generator, domain, tolerance and boundary
  checks, plus the memory-control numerical equivalence check.
- `protocol.tex`: confirmed English protocol; includes `executed_results.tex`
  when the completed result narrative has been generated.

After numerical completion, run these stages **sequentially** with one numerical
thread and worker:

1. `patient_analysis/plot_figures.py --stage core`.
2. `patient_analysis/validate_analysis.py` (numerical/data checks, not the old
   standalone patient workbook/PDF artifact mode).
3. `plot_sensitivity.py`.
4. `export_results.py`.
5. Compile `protocol.tex`, render all PDF pages and inspect layout.
6. Create the main XLSX using `build_workbook.mjs` and the bundled Artifact Tool;
   render/review the eight sheets and reconcile exported data and formulas.

The delivered patient subfolder already contains the frozen fit, inputs and
sources, so it does not need to be re-created with `prepare_patient_analysis.py`.
For a rerun in a relocated folder, install `requirements.txt`, use a C++17 compiler
for the patient simulator, and retain the relative directory structure for TeX
figure inputs. The XLSX builder additionally requires the bundled Artifact Tool
runtime documented by the spreadsheet skill. `compare_simulations.py` performs
the separate full-process comparison without refitting. `validate_delivery.py`
reconciles outputs against their owning scenario sources. Original absolute paths
in manifests are provenance; calculation scripts resolve their local folders.

The main results, report and workbook are not marked complete until all compatible
scenarios and export checks pass. The unused `native_rhs.cpp`/`native_kernel.py`
prototype is not used for inference in this run. Exact Moran simulations use the
separate patient simulator in `patient_analysis/`.

Last execution-policy update: 2026-11-05, America/New_York.
