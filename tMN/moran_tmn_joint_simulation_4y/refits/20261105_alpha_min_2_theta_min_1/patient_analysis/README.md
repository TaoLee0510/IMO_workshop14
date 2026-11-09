# Patient-specific figures at the frozen latest joint fit

The current point is alpha=15.8479717593, theta=1 day. The Gamma PDF
integrates to 1 over [0,infinity). Inputs and inference source are frozen
by SHA256 in analysis_manifest.json. Every patient has their own measured
baseline and initial C fraction = 2 x baseline VAF, with D=0 at baseline.
Negative baseline dates are evolved neutrally to CAR infusion before
post-CAR selection begins.

Execution order: write code, write and compile report.tex, then run:

    python compute.py --stage all --workers 3
    python plot_figures.py
    python build_report.py
    python validate_analysis.py
    node build_workbook.mjs
    pdflatex -interaction=nonstopmode -halt-on-error -output-directory output/pdf -jobname Patient_specific_Moran_analysis report.tex
    python validate_analysis.py --artifacts

The workbook requires the bundled @oai/artifact-tool Node package and its
operation marker. The local scientific Python uses
PYTHONPATH=/private/tmp/tmn_scientific_runtime. Required Python versions
are in requirements.txt; the native compiler is clang++ with C++17/O3.

Patient risks are forward-model first-passage probabilities based only
on baseline and event-free CAR day 0. Follow-up VAF and eventual outcome
are not used to set individual risk values. VAF diagnostics use the
shared-state sequential filter and are checked against the saved fit.

Figure files are in figures/, with four separate figure types per patient
(VAF, C/D trajectories, tMN risk, combined dashboard). Overview figures,
Gamma/fitness curves, likelihood diagnostics and the heatmap accompany them.
All figures use English labels, and all 200 paths contribute to means.
Daily integer W/C/D counts and seeds are in results/trajectories.
Result tables include CSV and an English workbook. The PDF is in output/pdf.

Prediction bands are conditional on the frozen parameters. They are not
parameter confidence intervals. These are in-sample model diagnostics and
conditional threshold predictions, not calibrated population incidence.
