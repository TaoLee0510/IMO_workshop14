# Reference joint refit with alpha lower bound 2

This isolated refit changes alpha's search interval from [1.01, 50] to
[2, 50]. Theta remains in [0.25, 1460] days. The integrated additional
fitness advantage remains fixed at 1 day; peak height is derived.
All reference biological and observation parameters and all cohort inputs
are unchanged. The two exact tMN events and one censor record remain in the
same shared-state joint likelihood as longitudinal VAF.

Run `python run_refit.py` with the scientific requirements installed.
The local execution uses the bundled Python runtime and
`PYTHONPATH=/private/tmp/tmn_scientific_runtime`.

`fitting.py` reads the physical bounds from this directory's `config.json`.
Starting points are clipped to those bounds. Checkpoint reuse requires the
same bounds, fixed parameters and cohort hash. The runner checks every
recorded evaluation against the requested bounds, optimizer convergence,
clinical inclusion and the likelihood decomposition. Numerical state-cap,
solver-tolerance and full-population-chain checks are retained.

Results are in `results/fits/reference.json`, `results/joint_fits.csv`,
`results/comparison.csv` and `results/refit_summary.json`. The preceding
reference result is retained in `data/previous_reference.json`, with original
source paths and hashes in `provenance.json`. The completed study's original
PDF, workbook, simulations and sensitivity results remain the preceding run;
this directory contains a reference fit only.
