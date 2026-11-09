# Reference joint refit: alpha >= 2 and theta >= 1 day

This isolated refit changes only theta's lower search bound from 0.25 to
1 day relative to the immediately preceding alpha >= 2 fit. The complete
search intervals are alpha in [2, 50] and theta in [1, 1460] days.
The standard Gamma density integrates to 1 over [0, infinity), independently
of alpha and theta. Its multiplicative coefficient remains numerically 1.
The peak is derived and is not constrained to a fixed height.

All reference biological and observation parameters, longitudinal VAF inputs,
two exact tMN events and one censor record remain unchanged. The same
shared-state joint likelihood is used. No simulations or sensitivity fits
are rerun in this directory.

Run `python run_refit.py` with the scientific requirements installed.
Local execution uses the bundled Python runtime and
`PYTHONPATH=/private/tmp/tmn_scientific_runtime`.

The optimizer reads bounds from `config.json` and clips all starting points
to the allowed intervals. Checkpoint reuse requires matching bounds, fixed
parameters and cohort hash. The runner checks every recorded evaluation,
convergence, outcome inclusion and likelihood decomposition. Numerical cap,
tolerance and full-population-chain comparisons are retained.

Outputs: `results/fits/reference.json`, `results/joint_fits.csv`,
`results/comparison.csv`, `results/refit_summary.json`, patient-level
likelihoods and VAF predictions. The comparison includes the original fit
and the immediately preceding alpha >= 2 fit. Earlier runs remain unchanged.
