# GreenLight · MarrowTwin (Model3) · Results bundle and data contract (schema `glm3-1.0`)

GreenLight is static (HTML/JS, runs from `file://`, no CDN). Everything it shows comes from a
**results bundle** written by step G (`Model3/Code/R/an_G_export.R`, `g_export(cfg)`):

```
data/model3_index.js          window.GLM3_INDEX = {...}         public; loaded at start (no real-patient data)
data/patients/DP<nn>.js       GLM3.addPatient({...})            public digital patients, loaded on demand
data/sealed/<fid>.js          GLM3.addSealed(fid, iv, ct)       real patients, AES-256-GCM (tools/seal.js)
data/users.js                 window.GLM3_USERS = {...}         accounts and wrapped keys (tools/users.js)
tools/ref_model3.json         reference trajectories (digital designs; 10 cases: the 6 Model2 cases + neutrophil chains 8/6, 5/1, 3/2, th_cp > 0, d_zpl != 0) for tools/test_model3_ode.js
tools/digital_patients.rds    the digital patients with their truth (R object)
```
Step G writes real patients only to `Model3/Results/real/greenlight_stage/`
(`clinical/<id>.json` = all fields, `research/<id>.json` = minimised, `cohort_research.json`); the
pseudonym key and the keystore stay in `Model3/Results/real/greenlight_admin/` (copied from Model2, which copied them from Model1: same pseudonyms). `tools/seal.js` encrypts
them into `data/sealed/`; the browser opens a payload only with a key wrapped for the signed-in user
(`docs/PRIVACY.md`). Clinical payloads have the patient-file layout below; research payloads keep
`id` (pseudonym), `summary` (bands), `base`, `obs`, `label` and per as-of day `L, n, map, chol, names,
lo, hi, p, rec, nadir, band, when, who, f` (no plans, top3, driver or kernel).

## Index (`GLM3_INDEX`)

| field | content |
|---|---|
| `schema`, `generated`, `model_name`, `model_version`, `model_label` | `glm3-1.0`, time stamp, `MarrowTwin`, `Model3 v1.0`, `MarrowTwin (Model3 v1.0)` (shown in the UI) |
| `population` | `source` (`prior` / `fit`), `model` (final model id), `struct` (natural scale), `mu`, `omega`, `beta`, `sigma`, `tau`, `base_s`, `re`, `clone_types` |
| `asof_days`, `grid_days` | as-of days (D−5..D+30 daily, then weekly to D+86); days of the bands |
| `plans`, `utility` | plan library and utility weights |
| `phenotype` | primary rule and the phenotype term (`w`, `s`, `eps`) |
| `params` | parameter table from `config.R` (struct: value, prior, estimated/fixed, prior SD, box; random effects: centre, ω, box, covariates) |
| `k1` | K1 tables (unit tests, calibrations, Aim4 comparison, stability) |
| `results` | assembled cluster fits found under `Results/real/hpc_runs` (BIC, rounds, converged; group `synthetic` for the S1–S4 fits, never ranked), V1 metrics, pre-registration flag, synthetic cohort phenotype mixes, `s`: step S tables (S1–S4 summaries and convergence, S5 visual predictive check; no patient-level rows) |
| `logistic` | coefficients of the example ML module |
| `real_cohort` | `{ n, sealed: true }` |
| `roster` | digital patients only: `id`, `type` (`digital`), `status` (`completed` / `active`), `horizon` (data horizon of an in-treatment patient), `summary`, `label`, `mechanism` (digital), `evolution` ([L, P quick, P second wave, P delayed] per as-of day) |

## Patient file

| field | content |
|---|---|
| `id`, `type`, `status`, `horizon` | as in the roster |
| `summary`, `base`, `design`, `cov` | baseline information, set points, inputs for make_segments, covariates |
| `obs` | per count: `[day, value, cens, keep, setpoint]` |
| `ferritin` | `[day, value]` |
| `label` | observed phenotype under the primary rule and the five variants |
| `truth` | digital patients only: individual parameters, ferritin kernel, latent phenotype |
| `asof[]` | per as-of day `L`: `map` (posterior mode of ψ), `chol` (lower Cholesky factor of the Laplace covariance, column-major), `names`, `lo` / `hi` (search box), `kern` (ferritin kernel median), `drv` (driver on a 1-day grid), `p` (P quick / second wave / delayed), `rec` (recovery day median, 10%, 90%, P(by D+30)), `nadir`, `band` (ANC/PLT/AMC/ALC: 10/50/90% on `grid_days`), `when`, `who`, `f` (stem-cell loss a / inflammation / CAR T signal, clone fraction D+30), `plans` (ranked), `top3` |

As-of results use only the data up to `L`, no G-CSF after `L`, ferritin after `L` from the kernel
population conditioned on the ferritin up to `L`, and **no phenotype term** (the phenotype is
what is being predicted).

## Regenerating

```r
# from Model3/
source("Code/config.R"); for (f in list.files("Code/R", full.names = TRUE)) source(f)
cfg <- make_config("real"); load_model3(cfg)
g_export(cfg, cores = 8, n_draws = 60L, n_cf = 24L)      # ~4 min on 8 cores (51 patients x 44 as-of days)
```
then `node GreenLight/tools/seal.js` (existing keys are kept, so accounts stay valid).
When Model3's final fit is copied back (`Model3/Results/real/hpc_runs/_chain/final` and its `fit.rds`), the same
call uses it and the top bar switches from "prior centres" to "fitted".
