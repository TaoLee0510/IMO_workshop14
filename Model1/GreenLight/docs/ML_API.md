# GreenLight · Model1 · ML module contract

An ML module predicts the recovery phenotype from what is known **on a given day**. It is
registered like any module (`docs/MODULE_API.md`) with `type: 'ml'`.

## Input: patient snapshot as of day *L*

```js
{
  patient: 'P07',
  asOf: 14,                                   // day relative to the CAR T infusion
  obs: { anc: [[day, value, cens]], plt: [...], amc: [...], alc: [...] },  // keep == 1 rows, day <= asOf
  summary: { age, female, prior_lines, prior_asct, bridging, egfr, ferritin_ld, crp_ld, ldh_ld,
             hgb_ld, plt_ld, anc_ld, max_crs, ... },  // baseline information (no outcome)
  design: { t0, egfr, bridging, apheresis_day, ldh_c, gcsf: [[start, end]] }
}
```
Values in ×10⁹/L; `cens = 1` marks a left-censored value (below the detection limit). Pre-treatment
set-point rows are included (day < −5). No value after `asOf` is ever passed.

## Output

```js
{ p: { quick, second_wave, delayed } }        // preferred: the three phenotypes (sum 1)
{ p_nonquick, landmark? }                     // or: probability of a non-quick recovery
{ forecast: { days: [...], anc: { lo, med, hi } } }   // optional count forecast with an 80% band
```
The clinician dashboard shows P(not quick) of every module next to the digital twin; the
researcher portal lists the module and its card.

## Evaluation

Out-of-sample comparison uses the pre-registered V1 folds (`Results/real/prereg/V1_folds.csv`,
5 folds × 7 patients, stratified by phenotype): train on four folds, predict the fifth at the
landmarks D0, D7, D10, D14, D21; report AUROC, calibration and Brier score next to the digital twin
and CAR-HEMATOTOX. The example module `counts_logistic` is fitted in sample (interface
demonstration only); its V1 version is part of step F3 of the Model1 pipeline.

## Example

`js/modules/logistic.js`: ridge logistic regression of non-quick recovery on log ANC and log PLT
(last value in the 3 days before the landmark), one model per landmark; coefficients written by
`g_logistic_module()` in `Model1/Code/R/an_G_export.R` into `data/model1_index.js` (`logistic`).
