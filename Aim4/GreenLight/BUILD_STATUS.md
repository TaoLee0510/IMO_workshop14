# GreenLight · Aim 4 copy — build status

Updated 2026-11-03 · copy of `Aim2/GreenLight` v0.5 → v0.6 · Aim 4 → **v0.7 · Aim 4 only** (Aim 2 copy untouched).

## v0.7 (2026-11-03, evening)
- **Aim 2 removed**: `js/modules/marrow.js`, the Aim 2 / Aim 4 switch (Integration page now shows a read-only
  *Marrow model* card), the Aim 2-only tools (`test_pipeline.js`, `test_demo_fit.js`, `phenomap.js`) and every
  visible Aim 2 string. `GL.plans.marrowOut` reads `marrow4` only.
- **Model aligned with Aim 4 v1.2**: `aim4_core.js` ports the 24-state lineage tree and the age / generation
  moments (`opt.track`); `marrow4.js` forecasts ALC (lymphoid branch, not fitted to ALC) and reports the mean
  time and divisions since MPP entry of circulating ANC, PLT, AMC, ALC (table in the "when · why" panel, tile
  for neutrophils); help texts of the model in 7 languages describe v1.2; intervals labelled as the
  10th–90th percentile (what `GL.smc.val` returns).
- **Demo cohort regenerated with the Aim 4 model** (`tools/build_demo_data.js`): 16 fully synthetic patients,
  one mechanism each; recorded phenotypes 5 quick / 5 second wave / 4 delayed / 2 persistent.
- **Theme chrome**: light theme with a light sidebar and top bar (green accent line), dark theme unchanged.
- Stored versions are cleared automatically when module version, population file or demo cohort change.
- Tests: `test_aim4_ode.js` PASS (max relative difference 2e-11 incl. time / divisions since MPP entry);
  `browser_test.js` 0 console errors/warnings (modules: tumour_immune, marrow4).

## Done
- `js/modules/aim4_core.js`: line-by-line port of `Aim4/Code/src/aim4_model.cpp` (21 states, parameter layout of
  `a4_param_names()`, same Dormand–Prince RK45 and step control), inputs as `make_segments()`, ferritin driver
  (Whittaker smoother of log ferritin/baseline + Gamma-kernel continuation), phenotype rule of `phenotype.R`
  (primary rule + persistent class), when/who classes of `attribution.R`.
- `js/modules/marrow4.js`: particle posterior (≤ 400 particles) over r_s, th_ap, th_as, th_bs, f0; covariate
  priors from `data/aim4_pop.js`; censored and G-CSF upper-bound counts enter as log Φ; forecasts with plan
  policies (LD dose before D−5, anti-inflammatory scaling of the driver after the decision day, reactive /
  prophylactic G-CSF, TPO-RA, CD34+ boost); outcomes p_quick/p_second_wave/p_delayed/p_persistent,
  p_a_driven/p_b_driven/p_ab_driven, p_clone/p_reserve, f_a, f_b, f_clone30, plus the utility inputs.
- Real cohort tab (35 records, `data/real_cohort.js`, de-identified banner), "When … · why" panel,
  marrow-model switch, 4-class phenotype bars everywhere, 7-language strings (`js/i18n/aim4_dict.js`).
- `Aim4/Code/R/an_G_export.R` (step G) + `tools/build_aim4_data.R` (same files without a fit).

## Tests (2026-11-03)
- `node tools/test_aim4_ode.js`: 4 C++ reference cases, max relative difference on ANC/PLT/AMC 1.0e-8
  (S and clone fraction also < 1e-8) — PASS.
- `node tools/test_demo_fit.js` (Aim 4) and `GL_MARROW=aim2 node tools/test_demo_fit.js` (Aim 2): both run.
- `node tools/test_aim4_fit.js`: real records at D0 / D+14 / D+28, ~1.2–2.3 s per version in node.
- `node tools/browser_test.js docs/screenshots`: full walkthrough incl. panel, real tab (35 rows), real record
  scrubbed D+28 → D+7, switch Aim 4 → Aim 2 → Aim 4: **0 console errors/warnings**.
- i18n: 7 languages key-complete (536 keys).

## Placeholders / assumptions
- `data/aim4_pop.js` holds **prior centres** until step G runs on the real fit (`source: "prior"` is shown
  on the panel and the Integration page).
- Anti-inflammatory plan factors on the driver (tocilizumab at grade 1 ×0.8, grade 3 ×1.1, steroids at
  grade 2 ×0.9, anakinra prophylaxis ×0.6) are illustrative.
- Bridging choices in plans map "chemo" to the oral-cyclophosphamide intensity unless the record is R-Gem-Ox;
  bridging is locked from D−21, so in practice the recorded bridging is used.
- Demo (synthetic) records lack eGFR, bridging regimen and AMC: defaults eGFR 90, bridging "chemo" → oral Cy,
  set points from the first counts; a CH panel, when present, centres the f0 prior on 2 × VAF.
- No CRS model in `marrow4`: the P(CRS ≥ 3) tile shows "–" when Aim 4 is active.
- Attribution in the browser uses 40 particles; the R engine's MCMC attribution is shown beside it for real
  records once step G has run with F2.
