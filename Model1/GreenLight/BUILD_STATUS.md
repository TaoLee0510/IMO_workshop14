# GreenLight · Model1 — build status

## 2026-11-04 17:46 · logo
- Logo replaced by `image.png` (user): `assets/logo_full.png` (top bar, sign-in, lock and home page), `assets/logo_full_dark.png`
  (dark theme: the dark "Green" letters recoloured light), `assets/logo_mark.png` (icon), `assets/favicon.png`. The old SVG
  beacon mark is no longer used. `node tools/browser_test.js`: 60/60, 0 console errors.


## Access control (2026-11-04, plan §13 v0.4, A1–A7)
- **Bundle:** real patients are no longer plaintext anywhere in GreenLight. Step G stages clinical and
  research payloads in `Model1/Results/real/greenlight_stage/` (research payload minimised at generation:
  HMAC pseudonym, age / lines / LD-lab bands, bridging class, no drugs, design, G-CSF windows, ferritin,
  plans or what-if inputs; aggregate cells 1–4 as "<5"); `tools/seal.js` encrypts 71 payloads
  (AES-256-GCM, per-payload keys, random file ids); keystore and pseudonym key in
  `Model1/Results/real/greenlight_admin/` (mode 600). Reference trajectories now use digital designs.
- **Accounts:** `tools/users.js` (add / role / assign / unassign / passwd / revoke / list; hidden password
  prompt); per-user RSA-OAEP key pair, private key under a PBKDF2-SHA256 (600k) KEK; removing an
  entitlement rotates the affected keys. No accounts exist yet (the steward creates them).
- **Front end:** sign-in before everything (+ guest = digital patients only); role routing (clinician:
  assigned real + digital patients, research tabs Model/Equations/Symbols/Interfaces; researcher: research
  portal incl. Data with the research view, clinician portal digital only; steward: `#/admin` user list +
  audit export only); WebCrypto decryption on `file://` (verified in Chrome); 15-min idle lock;
  sign-out/lock wipe decrypted data; watermark (user + time) over patient views; banner "Authorized use
  only – access is logged"; real-patient UI state and what-if scenarios memory-only; hash-chained audit log
  (IndexedDB) with export; legacy real-patient browser data purged; 65 new strings × 7 languages.
- **Tests:** `test_access.js` 18/18; `browser_test.js` 60/60 incl. guest, clinician A/B, researcher,
  steward, lock/unlock, logout, one screen at 1440×900 and 1920×1080, 0 console errors; temporary accounts
  removed afterwards. `docs/PRIVACY.md` (classification, role/field matrix, controls, limitations,
  deployment target).

Updated 2026-11-04 · built from a copy of `Aim4/GreenLight` v0.7 (Aim4 copy untouched) · plan
`Model1/Model1_plan.md` v0.3 section 9, decisions G1–G7 (G4: all content in 7 languages).

## Done
- **Home**: enlarged GreenLight beacon logo + wordmark, digital-twin claim, two entries (Clinician,
  Researcher), animated patient ⇄ digital-twin picture, cohort counts with type marks.
- **Clinician**: one-screen dashboard (1440×900 and 1920×1080 without vertical scroll); patient list
  with status colour + origin symbol and filters; summary and data available; four count charts with
  band, median, data points and fixed axes; ferritin driver; "how the twin learned" strip; forecast
  card, why (when / who), top-3 strategies, other predictors; pinned scrubber (drag, step, keys, play
  once / loop / ping-pong, speed); data horizon and "next data batch" for patients in treatment;
  What-if drawer (JS port of Model1, presets from the plan library, scenarios saved in IndexedDB).
- **Researcher**: Model (clickable diagram, equations and parameters per box), Equations (all nested
  models, switched-off terms grey, KaTeX vendored), Symbols & parameters (searchable, from config),
  Data (cohort, rules, coverage heat map), Fit & accuracy (K1 21/21, calibrations, Aim4 comparison,
  synthetic cohort mixes; BIC, V1, VPC, S1–S4, residuals pending), Twins vs truth (16 digital
  patients), Versions (two as-of days of one twin), Interfaces (module table, model-card JSON, remote
  module form, bundle info, docs).
- **Interfaces**: model-card registry (`js/core/registry.js`), ML-module contract (snapshot as of a
  day → P(phenotype)), example modules CAR-HEMATOTOX (score) and counts logistic (ML), remote REST
  adapter + stubs (`api/`), results bundle schema `glm1-1.0` (`docs/DATA_API.md`).
- **Precompute** (`Model1/Code/R/an_G_export.R`): 35 real + 16 digital patients × 44 as-of days
  (D−5..D+30 daily, weekly to D+86): MAP + Laplace without the phenotype term, 60 draws (bands,
  phenotype, recovery day), 24 draws for counterfactual when/who and the 7-plan library; 3.4 min on
  8 cores (4.6 min with staging); digital patients public, real patients sealed (see above).
- **i18n**: en, ja, de, fr, es, zh-CN, zh-TW; 621 keys each.

## Tests (2026-11-04)
- `node tools/test_model1_ode.js`: 4 reference cases on real designs (bridging, G-CSF, TPO-RA,
  CD34+ boost, synergy, route variants) — max relative difference 3.7e-11 — PASS.
- `node tools/check_i18n.js`: 7 languages 556/556 keys, 507 keys used by the code all present — PASS.
- `node tools/browser_test.js`: superseded by the access-control walkthrough above (60 checks).

## Placeholders / pending
- Population = **prior centres** (MabC_ddr) until the cluster fits are in; then rerun `g_export(cfg)`
  (it reads `Results/real/hpc_runs/_chain/final`). Twin forecasts of real patients are therefore
  population-driven early on and can disagree with the recorded phenotype.
- Fit & accuracy: BIC table, V1, VPC (S5), S1–S4 recovery and residuals show "pending".
- Plan library and utility weights are illustrative (not clinical guidance); anti-inflammatory
  effect = driver × 0.5; CD34+ boost = 30% of the stem-cell pool.
- The counts-logistic module is fitted in sample (interface demonstration only).
- Digital patients use synthetic covariates and ferritin kernels; the browser what-if uses the median
  ferritin kernel of the as-of day (R uses one kernel draw per posterior draw).
