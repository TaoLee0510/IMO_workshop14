# GreenLight · MarrowTwin (Model3) · Module API (model card contract)

Every model shown in GreenLight is a **module** that registers a **model card** with
`GL.registry.register(...)` (`js/core/registry.js`). The researcher portal (Model, Equations,
Symbols & parameters, Fit & accuracy, Interfaces) is generated from the cards, and the clinician
dashboard shows every module's prediction beside the digital twin. A new mechanistic model, an ML
predictor or a clinical score therefore plugs in without touching the UI.

Modules that ship with this GreenLight (MarrowTwin = Model3):

| id | file | type | what it gives |
|---|---|---|---|
| `marrowtwin` | `js/modules/model3.js` + `js/engine/model3_core.js` | mechanistic | the digital twin: precomputed as-of predictions (R engine, step G) and the in-browser supportive-care simulation ("what-if") (JS port of the 32-state Model3 engine, < 1e-8 vs C++; 7.7e-12 on the reference cases) |
| `car_hematotox` | `js/modules/hematotox.js` | score | CAR-HEMATOTOX at lymphodepletion (Rejeski 2021) |
| `counts_logistic` | `js/modules/logistic.js` | ml | example ML module (ridge logistic on log ANC / log PLT at the landmark; coefficients from step G) |
| remote | `js/modules/remote.js` | ml (or any) | a REST service registered on Researcher → Interfaces |

## 1. Module object

```js
GL.registry.register({
  id: 'my_model', version: '0.1', type: 'mechanistic' | 'ml' | 'score',
  nameKey: 'mod.my_model', descKey: 'mod.my_model.d',      // i18n keys (all 7 languages, tools/check_i18n.js)
  source: 'where the model lives (repo path, paper, URL)',
  card: {
    equations(spec) -> [{ id, titleKey, descKey, latex }],   // optional; LaTeX rendered with the vendored KaTeX;
                                                             // wrap switched-off terms in \textcolor{#a3aab0}{...}
    symbols: { param_id: ['LaTeX symbol', 'unit'] },         // optional; meanings are i18n keys sym.<param_id>
    outputs: [{ id, labelKey, unit, direction: 'higher' | 'lower' | 'neutral' }],
    inputs:  ['plt_ld', 'anc_ld', ...],                      // optional, documentation
    validation: () => ({ ... })                              // optional; metrics shown under Fit & accuracy
  },
  predict(snapshot, patient) -> { p: { quick, second_wave, delayed } } | { p_nonquick } | { score, group } | null,
  whatIf(patient, asofEntry, plan, opt) -> { p, plt20, anc, plt }   // optional (mechanistic models)
});
```

`snapshot` (built by `GL.registry.snapshot(patient, asOf)`) contains only the data up to the
as-of day: `{ patient, asOf, obs: { anc: [[day, value, cens]], plt, amc, alc }, summary, design }`.
Modules must not look at later data; the clinician view relies on this for the day scrubber.

Add the script tag to `index.html` after `js/core/registry.js` (modules section) and the i18n keys
`mod.<id>`, `mod.<id>.d` and the output labels to all seven `js/i18n/*.js` files.

## 2. Remote modules (REST)

Researcher → Interfaces → *Add a remote module* stores a base URL (`localStorage glm3.remotes`). The
service implements:

| method | path | body / response |
|---|---|---|
| `GET` | `/card` | `{ id, version, type, outputs: [...], symbols?: {...} }` |
| `POST` | `/predict` | body = snapshot (above) → `{ p: {quick, second_wave, delayed} }` or `{ p_nonquick }` or `{ score }` |

Stubs: `api/fastapi_stub.py` (Python) and `api/plumber_stub.R` (R). Remote calls are asynchronous
and never block the dashboard.

## 3. Mechanistic models: what the R side provides

The displayed as-of results of a mechanistic model are precomputed by its R engine and written as a
results bundle (`docs/DATA_API.md`). For MarrowTwin this is step G of Model3, `Model3/Code/R/an_G_export.R`:
`g_export(cfg)` uses the final fit (`Results/real/hpc_runs/_chain/final`) when it exists, otherwise
the prior centres (the top bar shows which).
