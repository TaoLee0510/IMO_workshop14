# GreenLight · Module API (model plug-in contract)

GreenLight's twin is an **orchestrator over model modules**. A module is any
object (local JavaScript) or REST service (R plumber, Python FastAPI, Julia…)
that implements the contract below. The UI never knows what is inside a module:
it renders whatever variables, parameters, outcomes and signals a module declares.

Two modules ship with the demo:

| id | file | status |
|---|---|---|
| `tumor_immune` | `js/modules/tumor_immune.js` | **Placeholder** for the tumour–immune ODE of the other Aim. Variable names are provisional. |
| `marrow4` | `js/modules/marrow4.js` (+ `aim4_core.js`) | **Aim 4 parent model v1.2** (who & when; time and divisions since MPP entry). The only marrow model. Population constants in `data/aim4_pop.js`. |

### Aim 4 population file (`data/aim4_pop.js`)
```js
window.GL_AIM4_POP = {
  source: 'prior' | 'fitted', data_source: 'real' | null, fitted_at, model: 'MabC_ddr',   // any Aim 4 model id
  struct: { alpha, gmax, gamma_m, k_mpp, sr_mpp, k_cmp, sr_cmp, k_gmp, sr_gmp, k_mep, sr_mep,
            k_clp, sr_clp, k_ly, sr_ly, fb_mpp, fb_cmp, fb_gmp, fb_mep, fb_clp, fb_ly, r_a, r_b, delta,
            mtt_n, mtt_mo, mtt_pl, c_n, c_mo, c_pl, c_l, gam_n, gam_mo, gam_pl, k_i, emax, ic50, hill,
            e_mo, e_pl, e_l, eta_syn, e_sr, th_al, k_e, beta_renal, u_br_rgo, u_br_cy, gcsf_fold, tpo_fold,
            d_n, d_mo, d_pl },                                          // natural scale (model v1.2; Aim4/Code/config.R)
  mu: { r_s, th_ap, th_as?, th_bs?, f0? },  omega: { ... },          // eta scale (log; logit for f0)
  beta: { r_s: { prior_lines_c, prior_asct, age_c }, f0: { age_c, prior_lines_c, prior_asct, log_rdw_c, mcv_c } },
  sigma: { anc, plt, amc },
  covariates: { age_ref, age_scale, prior_lines_median, log_rdw_mean, mcv_mean },
  kernel_pop: { mu: [log A, log tp], cov: [[..],[..]], h },            // ferritin Gamma-kernel population
  cd34_add }
```
Written by `an_G_export()` (step G, `Aim4/Code/R/an_G_export.R`) from `F1$fits[[F1$selected$final]]`.

## 1. Module object

```js
{
  meta:     { id, nameKey, owner, version, kind: 'local' | 'remote', placeholder?: bool, engine?: string },
  consumes: { variables: ['anc','plt','crp'],          // observed variables it fits
              planKeys:  ['ld.flu_frac','cytopenia.*'], // plan parameters it reads
              signals:   ['car_expansion'] },           // signals from upstream modules
  provides: { signals:   ['inflammation'],
              variables: ['anc','plt','crp'],           // variables it forecasts
              params:    [{ id, labelKey, unit?, log?, lo, hi, descKey? }],
              outcomes:  [{ id, labelKey, kind: 'prob'|'value', direction: 'higher'|'lower', unit? }] },
  async fit(snapshot, ctx)          -> state,
  async forecast(state, plan, ctx)  -> forecast
}
```

### `fit(snapshot, ctx) -> state`
* `snapshot` — the patient as knowable at `ctx.asOfDay` (see DATA_API.md: observations/events filtered, CRS known only after onset).
* `ctx` — `{ asOfDay, horizon, nParticles, seed, signals: { id: { grid, values } } }`.
* `state` — serialisable. Must contain `nObs`, `ess`, `params: { id: { median, lo, hi } }`. Anything else (particles, weights, MCMC draws) is opaque and passed back to `forecast`.

### `forecast(state, plan, ctx) -> forecast`
* `plan` — a treatment plan (section 3). Stages before `ctx.asOfDay` are locked: use the recorded history in `ctx.snapshot`, apply the plan's policies only to the future.
* returns
```js
{
  grid: [-6, -5, …, 120],                                  // daily grid
  variables: { anc: { q05, q25, q50, q75, q95 }, … },      // arrays aligned to grid (predictive quantiles)
  outcomes:  { p_persistent: { value }, sev_neut_days: { value, lo, hi }, … },
  signals:   { inflammation: { grid, values } },           // for downstream modules
  extra?:    {}
}
```

### Signals (module coupling)
Modules are run in dependency order (Kahn topological sort on `provides.signals` → `consumes.signals`).
The demo passes the tumour–immune module's posterior-median `car_expansion(t)` into the marrow
module as its inflammation driver `D(t)`; a module must fall back to a default when a consumed
signal is absent.

## 2. REST transport (`js/modules/remote.js`)

```
GET  {base}/meta       -> { meta, consumes, provides }
POST {base}/fit        { snapshot, ctx }          -> state
POST {base}/forecast   { state, plan, ctx }       -> forecast      (ctx includes snapshot)
```
Register from the Integration page or in code:
```js
const m = await GL.modules.makeRemote({ baseUrl: 'http://localhost:8788/module', token: '…' }).load();
GL.registry.register(m);
```
Stubs: `api/plumber_stub.R` (R), `api/fastapi_stub.py` (Python). Both return the
contract shapes with dummy numbers so the wiring can be tested before a real model exists.

## 3. Treatment plan

```js
{
  id, nameKey | name, custom: bool,
  stages: {
    bridging:  { type: 'as_recorded' | 'none' | 'chemo' | 'radiation' },
    ld:        { flu_frac: 1.0, cy_frac: 1.0 },                     // fraction of Flu 30 mg/m² ×3, Cy 500 mg/m² ×3
    infusion:  { car_dose_frac: 1.0 },
    crs:       { toci_grade: 2, steroid_grade: 3, anakinra_prophylaxis: false },
    cytopenia: { gcsf: 'none'|'reactive'|'prophylactic', gcsf_start_day, gcsf_anc_thr, gcsf_days,
                 cd34_backup: bool, cd34_boost: 'none'|'if_persistent', cd34_boost_day,
                 tpo_ra: bool, tpo_plt_thr, plt_tx_thr },
    followup:  { cbc_interval_days }
  }
}
```
Stage windows and locking rules: `GL.schema.stageWindows`, `GL.schema.stageLocked(stage, asOfDay)`.
The candidate library and the multi-objective utility live in `js/core/plans.js`:

```
score = w_eff · P(response)                                   (tumor_immune.p_response)
      − w_rec · [1.0 · P(persistent) + 0.6 · P(second wave)]  (marrow recovery phenotype)
      − w_tox · [0.5 · P(CRS≥3) + 0.3 · min(1, sev_neut_days/30) + 0.25 · P(PLT<20)]
      − burden                                                (small per-intervention penalties)
```
`GL.plans.utility(outcomes, plan, weights)` returns `{ score, components: { efficacy, recovery, toxicity, burden } }`;
`weights = { eff, rec, tox }` (defaults 1.0 / 1.2 / 1.0, `GL.plans.normWeights` fills missing keys).
Any module that provides `p_persistent` / `p_second_wave` feeds the recovery term; a module without them
simply contributes nothing to it. `GL.plans.phenoBest(rankedList)` reports which candidate is best on each
recovery pattern (used by the "Recovery phenotype by treatment plan" panel).

## 4. Version object (what the store keeps)

`GL.orchestrator.buildVersion(patient, asOfDay)` returns an immutable `TwinVersion`:
`{ id, patient_id, version, asOfDay, createdAt, dataSummary, modules: { id: { meta, coverage, nObs, ess, params, fitMs } },
   states, plans: { planId: { plan, forecasts: { moduleId: forecast }, outcomes, score, components } }, ranking, activePlanId, weights, signals }`.
`GL.orchestrator.whatIf(patient, version, plan)` evaluates a custom plan on a stored posterior;
`GL.orchestrator.evaluateForecast(prev, patient, newAsOfDay)` scores an older version against data that arrived later.
