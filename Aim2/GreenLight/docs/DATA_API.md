# GreenLight · Data API (data adapter contract)

A **DataAdapter** turns any source (embedded JSON, the Aim 2 CSV files, a REST service, an EHR
feed) into the patient object below. Three adapters ship with the demo (`js/adapters/`):

| id | source | notes |
|---|---|---|
| `demo` | `data/demo_patients.js` | 8 synthetic patients built from `Aim2/TestData` by `tools/build_demo_data.js`; runs from `file://` |
| `csv` | Aim 2 CSV contract | `patients.csv, labs_long.csv, inflam_long.csv, events.csv, mutation_tested.csv, mutations_pre.csv`; needs an HTTP server |
| `rest` | REST service | endpoints below |

## Adapter interface
```js
{ id, labelKey, kind,
  async connect(...args)            -> { ok, patients, source },
  async listPatients()              -> [{ patient_id, label, demographics, product, demo }],
  async getPatient(id)              -> patient,
  async getObservations(id, uptoDay), async getEvents(id, uptoDay), async getMutations(id) }
```

## Patient object
```js
{
  patient_id, label,
  demographics: { age, sex: 'M'|'F' }, product, costim: 'CD28'|'4-1BB',
  history: { prior_lines, prior_asct, prior_bendamustine, months_since_last_chemo, marrow_involvement },
  bridging: 'none'|'chemo'|'radiation',
  ld: { flu_total, cy_total, dose_reduced },        // mg/m² totals; standard = 90 / 1500
  car_dose_frac: 1,
  crs: { grade, onset_day, end_day } | null,        // revealed to the twin only once onset_day <= asOfDay
  mutations: { tested, sample_day, sample_timepoint, panel, variants: [{ gene, protein_change, vaf, depth }] },
  observations: [{ day, variable, value, unit?, source? }],   // LONG format — any variable
  events: [{ type, day_start, day_end, dose }],
  recorded: { relapse_day, death_day, last_followup_day, phenotype_clinical },  // never shown to the twin
  demo: { checkpoints: [-5, 0, 7, 14, 21, 28, 42, 60, 90], initialCheckpoint, tag, synthetic_variables }
}
```
Days are relative to infusion (day 0); lymphodepletion is days −5…−3. No calendar dates are stored.

### Variable dictionary (`GL.schema.variables`)
`anc, plt, hgb, wbc, alc` (×10⁹/L, g/dL) · `crp` (mg/L) · `ferritin` (ng/mL) · `car_copies` (×10³ copies/µg) ·
`tumor_burden` (mL MTV) · `ldh` (U/L). Add a variable by adding a dictionary entry and a translation key `var.<id>`.

### Event types
`tocilizumab, steroid, anakinra, gcsf, tpo_ra, plt_tx, rbc_tx, cd34_boost, cd34_backup, bridging, lymphodepletion, infusion`.

## Snapshot (what modules see)
`GL.schema.snapshot(patient, asOfDay)` keeps only observations with `day <= asOfDay`, events started by then
(ongoing windows clipped), and the CRS record only if its onset is in the past. The twin never sees
`recorded.*`.

## REST endpoints (`js/adapters/data_rest.js`)
```
GET {base}/health
GET {base}/patients
GET {base}/patients/{id}
GET {base}/patients/{id}/observations?upto={day}
GET {base}/patients/{id}/events?upto={day}
GET {base}/patients/{id}/mutations
```
Optional bearer token via the `Authorization` header.
