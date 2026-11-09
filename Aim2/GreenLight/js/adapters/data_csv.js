/* GreenLight · adapters/data_csv.js — DataAdapter over the Aim 2 CSV contract
 * (patients.csv, labs_long.csv, inflam_long.csv, events.csv, mutation_tested.csv, mutations_pre.csv).
 * Needs an HTTP server (fetch is blocked on file://). baseUrl e.g. '../TestData/'.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.adapters = GL.adapters || {};
  function parseCsv(text) {
    const lines = text.trim().split(/\r?\n/); const h = lines[0].split(',');
    return lines.slice(1).map(l => { const c = l.split(','); const o = {}; h.forEach((k, i) => { o[k] = c[i] === undefined ? '' : c[i]; }); return o; });
  }
  const num = (x) => (x === '' || x === undefined) ? null : Number(x);
  GL.adapters.CsvDataAdapter = {
    id: 'csv', labelKey: 'adapter.csv', kind: 'csv', baseUrl: '../TestData/',
    async connect(baseUrl) {
      if (baseUrl) this.baseUrl = baseUrl;
      const get = async (f) => parseCsv(await (await fetch(this.baseUrl + f)).text());
      this.patients = await get('patients.csv'); this.labs = await get('labs_long.csv'); this.inflam = await get('inflam_long.csv');
      this.events = await get('events.csv'); this.tested = await get('mutation_tested.csv'); this.muts = await get('mutations_pre.csv');
      return { ok: true, patients: this.patients.length, source: this.baseUrl };
    },
    async listPatients() { return this.patients.map(p => ({ patient_id: p.patient_id, label: p.patient_id, demographics: { age: Number(p.age), sex: p.sex }, product: p.product, demo: { status: 'completed', checkpoints: [-5, 0, 7, 14, 21, 28, 42, 60, 90], initialCheckpoint: 8 } })); },
    async getPatient(id) {
      const p = this.patients.find(x => x.patient_id === id); if (!p) return null;
      const obs = [];
      for (const r of this.labs.filter(x => x.patient_id === id)) for (const v of ['anc', 'plt', 'wbc', 'hgb', 'alc']) if (r[v] !== '') obs.push({ day: Number(r.day), variable: v, value: Number(r[v]), source: 'cbc' });
      for (const r of this.inflam.filter(x => x.patient_id === id)) for (const v of ['crp', 'ferritin']) if (r[v] !== '') obs.push({ day: Number(r.day), variable: v, value: Number(r[v]), source: 'chem' });
      const t = this.tested.find(x => x.patient_id === id);
      return {
        patient_id: id, label: id, demographics: { age: Number(p.age), sex: p.sex }, product: p.product, costim: p.costim,
        history: { prior_lines: Number(p.prior_lines), prior_asct: Number(p.prior_asct), prior_bendamustine: Number(p.prior_bendamustine), months_since_last_chemo: Number(p.months_since_last_chemo), marrow_involvement: Number(p.marrow_involvement) },
        bridging: p.bridging, ld: { flu_total: Number(p.ld_flu_total), cy_total: Number(p.ld_cy_total), dose_reduced: Number(p.ld_dose_reduced) }, car_dose_frac: 1,
        crs: num(p.crs_grade) > 0 ? { grade: Number(p.crs_grade), onset_day: num(p.crs_onset_day), end_day: num(p.crs_end_day) } : null,
        mutations: { tested: !!t, sample_day: t ? Number(t.sample_day) : null, variants: this.muts.filter(x => x.patient_id === id).map(m => ({ gene: m.gene, protein_change: m.protein_change, vaf: Number(m.vaf), depth: Number(m.depth) })) },
        observations: obs.sort((a, b) => a.day - b.day),
        events: this.events.filter(x => x.patient_id === id).map(e => ({ type: e.event, day_start: Number(e.start_day), day_end: Number(e.end_day), dose: num(e.dose) })),
        recorded: { relapse_day: num(p.relapse_day), death_day: num(p.death_day), last_followup_day: num(p.last_followup_day), phenotype_clinical: p.phenotype_clinical || null },
        demo: { status: 'completed', checkpoints: [-5, 0, 7, 14, 21, 28, 42, 60, 90], initialCheckpoint: 8, tag: 'csv', synthetic_variables: [] }
      };
    },
    async getObservations(id, uptoDay) { const p = await this.getPatient(id); return p ? p.observations.filter(o => uptoDay === undefined || o.day <= uptoDay) : []; },
    async getEvents(id, uptoDay) { const p = await this.getPatient(id); return p ? p.events.filter(e => uptoDay === undefined || e.day_start <= uptoDay) : []; },
    async getMutations(id) { const p = await this.getPatient(id); return p ? p.mutations : null; }
  };
})(typeof window !== 'undefined' ? window : global);
