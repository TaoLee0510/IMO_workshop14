/* GreenLight · adapters/data_demo.js — DataAdapter over the embedded demo cohort (works from file://). */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.adapters = GL.adapters || {};
  GL.adapters.DemoDataAdapter = {
    id: 'demo', labelKey: 'adapter.demo', kind: 'embedded',
    async connect() {
      const demo = root.GL_DEMO_DATA; if (!demo) throw new Error('demo data not loaded');
      const real = root.GL_REAL_DATA ? root.GL_REAL_DATA.patients : [];   // de-identified real cohort, local use only
      this.data = Object.assign({}, demo, { patients: demo.patients.concat(real) });
      return { ok: true, patients: this.data.patients.length, generated: demo.generated, source: demo.source + (real.length ? ' + real cohort (' + real.length + ', de-identified)' : '') };
    },
    async listPatients() { return this.data.patients.map(p => ({ patient_id: p.patient_id, label: p.label, demographics: p.demographics, product: p.product, demo: p.demo, nObs: p.observations.length })); },
    async getPatient(id) { const p = this.data.patients.find(x => x.patient_id === id); return p ? GL.util.deepClone(p) : null; },
    async getObservations(id, uptoDay) { const p = await this.getPatient(id); return p ? p.observations.filter(o => uptoDay === undefined || o.day <= uptoDay) : []; },
    async getEvents(id, uptoDay) { const p = await this.getPatient(id); return p ? p.events.filter(e => uptoDay === undefined || e.day_start <= uptoDay) : []; },
    async getMutations(id) { const p = await this.getPatient(id); return p ? p.mutations : null; }
  };
})(typeof window !== 'undefined' ? window : global);
