/* GreenLight · adapters/data_demo.js — DataAdapter over the embedded demo cohort (works from file://). */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.adapters = GL.adapters || {};
  GL.adapters.DemoDataAdapter = {
    id: 'demo', labelKey: 'adapter.demo', kind: 'embedded',
    async connect() { this.data = root.GL_DEMO_DATA; if (!this.data) throw new Error('demo data not loaded'); return { ok: true, patients: this.data.patients.length, generated: this.data.generated, source: this.data.source }; },
    async listPatients() { return this.data.patients.map(p => ({ patient_id: p.patient_id, label: p.label, demographics: p.demographics, product: p.product, demo: p.demo, nObs: p.observations.length })); },
    async getPatient(id) { const p = this.data.patients.find(x => x.patient_id === id); return p ? GL.util.deepClone(p) : null; },
    async getObservations(id, uptoDay) { const p = await this.getPatient(id); return p ? p.observations.filter(o => uptoDay === undefined || o.day <= uptoDay) : []; },
    async getEvents(id, uptoDay) { const p = await this.getPatient(id); return p ? p.events.filter(e => uptoDay === undefined || e.day_start <= uptoDay) : []; },
    async getMutations(id) { const p = await this.getPatient(id); return p ? p.mutations : null; }
  };
})(typeof window !== 'undefined' ? window : global);
