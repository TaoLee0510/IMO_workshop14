/* GreenLight · adapters/data_rest.js — DataAdapter over a REST endpoint (stub).
 * Expected endpoints (see docs/DATA_API.md):
 *   GET {base}/patients                      -> [{patient_id, label, demographics, product, demo}]
 *   GET {base}/patients/{id}                 -> patient object (same shape as demo_patients.json entries)
 *   GET {base}/patients/{id}/observations?upto={day}
 *   GET {base}/patients/{id}/events?upto={day}
 *   GET {base}/patients/{id}/mutations
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  GL.adapters = GL.adapters || {};
  GL.adapters.RestDataAdapter = {
    id: 'rest', labelKey: 'adapter.rest', kind: 'rest', baseUrl: 'http://localhost:8787/api', token: '',
    headers() { const h = { 'Content-Type': 'application/json' }; if (this.token) h.Authorization = 'Bearer ' + this.token; return h; },
    async connect(baseUrl, token) {
      if (baseUrl) this.baseUrl = baseUrl; if (token !== undefined) this.token = token;
      const r = await fetch(this.baseUrl + '/health', { headers: this.headers() });
      if (!r.ok) throw new Error('HTTP ' + r.status);
      return { ok: true, ...(await r.json()) };
    },
    async listPatients() { return (await fetch(this.baseUrl + '/patients', { headers: this.headers() })).json(); },
    async getPatient(id) { return (await fetch(this.baseUrl + '/patients/' + encodeURIComponent(id), { headers: this.headers() })).json(); },
    async getObservations(id, uptoDay) { return (await fetch(this.baseUrl + '/patients/' + encodeURIComponent(id) + '/observations' + (uptoDay !== undefined ? '?upto=' + uptoDay : ''), { headers: this.headers() })).json(); },
    async getEvents(id, uptoDay) { return (await fetch(this.baseUrl + '/patients/' + encodeURIComponent(id) + '/events' + (uptoDay !== undefined ? '?upto=' + uptoDay : ''), { headers: this.headers() })).json(); },
    async getMutations(id) { return (await fetch(this.baseUrl + '/patients/' + encodeURIComponent(id) + '/mutations', { headers: this.headers() })).json(); }
  };
})(typeof window !== 'undefined' ? window : global);
