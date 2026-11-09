/* GreenLight · modules/hematotox.js
 * Clinical score module: CAR-HEMATOTOX (Rejeski et al., Blood 2021) at lymphodepletion from
 * platelets, ANC, haemoglobin, CRP (mg/dL) and ferritin; high risk at a score >= 2.
 * Same function as car_hematotox() in Model2/Code/R/an_F3_landmark.R.
 */
(function (root) {
  const GL = root.GL = root.GL || {};
  function score(s) {
    if ([s.plt_ld, s.anc_ld, s.hgb_ld, s.crp_ld, s.ferritin_ld].some(x => x === null || x === undefined)) return null;
    return (s.plt_ld < 75 ? 2 : s.plt_ld <= 175 ? 1 : 0) + (s.anc_ld < 1.2 ? 1 : 0) + (s.hgb_ld < 9 ? 1 : 0) +
      (s.crp_ld >= 3 ? 1 : 0) + (s.ferritin_ld > 2000 ? 2 : s.ferritin_ld >= 650 ? 1 : 0);
  }
  GL.registry.register({
    id: 'car_hematotox', version: 'Rejeski 2021', type: 'score', nameKey: 'mod.hematotox', descKey: 'mod.hematotox.d',
    source: 'Rejeski K et al. Blood 2021;138:2499', card: { outputs: [{ id: 'score', labelKey: 'out.hx_score', unit: '0-7', direction: 'lower' }],
      inputs: ['plt_ld', 'anc_ld', 'hgb_ld', 'crp_ld', 'ferritin_ld'] },
    predict(snap) { const x = score(snap.summary || {}); return x === null ? null : { score: x, group: x >= 2 ? 'high' : 'low' }; },
    score
  });
})(typeof window !== 'undefined' ? window : global);
