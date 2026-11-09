// Fit real-cohort patients with the Aim 4 module at several checkpoints (node); print the
// posterior, phenotype probabilities against the recorded phenotype, attribution and timing.
// Usage: node tools/test_aim4_fit.js [patient_id ...]   (default: a few real-cohort records)
const GL = require('./load_core');
const ids = process.argv.slice(2);
const pts = (global.GL_REAL_DATA ? global.GL_REAL_DATA.patients : []);
const pick = ids.length ? pts.filter(p => ids.includes(p.patient_id) || ids.includes(p.label)) : pts.filter(p => ['P05', 'P06', 'P10', 'P22', 'P29'].includes(p.label));
(async () => {
  console.log('population:', global.GL_AIM4_POP.source, global.GL_AIM4_POP.model);
  let hit = 0, n = 0;
  for (const pt of pick) {
    for (const asOf of [0, 14, 28]) {
      const v = await GL.orchestrator.buildVersion(pt, asOf, {});
      const m = v.modules.marrow4, soc = v.plans.soc.outcomes.marrow4, f = (x, d) => x.toFixed(d === undefined ? 2 : d);
      const pr = { quick: soc.p_quick.value, second_wave: soc.p_second_wave.value, delayed: soc.p_delayed.value, persistent: soc.p_persistent.value };
      const best = Object.keys(pr).sort((a, b) => pr[b] - pr[a])[0];
      if (asOf === 28) { n++; if (best === pt.recorded.phenotype_clinical) hit++; }
      console.log(`${pt.label} [${pt.recorded.phenotype_clinical}] D${asOf >= 0 ? '+' : ''}${asOf}: ${v.computeMs}ms nObs=${m.nObs} ess=${f(m.ess, 0)} r_s=${f(m.params.r_s.median, 3)} th_ap=${f(m.params.th_ap.median)} th_as=${f(m.params.th_as.median, 3)} th_bs=${f(m.params.th_bs.median, 4)} f0=${f(m.params.f0.median, 3)} | Q ${f(pr.quick)} SW ${f(pr.second_wave)} DEL ${f(pr.delayed)} PER ${f(pr.persistent)} | a ${f(soc.p_a_driven.value)} b ${f(soc.p_b_driven.value)} ab ${f(soc.p_ab_driven.value)} clone ${f(soc.p_clone.value)} reserve ${f(soc.p_reserve.value)} | top ${v.ranking[0]}`);
    }
  }
  console.log(`modal phenotype at D+28 matches the recorded one in ${hit}/${n}`);
})().catch(e => { console.error(e); process.exit(1); });
