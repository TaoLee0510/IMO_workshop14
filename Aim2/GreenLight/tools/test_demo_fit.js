// Fit each demo patient at its initial checkpoint and a later one; compare with recorded phenotype.
const GL = require('./load_core'); const fs = require('fs'), path = require('path');
const data = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'data', 'demo_patients.json'), 'utf8'));
const which = process.argv[2];
(async () => {
  for (const pt of data.patients) {
    if (which && pt.patient_id !== which) continue;
    const cps = pt.demo.checkpoints;
    const days = which ? cps : [cps[pt.demo.initialCheckpoint], 90];
    let prev = null;
    for (const asOf of days) {
      const v = await GL.orchestrator.buildVersion(pt, asOf, {});
      const m = v.modules.marrow, soc = v.plans.soc.outcomes.marrow, ti = v.plans.soc.outcomes.tumor_immune;
      let skill = '';
      if (prev) { const e = GL.orchestrator.evaluateForecast(prev, pt, asOf); skill = ' | prev-forecast skill: ' + Object.entries(e).map(([k, s]) => `${k} rmsle=${s.rmsle.toFixed(2)} cov90=${s.coverage90.toFixed(2)}`).join(', '); }
      console.log(`${pt.patient_id} [${pt.recorded.phenotype_clinical}] asOf ${String(asOf).padStart(3)}: ${v.computeMs}ms nObs=${m.nObs} ess=${m.ess} R=${m.params.R.median.toFixed(2)} g=${m.params.g.median.toFixed(0)} q=1/${(1 / m.params.q.median).toFixed(0)} f0=${m.params.f0.median.toFixed(3)} N0=${m.params.N0.median.toFixed(1)} | quick=${soc.p_quick.value.toFixed(2)} sw=${soc.p_second_wave.value.toFixed(2)} pers=${soc.p_persistent.value.toFixed(2)} crs3=${soc.p_crs3.value.toFixed(2)} | resp=${ti.p_response.value.toFixed(2)} peak=${ti.peak_expansion.value.toFixed(1)} | top: ${v.ranking.slice(0, 3).join(', ')}${skill}`);
      prev = v;
    }
  }
})();
