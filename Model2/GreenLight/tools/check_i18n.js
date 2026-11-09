/* GreenLight · tools/check_i18n.js
 * i18n completeness: (1) every language has exactly the keys of en (no missing, no extra, no
 * empty values, same {placeholders}); (2) every literal t('key') in js/ and every dynamic key
 * family used by the UI exists in en.
 * Usage: node tools/check_i18n.js
 */
const fs = require('fs'), path = require('path');
global.window = global;
const dir = path.join(__dirname, '..', 'js', 'i18n');
const langs = ['en', 'ja', 'de', 'fr', 'es', 'zh-CN', 'zh-TW'];
langs.forEach(l => require(path.join(dir, l + '.js')));
const D = global.GL.I18N, en = D.en, keys = Object.keys(en);
let bad = 0;
const ph = (s) => (String(s).match(/\{[a-z0-9_]+\}/g) || []).sort().join(',');
for (const l of langs.slice(1)) {
  const d = D[l] || {}, miss = keys.filter(k => !(k in d)), extra = Object.keys(d).filter(k => !(k in en));
  const empty = keys.filter(k => k in d && !String(d[k]).trim()), phd = keys.filter(k => k in d && ph(d[k]) !== ph(en[k]));
  console.log(l.padEnd(6), `keys ${Object.keys(d).length}/${keys.length}`, miss.length ? `MISSING ${miss.length}: ${miss.slice(0, 8).join(', ')}` : '', extra.length ? `EXTRA ${extra.length}: ${extra.slice(0, 5).join(', ')}` : '', empty.length ? `EMPTY ${empty.length}` : '', phd.length ? `PLACEHOLDER ${phd.slice(0, 5).join(', ')}` : '');
  bad += miss.length + extra.length + empty.length + phd.length;
}
// literal keys in the code
const src = [];
(function walk(p) { for (const f of fs.readdirSync(p)) { const q = path.join(p, f); if (fs.statSync(q).isDirectory()) { if (!/i18n/.test(q)) walk(q); } else if (q.endsWith('.js')) src.push(fs.readFileSync(q, 'utf8')); } })(path.join(__dirname, '..', 'js'));
const lit = new Set(); src.forEach(s => { for (const m of s.matchAll(/\bt\('([a-zA-Z0-9_.\-]+)'/g)) lit.add(m[1]); for (const m of s.matchAll(/(?:titleKey|descKey|nameKey|labelKey): '([a-zA-Z0-9_.\-]+)'/g)) lit.add(m[1]); });
// dynamic key families
const fam = { 'ph.': ['quick', 'second_wave', 'delayed'], 'plan.': ['standard', 'gcsf_early', 'gcsf_reactive', 'anti_inflam', 'tpo_ra', 'cd34_boost', 'combo'],
  'when.': ['a', 'b', 'ab_joint', 'ab_redundant', 'other', 'none'], 'who.': ['reserve', 'clone', 'both', 'other'], 'type.': ['active', 'completed', 'real', 'digital'],
  'filter.': ['all', 'active', 'completed', 'digital'], 'lin.': ['anc', 'plt', 'amc', 'alc', 'anc.u', 'plt.u', 'amc.u', 'alc.u'], 'br.': ['none', 'cy', 'rgo'], 'hx.': ['high', 'low'],
  'mtype.': ['mechanistic', 'score', 'ml'], 'kind.': ['individual', 'estimated', 'fixed'], 'src.': ['calib', 'aim4', 'lit', 'scenario', 'fixed', 'prior', 'mstep', 'plan', 'mpp_equal'],
  'rt.': ['model', 'equations', 'symbols', 'data', 'accuracy', 'twins', 'versions', 'interfaces'], 'pop.': ['prior', 'fit'], 'acc.gen.': ['model1', 'aim4'],
  'rd.rule.': ['blank', 'manual', 'gcsf', 'plt', 'follow', 'setpoint', 'phenotype', 'nottarget', 'research'], 'wi.k.': ['gcsf', 'reactive', 'anti', 'tpo', 'cd34'],
  'eqsel.q2.': ['M0', 'Ma', 'Mb', 'Mab', 'Mabx'], 'eqsel.clone.': ['none', 'ddr', 'dta'], 'eqsel.route.': ['both', 'hsc', 'mpp'],
  'mech.': ['typical', 'car_t_phase_second_wave', 'low_reserve', 'clone_ddr', 'chemo_hsc_injury', 'car_t_target_route', 'typical_strong_reserve', 'inflammation_mixed', 'low_reserve_delayed', 'clone_ddr_chemo', 'chemo_and_car_t', 'strong_inflammation', 'quick_young'],
  'rm.': ['k1', 'k1.t', 'k2', 'k2.t', 'k3', 'k3.t'], 'role.': ['clinician', 'researcher', 'steward', 'guest'],
  'rd.lvl.': ['chemo', 'none', '1-2', '3-4', '>=5'], 'rd.xt.': ['bridging_class', 'lines_band'], 'land.scope.': ['assigned', 'digital', 'research', 'model'] };
const nodes = ['in_c', 'in_i', 'in_z', 'in_g', 'hsc', 'mpp_n', 'mpp_mo', 'mpp_pl', 'mpp_l', 'ch_n', 'ch_mo', 'ch_pl', 'b_n', 'b_mo', 'b_pl', 'b_l', 'lik'];
fam['node.'] = nodes.flatMap(n => [n, n + '.s', n + '.d']);
const syms = ['r_s', 'f0', 'th_as', 'th_ap', 'th_bs', 'rho', 'r_a', 'r_b', 'delta', 'eta_syn', 'th_zs', 'k_p', 'sr_p', 'gmax', 'gamma', 'kappa_c', 'th_am', 'mtt_n', 'mtt_mo', 'mtt_pl', 'c_n', 'c_mo', 'c_pl', 'c_l', 'th_al', 'k_i', 'emax', 'ic50', 'hill', 'k_e', 'beta_renal', 'u_br_rgo', 'u_br_cy', 'gcsf_fold', 'gcsf_rel', 'k_g', 'tpo_fold', 'car_t0', 'k_ag', 'kappa', 'k_t', 'beta_a', 'k_z', 'z50', 'e_z', 'sigma', 'tau'];
fam['sym.'] = syms;
const help = ['population', 'asof', 'data_points', 'twin_forecast', 'nonquick', 'rec_day', 'why', 'plans', 'modules', 'ferritin', 'evolution', 'coverage', 'bic', 'lin.anc', 'lin.plt', 'lin.amc', 'lin.alc',
  'small_cells', 'research_view', 'plan.standard', 'plan.gcsf_early', 'plan.gcsf_reactive', 'plan.anti_inflam', 'plan.tpo_ra', 'plan.cd34_boost', 'plan.combo', 'wi.wGcsf', 'wi.wReactive', 'wi.wAnti', 'wi.wTpo', 'wi.wCd34'];
fam['help.'] = help.flatMap(h => [h, h + '.t']);
const need = new Set([...lit].filter(k => !k.endsWith('.')));
for (const [p, ks] of Object.entries(fam)) ks.forEach(k => need.add(p + k));
// every parameter of the public index (config of the exported model) has a meaning (sym.<id>) and a symbol / unit
global.GL.registry = { register: (m) => { global.GL.mod = m; } }; global.GL.loader = { index: () => global.GLM2_INDEX };
require(path.join(__dirname, '..', 'data', 'model2_index.js')); require(path.join(__dirname, '..', 'js', 'engine', 'model2_core.js')); require(path.join(__dirname, '..', 'js', 'modules', 'model2.js'));
const pids = global.GLM2_INDEX.params.struct.map(r => r.id).concat(global.GLM2_INDEX.params.re.map(r => r.id));
pids.forEach(k => need.add('sym.' + k));
const noSym = pids.filter(k => !global.GL.mod.card.symbols[k]);
console.log(`parameters in the index: ${pids.length}; without symbol/unit in the model card: ${noSym.length}${noSym.length ? ' -> ' + noSym.join(', ') : ''}`);
bad += noSym.length;
const missEn = [...need].filter(k => !(k in en));
console.log(`en: ${keys.length} keys; literal + dynamic keys used: ${need.size}; missing in en: ${missEn.length}${missEn.length ? ' -> ' + missEn.join(', ') : ''}`);
bad += missEn.length;
console.log(bad ? 'FAIL' : 'PASS: all 7 languages key-complete');
process.exit(bad ? 1 : 0);
