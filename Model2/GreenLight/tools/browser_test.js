/* GreenLight · tools/browser_test.js
 * Headless-Chrome walkthrough over the DevTools protocol, with access control (Model1 plan 13.6 item 11, carried over by Model2 plan 6.5).
 * Temporary accounts (random names, random passwords kept in memory and given to tools/users.js on stdin)
 * are created at the start and always removed at the end. Flow: login screen -> wrong password -> guest
 * (both portals, digital patients only: filters, scrub, ping-pong, ingest, what-if, every allowed research tab,
 * languages, dark theme) -> clinician A (assigned real patients only, one screen at 1440x900 and 1920x1080,
 * watermark, real what-if kept in memory, audit entries, idle lock + unlock, logout clears data) ->
 * clinician B (only B's patient) -> researcher (research view: pseudonyms, bands, no recommendations or
 * what-if) -> data steward (user list, audit chain, export). Checks 0 console errors / warnings.
 * MarrowTwin (Model2) checks: public name and version label, no "Model1" in any view, the two chemotherapy kill
 * terms in the equations of every nested model and in the diagram, kappa_c / th_aM in the symbol table, glm2.*
 * browser storage. Temporary accounts left by an aborted earlier run (tb_<role>_<6 hex>) are removed first;
 * the existing accounts must be unchanged at the end.
 * Screenshots without real-patient data go to docs/screenshots; screenshots of real-patient views go to
 * Results/real/greenlight_stage/test_screenshots (outside GreenLight). Prints no patient-level data.
 * Usage: node tools/browser_test.js [chromePath]   (kills a stray Chrome on port 9333 first)
 */
const { spawn, execSync, execFileSync } = require('child_process'); const fs = require('fs'); const path = require('path'); const http = require('http'); const crypto = require('crypto');
const A = require('./lib_access.js');
const PUB = path.join(__dirname, '..', 'docs', 'screenshots'), PRIV = path.join(A.P.stage, 'test_screenshots');
fs.mkdirSync(PUB, { recursive: true }); fs.mkdirSync(PRIV, { recursive: true });
const chrome = process.argv[2] || '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const port = 9333, profile = '/tmp/glm2-chrome-cdp-' + Date.now();
try { execSync(`pkill -f "remote-debugging-port=${port}"`); } catch (e) { /* none */ }
const url = 'file://' + path.resolve(__dirname, '..', 'index.html');
const sleep = (ms) => new Promise(r => setTimeout(r, ms));
const getJson = (u) => new Promise((res, rej) => http.get(u, r => { let d = ''; r.on('data', c => d += c); r.on('end', () => res(JSON.parse(d))); }).on('error', rej));
const results = []; const check = (name, ok, info) => { results.push({ name, ok }); console.log((ok ? 'OK  ' : 'FAIL') + ' ' + name + (info !== undefined ? ' — ' + info : '')); };

// ---- temporary accounts ----------------------------------------------------------------------------------
const tag = crypto.randomBytes(3).toString('hex'), pw = () => crypto.randomBytes(18).toString('base64');
const U = { clinA: { u: 'tb_clinA_' + tag, pw: pw() }, clinB: { u: 'tb_clinB_' + tag, pw: pw() }, res: { u: 'tb_res_' + tag, pw: pw() }, stew: { u: 'tb_stew_' + tag, pw: pw() } };
const cli = (args, input) => execFileSync(process.execPath, [path.join(__dirname, 'users.js')].concat(args), { input: input === undefined ? '' : input + '\n', stdio: ['pipe', 'pipe', 'pipe'] });
const clinIds = Object.keys(A.loadKeystore().payloads).filter(p => p.startsWith('clinical:')).map(p => p.slice(9)).sort();
const [pa1, pa2, pb] = clinIds;
const regimens = [...new Set(clinIds.map(id => JSON.parse(fs.readFileSync(path.join(A.P.stage, 'clinical', id + '.json'), 'utf8')).summary.bridging_regimen).filter(r => r && r.length > 3 && r !== 'none'))];
const created = [];
const TEMP = /^tb_(clinA|clinB|res|stew)_[0-9a-f]{6}$/;
// the two Model2 kill terms as the model card writes them (teal = NEW in js/modules/model2.js)
const TEX = { kc: String.raw`\kappa_c\,\theta_{aP}\,C`, am: String.raw`\theta_{aM}\,C\,B_{\mathrm{Mo}}` };
TEX.kcNew = String.raw`\textcolor{#0e9fb5}{-` + TEX.kc + '}'; TEX.amNew = String.raw`\textcolor{#0e9fb5}{-\mathbf{1}[\ell=\mathrm{Mo}]\,` + TEX.am + '}';
const stale = A.loadUsers().users.map(x => x.u).filter(u => TEMP.test(u));
stale.forEach(u => cli(['revoke', u]));
if (stale.length) console.log(`removed ${stale.length} temporary account(s) left by an aborted earlier run`);
const creds = () => JSON.stringify(A.loadUsers().users.filter(x => !TEMP.test(x.u)).map(x => [x.u, x.role, x.salt, x.iter, x.ekey, x.pub]));
const credsBefore = creds();
// interrupted (Ctrl-C, timeout): still remove the temporary accounts
for (const sig of ['SIGINT', 'SIGTERM']) process.on(sig, () => { if (proc) proc.kill(); for (const u of created) { try { cli(['revoke', u]); } catch (e) { /* listed by the next run */ } } process.exit(130); });

let proc = null;
(async () => {
  cli(['add', U.clinA.u, '--role', 'clinician', '--patients', pa1 + ',' + pa2, '--password-stdin'], U.clinA.pw); created.push(U.clinA.u);
  cli(['add', U.clinB.u, '--role', 'clinician', '--patients', pb, '--password-stdin'], U.clinB.pw); created.push(U.clinB.u);
  cli(['add', U.res.u, '--role', 'researcher', '--password-stdin'], U.res.pw); created.push(U.res.u);
  cli(['add', U.stew.u, '--role', 'steward', '--password-stdin'], U.stew.pw); created.push(U.stew.u);

  proc = spawn(chrome, ['--headless=new', '--disable-gpu', '--no-sandbox', '--hide-scrollbars', '--window-size=1440,900', '--remote-debugging-port=' + port, '--user-data-dir=' + profile, 'about:blank'], { stdio: 'ignore' });
  let targets = null;
  for (let i = 0; i < 150 && !(targets && targets.length); i++) { try { targets = await getJson(`http://localhost:${port}/json/list`); } catch (e) { targets = null; } if (!(targets && targets.length)) await sleep(300); }
  const page = targets.find(t => t.type === 'page');
  const ws = new WebSocket(page.webSocketDebuggerUrl);
  await new Promise(r => { ws.onopen = r; });
  let id = 0; const pending = {}; const errors = [];
  ws.onmessage = (m) => { const msg = JSON.parse(m.data); if (msg.id && pending[msg.id]) { pending[msg.id](msg); delete pending[msg.id]; }
    if (msg.method === 'Runtime.exceptionThrown') errors.push('EXC: ' + (msg.params.exceptionDetails.exception && msg.params.exceptionDetails.exception.description || msg.params.exceptionDetails.text));
    if (msg.method === 'Runtime.consoleAPICalled' && (msg.params.type === 'error' || msg.params.type === 'warning')) errors.push(msg.params.type.toUpperCase() + ': ' + msg.params.args.map(a => a.value || a.description).join(' '));
    if (msg.method === 'Log.entryAdded' && msg.params.entry.level === 'error') errors.push('LOG: ' + msg.params.entry.text); };
  const send = (method, params) => new Promise(res => { const i = ++id; pending[i] = res; ws.send(JSON.stringify({ id: i, method, params: params || {} })); });
  const evalJs = async (expr) => { const r = await send('Runtime.evaluate', { expression: expr, awaitPromise: true, returnByValue: true }); if (r.result && r.result.exceptionDetails) errors.push('EVAL: ' + JSON.stringify(r.result.exceptionDetails.exception && r.result.exceptionDetails.exception.description)); return r.result && r.result.result ? r.result.result.value : undefined; };
  const viewport = async (w, h) => { await send('Emulation.setDeviceMetricsOverride', { width: w, height: h, deviceScaleFactor: 1, mobile: false }); await sleep(250); };
  const shot = async (dir, name, full) => { const r = await send('Page.captureScreenshot', { format: 'png', captureBeyondViewport: !!full }); fs.writeFileSync(path.join(dir, name + '.png'), Buffer.from(r.result.data, 'base64')); };
  const waitFor = async (expr, ms) => { const t0 = Date.now(); while (Date.now() - t0 < (ms || 20000)) { if (await evalJs(expr)) return true; await sleep(150); } console.log('TIMEOUT waiting for', expr.slice(0, 80)); return false; };
  const go = async (hash, ready) => { await evalJs(`location.hash = ${JSON.stringify(hash)}; true`); await sleep(300); if (ready) await waitFor(ready); await sleep(250); };
  const noScroll = async (label) => { const r = await evalJs('[document.documentElement.scrollHeight, window.innerHeight, document.body.scrollHeight]'); check(`no vertical scroll ${label} (${r[0]} <= ${r[1]})`, r[0] <= r[1] && r[2] <= r[1]); };
  const noOldName = async (label) => { const r = await evalJs('(() => { const h = document.title + document.body.innerHTML; const m = h.match(/Model ?1\\b|model1/i); return m ? m[0] : ""; })()'); check(`no "Model1" in the page ${label}`, !r, r || undefined); };
  const rawKeys = async (label) => { const n = await evalJs(`(() => { const m = document.body.innerText.match(/\\b(?:est|ph|plan|when|who|lin|scrub|roster|rt|rd|acc|sym|node|wi|help|mods|land|why|evo|fer|eq|eqd|ri|rv|rtw|mech|type|filter|pop|auth|role|admin)\\.[a-z0-9_]+(?:\\.[a-z0-9_]+)?\\b/g); return m ? m.slice(0, 5).join(', ') : ''; })()`); check(`no raw i18n keys ${label}`, !n, n || undefined); };
  const clinicReady = '!!document.querySelector("#cright .card") && !!document.querySelector("#ch_anc svg path.med")';
  const login = async (u) => {
    await evalJs(`(() => { const a = document.getElementById('lgUser'), b = document.getElementById('lgPw'); a.value = ${JSON.stringify(u.u)}; b.value = ${JSON.stringify(u.pw)}; document.getElementById('loginForm').requestSubmit(); return true; })()`);
    return waitFor(`GL.session.user === ${JSON.stringify(u.u)} && !!document.querySelector('#view > *')`, 30000);
  };
  const logout = async () => { await evalJs('GL.session.logout(); true'); await waitFor('!!document.getElementById("loginForm")'); };
  const realIdsInRoster = () => evalJs('Array.from(document.querySelectorAll(".ritem")).map(r => r.dataset.pid).filter(x => !/^DP/.test(x))');
  await send('Runtime.enable'); await send('Page.enable'); await send('Log.enable');
  await viewport(1440, 900);
  await send('Page.navigate', { url: url + '#/' }); await sleep(600);
  await waitFor('!!document.getElementById("loginForm")');

  // 1. login screen, wrong password
  check('login screen first: banner, sign-in form, guest button', await evalJs('!!document.querySelector(".authbar") && !!document.getElementById("lgPw") && document.getElementById("lgPw").type === "password" && !!document.getElementById("lgGuest")'));
  await shot(PUB, '00_login_1440');
  await evalJs(`(() => { document.getElementById('lgUser').value = ${JSON.stringify(U.clinA.u)}; document.getElementById('lgPw').value = ${JSON.stringify(U.clinA.pw + 'x')}; document.getElementById('loginForm').requestSubmit(); return true; })()`);
  await waitFor('document.getElementById("lgErr").textContent.length > 0 && !document.getElementById("lgGo").disabled', 20000);
  check('wrong password fails, nothing decrypted', await evalJs('!GL.session.user && Object.keys(GL.session.real).length === 0 && !!document.getElementById("loginForm")'));

  // 2. guest: digital patients only
  await evalJs('document.getElementById("lgGuest").click(); true'); await waitFor('!!document.querySelector(".entries")');
  check('guest: landing with both portals', await evalJs('document.querySelectorAll(".entry").length === 2'));
  check('MarrowTwin: landing lead, version chip "MarrowTwin (Model2 v1.0)", page title', await evalJs('document.body.innerText.includes(GL.i18n.t("land.lead")) && GL.i18n.t("land.lead").includes("MarrowTwin") && Array.from(document.querySelectorAll(".land-foot .chip")).some(c => c.innerText.includes("MarrowTwin (Model2 v1.0)")) && document.title.includes("MarrowTwin")'));
  await noOldName('landing');
  await shot(PUB, '01_landing_guest_1440');
  await go('#/clinic', clinicReady);
  check('guest: no real patient in the clinician roster (16 digital)', (await realIdsInRoster()).length === 0 && await evalJs('document.querySelectorAll(".ritem").length === 16'));
  const counts = {};
  for (const f of ['digital', 'active', 'completed', 'all']) { await evalJs(`document.querySelector('.ftabs [data-f="${f}"]').click(); true`); await sleep(120); counts[f] = await evalJs('document.querySelectorAll(".ritem").length'); }
  check('guest roster filters (16 / 8 / 8 / 16)', counts.digital === 16 && counts.active === 8 && counts.completed === 8 && counts.all === 16, JSON.stringify(counts));
  for (const pid of ['DP02', 'DP10']) {
    await go('#/clinic/' + pid, clinicReady); await noScroll('guest ' + pid + ' 1440x900');
    check(`${pid}: 4 charts, ferritin, evolution, scrubber, top-3`, await evalJs('document.querySelectorAll(".charts4 svg").length === 4 && !!document.querySelector("#ch_evo svg") && !!document.querySelector("#scSvg circle.tl-handle") && document.querySelectorAll(".plan-row").length === 3'));
    await shot(PUB, `02_clinic_${pid}_guest_1440`);
  }
  const day0 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.30, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`); await sleep(300);
  check('scrubber click changes the as-of day', (await evalJs('GL.views.clinic.ui.asOf')) !== day0);
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.62, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`);
  await evalJs('document.querySelector(\'[data-mode="pingpong"]\').click(); document.getElementById("scSpeed").value = "4"; document.getElementById("scSpeed").dispatchEvent(new Event("change")); window.__d = []; window.__iv = setInterval(() => window.__d.push(GL.views.clinic.ui.asOf), 100); document.getElementById("scPlay").click(); true');
  await sleep(4200); await evalJs('document.getElementById("scPlay").click(); clearInterval(window.__iv); true');
  const seq = await evalJs('window.__d');
  check('ping-pong playback runs back and forth', seq.some((d, i) => i && d < seq[i - 1]) && seq.some((d, i) => i && d > seq[i - 1]), seq.length + ' samples');
  await go('#/clinic/DP03', clinicReady);
  const h0 = await evalJs('GL.views.clinic.ui.horizon');
  await evalJs('document.getElementById("scIngest").click(); true'); await waitFor(`GL.views.clinic.ui.horizon > ${h0}`);
  check('ingest next data batch moves the data horizon', (await evalJs('GL.views.clinic.ui.horizon')) > h0);
  await evalJs('document.getElementById("btnWhatIf").click(); true'); await sleep(250);
  await evalJs('document.querySelector(\'[data-preset="combo"]\').click(); true'); await sleep(200);
  await evalJs('document.getElementById("wiRun").click(); true'); await waitFor('!!document.querySelector("#wiOut .stack")', 30000);
  await evalJs('document.getElementById("wiSave").click(); true'); await waitFor('!!document.querySelector("#wiSaved h3")', 5000);
  check('what-if on a digital patient: result, overlay, scenario saved', await evalJs('!!document.querySelector("#ch_anc path.whatif") && !!document.querySelector("#wiSaved h3")'));
  await shot(PUB, '03_clinic_DP03_whatif_guest_1440');
  await evalJs('document.getElementById("wiClose").click(); true');
  check('guest research tabs: no Data tab', await evalJs('(async () => { location.hash = "#/research/model"; await new Promise(r => setTimeout(r, 600)); return !document.querySelector(\'.rnav [data-tab="data"]\'); })()'));
  for (const tab of ['model', 'equations', 'symbols', 'accuracy', 'twins', 'versions', 'interfaces']) {
    await go('#/research/' + tab, '!!document.querySelector(".rmain h1")');
    if (tab === 'twins') await waitFor('!!document.querySelector("#twinBody table")');
    if (tab === 'versions') await waitFor('!!document.querySelector("#verBody table")');
    if (tab === 'model') { await evalJs('document.querySelector(\'[data-node="mpp_n"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200); check('model diagram node shows KaTeX equation', await evalJs('!!document.querySelector("#nodeDetail .katex")'));
      check('research pages: "MarrowTwin (Model2 v1.0)" in the version chip and the navigation', await evalJs('document.querySelector(".rhead .chip-twin").innerText.includes("MarrowTwin (Model2 v1.0)") && document.querySelector(".rnav .foot").innerText.includes("MarrowTwin (Model2 v1.0)")'));
      await evalJs('document.querySelector(\'[data-node="b_mo"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200);
      check('diagram: monocyte node shows both kill terms and θaM', await evalJs(`(() => { const tx = Array.from(document.querySelectorAll('#nodeDetail annotation')).map(a => a.textContent).join(' ');
        return tx.includes(${JSON.stringify(TEX.kc)}) && tx.includes(${JSON.stringify(TEX.am)}) && document.getElementById('nodeDetail').innerText.includes(GL.i18n.t('sym.th_am')); })()`));
      await evalJs('document.querySelector(\'[data-node="hsc"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(150); }
    if (tab === 'equations') {
      check('equations rendered (KaTeX)', await evalJs('document.querySelectorAll(".eq .katex").length >= 10'));
      // every nested model (5 x 3 clone x 3 route = 45): chain equation carries both kill terms, never greyed; MF unchanged
      const nested = await evalJs(`(() => { const c = GL.registry.get('marrowtwin').card, ids = [], KC = ${JSON.stringify(TEX.kcNew)}, AM = ${JSON.stringify(TEX.amNew)};
        ['M0', 'Ma', 'Mb', 'Mab', 'Mabx'].forEach(q => ['', 'C_ddr', 'C_dta'].forEach(cl => ['', '_Zh', '_Zm'].forEach(r => ids.push(q + cl + r))));
        const ok = ids.filter(id => { const e = c.equations(GL.model2.spec(id)).find(x => x.id === 'chain'); return e && e.latex.split(KC).length === 3 && e.latex.includes(AM); });
        const mf = c.equations(GL.model2.spec('MF')).every(e => !e.latex.includes(${JSON.stringify(TEX.kc)}) && !e.latex.includes(${JSON.stringify(TEX.am)}));
        return { n: ids.length, ok: ok.length, mf, shown: Array.from(document.querySelectorAll('.eq annotation')).some(a => a.textContent.includes(KC)) && !!document.querySelector('.eq-new-key') }; })()`);
      check('equations: kill terms κc·θaP·C (T1, T2) and θaM·C (B_Mo) in all 45 nested models, rendered with the legend; MF without them', nested.ok === nested.n && nested.n === 45 && nested.mf && nested.shown, JSON.stringify(nested));
    }
    if (tab === 'symbols') {
      check('symbol table lists all parameters', await evalJs('document.querySelectorAll(".tbl tbody tr").length >= 47'));
      const sy = await evalJs(`(() => { const row = (id) => Array.from(document.querySelectorAll('.tbl tbody tr')).find(r => r.cells[1] && r.cells[1].querySelector('.tiny') && r.cells[1].querySelector('.tiny').innerText.trim() === id);
        const st = GL.loader.index().population.struct, out = {};
        for (const [id, unit] of [['kappa_c', '–'], ['th_am', '1/(d·C)']]) { const r = row(id); out[id] = !!r && r.cells[1].innerText.includes(GL.i18n.t('sym.' + id)) && r.cells[2].innerText.trim() === unit &&
          r.cells[4].innerText.includes(GL.i18n.t('kind.estimated')) && r.cells[5].innerText.includes(GL.i18n.t('sym.centre') + ' ' + GL.util.fmt(st[id], 3)) && !!r.querySelector('.katex'); }
        return out; })()`);
      check('symbols: kappa_c and th_aM with meaning, unit, centre, estimated', sy.kappa_c && sy.th_am, JSON.stringify(sy));
    }
    if (tab === 'twins') check('twins vs truth: 16 digital patients', await evalJs('document.querySelectorAll("#twinBody tbody tr").length === 16'));
    if (tab === 'versions') check('guest versions: digital patients only', await evalJs('Array.from(document.querySelectorAll("#verP option")).every(o => /^DP/.test(o.value))'));
    await shot(PUB, `04_research_${tab}_guest_1440`, true);
  }
  await go('#/research/data', '!!document.querySelector(".rmain h1")');
  check('guest: Data tab not reachable', await evalJs('!location.hash.includes("/data")'));
  await rawKeys('guest en');
  check('browser engine is the Model2 port (41 parameters incl. kappa_c, th_am; population carries both)', await evalJs('GL.model2.PAR.length === 41 && GL.model2.PAR.indexOf("kappa_c") === 16 && GL.model2.PAR.indexOf("th_am") === 17 && "kappa_c" in GL.loader.index().population.struct && "th_am" in GL.loader.index().population.struct'));
  for (const lang of ['ja', 'de', 'fr', 'es', 'zh-CN', 'zh-TW']) {
    await evalJs(`GL.i18n.set('${lang}'); true`);
    await go('#/clinic/DP02', clinicReady); await sleep(200); await evalJs('GL.app.render(); true'); await waitFor(clinicReady); await noScroll('guest ' + lang); await rawKeys('clinic ' + lang);
    if (['ja', 'zh-TW'].includes(lang)) { for (const tab of ['model', 'equations', 'symbols', 'accuracy']) { await go('#/research/' + tab, '!!document.querySelector(".rmain h1")'); await noOldName(`research ${tab} ${lang}`); } await go('#/clinic/DP02', clinicReady); }
    if (['ja', 'zh-CN'].includes(lang)) await shot(PUB, `05_clinic_DP02_${lang}_1440`);
  }
  await evalJs(`GL.i18n.set('en'); GL.app.theme('dark'); GL.app.render(); true`); await waitFor(clinicReady); await noScroll('guest dark'); await shot(PUB, '06_clinic_DP02_dark_1440');
  await evalJs(`GL.app.theme('light'); true`);
  const stor = await evalJs('(async () => ({ ls: Object.keys(localStorage), idb: indexedDB.databases ? (await indexedDB.databases()).map(d => d.name) : [] }))()');
  check('browser storage: localStorage glm2.* only, IndexedDB greenlight_model2', stor.ls.length > 0 && stor.ls.every(k => k.startsWith('glm2.')) && stor.idb.includes('greenlight_model2') && !stor.idb.includes('greenlight_model1'), JSON.stringify({ ls: stor.ls.length, idb: stor.idb }));
  await logout();

  // 3. clinician A
  check('clinician A signs in', await login(U.clinA));
  await go('#/clinic', clinicReady);
  const ra = (await realIdsInRoster()).sort();
  check('clinician A: exactly the 2 assigned real patients + 16 digital', JSON.stringify(ra) === JSON.stringify([pa1, pa2].sort()) && await evalJs('document.querySelectorAll(".ritem").length === 18'));
  await go('#/clinic/' + pb, clinicReady);
  check('clinician A: B\'s patient is not reachable', await evalJs(`!location.hash.endsWith(${JSON.stringify('/' + pb)}) && !new RegExp('\\\\b' + ${JSON.stringify(pb)} + '\\\\b').test(document.body.innerText)`));
  await go('#/clinic/' + pa1, clinicReady);
  await noScroll('clinician real 1440x900');
  check('watermark with user name and access banner over the patient view', await evalJs(`!!document.getElementById('wm') && document.getElementById('wm').dataset.text.startsWith(${JSON.stringify(U.clinA.u)}) && !!document.querySelector('.authbar')`));
  await shot(PRIV, '10_clinic_real_clinicianA_1440');
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.25, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`); await sleep(300);
  await evalJs('document.querySelector(\'[data-mode="once"]\').click(); document.getElementById("scPlay").click(); true'); await sleep(900); await evalJs('GL.ui.scrubber.stopAll(); true');
  await evalJs('document.getElementById("btnWhatIf").click(); true'); await sleep(250);
  await evalJs('document.querySelector(\'[data-preset="gcsf_reactive"]\').click(); true'); await sleep(200);
  await evalJs('document.getElementById("wiRun").click(); true'); await waitFor('!!document.querySelector("#wiOut .stack")', 30000);
  await evalJs('document.getElementById("wiSave").click(); true'); await waitFor('!!document.querySelector("#wiSaved h3")', 5000);
  await shot(PRIV, '11_clinic_real_whatif_1440');
  const store = await evalJs(`(async () => { const d = await GL.store.db(); const sc = await new Promise(r => { const o = []; const c = d.transaction('scenarios').objectStore('scenarios').openCursor(); c.onsuccess = () => { const x = c.result; if (x) { o.push(x.value.patient); x.continue(); } else r(o); }; });
    return { mem: GL.session.scenarios.map(s => s.patient), idb: sc, ls: Object.keys(localStorage).filter(k => /\\bP\\d{2}\\b/.test(k) || /\\bP\\d{2}\\b/.test(localStorage.getItem(k) || '')) }; })()`);
  check('real what-if scenario only in memory; no real patient in IndexedDB or localStorage', store.mem.includes(pa1) && !store.idb.some(p => !/^DP/.test(p)) && store.ls.length === 0);
  const audit = await evalJs(`(async () => (await GL.audit.all()).filter(e => e.subject === ${JSON.stringify(pa1)}).map(e => e.action))()`);
  check('audit log records open, scrub, play and what-if', ['open', 'scrub', 'play', 'what-if'].every(a => audit.includes(a)), audit.join(','));
  const tabs = await evalJs('(async () => { location.hash = "#/research/model"; await new Promise(r => setTimeout(r, 600)); return Array.from(document.querySelectorAll(".rnav [data-tab]")).map(a => a.dataset.tab); })()');
  check('clinician research tabs: Model, Equations, Symbols, Interfaces only', JSON.stringify(tabs) === JSON.stringify(['model', 'equations', 'symbols', 'interfaces']), tabs.join(','));
  await viewport(1920, 1080); await go('#/clinic/' + pa1, clinicReady); await noScroll('clinician real 1920x1080'); await shot(PRIV, '12_clinic_real_1920'); await viewport(1440, 900);
  // idle lock
  await go('#/clinic/' + pa1, clinicReady);
  await evalJs('GL.session.setIdle(1200); true'); await waitFor('!!document.getElementById("ulPw")', 8000);
  check('idle lock: lock screen, decrypted data cleared, no watermark', await evalJs(`Object.keys(GL.session.real).length === 0 && GL.session.scenarios.length === 0 && !GL.views.clinic.ui.p && !new RegExp('\\\\b' + ${JSON.stringify(pa1)} + '\\\\b').test(document.body.innerText) && !document.getElementById('wm')`));
  await evalJs(`(() => { document.getElementById('ulPw').value = ${JSON.stringify(U.clinA.pw)}; document.getElementById('lockForm').requestSubmit(); return true; })()`);
  await waitFor(`!GL.session.locked && Object.keys(GL.session.real).length === 2`, 20000);
  await evalJs('GL.session.setIdle(15 * 60 * 1000); true');
  check('unlock with the password restores the data', await evalJs(`!GL.session.locked && !!GL.session.real[${JSON.stringify(pa1)}]`));
  await logout();
  check('logout clears decrypted data from memory and screen', await evalJs(`Object.keys(GL.session.real).length === 0 && Object.keys(GL.session.research).length === 0 && !GL.session.cohort && !GL.session.user && !GL.views.clinic.ui.p && !new RegExp('\\\\b' + ${JSON.stringify(pa1)} + '\\\\b').test(document.body.innerText)`));

  // 4. clinician B
  check('clinician B signs in', await login(U.clinB));
  await go('#/clinic', clinicReady);
  check('clinician B sees only B\'s patient', JSON.stringify(await realIdsInRoster()) === JSON.stringify([pb]));
  await logout();

  // 5. researcher
  check('researcher signs in', await login(U.res));
  await go('#/clinic', clinicReady);
  check('researcher: no real patient in the clinician portal', (await realIdsInRoster()).length === 0);
  await go('#/research/data', '!!document.getElementById("researchTable")');
  const rd = await evalJs(`(() => { const rows = Array.from(document.querySelectorAll('#researchTable tbody tr')); const txt = document.body.innerText;
    return { n: rows.length, ids: rows.every(r => /^R-[0-9A-F]{6}$/.test(r.cells[0].innerText.trim())), ages: rows.every(r => /^\\d{2}-\\d{2}$|^90\\+$/.test(r.cells[2].innerText.trim())),
      clinIds: /\\bP\\d{2}\\b/.test(txt), rec: !!document.querySelector('.plan-row, #btnWhatIf'), small: txt.includes('<5') }; })()`);
  check('research Data tab: 35 pseudonyms (no clinical id), age bands, small cells "<5", no recommendations', rd.n === clinIds.length && rd.ids && rd.ages && !rd.clinIds && !rd.rec && rd.small, JSON.stringify({ n: rd.n }));
  check('research DOM and payloads: no bridging regimen, exact age or plans', await evalJs(`(() => { const t = document.body.innerText, j = JSON.stringify(GL.session.research) + JSON.stringify(GL.session.cohort);
    return ${JSON.stringify(regimens)}.every(r => !t.includes(r) && !j.includes(r)) && !j.includes('"age"') && !/"plans"|"top3"|"drv"/.test(j); })()`));
  await shot(PRIV, '20_research_data_researcher_1440', true);
  const psd = await evalJs('Object.keys(GL.session.research).sort()[0]');
  await evalJs(`GL.views.research.ui.vp = ${JSON.stringify(psd)}; true`);
  await go('#/research/versions', '!!document.querySelector("#verBody table")');
  check('research Versions of a pseudonym: no top strategy, no what-if', await evalJs(`!document.body.innerText.includes(GL.i18n.t('rv.top')) && !document.getElementById('btnWhatIf') && document.getElementById('verP').value === ${JSON.stringify(psd)}`));
  await shot(PRIV, '21_research_versions_researcher_1440', true);
  await go('#/clinic/' + psd, clinicReady);
  check('researcher cannot open a pseudonym in the clinician portal', await evalJs(`!location.hash.endsWith(${JSON.stringify('/' + psd)})`));
  await logout();

  // 6. data steward
  check('data steward signs in', await login(U.stew));
  await waitFor('!!document.getElementById("auditExport")');
  check('steward: admin view only, user list, audit chain intact', await evalJs(`location.hash === '#/admin' && ${JSON.stringify(created)}.every(u => document.body.innerText.includes(u)) && +document.getElementById('auditN').textContent > 0 && document.getElementById('auditChain').classList.contains('ok-badge') && Object.keys(GL.session.real).length === 0`));
  await evalJs('document.getElementById("auditExport").click(); true'); await waitFor('!!GL.views.admin.lastExport', 8000);
  check('audit log export (hash chain verified)', await evalJs('GL.views.admin.lastExport.ok && GL.views.admin.lastExport.n > 0'));
  await shot(PRIV, '30_admin_steward_1440');
  await go('#/clinic', null); await sleep(500);
  check('steward cannot open a portal', await evalJs('location.hash === "#/admin"'));
  await logout();

  check('console errors / warnings: 0', errors.length === 0, errors.length ? errors.slice(0, 5).join(' | ') : undefined);
  ws.close();
})().catch(e => { check('walkthrough ran', false, String(e && e.message || e)); }).finally(async () => {
  if (proc) { proc.kill(); await sleep(800); }
  // the throw-away Chrome profile holds this run's audit log (study ids + actions): remove it
  try { fs.rmSync(profile, { recursive: true, force: true }); } catch (e) { console.log('note: Chrome profile not removed: ' + profile); }
  for (const u of created) { try { cli(['revoke', u]); } catch (e) { results.push({ ok: false }); console.log('FAIL cleanup of a temporary account'); } }
  const left = A.loadUsers().users.filter(x => x.u.endsWith('_' + tag));
  check('temporary accounts removed', left.length === 0);
  check('existing accounts unchanged (role, salt, encrypted key, public key)', creds() === credsBefore);
  const nf = results.filter(r => !r.ok).length;
  console.log(nf ? `FAIL: ${nf} of ${results.length} checks failed` : `PASS: ${results.length} checks, 0 console errors`);
  process.exit(nf ? 1 : 0);
});
