/* GreenLight · tools/browser_test.js
 * Headless-Chrome walkthrough over the DevTools protocol, with access control (Model1 plan 13.6 item 11, carried over by Model2 plan 6.5 and Model3).
 * Temporary accounts (random names, random passwords kept in memory and given to tools/users.js on stdin)
 * are created at the start and always removed at the end. Flow: login screen -> wrong password -> guest
 * (both portals, digital patients only: filters incl. real, roster scroll kept + green selection, prediction wording,
 * x-aligned time axes, scrub steps (no playback), ingest, supportive-care panel over the right column only
 * that follows the As-of day and closes on an outside click, every allowed research tab,
 * languages, dark theme) -> clinician A (assigned real patients only, one screen at 1440x900 and 1920x1080,
 * watermark, real what-if kept in memory, audit entries, idle lock + unlock, logout clears data) ->
 * clinician B (only B's patient) -> researcher (research view: pseudonyms, bands, no recommendations or
 * what-if) -> data steward (user list, audit chain, export). Checks 0 console errors / warnings.
 * MarrowTwin (Model3) checks: public name and version label "MarrowTwin (Model3 v1.0)", no "Model1" / "Model2" in any
 * view; the Model3 terms in the equations of all 45 nested models (teal neutrophil storage pool N, teal platelet e_Z,Pl
 * P2 where the model has the MPP route of Z, grey P1 platelet consumption; the chemotherapy kill terms of Model2 kept)
 * and in the diagram (7 chain cells, 3 storage); symbol table rows th_cP / n_N / k_N / d_ZPl and the stability limits
 * gamma <= 1.0, MTT_N <= 15; the 45-parameter / 32-state browser engine; glm3.* / greenlight_model3 browser storage.
 * Temporary accounts left by an aborted earlier run (tb_<role>_<6 hex>) are removed first;
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
const port = 9333, profile = '/tmp/glm3-chrome-cdp-' + Date.now();
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
// terms as the model card writes them (js/modules/model3.js): Model2 kill terms (plain), Model3 N and P2 (teal = NEW), P1 (grey = OFF)
const TEX = { kc: String.raw`\kappa_c\,\theta_{aP}\,C`, am: String.raw`\theta_{aM}\,C\,B_{\mathrm{Mo}}` };
const NEWC = String.raw`\textcolor{#0e9fb5}{`, OFFC = String.raw`\textcolor{#a3aab0}{`;
TEX.store = NEWC + String.raw`O_j=r_N(t)\,T_{Nj}\ (j>k_N)` + '}'; TEX.div = NEWC + String.raw`\mathbf{1}[j\le k_N]` + '}';
TEX.nk = NEWC + String.raw`n_N=7,\ k_N=4` + '}';
TEX.ezpl = String.raw`e_{Z,Pl}=\mathrm{logit}^{-1}\big(\mathrm{logit}\,e_Z+d_{Z,Pl}\big)`;
TEX.p1 = OFFC + String.raw`-\mathbf{1}[\ell=\mathrm{Pl}]\,\theta_{cP}\,D(t)\,B_{\mathrm{Pl}}` + '}';
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
  const atDay = (sel, day, type) => `(() => { const p = document.querySelector(${JSON.stringify(sel)}), r = p.getBoundingClientRect(), m = GL.charts.xMargin, x = r.left + GL.charts.timeScale(m.l, r.width - m.r)(${day}); p.dispatchEvent(new PointerEvent(${JSON.stringify(type || 'pointerdown')}, { bubbles: true, clientX: x, clientY: r.top + 20, buttons: 1, pointerId: 1 })); return true; })()`;
  const go = async (hash, ready) => { await evalJs(`location.hash = ${JSON.stringify(hash)}; true`); await sleep(300); if (ready) await waitFor(ready); await sleep(250); };
  const noScroll = async (label) => { const r = await evalJs('[document.documentElement.scrollHeight, window.innerHeight, document.body.scrollHeight]'); check(`no vertical scroll ${label} (${r[0]} <= ${r[1]})`, r[0] <= r[1] && r[2] <= r[1]); };
  const noOldName = async (label) => { const r = await evalJs('(() => { const h = document.title + document.body.innerHTML; const m = h.match(/Model ?[12]\\b|model[12]/i); return m ? m[0] : ""; })()'); check(`no "Model1" / "Model2" in the page ${label}`, !r, r || undefined); };
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
  check('MarrowTwin: landing lead, version chip "MarrowTwin (Model3 v1.0)", page title', await evalJs('document.body.innerText.includes(GL.i18n.t("land.lead")) && GL.i18n.t("land.lead").includes("MarrowTwin") && Array.from(document.querySelectorAll(".land-foot .chip")).some(c => c.innerText.includes("MarrowTwin (Model3 v1.0)")) && document.title.includes("MarrowTwin")'));
  await noOldName('landing');
  await shot(PUB, '01_landing_guest_1440');
  await go('#/clinic', clinicReady);
  check('guest: no real patient in the clinician roster (16 digital)', (await realIdsInRoster()).length === 0 && await evalJs('document.querySelectorAll(".ritem").length === 16'));
  const counts = {};
  for (const f of ['digital', 'real', 'active', 'completed', 'all']) { await evalJs(`document.querySelector('.ftabs [data-f="all"]').click(); document.querySelector('.ftabs [data-f="${f}"]').click(); true`); await sleep(120); counts[f] = await evalJs('document.querySelectorAll(".ritem").length'); }
  check('guest roster filters (digital 16 / real 0 / 8 / 8 / 16)', counts.digital === 16 && counts.real === 0 && counts.active === 8 && counts.completed === 8 && counts.all === 16, JSON.stringify(counts));
  await evalJs('document.querySelector(\'.ftabs [data-f="gcsf"]\').click(); true'); await sleep(120);
  const gl = await evalJs('({ ids: Array.from(document.querySelectorAll(".ritem")).map(r => r.dataset.pid), tip: document.querySelector(\'.ftabs [data-f="gcsf"]\').title, n: document.querySelector(\'.ftabs [data-f="gcsf"] .muted\').textContent })');
  check('guest "Received G-CSF" filter: the completed digital patients with G-CSF (DP10, DP11, DP12, DP13, DP15); none in treatment before their data horizon', JSON.stringify(gl.ids) === JSON.stringify(['DP10', 'DP11', 'DP12', 'DP13', 'DP15']) && gl.n === '5' && gl.tip.length > 20, JSON.stringify(gl));
  await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); true'); await sleep(120);
  const fc = async (keys) => { await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); true'); for (const k of keys) await evalJs(`document.querySelector('.ftabs [data-f="${k}"]').click(); true`); await sleep(120);
    return evalJs('({ ids: Array.from(document.querySelectorAll(".ritem")).map(r => r.dataset.pid), pressed: Array.from(document.querySelectorAll(".ftabs [aria-pressed=\\"true\\"]")).map(b => b.dataset.f) })'); };
  const c1 = await fc(['completed', 'gcsf']), c2 = await fc(['active', 'completed']), c3 = await fc(['digital', 'active']), c4 = await fc(['active', 'active']);
  check('filter buttons combine: completed + G-CSF = both (5); in treatment + completed = either (16); digital + in treatment (8); a second click turns one off', JSON.stringify(c1.ids) === '["DP10","DP11","DP12","DP13","DP15"]' && JSON.stringify(c1.pressed) === '["completed","gcsf"]' && c2.ids.length === 16 && c3.ids.length === 8 && c4.ids.length === 16 && JSON.stringify(c4.pressed) === '["all"]', JSON.stringify([c1, c2.ids.length, c3.ids.length, c4.pressed]));
  await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); true'); await sleep(120);
  // roster: opening a patient low in the list keeps the scroll position; the open patient is marked green
  await evalJs('(() => { const l = document.getElementById("rosterList"); l.scrollTop = 260; l.dispatchEvent(new Event("scroll")); return true; })()'); await sleep(100);
  const top0 = await evalJs('document.getElementById("rosterList").scrollTop');
  await evalJs('(() => { const l = document.getElementById("rosterList"), b = l.getBoundingClientRect(); const r = Array.from(l.querySelectorAll(".ritem")).find(x => { const q = x.getBoundingClientRect(); return q.top > b.top + 200 && q.bottom < b.bottom; }); window.__pick = r.dataset.pid; r.click(); return true; })()');
  await waitFor(`location.hash.endsWith("/" + window.__pick) && ${clinicReady}`); await sleep(200);
  const rs = await evalJs('(() => { const l = document.getElementById("rosterList"), a = l.querySelector(".ritem.active"), b = l.getBoundingClientRect(), q = a.getBoundingClientRect(); return { top: l.scrollTop, pid: a.dataset.pid, vis: q.top >= b.top && q.bottom <= b.bottom, bg: getComputedStyle(a).backgroundColor, bc: getComputedStyle(a).borderTopColor, cur: a.getAttribute("aria-current") }; })()');
  check('roster keeps its scroll position on a patient change; the open patient is green (brand-100 / brand-600) and in sight', top0 > 0 && Math.abs(rs.top - top0) <= 2 && rs.vis && rs.cur === 'true' && rs.bg === 'rgb(220, 252, 231)' && rs.bc === 'rgb(22, 163, 74)', JSON.stringify(Object.assign({ top0 }, rs)));
  for (const pid of ['DP02', 'DP10']) {
    await go('#/clinic/' + pid, clinicReady); await noScroll('guest ' + pid + ' 1440x900'); await noOldName('clinic ' + pid);
    check(`${pid}: 4 charts, ferritin, evolution, scrubber in the left column, top-3`, await evalJs('document.querySelectorAll(".charts4 svg").length === 4 && !!document.querySelector("#ch_evo svg") && !!document.querySelector(".cleft #scrub #scSvg circle.tl-handle") && document.querySelectorAll(".plan-row").length === 3'));
    check(`${pid}: "prediction", never "forecast"; predictions name their target`, await evalJs('!/forecast/i.test(document.body.innerText) && document.body.innerText.includes(GL.i18n.t("est.ph")) && document.body.innerText.includes(GL.i18n.t("est.nadir")) && document.body.innerText.includes(GL.i18n.t("evo.sub"))'));
    // inflammation driver, how the twin learned and the scrubber share one time axis: the as-of day sits at the same x
    const ax = await evalJs('(() => { const x = (el) => { const r = el.getBoundingClientRect(); return r.left + r.width / 2; }; return [x(document.querySelector("#ch_fer line.asof-line")), x(document.querySelector("#ch_evo line.asof-line")), x(document.querySelector("#scSvg circle.tl-handle")), x(document.querySelector("#ch_fer .plot")) - x(document.querySelector("#scSvg"))]; })()');
    check(`${pid}: ferritin, evolution and scrubber x-aligned (as-of x within 1 px)`, Math.max(...ax.slice(0, 3)) - Math.min(...ax.slice(0, 3)) <= 1 && Math.abs(ax[3]) <= 1, ax.map(v => Math.round(v * 10) / 10).join(', '));
    check(`${pid}: right column reaches the bottom of the page area`, await evalJs('(() => { const r = document.getElementById("cright").getBoundingClientRect(), l = document.querySelector(".cleft").getBoundingClientRect(); return Math.abs(r.bottom - l.bottom) <= 1 && Math.abs(r.top - l.top) <= 1; })()'));
    await shot(PUB, `02_clinic_${pid}_guest_1440`);
  }
  await viewport(1920, 1080); await go('#/clinic/DP02', clinicReady); await noScroll('guest DP02 1920x1080'); await viewport(1440, 900); await go('#/clinic/DP10', clinicReady);
  const day0 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.30, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`); await sleep(300);
  check('scrubber click changes the as-of day', (await evalJs('GL.views.clinic.ui.asOf')) !== day0);
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.62, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`);
  await sleep(200);
  const st0 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs('document.getElementById("scNext").click(); true'); await sleep(150); const st1 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs('document.getElementById("scPrev").click(); true'); await sleep(150); const st2 = await evalJs('GL.views.clinic.ui.asOf');
  check('scrubber steps forward and back; no playback controls', st1 > st0 && st2 === st0 && await evalJs('!document.getElementById("scPlay") && !document.querySelector("#scrub [data-mode]") && !document.getElementById("scSpeed")'), [st0, st1, st2].join(' → '));
  // G-CSF after the selected day is grey (scrubber and ANC chart); given by then it is blue
  const gq = '({ blue: document.querySelectorAll("#scSvg rect.tl-gcsf:not(.tl-gcsf-future)").length, grey: document.querySelectorAll("#scSvg rect.tl-gcsf-future").length, ancBlue: document.querySelectorAll("#ch_anc rect.gcsf-bar:not(.gcsf-future)").length, ancGrey: document.querySelectorAll("#ch_anc rect.gcsf-bar.gcsf-future").length, w: GL.views.clinic.ui.p.design.gcsf.length })';
  await evalJs('GL.views.clinic.ui.scrub.set(-5); true'); await sleep(200); const g0 = await evalJs(gq);
  await evalJs('GL.views.clinic.ui.scrub.set(14); true'); await sleep(250); await shot(PUB, '07_clinic_DP10_D14_gcsf_grey_guest_1440');
  await evalJs('GL.views.clinic.ui.scrub.set(86); true'); await sleep(200); const g1 = await evalJs(gq);
  check('G-CSF windows grey after the selected day, blue up to it (scrubber + ANC chart)', g0.w > 0 && g0.blue === 0 && g0.grey > 0 && g0.ancBlue === 0 && g0.ancGrey > 0 && g1.blue > 0 && g1.ancBlue > 0, JSON.stringify([g0, g1]));
  // a click or drag on any time chart sets the As-of day (like the scrubber); hover guide shows the day
  await evalJs(atDay('#ch_anc .plot', 14)); await evalJs(atDay('#ch_anc .plot', 14, 'pointerup')); await sleep(150);
  const k1 = await evalJs('({ d: GL.views.clinic.ui.asOf, s: document.getElementById("scrubDay").textContent, g: (document.querySelector("#ch_anc .tguide span") || {}).textContent })');
  await evalJs(atDay('#ch_evo .plot', 7)); await evalJs(atDay('#ch_evo .plot', 7, 'pointerup')); await sleep(150); const k2 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs(atDay('#ch_fer .plot', 21)); await evalJs(atDay('#ch_fer .plot', 21, 'pointerup')); await sleep(150); const k3 = await evalJs('GL.views.clinic.ui.asOf');
  await evalJs(atDay('#ch_plt .plot', 3)); await evalJs(atDay('#ch_plt .plot', 10, 'pointermove')); await evalJs(atDay('#ch_plt .plot', 10, 'pointerup')); await sleep(150); const k4 = await evalJs('GL.views.clinic.ui.asOf');
  check('single patient: a click on a chart sets the As-of day (ANC, learning, ferritin), a drag moves it, the guide names the day', k1.d === 14 && k1.s === 'D+14' && k1.g === 'D+14' && k2 === 7 && k3 === 21 && k4 === 10, JSON.stringify([k1, k2, k3, k4]));
  // data after the selected day are hollow in the four count charts, also for a patient in treatment
  await go('#/clinic/DP06', clinicReady);
  await evalJs('GL.views.clinic.ui.scrub.set(5); true'); await sleep(250);
  const hol = await evalJs('(() => { const p = GL.views.clinic.ui.p, L = GL.views.clinic.ui.asOf, hz = GL.views.clinic.ui.horizon; const n = (v) => p.obs[v].filter(o => o[0] > L && o[0] <= hz && o[0] <= 90).length; const c = ["anc", "plt", "amc", "alc"].map(v => document.querySelectorAll("#ch_" + v + " .obs-future").length); const bg = getComputedStyle(document.querySelector(".card")).backgroundColor, cs = Array.from(document.querySelectorAll(".charts4 .obs-future")).map(e => getComputedStyle(e)); return { L, hz, want: ["anc", "plt", "amc", "alc"].map(n), got: c, allHollow: cs.every(x => x.fill === "none"), outlined: cs.every(x => x.stroke !== bg && x.stroke !== "none" && +x.opacity >= 0.6) }; })()');
  await shot(PUB, '08_clinic_DP06_D5_hollow_guest_1440');
  check('count charts: data after the selected day hollow and visible (outlined in the series colour; patient in treatment, up to the data horizon)', hol.hz > hol.L && hol.got.every((g, i) => g === hol.want[i]) && hol.want[0] > 0 && hol.allHollow && hol.outlined, JSON.stringify(hol));
  await go('#/clinic/DP03', clinicReady);
  const h0 = await evalJs('GL.views.clinic.ui.horizon');
  await evalJs('document.getElementById("scIngest").click(); true'); await waitFor(`GL.views.clinic.ui.horizon > ${h0}`);
  check('ingest next data batch moves the data horizon', (await evalJs('GL.views.clinic.ui.horizon')) > h0);
  // learning chart with the adjusted plan: a plan without effect on ANC (TPO-RA) must leave it exactly unchanged
  await evalJs('document.getElementById("btnWhatIf").click(); true'); await sleep(250);
  await evalJs('document.querySelector(\'[data-preset="tpo_ra"]\').click(); true'); await sleep(150);
  await evalJs('document.getElementById("wiRun").click(); true'); await waitFor('GL.views.clinic.ui.evo && GL.views.clinic.ui.evo.d.size === GL.views.clinic.ui.p.asof.filter(x => x.L <= GL.views.clinic.ui.a.L && x.L <= GL.views.clinic.ui.horizon).length', 60000);
  const tpo = await evalJs('(() => { const d = Array.from(GL.views.clinic.ui.evo.d.values()); return { n: d.length, maxAbs: Math.max(...d.map(x => Math.max(Math.abs(x[0]), Math.abs(x[1])))), paths: document.querySelectorAll("#ch_evo path.evo-plan").length }; })()');
  check('learning chart: TPO-RA alone (no ANC effect) adds 0 (< 1e-12) on every as-of day (same draws with and without the plan)', tpo.n > 1 && tpo.maxAbs < 1e-12 && tpo.paths === 2, JSON.stringify(tpo));
  await evalJs('document.getElementById("wiClose").click(); true'); await sleep(100);
  await evalJs('document.getElementById("btnWhatIf").click(); true'); await sleep(250);
  await evalJs('document.querySelector(\'[data-preset="combo"]\').click(); true'); await sleep(200);
  await evalJs('document.getElementById("wiRun").click(); true'); await waitFor('!!document.querySelector("#wiOut .stack")', 30000);
  await evalJs('document.getElementById("wiSave").click(); true'); await waitFor('!!document.querySelector("#wiSaved h3")', 5000);
  await waitFor('GL.views.clinic.ui.evo && GL.views.clinic.ui.evo.d.size === GL.views.clinic.ui.p.asof.filter(x => x.L <= GL.views.clinic.ui.a.L && x.L <= GL.views.clinic.ui.horizon).length', 60000);
  const evo = await evalJs('(() => { const ui = GL.views.clinic.ui, d = Array.from(ui.evo.d.values()); return { n: d.length, days: ui.p.asof.filter(x => x.L <= ui.a.L && x.L <= ui.horizon).length, moved: d.some(x => Math.abs(x[0]) > 0 || Math.abs(x[1]) > 0), paths: document.querySelectorAll("#ch_evo path.evo-plan").length, lbl: document.querySelector("#ch_evo .card-head .sub").textContent }; })()');
  check('learning chart changes with the adjusted plan: dashed P boundaries on every as-of day up to the selected one, labelled', evo.n === evo.days && evo.days > 1 && evo.moved && evo.paths === 2 && evo.lbl.endsWith(await evalJs('GL.i18n.t("evo.plan")')), JSON.stringify(evo));
  check('supportive-care panel on a digital patient: result, overlay, scenario saved', await evalJs('!!document.querySelector("#ch_anc path.whatif") && !!document.querySelector("#wiSaved h3")'));
  check('supportive-care panel lies over the right column only (charts uncovered)', await evalJs('(() => { const d = document.getElementById("drawer").getBoundingClientRect(), r = document.getElementById("cright").getBoundingClientRect(), l = document.querySelector(".cleft").getBoundingClientRect(); return d.width > 0 && Math.abs(d.left - r.left) <= 1 && Math.abs(d.right - r.right) <= 1 && d.left >= l.right && d.bottom <= r.bottom + 1; })()'));
  const anim = await evalJs('(() => { const d = document.getElementById("drawer"), b = getComputedStyle(d, "::before"), a = getComputedStyle(d, "::after"), cards = Array.from(d.querySelectorAll(".card")).map(c => getComputedStyle(c).backgroundColor), surf = getComputedStyle(document.querySelector(".cright .card") || document.body).backgroundColor; return { before: b.animationName, after: a.animationName, pe: b.pointerEvents + "/" + a.pointerEvents, cardsOpaque: cards.length > 0 && cards.every(c => c === surf) }; })()');
  check('open panel is highlighted by a subtle animation (green flow behind, light round the border); cards stay opaque, no pointer capture', anim.before === 'glFlow' && anim.after === 'glSpin' && anim.pe === 'none/none' && anim.cardsOpaque, JSON.stringify(anim));
  check('button and panel title name the action (no "What-if")', await evalJs('!/what-if/i.test(document.getElementById("btnWhatIf").innerText + document.querySelector("#drawer .drawer-head").innerText) && document.querySelector("#drawer .drawer-head").innerText.includes(GL.i18n.t("plans.whatif"))'));
  await shot(PUB, '03_clinic_DP03_whatif_guest_1440');
  // the panel follows the As-of day: title moves to the new decision day, the plan is re-simulated there
  const wl0 = await evalJs('GL.views.clinic.ui.whatif.L');
  await evalJs('document.getElementById("scPrev").click(); true');
  await waitFor(`GL.views.clinic.ui.whatif && GL.views.clinic.ui.whatif.L !== ${wl0} && GL.views.clinic.ui.whatif.L === GL.views.clinic.ui.a.L && !!document.querySelector("#wiOut .stack")`, 30000);
  const fol = await evalJs('({ open: document.getElementById("drawer").classList.contains("open"), title: document.querySelector("#drawer .drawer-head").innerText, L: GL.views.clinic.ui.a.L, line: !!document.querySelector("#ch_anc path.whatif"), evo: document.querySelectorAll("#ch_evo path.evo-plan").length })');
  check('panel stays open and follows the As-of day (title, re-simulation, overlay)', fol.open && fol.line && fol.evo === 2 && fol.L < wl0 && fol.title.includes('D+' + (Math.floor(fol.L) + 1)), JSON.stringify(fol));
  await evalJs(atDay('#ch_anc .plot', 3)); await evalJs(atDay('#ch_anc .plot', 3, 'pointerup')); await sleep(200);
  check('panel open: a click on a chart sets the day and the panel stays open', await evalJs('GL.views.clinic.ui.asOf === 3 && document.getElementById("drawer").classList.contains("open")'));
  await evalJs('(() => { const p = document.querySelector(".psum .facts"), r = p.getBoundingClientRect(); p.dispatchEvent(new PointerEvent("pointerdown", { bubbles: true, clientX: r.left + 4, clientY: r.top + 4 })); return true; })()'); await sleep(120);
  check('a click on empty space outside the panel closes it', await evalJs('!document.getElementById("drawer").classList.contains("open")'));
  // multi-patient view: tick boxes in the same list, rows in one column under count tabs, shared day
  await go('#/clinic/DP10', clinicReady);
  await evalJs('document.querySelector(\'.modesw [data-mode="multi"]\').click(); true');
  await waitFor('location.hash === "#/clinic/multi" && !!document.querySelector(".mempty, .mrow")');
  const m0 = await evalJs('({ checked: Array.from(document.querySelectorAll(".ritem[aria-checked=\\"true\\"]")).map(r => r.dataset.pid), boxes: document.querySelectorAll(".ritem .rsel").length, n: document.querySelectorAll(".ritem").length, empty: !!document.querySelector(".mempty") })');
  check('multi-patient mode: same list with a tick box per patient, the open patient pre-ticked, empty drawing area until "Show selected"', m0.boxes === m0.n && m0.n === 16 && JSON.stringify(m0.checked) === '["DP10"]' && m0.empty, JSON.stringify(m0));
  await evalJs('["DP02", "DP11", "DP06"].forEach(id => document.querySelector(`.ritem[data-pid="${id}"]`).click()); true'); await sleep(100);
  check('multi-patient: the selection buttons fit the list width; small icons keep their size', await evalJs('(() => { const r = document.getElementById("roster").getBoundingClientRect(); return Array.from(document.querySelectorAll(".mselrow .btn")).every(b => b.getBoundingClientRect().right <= r.right + 0.5) && Array.from(document.querySelectorAll("#view svg.pmark")).every(m => m.getBoundingClientRect().width <= 14); })()'));
  check('ticking updates the "Show selected" count', await evalJs('document.getElementById("mShow").textContent === GL.i18n.t("multi.show", { n: 4 }) && !document.getElementById("mShow").disabled'));
  await evalJs('document.getElementById("mShow").click(); true');
  await waitFor('document.querySelectorAll(".mrow svg path.med").length === 4'); await sleep(300);
  await noScroll('guest multi-patient 1440x900');
  const mv = await evalJs('(() => { const rows = Array.from(document.querySelectorAll(".mrow")), rr = rows.map(r => r.getBoundingClientRect()); const x = (el) => { const b = el.getBoundingClientRect(); return b.left + b.width / 2; }; return { ids: rows.map(r => r.dataset.pid), oneCol: rr.every((b, i) => Math.abs(b.left - rr[0].left) < 1 && Math.abs(b.width - rr[0].width) < 1 && (!i || b.top >= rr[i - 1].bottom - 1)), tabs: document.querySelectorAll(".mtabs [data-tab]").length, table: document.querySelectorAll("#cright .mtbl tbody tr").length, axis: Math.abs(x(rows[0].querySelector("line.asof-line")) - x(document.querySelector("#scSvg circle.tl-handle"))), yAxes: new Set(rows.map(r => Array.from(r.querySelectorAll(".ax text")).filter(t => !/D|\\+|-/.test(t.textContent)).map(t => t.textContent).join(","))).size }; })()');
  check('multi-patient: 4 rows in one column (ticked order), 4 count tabs, a table line per patient, the same y axis in every row, as-of x aligned with the scrubber', JSON.stringify(mv.ids) === '["DP10","DP02","DP11","DP06"]' && mv.oneCol && mv.tabs === 4 && mv.table === 4 && mv.yAxes === 1 && mv.axis <= 1, JSON.stringify(mv));
  await evalJs('document.querySelector(\'.mtabs [data-tab="plt"]\').click(); true'); await sleep(200);
  check('multi-patient tabs switch the count in every row', await evalJs('Array.from(document.querySelectorAll(".mrow path.med")).length === 4 && Array.from(document.querySelectorAll(".mrow path.med")).every(p => p.getAttribute("stroke") === "var(--l-plt)") && document.querySelector(\'.mtabs [data-tab="plt"]\').classList.contains("on")'));
  await evalJs('document.querySelector(\'.mtabs [data-tab="anc"]\').click(); true'); await sleep(150);
  await evalJs(atDay('.mrow[data-pid="DP11"] .plot', 21)); await evalJs(atDay('.mrow[data-pid="DP11"] .plot', 21, 'pointerup')); await sleep(200);
  const mc = await evalJs('({ day: GL.views.clinic.ui.masOf, scrub: document.getElementById("scrubDay").textContent, dp02: GL.views.clinic.multi.pts.find(x => x.p.id === "DP02").eff, title: document.querySelector("#cright .card-head h3").textContent })');
  check('multi-patient: a click on a row sets the shared day; a patient in treatment stays at its data horizon', mc.day === 21 && mc.scrub === 'D+21' && mc.dp02 === 1 && mc.title.includes('D+21'), JSON.stringify(mc));
  check('count charts clip their curves to the plot area', await evalJs('Array.from(document.querySelectorAll(".mrow .plot > svg")).every(v => { const g = v.querySelector("g[clip-path]"); return g && g.querySelector("path.med"); })'));
  check('multi-patient rows: no element outside its row (icons stay 12 px)', await evalJs('Array.from(document.querySelectorAll(".mrow svg.pmark, .mrow-head svg")).every(m => m.getBoundingClientRect().width <= 14)'));
  await shot(PUB, '09_clinic_multi_guest_1440');
  await evalJs('document.querySelector(\'.modesw [data-mode="single"]\').click(); true'); await waitFor(clinicReady);
  check('back to the single-patient view of the last patient', await evalJs('location.hash === "#/clinic/DP10" && !document.querySelector(".ritem .rsel")'));
  check('guest research tabs: no Data tab', await evalJs('(async () => { location.hash = "#/research/model"; await new Promise(r => setTimeout(r, 600)); return !document.querySelector(\'.rnav [data-tab="data"]\'); })()'));
  for (const tab of ['model', 'equations', 'symbols', 'accuracy', 'twins', 'versions', 'interfaces']) {
    await go('#/research/' + tab, '!!document.querySelector(".rmain h1")');
    if (tab === 'twins') await waitFor('!!document.querySelector("#twinBody table")');
    if (tab === 'versions') await waitFor('!!document.querySelector("#verBody table")');
    if (tab === 'model') { await evalJs('document.querySelector(\'[data-node="mpp_n"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200); check('model diagram node shows KaTeX equation', await evalJs('!!document.querySelector("#nodeDetail .katex")'));
      check('research pages: "MarrowTwin (Model3 v1.0)" in the version chip and the navigation', await evalJs('document.querySelector(".rhead .chip-twin").innerText.includes("MarrowTwin (Model3 v1.0)") && document.querySelector(".rnav .foot").innerText.includes("MarrowTwin (Model3 v1.0)")'));
      await evalJs('document.querySelector(\'[data-node="b_mo"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200);
      check('diagram: monocyte node shows both kill terms and θaM', await evalJs(`(() => { const tx = Array.from(document.querySelectorAll('#nodeDetail annotation')).map(a => a.textContent).join(' ');
        return tx.includes(${JSON.stringify(TEX.kc)}) && tx.includes(${JSON.stringify(TEX.am)}) && document.getElementById('nodeDetail').innerText.includes(GL.i18n.t('sym.th_am')); })()`));
      await evalJs('document.querySelector(\'[data-node="ch_n"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200);
      const dg = await evalJs(`(() => { const tx = Array.from(document.querySelectorAll('#nodeDetail annotation')).map(a => a.textContent).join(' '), nd = document.getElementById('nodeDetail').innerText;
        const box = document.querySelector('[data-node="ch_n"]');
        return { cells: box.querySelectorAll('rect.dg-cell').length, store: box.querySelectorAll('rect.dg-cell-store').length, dots: document.querySelectorAll('.dg-node .dg-new').length,
          eq: tx.includes(${JSON.stringify(TEX.store)}) && tx.includes(${JSON.stringify(TEX.kc)}), par: ['n_tr_n', 'nk_n', 'mtt_n'].every(k => nd.includes(GL.i18n.t('sym.' + k))),
          sub: box.textContent.includes(GL.i18n.t('node.ch_n.s', { n: 7, k: 4, s: 3 })), card: !!document.querySelector('.rcard-new') && document.querySelector('.rcard-new').innerText.includes(GL.i18n.t('rm.k4.t')) }; })()`);
      check('diagram: neutrophil chain with 7 cells (3 storage, teal), storage-pool equation, n_N / k_N / MTT_N; "new" marks on the chain and the platelet MPP; card "New in this version"', dg.cells === 7 && dg.store === 3 && dg.dots === 2 && dg.eq && dg.par && dg.sub && dg.card, JSON.stringify(dg));
      await evalJs('document.querySelector(\'[data-node="mpp_pl"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(200);
      check('diagram: platelet MPP node shows e_Z,Pl (teal) and d_ZPl', await evalJs(`(() => { const tx = Array.from(document.querySelectorAll('#nodeDetail annotation')).map(a => a.textContent).join(' ');
        return tx.includes(${JSON.stringify(NEWC + TEX.ezpl + '}')}) && document.getElementById('nodeDetail').innerText.includes(GL.i18n.t('sym.d_zpl')); })()`));
      await evalJs('document.querySelector(\'[data-node="hsc"]\').dispatchEvent(new MouseEvent("click", { bubbles: true })); true'); await sleep(150); }
    if (tab === 'equations') {
      check('equations rendered (KaTeX)', await evalJs('document.querySelectorAll(".eq .katex").length >= 10'));
      // every nested model (5 x 3 clone x 3 route = 45): neutrophil chain with the storage pool (teal, N), Mo/Pl chains with the
      // kill terms (T1, T2) and θaM (plain, Model2), P1 platelet consumption grey, e_Z,Pl teal when the model has the MPP route
      // of Z (route both / _Zm) and grey without it (_Zh); MF has none of these terms
      const nested = await evalJs(`(() => { const c = GL.registry.get('marrowtwin').card, ids = [], T = ${JSON.stringify(TEX)}, NEWC = ${JSON.stringify(NEWC)}, OFFC = ${JSON.stringify(OFFC)};
        ['M0', 'Ma', 'Mb', 'Mab', 'Mabx'].forEach(q => ['', 'C_ddr', 'C_dta'].forEach(cl => ['', '_Zh', '_Zm'].forEach(r => ids.push(q + cl + r))));
        const ok = ids.filter(id => { const sp = GL.model3.spec(id), eq = c.equations(sp), g = (k) => (eq.find(x => x.id === k) || {}).latex || '';
          const n = g('chain_n'), ch = g('chain'), fb = g('feedback');
          return n.includes(T.store) && n.includes(T.div) && n.includes(T.nk) && n.includes(T.kc) && !n.includes(OFFC) &&
            ch.split('-' + T.kc).length === 3 && ch.includes(T.am) && !ch.includes(NEWC + '-') && ch.includes(T.p1) &&
            fb.includes((sp.z_mpp ? NEWC : OFFC) + T.ezpl + '}'); });
        const mf = c.equations(GL.model3.spec('MF')).every(e => !e.latex.includes(T.kc) && !e.latex.includes(T.am) && !e.latex.includes('k_N') && !e.latex.includes('e_{Z,Pl}'));
        const ann = Array.from(document.querySelectorAll('.eq annotation')).map(a => a.textContent);
        return { n: ids.length, ok: ok.length, mf, shown: ann.some(a => a.includes(T.store)) && ann.some(a => a.includes(T.p1)) && !!document.querySelector('.eq-new-key') && document.querySelector('.eq-new-key').innerText === GL.i18n.t('re.new') }; })()`);
      check('equations: all 45 nested models carry the neutrophil storage pool (teal), e_Z,Pl (teal with the MPP route of Z, grey without), P1 grey and the kill terms κc·θaP·C / θaM·C; rendered with the legend; MF without them', nested.ok === nested.n && nested.n === 45 && nested.mf && nested.shown, JSON.stringify(nested));
    }
    if (tab === 'symbols') {
      check('symbol table lists all parameters (49 + 3)', await evalJs('document.querySelectorAll(".tbl tbody tr").length >= 52'));
      const sy = await evalJs(`(() => { const row = (id) => Array.from(document.querySelectorAll('.tbl tbody tr')).find(r => r.cells[1] && r.cells[1].querySelector('.tiny') && r.cells[1].querySelector('.tiny').innerText.trim() === id);
        const P = GL.loader.index().params.struct, out = {};
        for (const [id, unit] of [['kappa_c', '–'], ['th_am', '1/(d·C)']]) { const r = row(id); out[id] = !!r && r.cells[1].innerText.includes(GL.i18n.t('sym.' + id)) && r.cells[2].innerText.trim() === unit &&
          r.cells[4].innerText.includes(GL.i18n.t('kind.estimated')) && r.cells[5].innerText.includes(GL.i18n.t('sym.centre') + ' ' + GL.util.fmt(P.find(x => x.id === id).prior, 3)) && !!r.querySelector('.katex'); }
        return out; })()`);
      check('symbols: kappa_c and th_aM with meaning, unit, centre, estimated', sy.kappa_c && sy.th_am, JSON.stringify(sy));
      const s3 = await evalJs(`(() => { const row = (id) => Array.from(document.querySelectorAll('.tbl tbody tr')).find(r => r.cells[1] && r.cells[1].querySelector('.tiny') && r.cells[1].querySelector('.tiny').innerText.trim() === id);
        const P = GL.loader.index().params.struct, p = (id) => P.find(x => x.id === id), cp = { n: 7, k: 4, s: 3 }, txt = (id) => row(id) ? row(id).innerText : '';
        const est = (id) => txt(id).includes(GL.i18n.t('kind.estimated')), fix = (id) => txt(id).includes(GL.i18n.t('kind.fixed'));
        return { th_cp: fix('th_cp') && +p('th_cp').value === 0 && txt('th_cp').includes(GL.i18n.t('sym.note.p1')) && !!row('th_cp').querySelector('.katex'),
          n_tr_n: fix('n_tr_n') && +p('n_tr_n').value === 7 && txt('n_tr_n').includes(GL.i18n.t('sym.note.chain')), nk_n: fix('nk_n') && +p('nk_n').value === 4,
          d_zpl: est('d_zpl') && p('d_zpl').scale === 'identity' && +p('d_zpl').prior === 0 && txt('d_zpl').includes(GL.i18n.t('sym.centre') + ' 0') && txt('d_zpl').includes(GL.i18n.t('sym.note.p2')),
          gamma: est('gamma') && +p('gamma').box[1] === 1 && txt('gamma').includes(GL.i18n.t('sym.note.stab', cp)),
          mtt_n: est('mtt_n') && +p('mtt_n').box[1] === 15 && txt('mtt_n').includes(GL.i18n.t('sym.note.stab', cp)) }; })()`);
      check('symbols: th_cP (fixed 0, P1 not included), n_N 7 / k_N 4 (design constants), d_ZPl (estimated, centre 0, identity), box limits γ ≤ 1.0 and MTT_N ≤ 15 with the stability note', Object.values(s3).every(Boolean), JSON.stringify(s3));
    }
    if (tab === 'twins') check('twins vs truth: 16 digital patients', await evalJs('document.querySelectorAll("#twinBody tbody tr").length === 16'));
    // K1 of this version: its unit tests (all passed) and calibrations when the bundle has them, else marked pending
    // V1 card: the cohort-level out-of-sample table when the bundle has V1 (marked preliminary while a fold is provisional)
    if (tab === 'accuracy') check('accuracy: V1 table by landmark when exported (preliminary flag as exported)', await evalJs(`(() => { const V = (GL.loader.index().results || {}).v1;
      if (!V) return true; const card = Array.from(document.querySelectorAll('.rcard')).find(c => c.querySelector('h2') && c.querySelector('h2').innerText.startsWith(GL.i18n.t('acc.v1')));
      return !!card && card.querySelectorAll('tbody tr').length === V.metrics.filter(m => m.group === 'all').length && !!card.querySelector('.pending') === !!V.provisional; })()`));
    // model comparison: only runs with the final model's likelihood and data are ranked; MF and refits as notes
    if (tab === 'accuracy') check('accuracy: BIC table ranks candidates and drop-one runs only, final model marked', await evalJs(`(() => { const R = GL.loader.index().results || {}, F = R.fits || [];
      if (!F.length || !F[0].group) return true;
      const card = Array.from(document.querySelectorAll('.rcard')).find(c => c.querySelector('h2') && c.querySelector('h2').innerText.startsWith(GL.i18n.t('acc.bic')));
      const rows = card ? Array.from(card.querySelectorAll('.tbl tbody tr')).filter(r => r.cells.length === 6) : [], ids = rows.map(r => r.cells[0].querySelector('b').innerText);
      return ids.length === F.filter(f => f.group === 'candidate' || f.group === 'drop_one').length && !ids.includes('MF') && !ids.some(i => /_fold|_s0|_w[0-9]/.test(i)) &&
        (!R.final || rows.some(r => r.cells[0].innerText.includes(R.final) && r.querySelector('.ok-badge'))); })()`));
    // step S: VPC table (S5), synthetic cohorts with convergence (S1, S4), result tables; synthetic fits never ranked by BIC
    const sCard = (k) => `Array.from(document.querySelectorAll('.rcard')).find(c => c.querySelector('h2') && c.querySelector('h2').innerText.startsWith(GL.i18n.t('${k}')))`;
    if (tab === 'accuracy') check('accuracy: S5 VPC table (observed, simulated median and 90% range, inside or not) when exported', await evalJs(`(() => { const S = (GL.loader.index().results || {}).s || {}, q = S.S5_vpc || [], card = ${sCard('acc.vpc')};
      if (!q.length) return !!card && !!card.querySelector('.pending'); const rows = card.querySelectorAll('tbody tr');
      return rows.length === q.length && !card.querySelector('.pending') && card.querySelectorAll('.bad-badge').length === q.filter(x => !(x.inside_90 === true || x.inside_90 === 'TRUE')).length; })()`));
    if (tab === 'accuracy') check('accuracy: S1-S4 card lists every synthetic cohort with its fit status, rounds and the result tables', await evalJs(`(() => { const R = GL.loader.index().results || {}, S = R.s || {}, card = ${sCard('acc.syn')}, co = R.cohorts || [];
      if (!card) return false; const conv = [...(S.S1_convergence || []), ...(S.S4_convergence || [])];
      const first = card.querySelector('.tbl'), rows = Array.from(first.querySelectorAll('tbody tr'));
      if (rows.length !== co.length) return false;
      const capOk = conv.every(v => { const r = rows.find(x => x.cells[0].innerText === v.cohort); return r && r.cells[5].innerText.includes(String(v.rounds)) && !!r.cells[5].querySelector('.ok-badge') === (v.converged === true || v.converged === 'TRUE'); });
      const nT = card.querySelectorAll('.tbl').length, want = 1 + ['S1_indiv', 'S1_struct', 'S2_indiv', 'S3_auroc'].filter(k => (S[k] || []).length).length + (S.S1_loss || S.S4_loss ? 1 : 0);
      return capOk && nT === want && !!card.querySelector('.pending') === (want === 1); })()`));
    if (tab === 'accuracy') check('accuracy: synthetic fits kept out of the BIC ranking', await evalJs(`(() => { const F = (GL.loader.index().results || {}).fits || [], card = ${sCard('acc.bic')};
      const syn = F.filter(f => f.group === 'synthetic').map(f => f.id); return !card || !Array.from(card.querySelectorAll('tbody tr b')).some(b => syn.includes(b.innerText)); })()`));
    if (tab === 'accuracy') check('accuracy: K1 of this version shown (all unit tests passed) or marked pending', await evalJs(`(() => { const K = GL.loader.index().k1 || {}, ut = K.unit_tests || [];
      if (!ut.length) return !!document.querySelector(".rmain .pending") && document.body.innerText.includes(GL.i18n.t("acc.k1.p")) && document.body.innerText.includes(GL.i18n.t("acc.cal.cfg"));
      return ut.every(x => x.passed === true || x.passed === 'TRUE') && document.body.innerText.includes(ut.length + '/' + ut.length) && !document.body.innerText.includes(GL.i18n.t("acc.cal.cfg")); })()`));
    await noOldName('research ' + tab + ' en');
    if (tab === 'versions') check('guest versions: digital patients only', await evalJs('Array.from(document.querySelectorAll("#verP option")).every(o => /^DP/.test(o.value))'));
    await shot(PUB, `04_research_${tab}_guest_1440`, true);
  }
  await go('#/research/data', '!!document.querySelector(".rmain h1")');
  check('guest: Data tab not reachable', await evalJs('!location.hash.includes("/data")'));
  await rawKeys('guest en');
  check('browser engine is the Model3 port (45 parameters in the m3_param_names layout, 32 states; population carries th_cp, d_zpl, n_tr_n 7, nk_n 4)', await evalJs(`(() => { const C = GL.model3, st = GL.loader.index().population.struct;
    return C.PAR.length === 45 && C.STATES.length === 32 && C.PAR.indexOf('kappa_c') === 16 && C.PAR.indexOf('th_am') === 17 && C.PAR.indexOf('th_cp') === 18 && C.PAR.slice(-4).join() === 'e_z,e_z_pl,n_tr_n,nk_n' &&
      +st.n_tr_n === 7 && +st.nk_n === 4 && +st.th_cp === 0 && 'd_zpl' in st && GL.loader.index().model_label === 'MarrowTwin (Model3 v1.0)'; })()`));
  for (const lang of ['ja', 'de', 'fr', 'es', 'zh-CN', 'zh-TW']) {
    await evalJs(`GL.i18n.set('${lang}'); true`);
    await go('#/clinic/DP02', clinicReady); await sleep(200); await evalJs('GL.app.render(); true'); await waitFor(clinicReady); await noScroll('guest ' + lang); await rawKeys('clinic ' + lang);
    if (['ja', 'zh-TW'].includes(lang)) { for (const tab of ['model', 'equations', 'symbols', 'accuracy']) { await go('#/research/' + tab, '!!document.querySelector(".rmain h1")'); await noOldName(`research ${tab} ${lang}`); } await go('#/clinic/DP02', clinicReady); }
    if (['ja', 'zh-CN'].includes(lang)) await shot(PUB, `05_clinic_DP02_${lang}_1440`);
  }
  await evalJs(`GL.i18n.set('en'); GL.app.theme('dark'); GL.app.render(); true`); await waitFor(clinicReady); await noScroll('guest dark'); await shot(PUB, '06_clinic_DP02_dark_1440');
  await evalJs(`GL.app.theme('light'); true`);
  const stor = await evalJs('(async () => ({ ls: Object.keys(localStorage), idb: indexedDB.databases ? (await indexedDB.databases()).map(d => d.name) : [] }))()');
  check('browser storage: localStorage glm3.* only, IndexedDB greenlight_model3', stor.ls.length > 0 && stor.ls.every(k => k.startsWith('glm3.')) && stor.idb.includes('greenlight_model3') && !stor.idb.includes('greenlight_model1') && !stor.idb.includes('greenlight_model2'), JSON.stringify({ ls: stor.ls.length, idb: stor.idb }));
  await logout();

  // 3. clinician A
  check('clinician A signs in', await login(U.clinA));
  await go('#/clinic', clinicReady);
  const ra = (await realIdsInRoster()).sort();
  check('clinician A: exactly the 2 assigned real patients + 16 digital', JSON.stringify(ra) === JSON.stringify([pa1, pa2].sort()) && await evalJs('document.querySelectorAll(".ritem").length === 18'));
  await evalJs('document.querySelector(\'.ftabs [data-f="real"]\').click(); true'); await sleep(120);
  check('clinician A: "Real patients" filter shows the 2 real patients', await evalJs('document.querySelectorAll(".ritem").length === 2 && Array.from(document.querySelectorAll(".ritem")).every(r => !/^DP/.test(r.dataset.pid))'));
  await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); document.querySelector(\'.ftabs [data-f="gcsf"]\').click(); true'); await sleep(120);
  const gr = await evalJs('(() => { const real = Object.values(GL.session.real).filter(r => (r.design.gcsf || []).length > 0 || r.summary.gcsf_after_d14).map(r => r.id); const dig = GL.loader.index().roster.filter(r => r.status === "completed" ? r.gcsf.length > 0 : r.gcsf.some(w => w[0] < GL.store.patient(r).horizon)).map(r => r.id); return { got: Array.from(document.querySelectorAll(".ritem")).map(r => r.dataset.pid).sort(), want: real.concat(dig).sort() }; })()');
  check('clinician A: "Received G-CSF" lists exactly the assigned real and digital patients with G-CSF', JSON.stringify(gr.got) === JSON.stringify(gr.want), `${gr.got.length} listed`);
  await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); true'); await sleep(120);
  await evalJs('document.querySelector(\'.ftabs [data-f="all"]\').click(); true'); await sleep(120);
  await go('#/clinic/multi', '!!document.querySelector(".ritem-m")');
  await evalJs('document.getElementById("mAll").click(); true'); await sleep(100);
  await evalJs('document.getElementById("mShow").click(); true'); await waitFor('document.querySelectorAll(".mrow svg path.med").length === 18', 30000);
  const mreal = await evalJs('({ rows: document.querySelectorAll(".mrow").length, real: Array.from(document.querySelectorAll(".mrow")).filter(r => !/^DP/.test(r.dataset.pid)).length, ls: Object.keys(localStorage).filter(k => /\\bP\\d{2}\\b/.test(k) || /\\bP\\d{2}\\b/.test(localStorage.getItem(k) || "")).length })');
  check('clinician A multi-patient: all 18 listed patients incl. the 2 real ones; the selection is not in browser storage', mreal.rows === 18 && mreal.real === 2 && mreal.ls === 0, JSON.stringify(mreal));
  await noScroll('clinician multi-patient 1440x900');
  await shot(PRIV, '13_clinic_multi_clinicianA_1440');
  await go('#/clinic/' + pb, clinicReady);
  check('clinician A: B\'s patient is not reachable', await evalJs(`!location.hash.endsWith(${JSON.stringify('/' + pb)}) && !new RegExp('\\\\b' + ${JSON.stringify(pb)} + '\\\\b').test(document.body.innerText)`));
  await go('#/clinic/' + pa1, clinicReady);
  await noScroll('clinician real 1440x900');
  check('watermark with user name and access banner over the patient view', await evalJs(`!!document.getElementById('wm') && document.getElementById('wm').dataset.text.startsWith(${JSON.stringify(U.clinA.u)}) && !!document.querySelector('.authbar')`));
  await shot(PRIV, '10_clinic_real_clinicianA_1440');
  await evalJs(`(() => { const s = document.getElementById('scSvg'), r = s.getBoundingClientRect(); s.dispatchEvent(new PointerEvent('pointerdown', { bubbles: true, clientX: r.left + r.width * 0.25, clientY: r.top + 30 })); s.dispatchEvent(new PointerEvent('pointerup', { bubbles: true })); return true; })()`); await sleep(300);
  await evalJs('document.getElementById("scNext").click(); true'); await sleep(200);
  await evalJs('document.getElementById("btnWhatIf").click(); true'); await sleep(250);
  await evalJs('document.querySelector(\'[data-preset="gcsf_reactive"]\').click(); true'); await sleep(200);
  await evalJs('document.getElementById("wiRun").click(); true'); await waitFor('!!document.querySelector("#wiOut .stack")', 30000);
  await evalJs('document.getElementById("wiSave").click(); true'); await waitFor('!!document.querySelector("#wiSaved h3")', 5000);
  await shot(PRIV, '11_clinic_real_whatif_1440');
  const store = await evalJs(`(async () => { const d = await GL.store.db(); const sc = await new Promise(r => { const o = []; const c = d.transaction('scenarios').objectStore('scenarios').openCursor(); c.onsuccess = () => { const x = c.result; if (x) { o.push(x.value.patient); x.continue(); } else r(o); }; });
    return { mem: GL.session.scenarios.map(s => s.patient), idb: sc, ls: Object.keys(localStorage).filter(k => /\\bP\\d{2}\\b/.test(k) || /\\bP\\d{2}\\b/.test(localStorage.getItem(k) || '')) }; })()`);
  check('real what-if scenario only in memory; no real patient in IndexedDB or localStorage', store.mem.includes(pa1) && !store.idb.some(p => !/^DP/.test(p)) && store.ls.length === 0);
  const audit = await evalJs(`(async () => (await GL.audit.all()).filter(e => e.subject === ${JSON.stringify(pa1)}).map(e => e.action))()`);
  check('audit log records open, scrub and the supportive-care simulation (what-if)', ['open', 'scrub', 'what-if'].every(a => audit.includes(a)), audit.join(','));
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
