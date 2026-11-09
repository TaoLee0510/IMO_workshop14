/* GreenLight · tools/test_access.js — access-control test of the sealed bundle (Model1 plan 13.6 item 11, carried over by Model2 plan 6.5).
 * Creates temporary accounts (random names and random passwords held in memory only, given to
 * tools/users.js on stdin), checks, and removes them again (always, also on failure):
 *   - sealed files and users.js contain no plaintext patient id and no payload field name
 *   - a wrong password does not open an account
 *   - clinician A opens exactly the clinical payloads of A's patients, never B's
 *   - the researcher opens only research payloads: pseudonyms (not clinical ids), no exact age, no bridging
 *     regimen, no drugs, no plans / what-if inputs; aggregate cells 1-4 shown as "<5"
 *   - unassigning a patient rotates its key: A's old wrapped key no longer opens the new file
 *   - research pseudonyms are HMAC-SHA256(steward key, clinical id) (the key carried over from Model1, so the
 *     pseudonyms are the same in both bundles)
 *   - the existing (non-temporary) accounts keep their role and credentials
 * Temporary accounts left by an aborted earlier run (names t_<role>_<6 hex>) are removed first.
 * Prints no patient-level data. Usage: node tools/test_access.js
 */
const fs = require('fs'), path = require('path'), crypto = require('crypto'), { execFileSync } = require('child_process');
const A = require('./lib_access.js');
const results = []; const check = (name, ok, info) => { results.push(ok); console.log((ok ? 'OK  ' : 'FAIL') + ' ' + name + (info ? ' — ' + info : '')); };
const tag = crypto.randomBytes(3).toString('hex');
const pw = () => crypto.randomBytes(18).toString('base64');
const users = { clinA: { u: 't_clinA_' + tag, pw: pw() }, clinB: { u: 't_clinB_' + tag, pw: pw() }, res: { u: 't_res_' + tag, pw: pw() }, stew: { u: 't_stew_' + tag, pw: pw() },
                clinAll: { u: 't_clinAll_' + tag, pw: pw() } };
const cli = (args, input) => execFileSync(process.execPath, [path.join(__dirname, 'users.js')].concat(args), { input: input === undefined ? '' : input + '\n', stdio: ['pipe', 'pipe', 'pipe'] }).toString();
const findUser = (u) => A.loadUsers().users.find(x => x.u === u);
const TEMP = /^t_(clinA|clinB|res|stew|clinAll)_[0-9a-f]{6}$/;
const stale = A.loadUsers().users.map(x => x.u).filter(u => TEMP.test(u));
stale.forEach(u => cli(['revoke', u]));
if (stale.length) console.log(`removed ${stale.length} temporary account(s) left by an aborted earlier run`);
// the user's own accounts: role and credentials must be unchanged by this test
const creds = () => JSON.stringify(A.loadUsers().users.filter(x => !TEMP.test(x.u)).map(x => [x.u, x.role, x.salt, x.iter, x.ekey, x.pub]));
const credsBefore = creds();

let created = [];
try {
  const ks = A.loadKeystore();
  const clinIds = Object.keys(ks.payloads).filter(p => p.startsWith('clinical:')).map(p => p.slice(9)).sort();
  if (clinIds.length < 3) throw new Error('need sealed patients: run step G and tools/seal.js');
  const [pa1, pa2, pb] = clinIds;

  // 1. no plaintext ids / field names in what is distributed
  const fieldNames = ['summary', 'asof', 'obs', 'bridging_regimen', 'ferritin', 'label', 'design', 'patient_id', 'age_band', 'pseudo', 'clinical', 'research'];
  const dist = fs.readdirSync(A.P.sealed).map(f => fs.readFileSync(path.join(A.P.sealed, f), 'utf8'));
  const idRe = /\bP\d{2}\b/;
  // a sealed file is exactly GLM2.addSealed("<hex fid>","<base64 iv>","<base64 ciphertext>"); nothing else
  check('sealed files: ciphertext only, no plaintext patient id', dist.every(s => /^GLM2\.addSealed\("[0-9a-f]{24}","[A-Za-z0-9+/=]{16}","[A-Za-z0-9+/=]+"\);\n$/.test(s) && !/"P\d{2}"/.test(s)), dist.length + ' files');
  check('sealed files: no payload field names', dist.every(s => fieldNames.every(f => !s.includes('"' + f + '"') && !s.includes(f + ':'))));
  const plainLeft = fs.readdirSync(A.P.patients).filter(f => /^P\d+\.js$/.test(f));
  check('no real-patient plaintext in GreenLight/data/patients', plainLeft.length === 0);
  // pseudonyms: HMAC-SHA256 with the steward-held key (Results/real/greenlight_admin), first 6 hex digits
  const pkey = Buffer.from(fs.readFileSync(path.join(A.P.admin, 'pseudonym_key.hex'), 'utf8').trim(), 'hex');
  const pseud = (id) => 'R-' + crypto.createHmac('sha256', pkey).update(id).digest('hex').slice(0, 6).toUpperCase();
  check('research pseudonyms = HMAC(steward key, clinical id) for every patient', clinIds.every(id =>
    JSON.parse(fs.readFileSync(path.join(A.P.stage, 'research', id + '.json'), 'utf8')).id === pseud(id)), clinIds.length + ' patients');
  const idxTxt = fs.readFileSync(path.join(A.P.gl, 'data', 'model2_index.js'), 'utf8');
  check('public index: no real patient id', !idRe.test(idxTxt));

  // 2. temporary accounts
  cli(['add', users.clinA.u, '--role', 'clinician', '--patients', pa1 + ',' + pa2, '--password-stdin'], users.clinA.pw); created.push(users.clinA.u);
  cli(['add', users.clinB.u, '--role', 'clinician', '--patients', pb, '--password-stdin'], users.clinB.pw); created.push(users.clinB.u);
  cli(['add', users.res.u, '--role', 'researcher', '--password-stdin'], users.res.pw); created.push(users.res.u);
  cli(['add', users.stew.u, '--role', 'steward', '--password-stdin'], users.stew.pw); created.push(users.stew.u);
  const utxt = fs.readFileSync(A.P.users, 'utf8');
  // patient ids could only appear as plain strings; base64 blobs (salts, keys, ciphertexts) may contain a "P<nn>" pattern by chance
  const ustr = []; (function walk(o) { if (typeof o === 'string') ustr.push(o); else if (o && typeof o === 'object') for (const [k, v] of Object.entries(o)) { ustr.push(k); walk(v); } })(A.loadUsers());
  const plainStr = ustr.filter(x => !(x.length >= 16 && /^[A-Za-z0-9+/=]+$/.test(x)));
  check('users.js: no passwords and no patient ids', Object.values(users).every(x => !utxt.includes(x.pw)) && !/"P\d{2}"/.test(utxt) && !plainStr.some(x => idRe.test(x)), plainStr.length + ' plain strings');

  // 3. wrong password
  let wrongOpened = true; try { A.openPrivateKey(findUser(users.clinA.u), users.clinA.pw + 'x'); } catch (e) { wrongOpened = false; }
  check('wrong password does not open the account', !wrongOpened);

  // 4. clinicians
  const pA = A.openPayloads(findUser(users.clinA.u), users.clinA.pw), pB = A.openPayloads(findUser(users.clinB.u), users.clinB.pw);
  const idsA = pA.map(p => p.id).sort(), idsB = pB.map(p => p.id);
  check('clinician A opens exactly the assigned clinical payloads', pA.every(p => p.kind === 'clinical') && JSON.stringify(idsA) === JSON.stringify([pa1, pa2].sort()));
  check('clinician A cannot see clinician B\'s patient (and vice versa)', !idsA.includes(pb) && idsB.length === 1 && idsB[0] === pb &&
    !findUser(users.clinA.u).files.some(e => findUser(users.clinB.u).files.some(f => f.f === e.f)));

  // 5. researcher
  const pR = A.openPayloads(findUser(users.res.u), users.res.pw);
  const rec = pR.filter(p => p.kind === 'research'), coh = pR.find(p => p.kind === 'cohort');
  const rtxt = JSON.stringify(pR);
  check('researcher opens only research payloads (+ cohort aggregate)', rec.length === clinIds.length && !!coh && pR.every(p => p.kind === 'research' || p.kind === 'cohort'));
  check('research ids are pseudonyms, none equals a clinical id', rec.every(p => /^R-[0-9A-F]{6}$/.test(p.id)) && !idRe.test(rtxt) && new Set(rec.map(p => p.id)).size === rec.length);
  const clinTxt = clinIds.map(id => fs.readFileSync(path.join(A.P.stage, 'clinical', id + '.json'), 'utf8')).map(JSON.parse);
  const regimens = [...new Set(clinTxt.map(c => c.summary.bridging_regimen).filter(r => r && r !== 'none'))];
  check('research payload: no exact age, regimen, drugs, eGFR value, design or G-CSF windows', !rtxt.includes('"age"') && regimens.every(r => !rtxt.includes(r)) &&
    !/"tocilizumab"|"steroids"|"design"|"egfr":\d|"ferritin_ld"|"ldh_ld"|"gcsf"/.test(rtxt));
  check('research payload: ages in 10-year bands (90+)', rec.every(p => /^\d{2}-\d{2}$|^90\+$/.test(p.summary.age_band)));
  check('research payload: no treatment recommendations or what-if inputs', !/"plans"|"top3"|"drv"|"kern"/.test(rtxt));
  const cells = Object.values(coh.crosstab).flatMap(t => Object.values(t).flatMap(r => Object.values(r)));
  check('aggregate cells 1-4 shown as "<5"', cells.every(c => c === '<5' || c === '0' || +c >= 5) && cells.includes('<5'));
  check('steward opens no patient data', A.openPayloads(findUser(users.stew.u), users.stew.pw).length === 0);

  // 6. revocation by key rotation
  const before = findUser(users.clinA.u);
  const oldMap = Object.fromEntries(A.openPayloads(before, users.clinA.pw).map((p, i) => [p.id, before.files[i]]));
  cli(['unassign', users.clinA.u, pa2]);
  const ks2 = A.loadKeystore(), newFid = ks2.payloads['clinical:' + pa2].fid;
  let oldKeyWorks = false;
  try {
    const pk = A.openPrivateKey(findUser(users.clinA.u), users.clinA.pw);
    const k = crypto.privateDecrypt({ key: pk, padding: crypto.constants.RSA_PKCS1_OAEP_PADDING, oaepHash: 'sha256' }, A.unb64(oldMap[pa2].k));
    const s = A.readSealed(newFid); A.aesDecrypt(k, s.iv, s.ct); oldKeyWorks = true;
  } catch (e) { oldKeyWorks = false; }
  check('unassign rotates the key: the old key does not open the new file', newFid !== oldMap[pa2].f && !fs.existsSync(A.sealedFile(oldMap[pa2].f)) && !oldKeyWorks &&
    A.openPayloads(findUser(users.clinA.u), users.clinA.pw).map(p => p.id).join() === pa1);
  check('clinician B unaffected by the rotation', A.openPayloads(findUser(users.clinB.u), users.clinB.pw).map(p => p.id).join() === pb);

  // 7. clinician with access to all real patients (--all-patients): every clinical payload, no research payload,
  //    also patients sealed later (the grant is '*', not a list); unassigning one turns it into the explicit rest
  cli(['add', users.clinAll.u, '--role', 'clinician', '--all-patients', '--password-stdin'], users.clinAll.pw); created.push(users.clinAll.u);
  const pAll = A.openPayloads(findUser(users.clinAll.u), users.clinAll.pw);
  check('all-patients clinician opens every clinical payload and nothing else', pAll.every(p => p.kind === 'clinical') &&
    JSON.stringify(pAll.map(p => p.id).sort()) === JSON.stringify(clinIds));
  check('all-patients grant is stored as "*" (covers patients sealed later)', JSON.stringify(A.loadAssign()[users.clinAll.u]) === '["*"]' &&
    /ALL \(current and future\)/.test(cli(['list'])));
  cli(['unassign', users.clinAll.u, pb]);
  const restIds = A.openPayloads(findUser(users.clinAll.u), users.clinAll.pw).map(p => p.id).sort();
  check('unassigning one patient from an all-patients clinician leaves the explicit rest', !restIds.includes(pb) &&
    restIds.length === clinIds.length - 1 && !A.loadAssign()[users.clinAll.u].includes('*'));
  check('clinician B keeps access after the rotation caused by the all-patients clinician', A.openPayloads(findUser(users.clinB.u), users.clinB.pw).map(p => p.id).join() === pb);
} catch (e) {
  check('test ran', false, String(e && e.message || e));
} finally {
  for (const u of created) { try { cli(['revoke', u]); } catch (e) { console.log('cleanup failed for a temporary account'); results.push(false); } }
  const left = A.loadUsers().users.filter(x => x.u.endsWith('_' + tag));
  check('temporary accounts removed', left.length === 0);
  check('existing accounts unchanged (role, salt, encrypted key, public key)', creds() === credsBefore);
}
const nf = results.filter(x => !x).length;
console.log(nf ? `FAIL: ${nf} of ${results.length} checks failed` : `PASS: ${results.length} checks`);
process.exit(nf ? 1 : 0);
