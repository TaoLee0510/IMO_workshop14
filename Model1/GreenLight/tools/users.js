/* GreenLight · tools/users.js — user administration for the data steward (Model1 plan 13.6).
 * Passwords are typed interactively (hidden) and never stored: only a salt, an encrypted private key and
 * the wrapped data keys go into GreenLight/data/users.js. Removing an entitlement (unassign, role change,
 * revoke) rotates the affected data keys and re-seals those files, so old copies of a user's keys stop
 * working for new bundles (files already copied elsewhere cannot be recalled; docs/PRIVACY.md).
 *
 *   node tools/users.js add <user> --role clinician|researcher|steward [--patients P01,P02 | --all-patients]
 *   node tools/users.js role <user> <role>
 *   node tools/users.js assign <user> P01,P02 | all
 *   node tools/users.js unassign <user> P03 | all
 * --all-patients (or the list "all"): the clinician sees every real patient, also patients sealed later;
 * unassigning single patients from such a clinician turns the grant into the explicit list of the others.
 *   node tools/users.js passwd <user>            (new password = new key pair; entitlements kept)
 *   node tools/users.js revoke <user>            (removes the user, rotates every key the user held)
 *   node tools/users.js list
 * --password-stdin reads the password from standard input instead of the prompt (automated tests only).
 * Run tools/seal.js first (the keystore must exist before users get keys).
 */
const A = require('./lib_access.js');

const argv = process.argv.slice(2);
const cmd = argv[0], name = argv[1];
const opt = (f) => { const i = argv.indexOf(f); return i >= 0 ? argv[i + 1] : null; };
const fromStdin = argv.includes('--password-stdin');
const die = (m) => { console.error(m); process.exit(1); };

function readHidden(prompt) {
  return new Promise((resolve) => {
    const stdin = process.stdin;
    process.stdout.write(prompt);
    stdin.setRawMode(true); stdin.resume(); stdin.setEncoding('utf8');
    let s = '';
    const on = (ch) => {
      for (const c of ch) {
        if (c === '\r' || c === '\n') { stdin.setRawMode(false); stdin.pause(); stdin.removeListener('data', on); process.stdout.write('\n'); resolve(s); return; }
        if (c === '\u0003') { process.stdout.write('\n'); process.exit(130); }
        if (c === '\u007f' || c === '\b') s = s.slice(0, -1); else s += c;
      }
    };
    stdin.on('data', on);
  });
}
async function askPassword() {
  if (fromStdin || !process.stdin.isTTY) {
    const chunks = []; for await (const c of process.stdin) chunks.push(c);
    return Buffer.concat(chunks).toString('utf8').split(/\r?\n/)[0];
  }
  const a = await readHidden('Password (min 12 characters): '), b = await readHidden('Repeat password: ');
  if (a !== b) die('passwords do not match');
  return a;
}
const checkPw = (pw) => { if (!pw || pw.length < 12) die('password must have at least 12 characters'); };

function save(users, assign, ks, rotated) {
  if (rotated && rotated.length) A.rotate(ks, rotated);
  A.saveKeystore(ks); A.saveAssign(assign);
  A.rewrapAll(users, assign, ks);
  A.saveUsers(users);
}
const allIds = (ks) => Object.keys(ks.payloads).filter(p => p.startsWith('clinical:')).map(p => p.slice(9)).sort();
const parseIds = (s, ks) => {
  if ((s || '').trim().toLowerCase() === 'all') return ['*'];
  const ids = (s || '').split(',').map(x => x.trim()).filter(Boolean);
  const bad = ids.filter(id => !ks.payloads['clinical:' + id]);
  if (bad.length) die('unknown patient id(s): ' + bad.join(', ') + ' (sealed patients only)');
  return ids;
};

(async () => {
  const users = A.loadUsers(), assign = A.loadAssign(), ks = A.loadKeystore();
  const find = (u) => users.users.find(x => x.u === u);
  if (!cmd || cmd === 'help') { console.log(require('fs').readFileSync(__filename, 'utf8').split('*/')[0]); return; }
  if (cmd === 'list') {
    if (!users.users.length) { console.log('no users'); return; }
    for (const u of users.users) console.log(`${u.u.padEnd(20)} ${u.role.padEnd(11)} ${String(u.files.length).padStart(3)} file key(s)` + (u.role === 'clinician' ? `  patients: ${(assign[u.u] || []).includes('*') ? 'ALL (current and future)' : (assign[u.u] || []).join(',') || '-'}` : ''));
    return;
  }
  if (!name) die('user name missing');
  if (!Object.keys(ks.payloads).length) die('no keystore yet: run node tools/seal.js first');
  if (cmd === 'add') {
    if (!A.validName(name)) die('user name: 3-32 characters A-Z a-z 0-9 . _ -');
    if (find(name)) die('user exists: ' + name);
    const role = opt('--role');
    if (!A.ROLES.includes(role)) die('--role must be one of ' + A.ROLES.join(', '));
    const pw = await askPassword(); checkPw(pw);
    users.users.push(Object.assign({ u: name, role, files: [] }, A.makeCredentials(pw)));
    if (role === 'clinician') assign[name] = argv.includes('--all-patients') ? ['*'] : parseIds(opt('--patients'), ks);
    save(users, assign, ks);
    console.log(`added ${name} (${role})`);
    return;
  }
  const u = find(name);
  if (!u) die('no such user: ' + name);
  const before = new Set(A.entitlements(u, assign, ks));
  const lost = () => [...before].filter(pid => !new Set(A.entitlements(u, assign, ks)).has(pid));
  if (cmd === 'role') {
    const role = argv[2];
    if (!A.ROLES.includes(role)) die('role must be one of ' + A.ROLES.join(', '));
    u.role = role;
    if (role !== 'clinician') delete assign[name]; else assign[name] = assign[name] || [];
    save(users, assign, ks, lost());
    console.log(`${name}: role ${role}`);
  } else if (cmd === 'assign' || cmd === 'unassign') {
    if (u.role !== 'clinician') die('patients are assigned to clinicians only');
    const ids = parseIds(argv[2], ks);
    const had_all = (assign[name] || []).includes('*');
    if (ids.includes('*')) assign[name] = cmd === 'assign' ? ['*'] : [];
    else {
      const cur = new Set(had_all ? allIds(ks) : (assign[name] || []));
      ids.forEach(id => cmd === 'assign' ? cur.add(id) : cur.delete(id));
      assign[name] = had_all && cmd === 'assign' ? ['*'] : [...cur].sort();
    }
    save(users, assign, ks, lost());
    console.log(`${name}: ` + (assign[name].includes('*') ? 'ALL patients (current and future)' : `${assign[name].length} assigned patient(s)`) +
      (cmd === 'unassign' ? '; keys of the removed patients rotated' : ''));
  } else if (cmd === 'passwd') {
    const pw = await askPassword(); checkPw(pw);
    Object.assign(u, A.makeCredentials(pw));
    save(users, assign, ks);
    console.log(`${name}: password reset (new key pair)`);
  } else if (cmd === 'revoke') {
    const held = [...before];
    users.users = users.users.filter(x => x.u !== name);
    delete assign[name];
    save(users, assign, ks, held);
    console.log(`${name}: removed; ${held.length} data key(s) rotated`);
  } else die('unknown command: ' + cmd);
})().catch(e => die(String(e && e.message || e)));
