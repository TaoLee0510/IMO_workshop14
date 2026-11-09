# GreenLight · MarrowTwin (Model2) — privacy and access control

Model1 plan section 13 (v0.4), decisions A1–A7 (confirmed 2026-11-04 16:35), carried over unchanged to
the Model2 bundle (Model2 plan 6.5, decision C3): `data/users.js` and the steward's keystore and
pseudonym HMAC key were copied from Model1, so the accounts and the research pseudonyms are the same in
both bundles. The two bundles are separate: own sealed files, own browser storage (`glm2.*`,
IndexedDB `greenlight_model2`, own audit log). This page describes
the technical safeguards. **Whether a deployment is compliant is decided by the institution's
privacy / compliance office and the IRB; they must review it before real use.**

## 1. Data classification

| class | content | where | protection |
|---|---|---|---|
| Public | model description, equations, parameters, aggregate K1 / fit results, the 16 digital patients (simulated, incl. their truth) | `data/model2_index.js`, `data/patients/DP*.js` | none needed |
| Clinical (per real patient) | all fields of the de-identified research record: study ID (P01…), exact age, sex, prior therapy, bridging regimen, LD labs, CRS grade, drugs, daily counts, ferritin, G-CSF windows, model forecasts, strategies, what-if inputs | `data/sealed/<random id>.js` | AES-256-GCM, one key per patient |
| Research (per real patient) | minimised at generation time (table 3): pseudonym, bands, daily counts, phenotype, full model parameters and forecasts | `data/sealed/<random id>.js` | AES-256-GCM, one key per patient |
| Research cohort aggregate | counts, cross-tables (cells 1–4 shown as "<5"), sampling days by pseudonym | `data/sealed/<random id>.js` | AES-256-GCM |
| Steward secrets | data keys, file-id map, pseudonym HMAC key, clinician assignments; staged plaintext payloads | `Model2/Results/real/greenlight_admin/`, `Model2/Results/real/greenlight_stage/` (outside GreenLight, mode 600) | never distributed with GreenLight |

The source data are already de-identified (study codes, days relative to infusion, no dates,
names or record numbers). They remain patient-level health data and are handled on a
minimum-necessary basis.

## 2. Roles

| role | real patients | portals |
|---|---|---|
| Clinician | assigned patients only (care-team list), clinical view; or, when the steward grants it with `--all-patients`, every real patient (also patients sealed later), shown as `ALL` in `users.js list` | clinician portal (assigned real + all digital patients); research tabs without patient data: Model, Equations, Symbols & parameters, Interfaces |
| Researcher | all 35, research view only | research portal (all tabs incl. Data, Fit & accuracy, Twins vs truth, Versions); clinician portal with digital patients only |
| Data steward | none | user list and audit-log export (`#/admin`); users, roles and assignments via `tools/users.js` |
| Guest | none, no login | both portals, digital patients only |

One person with two roles uses two accounts; a session cannot switch roles (A5).

## 3. Field matrix (same patient, two views)

| information | clinician (assigned patient) | researcher |
|---|---|---|
| patient ID | study ID (P01) | pseudonym R-XXXXXX = HMAC-SHA256(steward key, study ID), first 6 hex; not reversible without the key; list sorted by pseudonym |
| age | exact | 10-year band; ≥ 90 → "90+" |
| sex | shown | shown |
| prior lines, ASCT, bridging | exact; regimen text | lines 1–2 / 3–4 / ≥ 5; ASCT yes/no; bridging class (chemotherapy / none) only |
| LD labs (ferritin, LDH, CRP, Hb, PLT, ANC, eGFR) | exact | eGFR < 60 / 60–89 / ≥ 90; the others by cohort tertiles (T1–T3) |
| CRS grade · tocilizumab, steroids | shown · shown | shown · not included |
| daily counts (days from infusion) | all | all |
| G-CSF windows, ferritin series, CAR-HEMATOTOX | shown | not included |
| recovery phenotype | shown | shown |
| forecasts (P phenotype, recovery day, bands) | shown | shown |
| individual parameters ψ̂ | clinical interpretation only | full values (posterior mode, Cholesky factor) |
| attribution (when / who) | clinical wording | shown |
| treatment strategies, what-if | shown | **not included** (no plans, no what-if inputs in the payload) |
| research aggregate tables | — | cells 1–4 shown as "<5" |

Minimisation happens when the payload is generated (`Model2/Code/R/an_G_export.R`,
`g_research_record()`, `g_cohort_research()`), not by hiding fields in the browser: the research
payload never contains the exact fields.

## 4. Technical controls

- **Encryption at rest.** Every real-patient payload has its own random AES-256-GCM key; sealed files
  are named by random ids and contain only ciphertext (`tools/seal.js`).
- **Key hierarchy.** Each user has an RSA-OAEP-2048 key pair. The private key is encrypted with a
  key derived from the password by PBKDF2-SHA256 (600 000 iterations, 16-byte random salt). Each data
  key the user is entitled to is wrapped with the user's public key and listed under its random file
  id. `data/users.js` holds user names, roles, salts, iteration counts, encrypted private keys, public
  keys and wrapped keys — no passwords, no patient ids. Because wrapping uses the public key, the
  steward can assign, revoke and re-seal without knowing any password.
- **Passwords** are typed by the user (hidden prompt in `tools/users.js`, at least 12 characters) and
  never written anywhere. A wrong password fails the AES-GCM authentication of the private key.
- **In the browser.** WebCrypto opens only the payloads wrapped for the account (works from
  `file://` in Chrome, a secure context; if a browser refuses, serve the folder via `http://localhost`).
  Decrypted data live only in memory (`GL.session`).
- **Session.** Idle lock after 15 minutes without input (password required); lock and sign-out drop all
  decrypted data, real-patient UI state and real-patient what-if scenarios from memory.
- **Browser storage.** localStorage: UI preferences (language, theme, filters, playback) and the demo
  state of digital patients only. IndexedDB: what-if scenarios of digital patients and the audit log.
  Anything older versions stored about real patients is purged at start.
- **Audit.** Hash-chained log (SHA-256 over the previous hash and the entry): user, role, patient id or
  pseudonym, action (login, logout, lock, unlock, open, scrub, play, what-if, export), timestamp. The
  steward verifies the chain and exports it as JSON (`#/admin`).
- **Deterrents.** Persistent banner "Authorized use only – access is logged"; semi-transparent watermark
  (user name + time) over patient views.
- **Revocation.** Removing an entitlement (unassign, role change, revoke) gives the affected payloads a
  new key and a new file id and re-wraps them for the remaining users; the old files are deleted.

## 5. Limitations

- **No instant revocation without a server.** A bundle that was already copied keeps working with the
  old keys for the person who had them. Revocation takes effect for newly distributed bundles only.
- **The audit log is local.** It can be exported or deleted by whoever controls the browser profile;
  the hash chain detects alteration inside the log but cannot prevent deletion (deleting the newest
  entries is not detectable).
- **Screenshots and photos cannot be prevented**; the watermark only deters.
- **Small-cell masking** is applied to cells, not to marginal totals; combinations of a masked cell
  with row or column totals can narrow it down. Further suppression should be agreed with the privacy
  office before publication of research tables.
- **Pseudonyms** are unlinkable only while the HMAC key in `greenlight_admin/` stays with the steward.
- **Strength of a password** is the strength of that account's access; PBKDF2 slows guessing but does
  not replace strong passwords.
- **The staged plaintext** (`Results/real/greenlight_stage/`) and the admin keystore are protected only
  by file permissions on the steward's machine.

## 6. Deployment target

For real use: a server behind the institution's single sign-on (SSO) that returns each role's payload
over HTTPS and keeps the audit log server-side. The REST stubs (`api/plumber_stub.R`,
`api/fastapi_stub.py`) and the data adapter / remote-module interfaces already define the calls; the
front end changes little (the session layer would request payloads instead of decrypting them).
Compliance (HIPAA minimum necessary §164.502(b) / §164.514, security rule access control, audit and
encryption §164.312, or the applicable local rules) must be confirmed by the institution's privacy
office and the IRB.

## 7. Operating it

```bash
# 1. regenerate the staged payloads (steward machine, from Model2/)
Rscript -e 'source("Code/config.R"); for (f in list.files("Code/R", full.names=TRUE)) source(f); cfg <- make_config("real"); load_model2(cfg); g_export(cfg, cores = 8, n_draws = 60L, n_cf = 24L)'
# 2. seal (keeps existing keys; --rotate gives every payload a new key)
node GreenLight/tools/seal.js
# 3. accounts (passwords typed by the user at the prompt)
node GreenLight/tools/users.js add <user> --role steward
node GreenLight/tools/users.js add <user> --role clinician --patients P01,P02
node GreenLight/tools/users.js add <user> --role clinician --all-patients      # every real patient, also future ones
node GreenLight/tools/users.js assign <user> all                              # grant all to an existing clinician
node GreenLight/tools/users.js unassign <user> all                            # withdraw all (keys rotated)
node GreenLight/tools/users.js add <user> --role researcher
node GreenLight/tools/users.js list
```
Distribute the `Model2/GreenLight/` folder only (never `Results/real/greenlight_*`). An account change made
with `Model2/GreenLight/tools/users.js` applies to the Model2 bundle only (the frozen Model1 bundle keeps its own copy).
