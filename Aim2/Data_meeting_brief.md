# Aim 2 (post-CAR T cytopenia): data request and phenotype sign-off

For the Tuesday 11/03 data meeting. Part A is for the data owner; Part B is for the clinicians.

**Timing.** If Part A1 arrives by **Wednesday 11/04, 12:00**, the results can go into the Friday slides. The code is already built and tested on a synthetic dataset in exactly this format.

---

## Part A: data request

### A1. Needed for this week, in priority order

Please send one CSV per table. Use study IDs and **days relative to CAR T infusion** (day 0). Do not include dates or names.

| # | Table | One row per | Key columns |
|---|---|---|---|
| 1 | Blood counts | draw | patient_id, day, **ANC, PLT**, Hb, WBC, ALC, **MCV, RDW** (if available). From the last pre-LD value to day +90 or later |
| 2 | Inflammation | draw | patient_id, day, CRP, ferritin (IL-6 if available), days 0–30 at least |
| 3 | Interventions | event | patient_id, event, start_day, end_day, dose. Events: G-CSF, PLT and RBC transfusion, TPO-RA, **CD34+ boost**, tocilizumab, steroids, anakinra |
| 4 | Patients | patient | age, sex, product, prior lines, prior ASCT, prior bendamustine, months since last chemo, bridging type, LD dose reduced, CRS grade and onset/end day, marrow lymphoma involvement, IEC-HS/HLH, viral infection, relapse day, next-therapy day, death day, last follow-up day |
| 5 | Pre-treatment sequencing (the ~15) | tested patient + variant | sample day, timepoint, sample type, panel, limit of detection; for each variant: gene, protein change, VAF, depth. **Please list tested-negative patients too.** "Tested, nothing found" is information; "not tested" is not |

The full column list, with a synthetic example in the same format, is in `Aim2/TestData/` (`README.md`).

We will create the 70/30 train/hold-out split ourselves. It is made from the patient list at receipt, with a recorded random seed, before anyone looks at outcomes.

### A2. Three questions that do not need new sequencing

1. **Are MCV and RDW in the CBC?** If so, for how many patients is there a pre-lymphodepletion value?
   - Why: RDW is the blood index most consistently linked to CH, and CH also shifts other red-cell indices (Br J Haematol 2025, PMID 40754673). The MCV link is weaker.
   - These two numbers are the only clone-burden proxy available for the ~335 unsequenced patients, and the model uses them automatically if present.
2. **Did any patient get a marrow work-up for persistent cytopenia that included a clinical myeloid NGS panel?** If so: how many patients, the sample day, and the gene/VAF results.
   - These are post-treatment and selected, so we use them only to test predictions, not to fit the model. See A3.
3. **For the ~15 sequenced patients:**
   - Which timepoint was sampled: apheresis, after bridging, or before lymphodepletion?
   - Is leftover apheresis material or DNA stored for other patients?
   - Why it matters: if stored material exists, pre-treatment CH status could be added for the whole cohort. That changes what Aim 2 can claim.

### A3. Post-treatment mutation results: please hold

For the ~3 patients with post-treatment sequencing (and any found under A2.2):

- Send **only the patient IDs and sample days** now.
- Keep the VAF results until we confirm that our predictions for these patients are frozen. We will send a SHA-256 fingerprint of the prediction file.
- After that, send the results. The comparison is then a genuine test rather than a fit after the fact.
- Reason: these patients were most likely sampled *because* their counts stayed low. They cannot estimate how common anything is. They can only show whether the model predicted the direction and size of clonal change correctly.

### A4. Quick checks we will run on arrival

- How many patients have ≥ 4 ANC and PLT values between day +45 and +90. This decides how the late-plateau test (P1) is run.
- Censoring: we stop using counts after a CD34+ boost, relapse, next therapy or death. We stop using CRP after tocilizumab.

---

## Part B: phenotype definitions for clinician sign-off

The primary endpoint is a three-class ANC recovery phenotype. The thresholds below are **placeholders**.

- Each one can be changed in one place, and the full analysis reruns the same day.
- They must be **locked before the hold-out set is evaluated (Thursday)**.

### Current rule

ANC is in ×10⁹/L, using values from day 0 onward. The rules are applied in this order:

| Order | Class | Current definition |
|---|---|---|
| 1 | Not classifiable | Follow-up ends before day +45, or fewer than 3 ANC values in days +45 to +75 |
| 2 | **Persistent** | Median ANC over days +45 to +75 **< 0.5** |
| 3 | **Second wave** | Early recovery: ≥ 3 consecutive values **≥ 1.0**, with the run starting in days 0–28.<br>Then a dip: ≥ 2 consecutive values **< 1.0**, starting after that recovery and within days +21 to +60 |
| 4 | **Quick** | Everything else |

How this compares with Rejeski et al. (Blood 2021):

- Their "intermittent" recovery is ANC recovering above **1.5** and then falling below 1.0 after day +21. Our recovery threshold (1.0) is lower; 1.5 would match their definition exactly.
- Their "aplastic" phenotype is ANC < 0.5 continuously for ≥ 14 days. That is an early-duration criterion. Our "persistent" class is a late plateau (days +45 to +75), which is the pattern the clonal hypothesis is about.

Please confirm, or tell us if you prefer the EHA/EBMT ICAHT grading.

### Thresholds to choose

| Item | Current | Alternatives to consider |
|---|---|---|
| Persistent: ANC threshold | < 0.5 | < 1.0 (broader) |
| Persistent: window | days +45 to +75 | days +30 to +90; or a single landmark at day +60 ± 7 |
| Recovery threshold | ≥ 1.0 | ≥ 1.5 (matches Rejeski); ≥ 0.5 |
| Dip threshold | < 1.0 | < 0.5 (only clinically significant dips) |
| Run lengths (recovery / dip) | 3 / 2 values | Depends on how often CBCs are drawn after day +14 |

### Handling questions

1. **G-CSF.** Does recovery *under* G-CSF count as recovery? Currently yes, and the model includes G-CSF as an input.
2. **CD34+ boost.** Counts after a boost are dropped, so a boost before about day +60 can make a patient unclassifiable. Should a boost given for poor graft function count as persistent instead?
3. **Early relapse or death** (before day +45). Currently not classifiable. Agree?
4. **Exclusions.** Marrow lymphoma involvement and IEC-HS/HLH are currently excluded. Agree, and are there others (e.g., CMV/parvovirus B19 marrow suppression)?
5. **Platelets.** Should PLT enter the definition, e.g. persistent = ANC *or* PLT low? Currently the definition uses ANC only, and PLT is a secondary check.
