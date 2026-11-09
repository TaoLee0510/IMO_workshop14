# Specific Aims page: Aim 2 draft (v6, real data)

Draft of 2026-11-03 afternoon, for the team's Aims page (due Thu 11/05, 16:20). Revised for the data that exist:

- 35 axi-cel patients;
- daily counts and CRP/ferritin to about day 30;
- no mutation data;
- tocilizumab, steroids and late G-CSF recorded only as given / not given.

Numbers come from `Results/real/summary.md` (real data) and `Results/real/` D3d (simulation at this scale). Section 3 lists what the text deliberately does **not** claim.

---

## 1. Aim 2 text

### Full version (~200 words)

**Aim 2. Determine why neutrophil recovery after CAR T is delayed or interrupted, using a mechanistic marrow digital twin that separates low marrow reserve from inflammatory marrow suppression.**

In our discovery cohort of 35 axi-cel patients, 43% did not recover quickly:
- 29% had an intermittent course (recovery, then a second fall below 1.0 × 10⁹/L after day 21);
- 14% had not recovered by day 30.

Low marrow reserve, CAR T-driven inflammation and pre-existing myeloid clones have all been proposed. They have not been separated within the same patients.

We fit one mechanistic model of stem, progenitor, neutrophil and platelet dynamics jointly to each patient's daily blood counts and CRP/ferritin over days 0–30. Each patient has two parameters, marrow reserve and inflammatory gain. Lymphodepletion and maturation kinetics are estimated from the cohort.

**2a.** Test whether inflammatory marrow suppression explains the counts beyond reserve alone, and which patients' courses it explains.

**2b.** Predict the day-30 recovery phenotype from day 0–10 data, by cross-validation, against a CAR-HEMATOTOX-like score and an elastic-net baseline.

**2c.** Use the model's clonal layer in silico to specify and size the cohort that could test clonal takeover: CH sequencing before lymphodepletion in every patient, and counts to day 90.

*Expected outcome:* per-patient estimates of reserve and inflammatory gain with early recovery forecasts, and a sized, falsifiable design for the clonal hypothesis.

### Short version (~110 words, if space is tight)

**Aim 2. Separate low marrow reserve from inflammatory marrow suppression in delayed neutrophil recovery after CAR T, with a mechanistic marrow digital twin.**

One model with two patient-level parameters (marrow reserve, inflammatory gain) is fitted to 35 patients' daily counts and CRP/ferritin over days 0–30.

- **2a:** test whether inflammation explains the intermittent courses beyond reserve.
- **2b:** predict the day-30 phenotype (quick / intermittent / no recovery) from day 0–10 data, by cross-validation.
- **2c:** design and size the cohort (pre-lymphodepletion CH sequencing, counts to day 90) that can test clonal takeover.

*Outcome:* early, mechanism-based recovery forecasts and a sized test of the clonal hypothesis.

---

## 2. Preliminary results that support the text

Real data, 35 patients, fit to day 30:

- **Inflammation layer.** Adding inflammatory marrow suppression (M1) to the reserve-only model (M0) improved the fit by ΔBIC = 54.
  - Most of the gain came from intermittent patients: 21.8 of 30.3 log-likelihood units, from 10 of 35 patients.
- **Timing (not testable here).** The intermittent nadir followed the ferritin peak by a median of 18.5 days (IQR 16–22, n = 10), against 15 days implied by the model's literature lag.
  - This is not evidence: the dip must fall in days 21–30 and the counts end at day 30, so the nadirs all sit on days 25–29 and the lag is largely set by the window.
  - Peak and nadir days were uncorrelated (Spearman 0.17).
  - A real test needs counts to day 60.
- **Platelets.** Platelets were also low on days 21–30 in 9/10 intermittent and 3/5 non-recovering patients, against 7/20 quick ones.
- **Forward prediction.** Predicted from the day 0–30 fit, ANC after day 30 matched the measured values on average (median log ratio −0.04 for days 31–44, −0.18 after day 75), but only weakly patient by patient (Spearman 0.18–0.64).

- **Is 35 patients enough? (D3d, simulation at this size and sampling, 50 cohorts per arm).**
  - When the fitted inflammation effect is real, ΔBIC > 10 in 100% of cohorts (median 140).
  - When it is absent, ΔBIC > 10 in 2% (median −10.5).
  - Approximating the treatment days costs almost nothing.
  - The real ΔBIC of 54 lies well inside the first range.
- **Early prediction (2b, cross-validated).** At day 10, the twin ranked non-recovery by day 30 best:
  - AUC 0.89 (95% CI 0.71–0.99), against 0.59 for a CAR-HEMATOTOX-like score and 0.63 for an elastic net;
  - but it over-called non-recovery and was poorly calibrated (Brier 0.94 vs 0.62 and 0.53);
  - it could not predict intermittent courses (AUC 0.50; elastic net 0.77).

Suggested sentence for 2a/2b:

> "In 35 patients, an inflammatory marrow-suppression layer improved the fit beyond reserve alone (ΔBIC 54), mostly in intermittent courses; in simulated cohorts of the same size and sampling this criterion separated the two mechanisms. From day 0–10 data the twin ranked non-recovery by day 30 with a cross-validated AUC of 0.89."

---

## 3. What the text deliberately does not claim (for the team)

**Anything about clones from these data.**
- There are no mutation data, and no counts after day 30 that could show a low plateau. The clonal layer is not fitted to real data.
- 2c is a design aim. Its basis is simulation:
  - with universal pre-treatment sequencing, the evidence about the clone type roughly doubles and clone-fraction recovery rises from r ≈ 0.73 to 0.93 (D3b);
  - at ~350 patients, the clone type was identified in 10/10 simulated cohorts (D3c).

**Validated prediction.**
- 35 patients allow cross-validation only. There is no locked hold-out, and all confidence intervals are wide.
- The text says "predict … by cross-validation", not "validate".

**Mechanism of the intermittent dip as settled.**
- The inflammation layer improves the fit of intermittent patients, but the fitted model does not generate the recovery-then-dip shape itself.
  - Simulated from each patient's own estimates, 0/10 intermittent courses are reproduced.
  - The fitted population gives about 2% intermittent courses, against 29% observed.
- Early G-CSF is not recorded. Some patients' ANC overshoots to 8–9 around day 12–15 and then falls, which looks like a G-CSF effect.
- Treatment days are approximated. Moving the CRS treatments 2 days later changed the reserve and gain estimates by almost nothing (correlation > 0.99).

**The intermittent label as robust.**
- The rule counts a single ANC < 1.0 after day 21. In D3d it labelled about 27% of patients intermittent even when inflammation had no marrow effect.
- That is an upper bound (the simulated noise includes model misfit), but a two-value version should be agreed with the clinicians.

**Phenotype labels.**
- The rule (Rejeski, recovery counted from the nadir) agrees with the dataset's own `phenotype_from_anc` in 32/35 patients.
- The original cohort labels (5 "aplastic") disagree. By Rejeski's definition (ANC < 0.5 for ≥ 14 days) no patient qualifies; the longest run is 13 days.
- **To confirm with the clinicians.**

---

## 4. Hooks for integration with the other aims (to agree with their owners)

**In:** CAR T expansion.
- Aim 2 drives inflammation with one population-level CAR T kernel; only its amplitude is patient-specific.
- Patient-level expansion curves from another aim would replace the kernel directly. This is the same link as `tumor_immune → marrow` in the GreenLight demo.
- Patient-specific CRS timing is also the most likely fix for the loose CRP fit.

**Out:** for each patient, Aim 2 provides:
- marrow reserve and inflammatory gain;
- predicted ANC and PLT paths;
- phenotype probabilities under alternative supportive-care plans (G-CSF timing).

These feed the decision layer.

---

## 5. Slides

### Thursday (1–2 slides with the Aims page)

**Slide 1, "Aim 2 at a glance"**

| Left | Centre | Right |
|---|---|---|
| ANC by phenotype, real data (`S0_anc_by_phenotype.png`) | Model schematic, two patient dials (reserve, inflammatory gain) + clonal layer greyed out | 2a / 2b / 2c with their test and data |

**Slide 2 (optional), "What 35 patients can say, and what we ask for"**
- Inflammation layer gain by phenotype (`D4_P5_gain_by_phenotype.png`) and nadir timing (`D4_P2_timing.png`).
- D3d at this scale (`D3d_model_choice.png`).
- The request:
  - pre-lymphodepletion CH sequencing for all;
  - counts to day 90;
  - treatment dates.

### Friday (final, 4–5 slides)

1. Problem and phenotypes, real data: `S0_anc_by_phenotype.png`, `S0_phenotype_bar.png`.
2. One model; what inflammation adds and where: `D4_P5_gain_by_phenotype.png`, `D4_example_fits.png`.
3. Timing and forward check: `D4_P2_timing.png`, `D4_late_check.png`.
4. Early prediction (cross-validated) against the baselines: `H_metrics.png`.
5. The clone as a design aim: `D3c_clone_type_evidence.png` (TestData design), `D6_power.png`.

Real-data figures are in `Results/real/figures/`; design figures in `Results/test/figures/`.
