# =============================================================================
# an_PR_prereg.R
# Purpose : Step PR - pre-registration of Model1, written before any real-data
#           fit and archived with SHA-256 (the cluster refuses a real-data fit
#           unless the hashes verify): model, code and data hashes, the
#           phenotype term, round 0, DE-ITS rules, model set and selection, the
#           V1 folds, the synthetic validation S1-S5, readouts and tests.
# Inputs  : ctx (config, data).
# Outputs : list(prereg = list(files), rules, folds, tables)
# =============================================================================

# V1 folds: stratified by the primary phenotype; within each class the patients are shuffled
# (fixed seed) and dealt round-robin, continuing across classes, so the folds have 7 patients
# and the class mix of the cohort.
v1_folds <- function(phenotype, config) {
  v <- config$validation$V1
  set.seed(v$seed)
  ph <- data.table::copy(phenotype[, .(patient_id, primary)])
  ph[, primary := factor(primary, levels = config$phenotype$classes)]
  data.table::setorder(ph, primary)
  ph[, u := stats::runif(.N)]
  data.table::setorder(ph, primary, u)
  ph[, fold := ((seq_len(.N) - 1L) %% v$k) + 1L]
  ph[order(patient_id), .(patient_id, phenotype = as.character(primary), fold)]
}

an_PR_prereg <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  code_files <- c("src/model1.cpp", "config.R", "R/model.R", "R/likelihood.R", "R/phenotype.R", "R/fit_de.R",
                  "R/fit_mcmc13.R", "R/attribution13.R", "R/io_convert.R", "R/io_read.R", "R/sim_cohort.R",
                  "R/aim4_ref.R", "tests/ref/aim4_v13.cpp")
  data_files <- c("patients.csv", "obs.csv", "inflam.csv", "gcsf.csv")
  folds <- v1_folds(d$phenotype, cfg)
  w <- cfg$fit$wave1
  ph <- cfg$fit$pheno
  rules <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    data_source = cfg$data_source,
    model_version = paste("Model1 v1.0: Aim4 v1.3 with a Friberg-like haematopoietic layer (stem cells with a",
                          "pre-existing clone -> four MPP components N, Mo, Pl, L with shared k_P, sr_P -> Friberg",
                          "maturation chains; one count feedback per myeloid lineage on MPP self-renewal; marrow",
                          "inflammation and the CAR T - target signal suppress MPP self-renewal; lymphocytes",
                          "without chain or feedback); recovery phenotype in the likelihood; round 0 + DE-ITS;",
                          "plan Model1/Model1_plan.md v0.3; Methods Model1/Methods/methods.pdf"),
    deviations_from_plan = c(
      "lymphoid MPP component without count feedback (K1: a feedback on P_L made ALC overshoot to ~1.3x baseline by day 30)",
      "upper limit 2.2 of the gamma search box (K1: Hopf bifurcation of the steady state at gamma = 2.37)"),
    code_sha256 = as.list(stats::setNames(vapply(code_files, function(f)
      digest::digest(file = file.path(cfg$paths$code, f), algo = "sha256"), character(1)), code_files)),
    data_sha256 = as.list(stats::setNames(vapply(data_files, function(f)
      digest::digest(file = file.path(cfg$paths$derived, f), algo = "sha256"), character(1)), data_files)),
    phenotype_rule = cfg$phenotype$primary,
    phenotype_variants = cfg$phenotype$variants,
    phenotype_term = list(w = ph$w, s = ph$s, eps = ph$eps, sensitivity_s = ph$sens_s, sensitivity_w = ph$sens_w,
                          form = "margin-based soft classification of the model ANC on the label days; composite likelihood; not used at landmarks"),
    gcsf_rule = cfg$data$gcsf_rule,
    baseline = cfg$model$baseline,
    estimation = list(method = "round 0 (independent fits, prior SD x sd_mult) + DE-ITS", round0 = cfg$fit$round0,
                      de = cfg$fit$de, outer = cfg$fit$outer, is = cfg$fit$is,
                      indiv_limits = cfg$model$indiv_limits, struct_prior_sd = cfg$model$struct_prior_sd,
                      struct_prior_sd_k = cfg$model$struct_prior_sd_k, struct_limits = cfg$model$struct_limits,
                      e_step_cores = 22L, m_step_cores = 60L),
    models = list(wave1 = w, wave2 = "for the selected Q2 model q: qC_ddr, qC_dta, q_Zh, q_Zm, qC_ddr_Zh, qC_ddr_Zm",
                  reference = "MF (with the same phenotype term)"),
    selection = list(
      q2 = "lowest BIC among M0, Ma, Mb, Mab, Mabx; 'supported' when the runner-up is >= 2 BIC units worse",
      q1 = "selected Q2 model with the clone layer (DDR- or DTA-type, fixed traits; the better of the two, counted as one extra parameter) vs without; clone layer supported when BIC decreases by >= 2",
      route = "route of the CAR T - target signal Z: both / HSC only (_Zh) / MPP only (_Zm) on the Q1-selected model, lowest BIC; a route is 'supported' when dropping it raises BIC by >= 2",
      reference = "MF (Friberg 2002) is compared on the ANC/PLT/AMC part of the marginal likelihood (MF has no ALC); both include the phenotype term",
      bic = "BIC = -2 sum log L_i + n_par log(35) of the composite likelihood (counts + w x phenotype term); log L_i by importance sampling (ESS >= 30), else the Laplace value (counted and reported)",
      convergence = "a fit that hits the outer safety cap is reported as not converged and is not used for selection unless no converged alternative exists (stated explicitly)",
      final_model = "Q1/route-selected model: used for attribution (when, route of b, who), forecasts, V1 and GreenLight"),
    V1 = list(design = paste("5 folds stratified by the primary phenotype (assignment below, fixed now); the final model is",
                             "refitted on the 28 patients of each training set (warm start from the full-data fit, no round 0);",
                             "each held-out patient is predicted at the landmarks with its fold's population, without the",
                             "phenotype term; phenotype calibration, Brier score and AUROC are reported only from V1"),
              folds = folds),
    synthetic = list(
      S1 = "10 cohorts from Model1 (MabC_ddr, structural truth = prior centres x truth_mult) on the real designs: recovery of the structural and individual parameters, of the HSC loss f_a, f_bI, f_bZ and of the attribution",
      S2 = "the S1 cohorts refitted with w = 0: bias of the phenotype term and coverage of the 90% posterior intervals of the individual parameters",
      S3 = "V1 on 2 S1 cohorts: prediction ceiling when the model is correct",
      S4 = "5 cohorts from the Aim4 v1.3 lineage-tree engine (MabC_ddr, Aim4 prior centres) fitted with Model1: does the simplification keep the a/b and clone conclusions",
      S5 = "visual predictive check: 500 replicate cohorts from the fitted population on the real designs vs the observed 19/10/6, nadir day, recovery day and second-wave day",
      seeds = "S1 cohorts seeds 101-110, S4 cohorts seeds 401-405"),
    readouts = list(
      when = "HSC level: a / b / a+b joint / a+b redundant / other, b = ferritin route (th_bS) + CAR T - target route (th_zS)",
      b_route = "among b-driven patients: ferritin route (no_bI improves) vs CAR T - target route (no_bZ improves) vs both",
      who = "clone / reserve / both / other (non-quick)"),
    known_properties = paste("K1: in Model1 a chemotherapy hit alone gives quick or delayed recovery, not a second wave",
                             "(no pseudo second wave as in Aim4); a second wave needs the CAR T phase (th_bS / th_zS with",
                             "self-renewal suppression); S4 quantifies the consequence for the a/b readout"),
    tests = list(
      T1 = "second-wave nadir day vs ferritin-driver peak day (descriptive)",
      T2 = "Q2 model comparison (primary)",
      T3 = "Q1 model comparison (secondary)",
      T4 = "landmark prediction of non-quick recovery at days 0, 7, 10, 14, 21, out of sample through the V1 folds; AUROC vs CAR-HEMATOTOX and a counts-only logistic model (same folds)",
      T5 = "dose-response of th_aS (bridging) and th_bS (ferritin AUC, CRS grade); secondary",
      T6 = "f0 vs risk factors from the clone model fitted without covariates on f0; descriptive",
      T7 = "multi-lineage hit time; secondary",
      T8 = "CAR T compartment vs peak_expansion (16 patients with a value): Spearman correlation of the model peak with the measured peak; external check"),
    sensitivity_planned = c("phenotype term: s = 0.15 and 0.6; w = 0 and 5 (final model)",
                            "six phenotype rule variants (label and model probability with the same variant)",
                            "Aim4 v1.2 handling of G-CSF rows (upper bounds) instead of the imputed windows",
                            "no-G-CSF subgroup"),
    landmarks = cfg$steps$F3$landmarks,
    looked_at_before_PR = paste("all Aim4 real-data results (v1.2 complete; v1.3 fits in progress: model comparison",
                                "and fit quality of the first rounds), which motivated Model1; Model1 K1 (prior-centre",
                                "simulations, calibrations to Aim4 v1.3 and to the median observed ALC at",
                                "lymphodepletion); the phenotype mix of two synthetic test cohorts (Model1 27/5/3,",
                                "Aim4 19/11/5)"),
    not_looked_at = "any Model1 fit to real counts")
  fl <- folds[, .(patients = paste(patient_id, collapse = ", "), quick = sum(phenotype == "quick"),
                  second_wave = sum(phenotype == "second_wave"), delayed = sum(phenotype == "delayed")), by = fold]
  md <- c("# Model1 pre-registration - v1.0",
          "",
          paste0("Created ", rules$created, " (data source: ", cfg$data_source, "). Archived with SHA-256 in `PR.sha256`; ",
                 "the cluster refuses to start a real-data fit unless the hashes verify."),
          "",
          paste0("Model ", rules$model_version, "."),
          "",
          "Deviations from the plan found in K1 (before any fit):",
          paste0("- ", rules$deviations_from_plan),
          "",
          "Code SHA-256 at archiving:",
          paste0("- `", names(rules$code_sha256), "`: ", unlist(rules$code_sha256)),
          "",
          "Data SHA-256 (Model1/Data/derived):",
          paste0("- `", names(rules$data_sha256), "`: ", unlist(rules$data_sha256)),
          "",
          "## Phenotype rule (primary, unchanged)",
          "Day 0-30 ANC; censored values < 0.5; isolated G-CSF spikes are not a recovery; nadir = first minimum on days 0-14;",
          "recovery = first ANC >= 1.5 after the nadir; second wave = any ANC < 1.0 after recovery; delayed = no recovery by day 30.",
          "",
          "## Phenotype term in the likelihood",
          sprintf("- w = %g, s = %g (log scale), eps = %g; %s.", ph$w, ph$s, ph$eps, rules$phenotype_term$form),
          "- The day 0-30 ANC enter twice (counts and phenotype): the BIC is that of a composite likelihood and is compared only between models with the same term.",
          "",
          "## Estimation",
          sprintf("Round 0: every patient fitted on its own with the prior SDs x %g (same search box); mu, beta, omega, tau, sigma initialised from these fits by the M-step formulas. Then DE-ITS as Aim4 v1.3 (DEoptim per patient, closed-form M step + structural DEoptim, outer loop with damping, over-relaxation with reset, averaging phase, cap %d rounds). Warm-started refits (V1, wave 2) skip round 0.",
                  cfg$fit$round0$sd_mult, cfg$fit$outer$max_iter),
          "",
          "## Models and selection",
          paste0("- Wave 1: ", paste(w, collapse = ", "), "; wave 2: ", rules$models$wave2),
          paste0("- Q2: ", rules$selection$q2),
          paste0("- Q1: ", rules$selection$q1),
          paste0("- Route of Z: ", rules$selection$route),
          paste0("- Reference: ", rules$selection$reference),
          paste0("- BIC: ", rules$selection$bic),
          paste0("- Convergence: ", rules$selection$convergence),
          paste0("- Final model: ", rules$selection$final_model),
          "",
          "## V1: 5-fold cross-validation (real data)",
          rules$V1$design,
          "",
          "| fold | quick | second wave | delayed | patients |", "|---|---|---|---|---|",
          paste0("| ", fl$fold, " | ", fl$quick, " | ", fl$second_wave, " | ", fl$delayed, " | ", fl$patients, " |"),
          "",
          "## Synthetic validation",
          paste0("- **", names(rules$synthetic), "**: ", unlist(rules$synthetic)),
          "Results of S1-S5 are reported as 'if the model (S4: Aim4) is correct' and kept apart from V1.",
          "",
          "## Known property (K1)",
          rules$known_properties,
          "",
          "## Readouts",
          paste0("- ", names(rules$readouts), ": ", unlist(rules$readouts)),
          "",
          "## Tests",
          paste0("- **", names(rules$tests), "**: ", unlist(rules$tests)),
          "",
          "## Planned sensitivity analyses",
          paste0("- ", rules$sensitivity_planned),
          "",
          "## Looked at before this pre-registration",
          rules$looked_at_before_PR,
          "",
          "## Not looked at",
          rules$not_looked_at)
  ctx$log("  PR: Model1 pre-registration prepared (V1 fold sizes ", paste(table(folds$fold), collapse = "/"), ")")
  list(prereg = list(files = list(PR_plan.md = md, PR_rules.json = rules, V1_folds.csv = folds)),
       rules = rules, folds = folds, tables = list(V1_folds = folds))
}
