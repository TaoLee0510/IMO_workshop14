# =============================================================================
# an_PR_prereg.R
# Purpose : Step PR - pre-registration, written before any real-data fit and
#           archived with SHA-256 (run_all.R refuses F1 on real data unless
#           the hashes verify). Fixes the tests T1-T7, the phenotype rule, the
#           model-selection rules, the K3 reporting routes and which random
#           effects are interpreted per patient.
# Inputs  : ctx (config, results K2, K3).
# Outputs : list(prereg = list(files), rules, tables)
# =============================================================================

an_PR_prereg <- function(ctx) {
  cfg <- ctx$config
  k3 <- ctx$results$K3
  k2 <- ctx$results$K2
  routes <- k3$routes
  # random effects interpreted per patient: K3 recovery (truth correlation on the real design)
  att <- k3$summary$attribution
  rec_r_s <- if (!is.null(att)) suppressWarnings(stats::cor(att[design == "real35"]$eta_r_s,
                                                            log(att[design == "real35"]$true_r_s))) else NA
  code_files <- c("src/aim4_model.cpp", "config.R", "R/model.R", "R/likelihood.R", "R/fit_its.R",
                  "R/fit_mcmc.R", "R/attribution.R", "R/phenotype.R")
  rules <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    data_source = cfg$data_source,
    model_version = "v1.2: lineage tree of Aim4/image (2).png (HSC -> MPP -> CLP / CMP -> GMP / MEP), seven self-renewal rates, e_sr estimated; Methods: Aim4/Methods/methods.pdf",
    code_sha256 = as.list(stats::setNames(vapply(code_files, function(f)
      digest::digest(file = file.path(cfg$paths$code, f), algo = "sha256"), character(1)), code_files)),
    phenotype_rule = cfg$phenotype$primary,
    phenotype_variants = cfg$phenotype$variants,
    models = list(q2 = cfg$fit$q2_models, clone_types = names(cfg$model$clone_types), reference = "MF"),
    selection = list(
      q2 = "lowest BIC among M0, Ma, Mb, Mab, Mabx (no clone); 'supported' when the runner-up is >= 2 BIC units worse",
      q1 = "selected Q2 model with the clone layer (DDR- or DTA-type, fixed traits; the better of the two, counted as one extra parameter) vs without; clone layer supported when BIC decreases by >= 2",
      reference = "every Aim 4 model is compared with MF (Friberg 2002) by BIC and by landmark AUROC",
      final_model = "selected Q2 model + DDR-type clone layer: used for attribution (both readouts), forecasts and GreenLight"),
    routes = if (is.null(routes)) NULL else lapply(seq_len(nrow(routes)), function(i) as.list(routes[i])),
    route_rules = cfg$steps$K3$route,
    individual_interpretation = list(
      r_s_truth_correlation_K3 = rec_r_s,
      rule = "a random effect whose K3 truth correlation on the real design is < 0.3 is reported at cohort level only; per-patient readouts are counterfactual probabilities from the MCMC posterior, never point estimates"),
    tests = list(
      T1 = "second-wave nadir day (minimum ANC after the first recovery, days 0-30) regressed on the patient's ferritin-driver peak day; slope and 95% CI; a-dominant ~ 0, b-dominant ~ 1; descriptive (n ~ 10, most had G-CSF)",
      T2 = "Q2 model comparison (selection rule above), also vs MF; primary",
      T3 = "Q1 model comparison (clone layer vs none); secondary, reported along the K3 route for Q1",
      T4 = "leave-one-out landmark prediction of non-quick recovery (second wave or delayed) at days 0, 7, 10, 14, 21: AUROC of the final model vs CAR-HEMATOTOX (at LD) and vs a logistic model on log ANC, log PLT and their 3-day slopes at L; primary at L = 10 and 14; all 35 and the 23 without G-CSF",
      T5 = "dose-response: individual th_aS higher with bridging chemotherapy (Wilcoxon); th_bS increasing with ferritin-driver AUC (days 0-30) and with maximum CRS grade (Spearman); secondary",
      T6 = "posterior-mean f0 from the clone model fitted WITHOUT covariates on f0, vs age, prior lines, prior ASCT (Spearman / Wilcoxon); descriptive",
      T7 = "multi-lineage hit time: per patient, mean over ANC, PLT, AMC of (nadir day in days -5..14 minus the lineage MTT); distribution vs day -4 (LD) and vs the ferritin-driver peak; secondary"),
    sensitivity = c("six phenotype rule variants", "no-G-CSF subgroup (23)", "bridging term removed",
                    "beta_renal 0.3 and 1.2", "DTA-type clone instead of DDR-type",
                    "each fixed self-renewal rate x 0.5 and x 1.5 (capped at 0.95 k)",
                    "self-renewal feedback on the MPP only", "e_sr = 0 (inflammation does not act on self-renewal)"),
    landmarks = cfg$steps$F3$landmarks,
    looked_at_before_PR = "daily ANC trajectories (phenotype definition, Aim 3 data review); ferritin curves only as model input (driver shape), never against outcomes; ALC relative to baseline on days -2..0 (median per patient 0.042) to calibrate th_aL (ALC is not in the likelihood)",
    not_looked_at = "ferritin vs nadir timing, any dose-response, any model fit to real counts")
  md <- c("# Aim 4 pre-registration",
          "",
          paste0("Created ", rules$created, " (data source: ", cfg$data_source, "). Archived with SHA-256 in `PR.sha256`; ",
                 "`run_all.R` refuses to fit real counts unless the hashes verify."),
          "",
          paste0("Model ", rules$model_version, "."),
          "",
          "Code SHA-256 at archiving:",
          paste0("- `", names(rules$code_sha256), "`: ", unlist(rules$code_sha256)),
          "",
          "## Phenotype rule (primary)",
          "Day 0-30 ANC; censored values < 0.5; isolated G-CSF spikes are not a recovery; nadir = first minimum on days 0-14;",
          "recovery = first ANC >= 1.5 after the nadir; second wave = any ANC < 1.0 after recovery; delayed = no recovery by day 30;",
          "quick = the rest. Six variants as sensitivity (window 35, Rejeski day-21, two consecutive values, recovery 1.0, dip 0.5).",
          "",
          "## Models and selection",
          paste0("- Q2: ", rules$selection$q2),
          paste0("- Q1: ", rules$selection$q1),
          paste0("- Reference: ", rules$selection$reference),
          paste0("- Final model: ", rules$selection$final_model),
          "",
          "## Tests",
          paste0("- **", names(rules$tests), "**: ", unlist(rules$tests)),
          "",
          "## Reporting routes from K3 (real design)",
          if (is.null(routes)) "K3 not available." else
            paste0("- ", routes$question, ": cohort-level correct ", round(routes$cohort_correct, 2),
                   ", individual mean P(true class) ", round(routes$indiv_mean_p_true, 2), " -> route **", routes$route, "**"),
          "",
          "Route 1: model comparison + individual attribution + prediction. Route 2: model comparison + cohort shares; individual attribution as examples only. Route 3: description and prediction only; the question goes to the measurement request.",
          "",
          "## Individual interpretation",
          rules$individual_interpretation$rule,
          paste0("K3 truth correlation of r_S on the real design: ", round(rec_r_s, 2), "."),
          "",
          "## Sensitivity analyses",
          paste0("- ", rules$sensitivity),
          "",
          "## Looked at before this pre-registration",
          rules$looked_at_before_PR,
          "",
          "## Not looked at",
          rules$not_looked_at)
  ctx$log("  PR: pre-registration prepared (routes: ",
          if (is.null(routes)) "none" else paste(routes$route, collapse = ", "), ")")
  list(prereg = list(files = list(PR_plan.md = md, PR_rules.json = rules)), rules = rules,
       tables = list(routes = routes))
}
