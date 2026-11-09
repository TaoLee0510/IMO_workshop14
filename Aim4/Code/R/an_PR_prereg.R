# =============================================================================
# an_PR_prereg.R
# Purpose : Step PR - pre-registration, written before any real-data fit and
#           archived with SHA-256 (run_all.R refuses F1 on real data unless
#           the hashes verify). v1.3: model version, code and data hashes, data
#           handling (G-CSF windows, individual baselines, ALC), DE-ITS rules,
#           model set and selection, readouts, tests T1-T8.
# Inputs  : ctx (config).
# Outputs : list(prereg = list(files), rules, tables)
# =============================================================================

an_PR_prereg <- function(ctx) {
  cfg <- ctx$config
  code_files <- c("src/aim4_model.cpp", "config.R", "R/model.R", "R/likelihood.R", "R/fit_de.R", "R/fit_its.R",
                  "R/fit_mcmc.R", "R/attribution.R", "R/phenotype.R", "R/io_convert.R", "R/io_read.R")
  data_files <- c("patients.csv", "obs.csv", "inflam.csv", "gcsf.csv")
  w <- cfg$fit$wave1
  rules <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    data_source = cfg$data_source,
    model_version = paste("v1.3: v1.2 lineage tree + CAR T - target interaction module (CAR T cells, CD19 antigen,",
                          "two-stage signal Z acting on HSC and MPP self-renewal), G-CSF as input (imputed windows),",
                          "individual baselines, DE-ITS estimation; Methods: Aim4/Methods/methods.pdf"),
    previous_version = "v1.2 pre-registration and results are kept unchanged in Results/real (code frozen in Code_v12_frozen)",
    code_sha256 = as.list(stats::setNames(vapply(code_files, function(f)
      digest::digest(file = file.path(cfg$paths$code, f), algo = "sha256"), character(1)), code_files)),
    data_sha256 = as.list(stats::setNames(vapply(data_files, function(f)
      digest::digest(file = file.path(cfg$paths$derived, f), algo = "sha256"), character(1)), data_files)),
    phenotype_rule = cfg$phenotype$primary,
    phenotype_variants = cfg$phenotype$variants,
    gcsf_rule = cfg$data$gcsf_rule,
    baseline = cfg$model$baseline,
    estimation = list(method = "DE-ITS", de = cfg$fit$de, outer = cfg$fit$outer, is = cfg$fit$is,
                      indiv_limits = cfg$model$indiv_limits, struct_prior_sd = cfg$model$struct_prior_sd,
                      struct_prior_sd_k = cfg$model$struct_prior_sd_k, struct_limits = cfg$model$struct_limits,
                      e_step_cores = 22L, m_step_cores = 60L),
    models = list(wave1 = w, wave2 = "for the selected Q2 model q: qC_ddr, qC_dta, q_Zh, q_Zm, qC_ddr_Zh, qC_ddr_Zm",
                  reference = "MF"),
    selection = list(
      q2 = "lowest BIC among M0, Ma, Mb, Mab, Mabx; 'supported' when the runner-up is >= 2 BIC units worse",
      q1 = "selected Q2 model with the clone layer (DDR- or DTA-type, fixed traits; the better of the two, counted as one extra parameter) vs without; clone layer supported when BIC decreases by >= 2",
      route = "route of the CAR T - target signal Z: both / HSC only (_Zh) / MPP only (_Zm) on the Q1-selected model, lowest BIC; a route is 'supported' when dropping it raises BIC by >= 2",
      reference = "MF (Friberg 2002) is compared on the ANC/PLT/AMC part of the marginal likelihood (MF has no ALC)",
      bic = "BIC = -2 sum log L_i + n_par log(35); log L_i by importance sampling (ESS >= 30), else the Laplace value (counted and reported)",
      convergence = "a fit that hits the outer safety cap is reported as not converged and is not used for selection unless no converged alternative exists (stated explicitly)",
      final_model = "Q1/route-selected model: used for attribution (when, route of b, who), forecasts and GreenLight"),
    readouts = list(
      when = "HSC level: a / b / a+b joint / a+b redundant / other, b = ferritin route (th_bS) + CAR T - target route (th_zS)",
      b_route = "among b-driven patients: ferritin route (no_bI improves) vs CAR T - target route (no_bZ improves) vs both",
      who = "clone / reserve / both / other (non-quick)"),
    tests = list(
      T1 = "second-wave nadir day vs ferritin-driver peak day (descriptive)",
      T2 = "Q2 model comparison (primary)",
      T3 = "Q1 model comparison (secondary)",
      T4 = "landmark prediction of non-quick recovery at days 0, 7, 10, 14, 21, APPROXIMATE: population of the full-data fit, individual posterior from the data up to the landmark, no G-CSF after the landmark; AUROC vs CAR-HEMATOTOX and a counts-only logistic model",
      T5 = "dose-response of th_aS (bridging) and th_bS (ferritin AUC, CRS grade); secondary",
      T6 = "f0 vs risk factors from the clone model fitted without covariates on f0; descriptive",
      T7 = "multi-lineage hit time; secondary",
      T8 = "CAR T compartment vs peak_expansion (16 patients with a value): Spearman correlation of the model peak with the measured peak; external check"),
    recovery = "10 synthetic real35 cohorts from the final structure with structural values off the prior centres, fitted with DE-ITS (replaces K3)",
    sensitivity_planned = c("v1.2 handling of G-CSF rows (upper bounds) instead of the imputed windows",
                            "six phenotype rule variants", "no-G-CSF subgroup"),
    landmarks = cfg$steps$F3$landmarks,
    looked_at_before_PR = paste("all v1.2 real-data results (F1-F4, 2026-11-03), which motivated v1.3 (poor fit of second",
                                "waves, extreme individual estimates, G-CSF upper bounds); the median observed ALC / baseline",
                                "by day window (CAR T peak calibration); the availability of G-CSF, CAR T, cytokine and",
                                "LDH fields in the extract; the imputed G-CSF windows"),
    not_looked_at = "any v1.3 fit to real counts")
  md <- c("# Aim 4 pre-registration - model v1.3",
          "",
          paste0("Created ", rules$created, " (data source: ", cfg$data_source, "). Archived with SHA-256 in `PR.sha256`; ",
                 "the cluster refuses to start a real-data fit unless the hashes verify."),
          "",
          paste0("Model ", rules$model_version, "."),
          rules$previous_version,
          "",
          "Code SHA-256 at archiving:",
          paste0("- `", names(rules$code_sha256), "`: ", unlist(rules$code_sha256)),
          "",
          "Data SHA-256 (Data/derived_v13):",
          paste0("- `", names(rules$data_sha256), "`: ", unlist(rules$data_sha256)),
          "",
          "## Phenotype rule (primary, unchanged)",
          "Day 0-30 ANC; censored values < 0.5; isolated G-CSF spikes are not a recovery; nadir = first minimum on days 0-14;",
          "recovery = first ANC >= 1.5 after the nadir; second wave = any ANC < 1.0 after recovery; delayed = no recovery by day 30.",
          "",
          "## Data handling new in v1.3",
          "- G-CSF windows imputed: one dose the day before each isolated ANC spike; flagged patients: 5-day course from the first day >= 15 with ANC < 1.0, repeated (to day 100). G-CSF rows are ordinary observations.",
          "- Baselines are individual parameters: log b = log(set-point median) + d, d ~ N(0, s^2/n + tau^2); the set-point rows enter only through this prior.",
          "- ALC is in the likelihood as lymphocytes + CAR T cells.",
          "",
          "## Estimation (DE-ITS)",
          "Per patient: DEoptim (strategy 2, NP 220, F 0.8, CR 0.9) in the physiological box, chunks of 50 generations, converged when the best value improved by < 1e-6 (|f|+1) over 100 generations and the population objective spread < 1e-3; L-BFGS-B polish; Laplace covariance. M step: closed forms for mu, beta, omega, tau, sigma; structural parameters by DEoptim (spread tolerance 1e-3 per patient). Outer loop converged when the population vector changes < 0.01 and sum(nll) < 1e-4 relative on 2 consecutive rounds; oscillating components damped; safety cap 50 rounds.",
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
          "## Readouts",
          paste0("- ", names(rules$readouts), ": ", unlist(rules$readouts)),
          "",
          "## Tests",
          paste0("- **", names(rules$tests), "**: ", unlist(rules$tests)),
          "",
          "## Identifiability",
          rules$recovery,
          "",
          "## Planned sensitivity analyses",
          paste0("- ", rules$sensitivity_planned),
          "",
          "## Looked at before this pre-registration",
          rules$looked_at_before_PR,
          "",
          "## Not looked at",
          rules$not_looked_at)
  ctx$log("  PR: v1.3 pre-registration prepared")
  list(prereg = list(files = list(PR_plan.md = md, PR_rules.json = rules)), rules = rules, tables = list())
}
