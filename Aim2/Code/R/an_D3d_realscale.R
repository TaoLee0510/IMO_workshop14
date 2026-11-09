# =============================================================================
# an_D3d_realscale.R
# Purpose : D3d - what can a cohort of the real-data size and sampling show?
#           Parametric bootstrap: cohorts of config$steps$D3d$n patients are
#           simulated from the real-data M1 fit (population parameters only:
#           estimated kinetics and inflammation effects, means, between-patient
#           and residual SDs; no clones, as in the fit), observed like the real
#           extract (counts to day 30 plus one late draw for some, treatment
#           days hidden and approximated with R/events_approx.R), and fitted
#           with M0 and M1 as the real data. The simulated phenotype mix is
#           reported next to the observed one: it is also a check of the model.
#           Two variants share each cohort:
#             inflam    : truth as fitted (inflammation suppresses the marrow)
#             no_inflam : q = emax = 0 in the truth (CRP still rises, but the
#                         marrow does not feel it)
# Inputs  : ctx$config (real-data mode), ctx$results$D4 (population fit).
# Outputs : list(per_rep, summary, truth, settings, tables, key_numbers)
#
# Metrics per cohort x variant
#   dbic_m1         : BIC(M0) - BIC(M1); > 0 favours the inflammation layer
#   q_hat, emax_hat : M1 estimates (truth: the fitted values, or 0)
#   cor_log_rho/g   : estimated vs true individual reserve / gain (M1)
#   *_exact         : the same with the true treatment days (what the missing
#                     dates cost)
#   gain_wave       : M1 - M0 log-likelihood summed over intermittent patients
#   twin_*          : day +10 prediction of the day-30 phenotype; population
#                     fitted on the same cohort, so optimistic
#
# Parallel runs as in D3c: parts under Results/real/objects/D3d_parts/, one
# Slurm array task each (Code/hpc/d3d_part.R). The truth is written once by
# Code/hpc/d3d_prepare.R to objects/D3d_truth.rds, so the cluster never needs
# the patient data.
# =============================================================================

d3d_labels <- c(inflam = "Inflammation suppresses the marrow (as fitted)",
                no_inflam = "No marrow effect of inflammation")

# Truth for the simulations from the real-data D4 result.
d3d_truth <- function(cfg, d4) {
  f <- d4$fits$M1
  n <- table(factor(d4$inferred$phenotype, phenotype_levels))
  list(struct = f$pop$struct, mu_rho = f$pop$mu[["rho"]], mu_g = f$pop$mu[["g"]],
       sd_rho = f$pop$omega[["rho"]], sd_g = f$pop$omega[["g"]],
       sigma = unlist(f$pop$sigma), mix = stats::setNames(as.numeric(n / sum(n)), phenotype_levels))
}

d3d_truth_path <- function(cfg) file.path(cfg$paths$results, "objects", "D3d_truth.rds")

# Structural truth: fitted kinetics, design clone traits; no_inflam removes
# the marrow effects of inflammation.
d3d_struct <- function(truth, variant) {
  st <- truth$struct
  st[c("r_ld", "r_inf", "delta")] <- truth_spec("main")$struct_override[c("r_ld", "r_inf", "delta")]
  if (variant == "no_inflam") st[c("q", "emax")] <- 0
  st
}

d3d_spec <- function(cfg, truth, variant = "inflam") {
  s <- cfg$steps$D3d
  spec <- truth_spec("main")
  spec$n <- s$n
  spec$n_sequenced <- 0L
  spec$n_post <- 0L
  spec$train_frac <- 1
  spec$all_cd28 <- TRUE
  spec$struct_override <- d3d_struct(truth, variant)
  spec$pop$ch_logit[["int"]] <- -30          # no clones: the fitted model has none
  spec$pop$mu_rho <- truth$mu_rho
  spec$pop$mu_g <- truth$mu_g
  spec$pop$sd_rho <- truth$sd_rho
  spec$pop$sd_g <- truth$sd_g
  spec$noise[c("anc", "plt", "crp", "ferritin")] <- truth$sigma[c("anc", "plt", "crp", "ferritin")]
  spec$sampling <- list(lab_days = s$lab_days, missing = s$lab_missing, late_frac = s$late_frac,
                        late_days = s$late_days, inflam_days = s$inflam_days,
                        inflam_missing = s$inflam_missing, sparse_frac = 0, sparse_after = 30)
  spec$events$p_gcsf <- s$p_gcsf
  spec
}

# Population means of the generating model: the fitted ones (no calibration to
# the observed mix, which the fitted model does not reproduce; see summary).
d3d_calibrate <- function(cfg, truth) {
  list(mu_rho = truth$mu_rho, mu_g = truth$mu_g)
}

# Treatment events as the fit sees them in the real extract: tocilizumab,
# steroids and G-CSF after day 14 known only as given / not given.
d3d_fit_events <- function(raw, traj, cfg) {
  ev <- raw$events
  p <- raw$patients
  crs <- traj[day >= 0, .(crs_max_day = day[which.max(is)]), by = patient_id]
  flags <- p[, .(patient_id,
                 tocilizumab = as.integer(patient_id %in% ev[event == "tocilizumab", patient_id]),
                 steroid = as.integer(patient_id %in% ev[event == "steroid", patient_id]),
                 gcsf_after_d14 = as.integer(patient_id %in% ev[event == "gcsf" & start_day >= 14, patient_id]))]
  flags <- merge(flags, crs, by = "patient_id", all.x = TRUE)
  flags[patient_id %in% p[crs_grade == 0, patient_id], crs_max_day := NA]
  approx <- approx_treatment_events(flags, raw$labs[, .(patient_id, day, anc)], cfg$data$approx_events)
  rbind(approx, ev[!event %in% c("tocilizumab", "steroid", "gcsf")], fill = TRUE)
}

# Area under the ROC curve (Mann-Whitney), no package needed on the cluster.
auc_mw <- function(score, y) {
  n1 <- sum(y)
  n0 <- sum(!y)
  if (!n1 || !n0) return(NA_real_)
  (sum(rank(score)[y]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

d3d_part_key <- function(cfg, truth) {
  code <- cfg$paths$code
  files <- c(Sys.glob(file.path(code, "R", c("model_*.R", "sim_*.R"))),
             file.path(code, "R", c("phenotype.R", "likelihood.R", "fit_its.R", "io_read.R", "events_approx.R",
                                    "an_H_holdout.R", "an_D3d_realscale.R")),
             file.path(code, "src", "marrow_model.cpp"))
  txt <- c(vapply(sort(files, method = "radix"), function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1)),
           deparse(list(cfg$model, cfg$fit, cfg$phenotype, cfg$steps$D3d, cfg$steps$H$landmark_day,
                        cfg$steps$H$n_draws, cfg$data$approx_events, cfg$data$baseline_window,
                        lapply(truth, signif, 6), cfg$runtime$seed, cfg$quick)))
  digest::digest(paste(txt, collapse = "\n"), algo = "md5", serialize = FALSE)
}

d3d_parts_dir <- function(cfg) file.path(cfg$paths$results, "objects", "D3d_parts")
d3d_part_file <- function(cfg, r, variant) file.path(d3d_parts_dir(cfg), sprintf("rep%02d_%s.rds", r, variant))

# One cohort x variant.
d3d_part <- function(cfg, truth, r, variant, cal, log_fun = message) {
  quiet <- function(...) invisible(NULL)
  t0 <- Sys.time()
  spec <- d3d_spec(cfg, truth, variant)
  sim <- suppressMessages(simulate_cohort(cfg, "main", seed = cfg$runtime$seed + 2000L + r, cal = cal, spec = spec))
  raw <- sim$raw
  # marrow involvement and IEC-HS are not recorded in the real extract
  raw$patients[, `:=`(marrow_involvement = NA_integer_, iec_hs = NA_integer_)]
  raw_a <- raw
  raw_a$events <- d3d_fit_events(raw, sim$truth$traj, cfg)
  tr_a <- analysis_set(standardise(raw_a, cfg), "train")
  tr_e <- analysis_set(standardise(raw, cfg), "train")
  pd_a <- build_pdats(tr_a, cfg)
  pd_e <- build_pdats(tr_e, cfg)
  m0 <- fit_population(pd_a, "M0", cfg, log_fun = quiet)
  m1 <- fit_population(pd_a, "M1", cfg, init = m0$pop, log_fun = quiet)
  m0e <- fit_population(pd_e, "M0", cfg, init = m0$pop, log_fun = quiet)
  m1e <- fit_population(pd_e, "M1", cfg, init = m1$pop, log_fun = quiet)

  tp <- sim$truth$patient
  rec <- function(f) {
    x <- merge(f$indiv[, .(patient_id, eta_rho, eta_g)], tp[, .(patient_id, rho, g)], by = "patient_id")
    c(rho = stats::cor(x$eta_rho, log(x$rho)), g = stats::cor(x$eta_g, log(x$g)))
  }
  ph <- tr_a$pheno[!is.na(phenotype)]
  wave <- ph[phenotype == "second_wave", patient_id]
  gain <- merge(m1$indiv[, .(patient_id, l1 = logL)], m0$indiv[, .(patient_id, l0 = logL)], by = "patient_id")

  # day +10 prediction with the cohort's own M1 population
  pds_lm <- build_pdats(tr_a, cfg, landmark = cfg$steps$H$landmark_day)
  fl <- fit_given_population(pds_lm, m1, cfg)
  tw <- data.table::rbindlist(twin_probabilities(fl, pds_lm, m1, tr_a, cfg, cfg$steps$H$n_draws,
                                                 seed = cfg$runtime$seed + 3000L + r))
  tw <- merge(tw, ph[, .(patient_id, phenotype)], by = "patient_id")
  P <- as.matrix(tw[, .(p_quick, p_second_wave, p_persistent)])

  r_a <- rec(m1)
  r_e <- rec(m1e)
  truth_ph <- tp[match(ph$patient_id, patient_id), phenotype_true]
  row <- data.table::data.table(
    rep = r, variant = variant, n = nrow(ph),
    n_quick = sum(ph$phenotype == "quick"), n_wave = sum(ph$phenotype == "second_wave"),
    n_pers = sum(ph$phenotype == "persistent"),
    dbic_m1 = m0$BIC - m1$BIC, dbic_m1_exact = m0e$BIC - m1e$BIC,
    q_hat = m1$pop$struct[["q"]], emax_hat = m1$pop$struct[["emax"]],
    q_hat_exact = m1e$pop$struct[["q"]], emax_hat_exact = m1e$pop$struct[["emax"]],
    cor_log_rho = r_a[["rho"]], cor_log_g = r_a[["g"]],
    cor_log_rho_exact = r_e[["rho"]], cor_log_g_exact = r_e[["g"]],
    gain_wave = sum(gain[patient_id %in% wave, l1 - l0]),
    twin_brier = mean(multiclass_brier(P, tw$phenotype)),
    twin_auc_wave = auc_mw(tw$p_second_wave, tw$phenotype == "second_wave"),
    twin_auc_pers = auc_mw(tw$p_persistent, tw$phenotype == "persistent"),
    phenotype_obs_vs_truth = mean(ph$phenotype == truth_ph, na.rm = TRUE))
  log_fun(sprintf("  D3d: cohort %d, %-9s dBIC %.1f (exact days %.1f), q %.3f emax %.2f, cor rho %.2f g %.2f, twin AUC wave %.2f (%.1f min)",
                  r, variant, row$dbic_m1, row$dbic_m1_exact, row$q_hat, row$emax_hat, row$cor_log_rho,
                  row$cor_log_g, row$twin_auc_wave, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(key = d3d_part_key(cfg, truth), rep = r, variant = variant, row = row)
}

an_D3d_realscale <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$D3d
  truth <- d3d_truth(cfg, ctx$results$D4)
  todo <- data.table::CJ(rep = seq_len(s$n_rep), variant = s$variants, sorted = FALSE)
  files <- d3d_part_file(cfg, todo$rep, todo$variant)
  if (any(file.exists(files))) {
    if (!all(file.exists(files)))
      stop("D3d: ", sum(!file.exists(files)), " of ", length(files), " parts missing in ", d3d_parts_dir(cfg),
           " (finish the array job, or delete the folder to compute here)")
    parts <- lapply(files, readRDS)
    key <- d3d_part_key(cfg, truth)
    if (!all(vapply(parts, `[[`, character(1), "key") == key))
      stop("D3d: parts in ", d3d_parts_dir(cfg), " were computed with other code, settings or truth",
           " (rerun the array job, or delete the folder to compute here)")
    ctx$log("  D3d: combining ", length(parts), " parts from ", d3d_parts_dir(cfg))
  } else {
    cal <- d3d_calibrate(cfg, truth)
    parts <- lapply(seq_len(nrow(todo)), function(i)
      d3d_part(cfg, truth, todo$rep[i], todo$variant[i], cal, log_fun = ctx$log))
  }
  per_rep <- data.table::rbindlist(lapply(parts, `[[`, "row"))
  q10 <- function(x) stats::quantile(x, 0.1, names = FALSE)
  q90 <- function(x) stats::quantile(x, 0.9, names = FALSE)
  summ <- per_rep[, .(cohorts = .N, n = round(mean(n)),
                      quick = round(mean(n_quick), 1), intermittent = round(mean(n_wave), 1),
                      no_recovery = round(mean(n_pers), 1),
                      dbic_median = round(stats::median(dbic_m1), 1),
                      dbic_q10 = round(q10(dbic_m1), 1), dbic_q90 = round(q90(dbic_m1), 1),
                      pct_dbic_gt_0 = round(100 * mean(dbic_m1 > 0)),
                      pct_dbic_gt_10 = round(100 * mean(dbic_m1 > 10)),
                      dbic_exact_median = round(stats::median(dbic_m1_exact), 1),
                      pct_exact_gt_10 = round(100 * mean(dbic_m1_exact > 10)),
                      q_hat = signif(stats::median(q_hat), 3), emax_hat = signif(stats::median(emax_hat), 3),
                      cor_log_rho = round(mean(cor_log_rho), 2), cor_log_g = round(mean(cor_log_g), 2),
                      cor_log_rho_exact = round(mean(cor_log_rho_exact), 2),
                      cor_log_g_exact = round(mean(cor_log_g_exact), 2),
                      twin_brier = round(mean(twin_brier), 3),
                      twin_auc_wave = round(mean(twin_auc_wave, na.rm = TRUE), 2),
                      twin_auc_no_recovery = round(mean(twin_auc_pers, na.rm = TRUE), 2),
                      phenotype_obs_vs_truth = round(mean(phenotype_obs_vs_truth), 2)), by = variant]
  truth_tab <- data.table::data.table(parameter = c(names(truth$struct), "mu_rho", "mu_g", "sd_rho", "sd_g",
                                                    paste0("sigma_", names(truth$sigma)), paste0("mix_", names(truth$mix))),
                                      value = signif(c(truth$struct, truth$mu_rho, truth$mu_g, truth$sd_rho, truth$sd_g,
                                                       truth$sigma, truth$mix), 4))
  list(per_rep = per_rep, summary = summ, truth = truth, settings = s,
       tables = list(summary = summ, per_rep = per_rep, truth = truth_tab),
       key_numbers = list(summary = summ))
}
