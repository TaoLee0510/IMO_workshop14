# =============================================================================
# an_D4_popfit.R
# Purpose : D4 - population fit on the TRAINING split only. Fits the nested
#           models in config$fit$models (M0/M1/M2; M0/M1 for the real-data
#           extract), compares them (overall and per phenotype), and tests
#           predictions P1-P5. Also writes the hold-out plan (model choice
#           fixed before the hold-out is touched).
# Inputs  : ctx$data() (train analysis set), config.
# Outputs : list(fits, compare, per_patient, P1..P5, examples, map_points,
#                tables, prereg)
#
# Predictions
#   P1 shape      : persistent patients - late plateau vs slow climb (day 45-90)
#   P2 timing     : second-wave nadir lags the inflammatory peak
#   P3 lineages   : persistent patients have low platelets as well
#   P4 risk       : clone fraction (covariate-free fit) rises with age,
#                   prior lines, prior ASCT, and with RDW / MCV when the CBC
#                   has them
#   P5 comparison : M1 helps second-wave patients, M2 helps persistent ones
# Mechanism attribution for persistent patients: the MAP call (d4_inferred) and
# posterior probabilities P(H1) / P(H3) (mechanism_probability); the expected
# number of clonal cases is the sum of P(H3), a lower bound (see the note there).
# P1, P4, the clone-trait profile and the mechanism probabilities need M2.
#
# Without M2 (real-data extract) two further checks:
#   event sensitivity : M1 refitted without the G-CSF effect, and with
#                       tocilizumab / steroids 2 days later (their days are
#                       approximations, R/events_approx.R)
#   late check        : counts after the fitting window (config$fit$obs_max_day)
#                       against the forward prediction of the last model
# =============================================================================

d4_inferred <- function(fit, d, cfg) {
  st <- fit$pop$struct
  rbindlist_safe(lapply(seq_len(nrow(fit$indiv)), function(i) {
    r <- fit$indiv[i]
    eta <- unlist(r[, paste0("eta_", fit$re), with = FALSE])
    names(eta) <- fit$re
    indiv <- eta_named_to_indiv(eta)
    prow <- d$patients[patient_id == r$patient_id]
    m <- infer_mechanism(indiv, st, fit$model, prow$anc_base, prow$ld_scale, cfg)
    data.table::data.table(patient_id = r$patient_id, phenotype_model = m$phenotype,
                           mechanism = m$mechanism, rho = indiv[["rho"]], g = indiv[["g"]],
                           f0 = if ("f0" %in% names(indiv)) indiv[["f0"]] else NA_real_)
  }))
}

rbindlist_safe <- function(x) data.table::rbindlist(x, fill = TRUE)

d4_P1 <- function(tr, ph, inferred, cfg) {
  s <- cfg$steps$P1
  pers <- ph[phenotype == "persistent", patient_id]
  o <- tr$obs[keep == TRUE & var == "anc" & patient_id %in% pers & day >= s$window[1] & day <= s$window[2]]
  shape <- o[, if (.N >= s$min_points) {
    fit <- stats::lm(log(value) ~ day)
    list(n = .N, slope_log_per_day = unname(stats::coef(fit)[2]), level_median = stats::median(value))
  }, by = patient_id]
  if (!nrow(shape)) return(list(table = NULL, cross = NULL))
  shape <- merge(shape, tr$patients[, .(patient_id, anc_base)], by = "patient_id")
  shape[, level_rel := level_median / anc_base]
  shape[, shape := data.table::fcase(slope_log_per_day >= s$slope_flat, "rising",
                                     abs(slope_log_per_day) < s$slope_flat & level_rel < s$level_low, "plateau",
                                     default = "other")]
  shape <- merge(shape, inferred[, .(patient_id, mechanism)], by = "patient_id", all.x = TRUE)
  list(table = shape, cross = shape[, .N, by = .(shape, mechanism)])
}

d4_P2 <- function(tr, ph, fit_m1, cfg) {
  s <- cfg$steps$P2
  wave <- ph[phenotype == "second_wave", patient_id]
  nadir <- tr$obs[keep == TRUE & var == "anc" & patient_id %in% wave &
                    day >= s$nadir_window[1] & day <= s$nadir_window[2],
                  .(nadir_day = day[which.min(value)]), by = patient_id]
  peak <- tr$obs[keep == TRUE & var == "ferritin" & patient_id %in% wave &
                   day >= s$peak_window[1] & day <= s$peak_window[2],
                 .(peak_day = day[which.max(value)]), by = patient_id]
  x <- merge(nadir, peak, by = "patient_id")
  x[, lag := nadir_day - peak_day]
  # model-implied lag for the typical patient
  st <- fit_m1$pop$struct
  typ <- c(rho = exp(fit_m1$pop$mu[["rho"]]), g = exp(fit_m1$pop$mu[["g"]]) * 1.8)
  sim <- simulate_anc_nat(typ, st, "M1", 3, 1, cfg)
  w <- sim$day >= s$nadir_window[1] & sim$day <= s$nadir_window[2]
  model_lag <- sim$day[w][which.min(sim$anc[w])] - sim$day[which.max(sim$is)]
  ct <- if (nrow(x) >= 5) suppressWarnings(stats::cor.test(x$peak_day, x$nadir_day, method = "spearman")) else NULL
  list(table = x,
       summary = data.table::data.table(
         n = nrow(x), median_lag_days = stats::median(x$lag), iqr_lag = paste(stats::quantile(x$lag, c(0.25, 0.75)), collapse = "-"),
         model_implied_lag = model_lag,
         spearman_rho = if (is.null(ct)) NA_real_ else round(unname(ct$estimate), 3),
         p_value = if (is.null(ct)) NA_real_ else signif(ct$p.value, 3)))
}

d4_P3 <- function(tr, ph, cfg) {
  s <- cfg$steps$P3
  x <- tr$obs[keep == TRUE & var == "plt" & day >= s$window[1] & day <= s$window[2],
              .(plt_late = stats::median(value)), by = patient_id]
  x <- merge(x, ph[, .(patient_id, phenotype)], by = "patient_id")
  x <- merge(x, tr$patients[, .(patient_id, plt_base)], by = "patient_id")
  x[, plt_low := plt_late < s$plt_low | plt_late < 0.5 * plt_base]
  tab <- x[!is.na(phenotype), .(n = .N, n_plt_low = sum(plt_low), pct_plt_low = round(100 * mean(plt_low), 1)),
           by = phenotype]
  ft <- NULL
  xx <- x[phenotype %in% c("persistent", "second_wave")]
  if (length(unique(xx$phenotype)) == 2) {
    ft <- stats::fisher.test(table(xx$phenotype == "persistent", xx$plt_low))
  }
  list(table = tab, test = data.table::data.table(
    comparison = "persistent vs second_wave",
    odds_ratio = if (is.null(ft)) NA_real_ else round(unname(ft$estimate), 2),
    ci = if (is.null(ft)) NA_character_ else paste(round(ft$conf.int, 2), collapse = "-"),
    p_value = if (is.null(ft)) NA_real_ else signif(ft$p.value, 3)))
}

d4_P4 <- function(fit_nocov, tr, cfg) {
  v <- intersect(cfg$model$covariates$f0, names(tr$patients))
  x <- merge(fit_nocov$indiv[, .(patient_id, eta_f0, se_f0)],
             tr$patients[, c("patient_id", v), with = FALSE], by = "patient_id")
  m <- stats::lm(stats::reformulate(v, "eta_f0"), data = x, weights = 1 / x$se_f0^2)
  cf <- summary(m)$coefficients
  ci <- suppressMessages(stats::confint(m))
  data.table::data.table(term = rownames(cf), estimate = round(cf[, 1], 4),
                         ci_low = round(ci[, 1], 4), ci_high = round(ci[, 2], 4),
                         p_value = signif(cf[, 4], 3))
}

d4_P5 <- function(fits, ph) {
  models <- intersect(c("M0", "M1", "M2"), names(fits))
  pp <- data.table::rbindlist(lapply(fits[models], function(f)
    f$indiv[, .(patient_id, logL, model = f$model)]))
  w <- data.table::dcast(pp, patient_id ~ model, value.var = "logL")
  w <- merge(w, ph[, .(patient_id, phenotype)], by = "patient_id")
  # gain = log-likelihood improvement; "better" = drop in -2 logL greater than 2
  w[, gain_M1_over_M0 := M1 - M0]
  if ("M2" %in% models) w[, gain_M2_over_M1 := M2 - M1]
  tab <- w[!is.na(phenotype), c(list(n = .N,
                                     sum_gain_M1_over_M0 = round(sum(gain_M1_over_M0), 1),
                                     pct_M1_better_by_2 = round(100 * mean(2 * gain_M1_over_M0 > 2), 1)),
                                if ("M2" %in% models)
                                  list(sum_gain_M2_over_M1 = round(sum(gain_M2_over_M1), 1),
                                       pct_M2_better_by_2 = round(100 * mean(2 * gain_M2_over_M1 > 2), 1))),
           by = phenotype]
  list(per_patient = w, table = tab)
}

# M1 refitted with modified treatment inputs: estimates and fit compared with
# the main M1 fit.
d4_event_sensitivity <- function(tr, fit_m1, cfg, log_fun) {
  variants <- list(
    no_gcsf = function(ev) ev[event != "gcsf"],
    crs_treatment_2d_later = function(ev) data.table::copy(ev)[event %in% c("tocilizumab", "steroid"),
                                                                `:=`(start_day = start_day + 2, end_day = end_day + 2)])
  quiet <- function(...) invisible(NULL)
  data.table::rbindlist(lapply(names(variants), function(v) {
    log_fun("  D4: event sensitivity - ", v)
    tv <- tr
    tv$events <- variants[[v]](tr$events)
    f <- fit_population(build_pdats(tv, cfg), "M1", cfg, init = fit_m1$pop, log_fun = quiet)
    x <- merge(fit_m1$indiv[, .(patient_id, r0 = eta_rho, g0 = eta_g)], f$indiv[, .(patient_id, eta_rho, eta_g)],
               by = "patient_id")
    data.table::data.table(variant = v, logLik = round(f$logLik, 1), delta_logLik = round(f$logLik - fit_m1$logLik, 1),
                           q = signif(f$pop$struct[["q"]], 3), emax = signif(f$pop$struct[["emax"]], 3),
                           cor_log_rho = round(stats::cor(x$r0, x$eta_rho), 3),
                           cor_log_g = round(stats::cor(x$g0, x$eta_g), 3),
                           median_change_log_rho = round(stats::median(x$eta_rho - x$r0), 3),
                           median_change_log_g = round(stats::median(x$eta_g - x$g0), 3))
  }))
}

# ANC observed after the fitting window versus the forward prediction of `fit`
# (MAP parameters, all recorded treatments).
d4_late_check <- function(fit, tr, cfg) {
  last <- cfg$fit$obs_max_day
  o <- tr$obs[keep == TRUE & var == "anc" & day > last & cens == FALSE]
  if (!nrow(o)) return(list(table = NULL, summary = NULL))
  pred <- data.table::rbindlist(lapply(unique(o$patient_id), function(pid) {
    oo <- o[patient_id == pid]
    prow <- tr$patients[patient_id == pid]
    r <- fit$indiv[patient_id == pid]
    eta <- unlist(r[, paste0("eta_", fit$re), with = FALSE])
    names(eta) <- fit$re
    seg <- make_segments(cfg$model$time$t0, max(oo$day) + 1, prow$ld_scale, tr$events[patient_id == pid])
    m <- marrow_sim(build_par(eta_named_to_indiv(eta), fit$pop$struct, fit$model), seg, oo$day)
    oo[, .(patient_id, day, anc_obs = value, anc_pred = round(m[, "anc_rel"] * prow$anc_base, 2),
           anc_base = prow$anc_base)]
  }))
  pred[, window := cut(day, c(last, 44, 75, Inf), c(paste0(last + 1, "-44"), "45-75", ">75"))]
  summ <- pred[, .(n_values = .N, n_patients = data.table::uniqueN(patient_id),
                   median_log_ratio = round(stats::median(log(anc_obs / anc_pred)), 2),
                   spearman = if (.N >= 5) round(stats::cor(anc_obs, anc_pred, method = "spearman"), 2) else NA_real_,
                   obs_below_1.5 = sum(anc_obs < 1.5), pred_below_1.5 = sum(anc_pred < 1.5),
                   obs_rel_base = round(stats::median(anc_obs / anc_base), 2),
                   pred_rel_base = round(stats::median(anc_pred / anc_base), 2)), by = window]
  list(table = pred, summary = summ)
}

d4_examples <- function(fits, tr, ph, cfg) {
  k <- cfg$steps$D4$n_examples_per_phenotype
  dens <- tr$obs[keep == TRUE & var == "anc", .N, by = patient_id]
  ex <- merge(ph[!is.na(phenotype)], dens, by = "patient_id")[order(-N)][, head(.SD, k), by = phenotype]
  days <- cfg$model$time$sim_days
  curves <- data.table::rbindlist(lapply(ex$patient_id, function(pid) {
    prow <- tr$patients[patient_id == pid]
    seg <- make_segments(cfg$model$time$t0, max(days), prow$ld_scale, tr$events[patient_id == pid])
    data.table::rbindlist(lapply(fits, function(f) {
      r <- f$indiv[patient_id == pid]
      eta <- unlist(r[, paste0("eta_", f$re), with = FALSE])
      names(eta) <- f$re
      m <- marrow_sim(build_par(eta_named_to_indiv(eta), f$pop$struct, f$model), seg, days)
      data.table::data.table(patient_id = pid, model = f$model, day = days,
                             anc = m[, "anc_rel"] * prow$anc_base, plt = m[, "plt_rel"] * prow$plt_base,
                             crp = prow$crp_base + cfg$model$obs$crp_scale * m[, "is"],
                             f_clone = m[, "f_clone"])
    }))
  }))
  obs <- tr$obs[patient_id %in% ex$patient_id]
  list(patients = ex[, .(patient_id, phenotype)], curves = curves, obs = obs)
}

an_D4_popfit <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  tr <- analysis_set(d, "train")
  ph <- tr$pheno[, .(patient_id, phenotype)]
  pds <- build_pdats(tr, cfg)
  models <- cfg$fit$models
  has_m2 <- "M2" %in% models
  fits <- list()
  ctx$log("  D4: fitting M0 (", length(pds), " training patients)")
  fits$M0 <- fit_population(pds, "M0", cfg, log_fun = ctx$log)
  ctx$log("  D4: fitting M1")
  fits$M1 <- fit_population(pds, "M1", cfg, init = fits$M0$pop, log_fun = ctx$log)
  prof <- NULL
  if (has_m2) {
    ctx$log("  D4: fitting M2 (covariate prior on f0)")
    fits$M2 <- fit_population(pds, "M2", cfg, init = fits$M1$pop, log_fun = ctx$log)
    ctx$log("  D4: profile likelihood for clone traits (delta x r_inf), M2 refitted at every point")
    prof <- profile_clone_traits(pds, fits$M2, cfg, log_fun = ctx$log)
    ctx$log("  D4: M2 at the profile maximum: delta = ", prof$best[["delta"]], ", r_inf = ", prof$best[["r_inf"]])
    fits$M2 <- prof$fit_best
    fits$M2$profile <- prof$table
    ctx$log("  D4: fitting M2 without covariates on f0 (for P4)")
    fits$M2_nocov <- fit_population(pds, "M2", cfg, cov_f0 = FALSE, init = fits$M2$pop, log_fun = ctx$log)
  }

  compare <- data.table::rbindlist(lapply(fits[models], function(f)
    data.table::data.table(model = f$model, n_patients = f$n, n_pop_params = f$n_par,
                           logLik = round(f$logLik, 1), BIC = round(f$BIC, 1),
                           its_iterations = nrow(f$trace), runtime_s = round(f$runtime_s))))
  compare[, delta_BIC := round(BIC - min(BIC), 1)]
  chosen <- compare[which.min(BIC), model]

  # mechanism attribution uses the largest model (only M2 can express a clonal
  # mechanism; with M1 a non-recovering patient is attributed to low reserve)
  inferred <- d4_inferred(fits[[utils::tail(models, 1)]], tr, cfg)
  inferred <- merge(inferred, ph, by = "patient_id", all.x = TRUE)
  mech_prob <- mech_summary <- P1 <- P4 <- sens <- late <- NULL
  if (has_m2) {
    ctx$log("  D4: mechanism probabilities for persistent patients")
    mech_prob <- merge(inferred[phenotype == "persistent", .(patient_id, mechanism_map = mechanism)],
                       mechanism_probability(fits$M2, tr, inferred[phenotype == "persistent", patient_id], cfg),
                       by = "patient_id")
    mech_summary <- mech_prob[, .(n_persistent = .N, map_h3 = sum(mechanism_map == "H3", na.rm = TRUE),
                                  expected_h3 = round(sum(p_h3), 1), mean_p_h3 = round(mean(p_h3), 3),
                                  expected_h1 = round(sum(p_h1), 1))]
    P1 <- d4_P1(tr, ph, inferred, cfg)
    P4 <- d4_P4(fits$M2_nocov, tr, cfg)
  } else {
    sens <- d4_event_sensitivity(tr, fits$M1, cfg, ctx$log)
    late <- d4_late_check(fits[[utils::tail(models, 1)]], tr, cfg)
  }
  P2 <- d4_P2(tr, ph, fits$M1, cfg)
  P3 <- d4_P3(tr, ph, cfg)
  P5 <- d4_P5(fits[models], ph)
  ex <- d4_examples(fits[models], tr, ph, cfg)

  pop_tab <- data.table::rbindlist(lapply(names(fits), function(nm) {
    f <- fits[[nm]]
    est <- cfg$fit$struct_est[[f$model]]
    data.table::data.table(fit = nm, parameter = c(paste0("mu_", names(f$pop$mu)), paste0("omega_", names(f$pop$omega)),
                                                  paste0("beta_", names(unlist(f$pop$beta))), paste0("sigma_", names(f$pop$sigma)), est),
                           value = round(c(f$pop$mu, f$pop$omega, unlist(f$pop$beta), f$pop$sigma, f$pop$struct[est]), 4))
  }))

  plan <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    data_source = cfg$data_source,
    rule = "model with the lowest BIC on the training split",
    chosen_model = chosen,
    landmark_day = cfg$steps$H$landmark_day,
    outcome = paste0("three-class phenotype (", paste(phenotype_display(cfg$phenotype), collapse = " / "),
                     "), rule in config$phenotype"),
    evaluation = if (identical(cfg$steps$H$mode, "cv")) paste0(cfg$steps$H$cv_folds, "-fold stratified cross-validation")
                 else "locked hold-out split",
    comparators = c("CAR-HEMATOTOX-like baseline score (multinomial logistic)",
                    "elastic-net multinomial logistic on day 0-10 features"),
    metrics = c("multiclass Brier score", "one-vs-rest AUC for second_wave and persistent",
                "calibration by risk bins", "bootstrap CIs"),
    excluded_from_holdout = if (is.null(d$post_planned)) character(0) else d$post_planned$patient_id,
    phenotype_rule = cfg$phenotype)

  list(fits = fits, compare = compare, chosen = chosen, inferred = inferred,
       mech_prob = mech_prob, mech_summary = mech_summary,
       P1 = P1, P2 = P2, P3 = P3, P4 = P4, P5 = P5, examples = ex, event_sensitivity = sens, late = late,
       tables = list(model_comparison = compare, population_parameters = pop_tab,
                     per_patient = merge(inferred, P5$per_patient[, c("patient_id", models), with = FALSE],
                                         by = "patient_id"),
                     P1_shape = P1$table, P1_cross = P1$cross, P2_timing = P2$summary, P2_per_patient = P2$table,
                     P3_platelets = P3$table, P3_test = P3$test, P4_risk_factors = P4,
                     P5_by_phenotype = P5$table, clone_trait_profile = prof$table,
                     mechanism_probability = mech_prob, mechanism_probability_summary = mech_summary,
                     event_sensitivity = sens, late_check = late$table, late_check_summary = late$summary),
       prereg = list(files = list(holdout_plan.json = plan), hash = "holdout_plan.sha256"),
       key_numbers = list(compare = compare, chosen = chosen, P4 = P4, P5 = P5$table,
                          P2 = P2$summary, P3 = P3$test, P1 = P1$cross, mechanism = mech_summary,
                          event_sensitivity = sens, late = late$summary))
}
