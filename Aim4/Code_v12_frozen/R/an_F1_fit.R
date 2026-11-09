# =============================================================================
# an_F1_fit.R
# Purpose : Step F1 - population fits of every model on the data, model
#           comparison on both axes (T2: Q2 when, T3: Q1 who) and vs the
#           Friberg reference MF, plus the pre-registered descriptive tests
#           T1 (nadir vs ferritin peak), T5 (dose-response), T6 (f0 vs risk
#           factors) and T7 (multi-lineage hit time).
# Inputs  : ctx (config, data, results PR).
# Outputs : list(fits, comparison, selected, tests, tables)
# =============================================================================

f1_models <- function(cfg) {
  if (cfg$quick) return(list(q2 = cfg$fit$q2_models, clone = "MabC_ddr", ref = "MF"))
  list(q2 = cfg$fit$q2_models, clone = vapply(cfg$fit$q2_models, model_id, character(1), clone = "ddr"),
       ref = "MF")
}

an_F1_fit <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  pdats <- build_pdats(d, cfg)
  ms <- f1_models(cfg)
  fits <- list()
  fit_one <- function(m) {
    fits[[m]] <<- fit_nested(pdats, m, fits, cfg)
    ctx$log(sprintf("  F1 %-10s logLik %8.1f  BIC %8.1f  (%.0fs%s)", m, fits[[m]]$logLik, fits[[m]]$BIC,
                    fits[[m]]$runtime_s, if (!is.null(fits[[m]]$start)) paste0(", start ", fits[[m]]$start) else ""))
  }
  for (m in c(ms$q2, ms$ref)) fit_one(m)
  bic <- function(ids) vapply(fits[ids], `[[`, numeric(1), "BIC")
  q2_bic <- sort(bic(ms$q2))
  q2_sel <- names(q2_bic)[1]
  for (m in ms$clone) fit_one(m)
  dta <- model_id(q2_sel, "dta")
  fit_one(dta)
  c_ids <- c(model_id(q2_sel, "ddr"), dta)
  best_c <- c_ids[which.min(bic(c_ids))]
  # the better of two fixed clone types counts as one extra parameter
  q1_delta <- (fits[[best_c]]$BIC + log(length(pdats))) - fits[[q2_sel]]$BIC
  final <- model_id(q2_sel, "ddr")

  comparison <- data.table::rbindlist(lapply(names(fits), function(m) data.table::data.table(
    model = m, family = fits[[m]]$spec$family, q2 = fits[[m]]$spec$q2, clone = fits[[m]]$spec$clone,
    logLik = fits[[m]]$logLik, n_par = fits[[m]]$n_par, BIC = fits[[m]]$BIC)))
  comparison[, delta_BIC := BIC - min(BIC)]
  data.table::setorder(comparison, BIC)
  selected <- list(q2 = q2_sel, q2_runner_up = names(q2_bic)[2], q2_margin = unname(q2_bic[2] - q2_bic[1]),
                   q2_supported = unname(q2_bic[2] - q2_bic[1]) >= 2,
                   q1_best_clone = best_c, q1_delta_bic = q1_delta, q1_supported = q1_delta <= -2,
                   mf_delta_bic = fits$MF$BIC - fits[[q2_sel]]$BIC, final = final)
  ctx$log(sprintf("  F1 Q2: %s (runner-up %s, +%.1f BIC); Q1: %s dBIC %.1f; MF +%.1f BIC", q2_sel,
                  selected$q2_runner_up, selected$q2_margin, best_c, q1_delta, selected$mf_delta_bic))

  # T6 needs the clone model without covariates on f0
  cfg_nc <- cfg
  cfg_nc$model$covariates$f0 <- character(0)
  pd_nc <- build_pdats(d, cfg_nc)
  fit_nc <- fit_population(pd_nc, final, cfg_nc, init = fits[[final]]$pop, eta_init = fit_eta(fits[[final]]),
                           log_fun = function(...) NULL)
  tests <- f1_tests(d, fits, selected, fit_nc, cfg)
  fitted <- f1_fitted_curves(d, pdats, fits[c(final, "MF")], cfg)
  list(fits = fits, comparison = comparison, selected = selected, tests = tests, fit_nocov = fit_nc,
       fitted = fitted, alc_obs = d$obs[var == "alc"],
       tables = list(model_comparison = comparison,
                     selected = data.table::as.data.table(lapply(selected, function(x) if (is.numeric(x)) round(x, 2) else x)),
                     population = f1_population_table(fits[[final]]),
                     individual = fits[[final]]$indiv,
                     T1 = tests$T1$data, T1_fit = tests$T1$fit, T5 = tests$T5, T6 = tests$T6, T7 = tests$T7$summary))
}

f1_population_table <- function(fit) {
  p <- fit$pop
  data.table::rbindlist(list(
    data.table::data.table(parameter = paste0("mu_", names(p$mu)), value = unname(p$mu), scale = "eta"),
    data.table::data.table(parameter = paste0("omega_", names(p$omega)), value = unname(p$omega), scale = "eta"),
    data.table::data.table(parameter = names(unlist(p$beta)), value = unname(unlist(p$beta)), scale = "beta"),
    data.table::data.table(parameter = names(p$struct), value = unname(p$struct), scale = "natural"),
    data.table::data.table(parameter = paste0("sigma_", names(p$sigma)), value = unname(p$sigma), scale = "log")))
}

# ---- pre-registered descriptive tests ------------------------------------------------------------
driver_peak_day <- function(dr) {
  tt <- dr$t0 + (seq_along(dr$v) - 1) * dr$dt
  sel <- tt >= 0 & tt <= 30
  tt[sel][which.max(dr$v[sel])]
}
driver_auc <- function(dr) {
  tt <- dr$t0 + (seq_along(dr$v) - 1) * dr$dt
  sel <- tt >= 0 & tt <= 30
  trapz(tt[sel], dr$v[sel])
}

second_wave_nadir_day <- function(a, rule) {
  a <- a[day >= 0 & day <= rule$window_end][order(day)]
  x <- ifelse(a$cens, rule$cens_value, a$value)
  nw <- which(a$day <= 14)
  nad <- nw[which.min(x[nw])]
  rec <- first_run(x >= rule$rec_thr, 1, nad + 1L)
  if (is.na(rec) || rec >= length(x)) return(NA_real_)
  after <- (rec + 1):length(x)
  a$day[after][which.min(x[after])]
}

f1_tests <- function(d, fits, selected, fit_nc, cfg) {
  rule <- phenotype_rule(cfg)
  ph <- d$phenotype
  pts <- d$patients
  # T1
  sw <- ph[primary == "second_wave", patient_id]
  t1 <- data.table::data.table(patient_id = sw)
  t1[, nadir2_day := vapply(patient_id, function(p) second_wave_nadir_day(d$obs[patient_id == p & var == "anc"], rule), numeric(1))]
  t1[, fer_peak_day := vapply(patient_id, function(p) driver_peak_day(d$drivers[[p]]), numeric(1))]
  t1[pts, gcsf := i.gcsf_after_d14, on = "patient_id"]
  t1_fit <- if (sum(stats::complete.cases(t1)) >= 3) {
    lm1 <- stats::lm(nadir2_day ~ fer_peak_day, data = t1)
    ci <- suppressWarnings(stats::confint(lm1)["fer_peak_day", ])
    data.table::data.table(n = stats::nobs(lm1), slope = stats::coef(lm1)[["fer_peak_day"]], lo = ci[[1]], hi = ci[[2]])
  } else data.table::data.table(n = sum(stats::complete.cases(t1)), slope = NA_real_, lo = NA_real_, hi = NA_real_)
  # T5 (individual estimates from the selected Q2 model, or Mab when it lacks a term)
  f5 <- if (all(c("th_as", "th_bs") %in% fits[[selected$q2]]$re)) fits[[selected$q2]] else fits$Mab
  e <- f5$indiv[, .(patient_id, eta_th_as, eta_th_bs)]
  e[pts, `:=`(bridging = i.bridging, max_crs = i.max_crs), on = "patient_id"]
  e[, drv_auc := vapply(patient_id, function(p) driver_auc(d$drivers[[p]]), numeric(1))]
  sp <- function(x, y) { ok <- stats::complete.cases(x, y); if (sum(ok) < 5) return(c(NA, NA))
    r <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman")); c(r$estimate, r$p.value) }
  wil <- function(x, g) { if (length(unique(g)) < 2) return(c(NA, NA))
    r <- suppressWarnings(stats::wilcox.test(x ~ g)); c(stats::median(x[g]) - stats::median(x[!g]), r$p.value) }
  t5 <- data.table::data.table(
    test = c("th_aS: bridging chemo vs none", "th_bS vs ferritin-driver AUC", "th_bS vs max CRS grade"),
    model = f5$model,
    estimate = c(wil(e$eta_th_as, e$bridging != "none")[1], sp(e$eta_th_bs, e$drv_auc)[1], sp(e$eta_th_bs, e$max_crs)[1]),
    p_value = c(wil(e$eta_th_as, e$bridging != "none")[2], sp(e$eta_th_bs, e$drv_auc)[2], sp(e$eta_th_bs, e$max_crs)[2]),
    statistic = c("difference in median eta", "Spearman rho", "Spearman rho"))
  # T6
  f0 <- fit_nc$indiv[, .(patient_id, f0 = stats::plogis(eta_f0))]
  f0[pts, `:=`(age = i.age, prior_lines = i.prior_lines, prior_asct = i.prior_asct), on = "patient_id"]
  t6 <- data.table::data.table(
    covariate = c("age", "prior lines", "prior ASCT"),
    estimate = c(sp(f0$f0, f0$age)[1], sp(f0$f0, f0$prior_lines)[1], wil(f0$f0, f0$prior_asct == 1)[1]),
    p_value = c(sp(f0$f0, f0$age)[2], sp(f0$f0, f0$prior_lines)[2], wil(f0$f0, f0$prior_asct == 1)[2]),
    statistic = c("Spearman rho", "Spearman rho", "difference in median f0"))
  # T7
  st <- fits[[selected$q2]]$pop$struct
  mtt <- c(anc = st[["mtt_n"]], plt = st[["mtt_pl"]], amc = st[["mtt_mo"]])
  t7 <- d$obs[keep == TRUE & var %in% lik_vars & day >= -5 & day <= 14, {
    x <- ifelse(cens, 0, value)
    .(nadir_day = day[which.min(x)])
  }, by = .(patient_id, var)]
  t7[, hit := nadir_day - mtt[var]]
  t7p <- t7[, .(hit_time = mean(hit)), by = patient_id]
  t7p[, fer_peak_day := vapply(patient_id, function(p) driver_peak_day(d$drivers[[p]]), numeric(1))]
  t7s <- data.table::data.table(median_hit_time = stats::median(t7p$hit_time),
                                iqr_lo = stats::quantile(t7p$hit_time, 0.25), iqr_hi = stats::quantile(t7p$hit_time, 0.75),
                                median_dist_to_ld = stats::median(abs(t7p$hit_time - (-4))),
                                median_dist_to_fer_peak = stats::median(abs(t7p$hit_time - t7p$fer_peak_day)))
  list(T1 = list(data = t1, fit = t1_fit), T5 = t5, T6 = t6, T7 = list(data = t7p, summary = t7s))
}

# Fitted trajectories at the MAP for every patient (for figures).
f1_fitted_curves <- function(d, pdats, fits, cfg) {
  data.table::rbindlist(lapply(names(fits), function(m) {
    f <- fits[[m]]
    pds <- prep_segments(pdats, f$pop$struct, cfg)
    tmpl <- par_template(f$pop$struct, f$spec)
    e <- fit_eta(f)
    data.table::rbindlist(lapply(names(pds), function(pid) {
      pd <- pds[[pid]]
      days <- seq(ceiling(pd$design$t0) + 1, max(pd$times))
      pd2 <- pd
      pd2$times <- sort(unique(c(days, pd$times)))
      m2 <- sim_model(fill_par(tmpl, e[pid, ], f$spec), pd$seg, pd$drv, pd2$times, f$spec, cfg)
      data.table::data.table(model = m, patient_id = pid, day = m2[, "t"], anc = m2[, "anc_rel"] * pd$base$anc,
                             plt = m2[, "plt_rel"] * pd$base$plt, amc = m2[, "amc_rel"] * pd$base$amc,
                             alc = if ("alc_rel" %in% colnames(m2) && !is.null(pd$base$alc)) m2[, "alc_rel"] * pd$base$alc else NA_real_,
                             S = m2[, "S"], I = m2[, "I"])
    }))
  }))
}
