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

# The fits are computed on the cluster (round 0 + DE-ITS, hpc/de_*.sbatch); F1 reads the
# assembled fits from Results/<data>/hpc_runs/<model>/fit.rds (copied back), applies the
# pre-registered selection rules and runs the descriptive tests. Only runs named by a model
# id enter (not the synthetic cohorts, V1 folds or sensitivity refits).
f1_read_fits <- function(cfg) {
  runs <- file.path(cfg$paths$results, "hpc_runs")
  ids <- list.dirs(runs, recursive = FALSE, full.names = FALSE)
  ids <- ids[vapply(ids, function(m) !is.null(tryCatch(model_spec(m), error = function(e) NULL)), logical(1))]
  fits <- lapply(stats::setNames(ids, ids), function(m) {
    f <- file.path(runs, m, "fit.rds")
    if (file.exists(f)) readRDS(f) else NULL
  })
  Filter(Negate(is.null), fits)
}

# psi matrix (rows = patients) of a v1.3 fit: random effects (eta_*) and baseline deviations (b_*)
fit_psi13 <- function(fit) {
  nm <- c(fit$re, grep("^b_", names(fit$indiv), value = TRUE))
  nm <- nm[!grepl("^b_.*_", nm)]
  cols <- ifelse(nm %in% fit$re, paste0("eta_", nm), nm)
  e <- as.matrix(fit$indiv[, cols, with = FALSE])
  colnames(e) <- nm
  rownames(e) <- fit$indiv$patient_id
  e
}

f1_select13 <- function(fits, cfg, n = 35) {
  bic <- function(m) if (m %in% names(fits)) fits[[m]]$BIC else NA_real_
  q2_ids <- intersect(cfg$fit$q2_models, names(fits))
  conv <- vapply(fits[q2_ids], function(f) isTRUE(f$converged), logical(1))
  pool <- if (any(conv)) q2_ids[conv] else q2_ids
  q2_bic <- sort(vapply(fits[pool], `[[`, numeric(1), "BIC"))
  q2_sel <- names(q2_bic)[1]
  c_ids <- intersect(c(model_id(q2_sel, "ddr"), model_id(q2_sel, "dta")), names(fits))
  best_c <- if (length(c_ids)) c_ids[which.min(vapply(c_ids, bic, numeric(1)))] else NA_character_
  q1_delta <- if (is.na(best_c)) NA_real_ else bic(best_c) + log(n) - bic(q2_sel)
  clone_ok <- isTRUE(q1_delta <= -2)
  ctype <- if (clone_ok) model_spec(best_c)$clone else "none"
  base <- model_id(q2_sel, ctype)
  r_ids <- intersect(c(base, model_id(q2_sel, ctype, "hsc"), model_id(q2_sel, ctype, "mpp")), names(fits))
  r_bic <- vapply(r_ids, bic, numeric(1))
  final <- r_ids[which.min(r_bic)]
  route_delta <- c(drop_mpp_route = bic(model_id(q2_sel, ctype, "hsc")) - bic(base),
                   drop_hsc_route = bic(model_id(q2_sel, ctype, "mpp")) - bic(base))
  mf_delta <- if ("MF" %in% names(fits)) fits$MF$BIC_core - fits[[q2_sel]]$BIC_core else NA_real_
  list(q2 = q2_sel, q2_runner_up = names(q2_bic)[2], q2_margin = unname(q2_bic[2] - q2_bic[1]),
       q2_supported = isTRUE(unname(q2_bic[2] - q2_bic[1]) >= 2), q2_pool_converged = any(conv),
       q1_best_clone = best_c, q1_delta_bic = q1_delta, q1_supported = clone_ok, clone_type = ctype,
       route_models = r_ids, route_bic = r_bic, route_delta = route_delta, final = final,
       final_route = model_spec(final)$route, mf_delta_bic_core = mf_delta,
       final_converged = isTRUE(fits[[final]]$converged))
}

an_F1_fit <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  fits <- f1_read_fits(cfg)
  if (!length(fits)) stop("no assembled fits under ", file.path(cfg$paths$results, "hpc_runs"))
  pdats <- build_pdats(d, cfg)
  comparison <- data.table::rbindlist(lapply(fits, function(f) data.table::data.table(
    model = f$model, family = f$spec$family, clone = f$spec$clone, route = f$spec$route %||% "none",
    logLik = f$logLik, logLik_laplace = f$logLik_laplace, logLik_core = f$logLik_core, n_par = f$n_par,
    BIC = f$BIC, BIC_core = f$BIC_core, rounds = f$rounds, converged = f$converged,
    converged_by = if (!is.null(f$converged_by)) f$converged_by else if (isTRUE(f$converged)) "change" else NA_character_,
    laplace_fallback = sum(f$indiv$logL_source == "Laplace"), n_at_bound = sum(nzchar(f$indiv$at_bound)))))
  comparison[, delta_BIC := BIC - min(BIC[family == "m1"])]
  data.table::setorder(comparison, BIC)
  selected <- f1_select13(fits, cfg, n = length(pdats))
  final <- selected$final
  ctx$log(sprintf("  F1 Q2: %s (runner-up %s, +%.1f BIC%s); Q1: %s dBIC %.1f (%s); route: %s; final %s%s; MF +%.1f BIC (ANC/PLT/AMC)",
                  selected$q2, selected$q2_runner_up, selected$q2_margin,
                  if (selected$q2_pool_converged) "" else ", none converged",
                  selected$q1_best_clone, selected$q1_delta_bic, if (selected$q1_supported) "supported" else "not supported",
                  selected$final_route, final, if (selected$final_converged) "" else " (NOT converged)",
                  selected$mf_delta_bic_core))
  fit_nc <- fits[[paste0(model_id(selected$q2, "ddr"), "_nc")]]
  tests <- f1_tests(d, fits, selected, fit_nc, cfg)
  fitted <- f1_fitted_curves13(d, pdats, fits[intersect(c(final, "MF"), names(fits))], cfg)
  tests$T8 <- f1_t8(d, fitted[model == final])
  list(fits = fits, comparison = comparison, selected = selected, tests = tests, fit_nocov = fit_nc,
       fitted = fitted, alc_obs = d$obs[var == "alc"],
       tables = list(model_comparison = comparison,
                     selected = data.table::as.data.table(lapply(selected[c("q2", "q2_runner_up", "q2_margin", "q2_supported",
                                                                         "q1_best_clone", "q1_delta_bic", "q1_supported",
                                                                         "final", "final_route", "mf_delta_bic_core",
                                                                         "final_converged")],
                                                                function(x) if (is.numeric(x)) round(x, 2) else x)),
                     route = data.table::data.table(model = names(selected$route_bic), BIC = round(selected$route_bic, 2)),
                     population = f1_population_table(fits[[final]]),
                     individual = fits[[final]]$indiv,
                     T1 = tests$T1$data, T1_fit = tests$T1$fit, T5 = tests$T5, T6 = tests$T6, T7 = tests$T7$summary,
                     T8 = tests$T8))
}

# T8: model CAR T peak (final model at the MAP) vs the measured peak expansion
f1_t8 <- function(d, fitted_final) {
  pk <- fitted_final[day >= 0 & day <= 30, .(model_peak = max(cart), model_peak_day = day[which.max(cart)]), by = patient_id]
  pk[d$patients, measured := i.peak_expansion, on = "patient_id"]
  ok <- stats::complete.cases(pk$model_peak, pk$measured)
  r <- if (sum(ok) >= 5) suppressWarnings(stats::cor.test(pk$model_peak[ok], pk$measured[ok], method = "spearman")) else NULL
  data.table::data.table(n = sum(ok), spearman_rho = if (is.null(r)) NA_real_ else unname(r$estimate),
                         p_value = if (is.null(r)) NA_real_ else r$p.value,
                         median_model_peak = stats::median(pk$model_peak), median_peak_day = stats::median(pk$model_peak_day))
}

f1_population_table <- function(fit) {
  p <- fit$pop
  data.table::rbindlist(list(
    data.table::data.table(parameter = paste0("mu_", names(p$mu)), value = unname(p$mu), scale = "eta"),
    data.table::data.table(parameter = paste0("omega_", names(p$omega)), value = unname(p$omega), scale = "eta"),
    data.table::data.table(parameter = names(unlist(p$beta)), value = unname(unlist(p$beta)), scale = "beta"),
    data.table::data.table(parameter = names(p$struct), value = unname(p$struct), scale = "natural"),
    data.table::data.table(parameter = paste0("sigma_", names(p$sigma)), value = unname(p$sigma), scale = "natural"),
    if (length(p$tau)) data.table::data.table(parameter = paste0("tau_", names(p$tau)), value = unname(p$tau),
                                              scale = "natural")))
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
  if (is.null(f5)) f5 <- fits[[selected$final]]
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
  # T6 (needs the clone model fitted without covariates on f0)
  if (is.null(fit_nc)) {
    t6 <- data.table::data.table(covariate = c("age", "prior lines", "prior ASCT"), estimate = NA_real_,
                                 p_value = NA_real_, statistic = "not run (no fit without covariates on f0)")
  } else {
  f0 <- fit_nc$indiv[, .(patient_id, f0 = stats::plogis(eta_f0))]
  f0[pts, `:=`(age = i.age, prior_lines = i.prior_lines, prior_asct = i.prior_asct), on = "patient_id"]
  t6 <- data.table::data.table(
    covariate = c("age", "prior lines", "prior ASCT"),
    estimate = c(sp(f0$f0, f0$age)[1], sp(f0$f0, f0$prior_lines)[1], wil(f0$f0, f0$prior_asct == 1)[1]),
    p_value = c(sp(f0$f0, f0$age)[2], sp(f0$f0, f0$prior_lines)[2], wil(f0$f0, f0$prior_asct == 1)[2]),
    statistic = c("Spearman rho", "Spearman rho", "difference in median f0"))
  }
  # T7
  st <- fits[[selected$q2]]$pop$struct
  mtt <- c(anc = st[["mtt_n"]], plt = st[["mtt_pl"]], amc = st[["mtt_mo"]])
  t7 <- d$obs[keep == TRUE & var %in% core_vars & day >= -5 & day <= 14, {
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

# Fitted trajectories at the MAP for every patient (v1.3: individual baselines, ALC =
# lymphocytes + CAR T, CAR T compartment, signal Z, G-CSF effect)
f1_fitted_curves13 <- function(d, pdats, fits, cfg) {
  data.table::rbindlist(lapply(names(fits), function(m) {
    f <- fits[[m]]
    pds <- prep_segments(pdats, f$pop$struct, cfg)
    tmpl <- par_template(f$pop$struct, f$spec)
    e <- fit_psi13(f)
    data.table::rbindlist(lapply(names(pds), function(pid) {
      pd <- pds[[pid]]
      psi <- e[pid, ]
      days <- seq(ceiling(pd$design$t0) + 1, max(pd$times))
      tt <- sort(unique(c(days, pd$times)))
      m2 <- sim_model(fill_par(tmpl, psi, f$spec, pd), pd$seg, pd$drv, tt, f$spec, cfg)
      b <- pdat_baselines(psi, pd)
      has <- function(k) k %in% colnames(m2)
      data.table::data.table(model = m, patient_id = pid, day = m2[, "t"], anc = m2[, "anc_rel"] * b[["anc"]],
                             plt = m2[, "plt_rel"] * b[["plt"]], amc = m2[, "amc_rel"] * b[["amc"]],
                             alc = if (has("alc_rel")) m2[, "alc_rel"] * b[["alc"]] + m2[, "cart"] else NA_real_,
                             cart = if (has("cart")) m2[, "cart"] else NA_real_,
                             z = if (has("z")) m2[, "z"] else NA_real_, gcsf = if (has("gcsf")) m2[, "gcsf"] else NA_real_,
                             S = m2[, "S"], I = m2[, "I"])
    }))
  }))
}
