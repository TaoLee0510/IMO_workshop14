# =============================================================================
# an_S_validation.R
# Purpose : Step S - synthetic validation of Model1 (plan section 8), from the
#           cohorts in Results/<data>/recovery/cohorts and the copied-back fits
#           in Results/<data>/hpc_runs:
#   S1 Model1 truth (synNN): structural and individual parameters vs truth, HSC
#      loss f_a, f_bI, f_bZ at the MAP vs truth, dominant source (a vs b) agreement
#   S2 the same cohorts refitted with w = 0 (synNNw0): bias of the phenotype term,
#      coverage of the 90% Laplace intervals of the random effects
#   S3 V1 on synNN (synNN_foldK): landmark AUROC when the model is correct
#   S4 Aim4 v1.3 truth (a4synNN) fitted with Model1: HSC loss and dominant source
#      agreement, clone fraction (the a/b and clone conclusions)
#   S5 visual predictive check of the final real-data fit (F1): phenotype counts,
#      nadir, recovery and second-wave days of replicate cohorts on the real designs
#   Results are 'if the model (S4: Aim4) is correct' and are kept apart from V1.
# Inputs  : ctx (config, data, results F1).
# Outputs : list(S1, S2, S3, S4, S5, tables)
# =============================================================================

s_cohorts <- function(cfg, pattern) {
  cdir <- file.path(cfg$paths$results, "recovery", "cohorts")
  sub("\\.rds$", "", list.files(cdir, pattern = paste0(pattern, "\\.rds$")))
}
s_cohort <- function(cfg, cn) readRDS(file.path(cfg$paths$results, "recovery", "cohorts", paste0(cn, ".rds")))
s_fit <- function(cfg, run) {
  f <- file.path(cfg$paths$results, "hpc_runs", run, "fit.rds")
  if (file.exists(f)) readRDS(f) else NULL
}

# HSC loss at the MAP of every patient of a cohort (counterfactual_draw13), with the truth
s_loss <- function(co, fit, cfg) {
  psi <- fit_psi13(fit)
  pdats <- prep_segments(build_pdats(co, cfg), fit$pop$struct, cfg)
  out <- data.table::rbindlist(lapply(names(pdats), function(pid) {
    x <- tryCatch(counterfactual_draw13(psi[pid, ], pdats[[pid]], fit, cfg), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    data.table::data.table(patient_id = pid, est_f_a = as.numeric(x[["f_a"]]), est_f_bi = as.numeric(x[["f_bi"]]),
                           est_f_bz = as.numeric(x[["f_bz"]]), est_actual = x[["actual"]])
  }))
  out[co$truth, `:=`(f_a = i.f_a, f_bi = i.f_bi, f_bz = i.f_bz), on = "patient_id"]
  out[, `:=`(f_b = 1 - (1 - f_bi) * (1 - f_bz), est_f_b = 1 - (1 - est_f_bi) * (1 - est_f_bz))]
  # dominant source among patients with a material HSC loss (truth, > 0.2 in total)
  out[, material := 1 - (1 - f_a) * (1 - f_b) > 0.2]
  out[, `:=`(dom = ifelse(f_a >= f_b, "a", "b"), est_dom = ifelse(est_f_a >= est_f_b, "a", "b"))]
  out
}

s_recovery <- function(cfg, pattern, struct = TRUE) {
  cns <- s_cohorts(cfg, pattern)
  per <- lapply(cns, function(cn) {
    fit <- s_fit(cfg, cn)
    if (is.null(fit)) return(NULL)
    co <- s_cohort(cfg, cn)
    st <- NULL
    if (struct) {
      est <- struct_est_names(fit$spec, cfg)
      st <- data.table::data.table(cohort = cn, parameter = est, truth = unname(co$truth_struct[est]),
                                   estimate = unname(fit$pop$struct[est]),
                                   prior_centre = unname(default_struct(fit$spec, cfg)[est]))
    }
    psi <- fit_psi13(fit)
    ind <- data.table::data.table(cohort = cn, patient_id = rownames(psi))
    for (k in intersect(fit$re, sub("^true_", "", names(co$truth)))) {
      tv <- co$truth[match(ind$patient_id, patient_id)][[paste0("true_", k)]]
      ind[, (paste0("true_", k)) := if (k == "f0") stats::qlogis(pmin(pmax(tv, 1e-4), 1 - 1e-4)) else log(pmax(tv, 1e-8))]
      ind[, (paste0("est_", k)) := psi[, k]]
      ind[, (paste0("se_", k)) := fit$indiv[match(ind$patient_id, patient_id)][[paste0("se_", k)]]]
    }
    loss <- s_loss(co, fit, cfg)[, cohort := cn]
    conv <- data.table::data.table(cohort = cn, rounds = fit$rounds, converged = fit$converged, BIC = fit$BIC,
                                   laplace_fallback = sum(fit$indiv$logL_source == "Laplace"))
    list(struct = st, indiv = ind, loss = loss, conv = conv)
  })
  per <- Filter(Negate(is.null), per)
  if (!length(per)) return(NULL)
  struct_t <- data.table::rbindlist(lapply(per, `[[`, "struct"))
  indiv <- data.table::rbindlist(lapply(per, `[[`, "indiv"), fill = TRUE)
  loss <- data.table::rbindlist(lapply(per, `[[`, "loss"))
  conv <- data.table::rbindlist(lapply(per, `[[`, "conv"))
  ssum <- if (nrow(struct_t)) {
    struct_t[, `:=`(ratio = estimate / truth, prior_ratio = prior_centre / truth)]
    struct_t[, .(n = .N, median_ratio = stats::median(ratio), q10 = stats::quantile(ratio, 0.1),
                 q90 = stats::quantile(ratio, 0.9), prior_over_truth = prior_ratio[1],
                 moved_towards_truth = mean(abs(log(ratio)) < abs(log(prior_ratio)))), by = parameter]
  } else NULL
  keys <- sub("^true_", "", grep("^true_", names(indiv), value = TRUE))
  isum <- data.table::rbindlist(lapply(keys, function(k) {
    x <- indiv[[paste0("true_", k)]]; y <- indiv[[paste0("est_", k)]]; se <- indiv[[paste0("se_", k)]]
    ok <- is.finite(x) & is.finite(y)
    data.table::data.table(parameter = k, n = sum(ok), pearson_r = stats::cor(x[ok], y[ok]),
                           bias = mean(y[ok] - x[ok]), rmse = sqrt(mean((y[ok] - x[ok])^2)),
                           coverage90 = mean(abs(y[ok] - x[ok]) <= 1.645 * se[ok], na.rm = TRUE))
  }))
  lsum <- loss[, .(n = .N, r_f_a = stats::cor(f_a, est_f_a), r_f_bi = stats::cor(f_bi, est_f_bi),
                   r_f_bz = stats::cor(f_bz, est_f_bz), r_f_b = stats::cor(f_b, est_f_b),
                   n_material = sum(material), dominant_agreement = mean((dom == est_dom)[material]),
                   truth_share_b = mean((dom == "b")[material]), est_share_b = mean((est_dom == "b")[material]))]
  list(struct = struct_t, indiv = indiv, loss = loss, conv = conv,
       summary = list(struct = ssum, indiv = isum, loss = lsum))
}

# S2: w = 1 vs w = 0 on the same cohorts
s_phenoterm <- function(cfg) {
  cns <- s_cohorts(cfg, "^syn[0-9]+")
  rows <- lapply(cns, function(cn) {
    f1 <- s_fit(cfg, cn); f0 <- s_fit(cfg, paste0(cn, "w0"))
    if (is.null(f1) || is.null(f0)) return(NULL)
    co <- s_cohort(cfg, cn)
    est <- struct_est_names(f1$spec, cfg)
    st <- data.table::data.table(cohort = cn, parameter = est, truth = unname(co$truth_struct[est]),
                                 w1 = unname(f1$pop$struct[est]), w0 = unname(f0$pop$struct[est]))
    ind <- data.table::rbindlist(lapply(list(w1 = f1, w0 = f0), function(f) {
      psi <- fit_psi13(f)
      data.table::rbindlist(lapply(intersect(f$re, sub("^true_", "", names(co$truth))), function(k) {
        tv <- co$truth[match(rownames(psi), patient_id)][[paste0("true_", k)]]
        x <- if (k == "f0") stats::qlogis(pmin(pmax(tv, 1e-4), 1 - 1e-4)) else log(pmax(tv, 1e-8))
        se <- f$indiv[match(rownames(psi), patient_id)][[paste0("se_", k)]]
        data.table::data.table(parameter = k, err = psi[, k] - x, covered = abs(psi[, k] - x) <= 1.645 * se)
      }))
    }), idcol = "w")
    ind[, cohort := cn]
    list(struct = st, indiv = ind)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  st <- data.table::rbindlist(lapply(rows, `[[`, "struct"))
  ind <- data.table::rbindlist(lapply(rows, `[[`, "indiv"))
  list(struct = st[, .(n = .N, log_bias_w1 = mean(log(w1 / truth)), log_bias_w0 = mean(log(w0 / truth)),
                       rmse_w1 = sqrt(mean(log(w1 / truth)^2)), rmse_w0 = sqrt(mean(log(w0 / truth)^2))), by = parameter],
       indiv = ind[, .(n = .N, bias = mean(err), rmse = sqrt(mean(err^2)), coverage90 = mean(covered, na.rm = TRUE)),
                   by = .(parameter, w)])
}

# S3: V1 on synthetic cohorts (landmark AUROC of non-quick recovery)
s_ceiling <- function(cfg) {
  cns <- s_cohorts(cfg, "^syn[0-9]+")
  cns <- cns[as.integer(sub("syn", "", cns)) <= cfg$validation$S3$n_cohorts]
  rows <- lapply(cns, function(cn) {
    full <- s_fit(cfg, cn)
    ff <- lapply(seq_len(cfg$validation$V1$k), function(k) s_fit(cfg, paste0(cn, "_fold", k)))
    if (is.null(full) || any(vapply(ff, is.null, logical(1)))) return(NULL)
    co <- s_cohort(cfg, cn)
    pred <- landmark_predictions(co, full, ff, v1_folds(co$phenotype, cfg), cfg)
    pred[co$phenotype, observed := i.primary, on = "patient_id"]
    pred[, cohort := cn]
    pred
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  pred <- data.table::rbindlist(rows, fill = TRUE)
  list(pred = pred, metrics = pred[, .(n = .N, auroc = auroc(1 - p_quick, as.numeric(observed != "quick"))),
                                   by = .(cohort, L)])
}

# S5: visual predictive check on the real designs from the final fit's population
s_vpc <- function(d, fit, cfg) {
  n_rep <- cfg$validation$S5$n_rep
  pdats <- prep_segments(build_pdats(d, cfg), fit$pop$struct, cfg)
  tmpl <- par_template(fit$pop$struct, fit$spec)
  rule <- phenotype_rule(cfg)
  sg <- fit$pop$sigma[["anc"]]
  box <- NULL
  one <- function(r) {
    set.seed(cfg$runtime$seed + r)
    data.table::rbindlist(lapply(names(pdats), function(pid) {
      pd <- pdats[[pid]]
      prior <- patient_prior(pd, fit$pop, fit$re)
      bx <- indiv_box(prior, cfg)
      psi <- pmin(pmax(prior$mean + prior$sd * stats::rnorm(length(prior$mean)), bx$lo), bx$hi)
      lab <- pheno_label(d$obs[patient_id == pid & var == "anc" & !is.na(value)], cfg)
      dd <- if (is.null(lab)) 0:30 else lab$days
      pd$times <- sort(unique(c(pd$times, dd)))
      m <- simulate_pdat(psi, pd, tmpl, fit$spec, cfg)
      if (anyNA(m[, 2])) return(NULL)
      anc <- m[match(dd, m[, "t"]), "anc_rel"] * pdat_baselines(psi, pd)[["anc"]] * exp(stats::rnorm(length(dd), 0, sg))
      cens <- anc < cfg$model$design$cens$anc
      ph <- phenotype_series(dd, ifelse(cens, cfg$model$design$cens$anc_bound, anc), cens, rule)
      x <- anc[dd >= 0 & dd <= 14]
      data.table::data.table(rep = r, patient_id = pid, phenotype = ph,
                             nadir_day = if (length(x)) dd[dd >= 0 & dd <= 14][which.min(x)] else NA_real_)
    }))
  }
  sims <- data.table::rbindlist(par_lapply(seq_len(n_rep), one, cfg$runtime$cores))
  counts <- sims[, .(quick = sum(phenotype == "quick"), second_wave = sum(phenotype == "second_wave"),
                     delayed = sum(phenotype == "delayed"), median_nadir_day = stats::median(nadir_day, na.rm = TRUE)), by = rep]
  obs <- table(factor(d$phenotype$primary, levels = cfg$phenotype$classes))
  q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  summ <- data.table::data.table(
    quantity = c("quick", "second_wave", "delayed", "median_nadir_day"),
    observed = c(obs[["quick"]], obs[["second_wave"]], obs[["delayed"]],
                 stats::median(d$obs[var == "anc" & day >= 0 & day <= 14, .(nd = day[which.min(ifelse(cens, 0, value))]),
                                     by = patient_id]$nd)),
    sim_q05 = c(q(counts$quick)[1], q(counts$second_wave)[1], q(counts$delayed)[1], q(counts$median_nadir_day)[1]),
    sim_median = c(q(counts$quick)[2], q(counts$second_wave)[2], q(counts$delayed)[2], q(counts$median_nadir_day)[2]),
    sim_q95 = c(q(counts$quick)[3], q(counts$second_wave)[3], q(counts$delayed)[3], q(counts$median_nadir_day)[3]))
  summ[, inside_90 := observed >= sim_q05 & observed <= sim_q95]
  list(counts = counts, summary = summ)
}

an_S_validation <- function(ctx) {
  cfg <- ctx$config
  s1 <- s_recovery(cfg, "^syn[0-9]+")
  s4 <- s_recovery(cfg, "^a4syn[0-9]+", struct = FALSE)
  s2 <- s_phenoterm(cfg)
  s3 <- s_ceiling(cfg)
  f1 <- ctx$results$F1
  s5 <- if (!is.null(f1)) s_vpc(ctx$data(), f1$fits[[f1$selected$final]], cfg) else NULL
  ctx$log(sprintf("  S: S1 %s cohorts; S2 %s; S3 %s; S4 %s cohorts; S5 %s",
                  if (is.null(s1)) 0 else nrow(s1$conv), if (is.null(s2)) "none" else "done",
                  if (is.null(s3)) "none" else "done", if (is.null(s4)) 0 else nrow(s4$conv),
                  if (is.null(s5)) "none" else paste(s5$summary$inside_90, collapse = "/")))
  tb <- list()
  if (!is.null(s1)) tb <- c(tb, list(S1_struct = s1$summary$struct, S1_indiv = s1$summary$indiv, S1_loss = s1$summary$loss,
                                     S1_convergence = s1$conv))
  if (!is.null(s2)) tb <- c(tb, list(S2_struct = s2$struct, S2_indiv = s2$indiv))
  if (!is.null(s3)) tb <- c(tb, list(S3_auroc = s3$metrics))
  if (!is.null(s4)) tb <- c(tb, list(S4_indiv = s4$summary$indiv, S4_loss = s4$summary$loss, S4_convergence = s4$conv))
  if (!is.null(s5)) tb <- c(tb, list(S5_vpc = s5$summary))
  list(S1 = s1, S2 = s2, S3 = s3, S4 = s4, S5 = s5, tables = tb)
}
