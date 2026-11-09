# =============================================================================
# an_H_holdout.R
# Purpose : H - the single pre-registered hold-out evaluation. At the day +10
#           landmark, three methods predict the three-class phenotype for the
#           test split:
#             twin      : model chosen in D4, individual fit to data <= day 10,
#                         posterior draws simulated forward and classified
#             baseline  : CAR-HEMATOTOX-like score (simplified; thresholds to
#                         verify against Rejeski et al. 2021) -> multinomial
#             ml        : elastic-net multinomial on day 0-10 features
#           Comparators are trained on the training split only.
#           With config$steps$H$mode = "cv" (real-data extract, too small for a
#           locked hold-out) the same three methods are evaluated by stratified
#           k-fold cross-validation instead: in every fold the population (the
#           model chosen in D4, started from the D4 estimates) and both
#           comparators are refitted on the other folds.
# Inputs  : ctx$data(), ctx$results$D4, ctx$holdout_plan (verified).
# Outputs : list(pred, metrics, calibration, tables)
# =============================================================================

hematotox_like_score <- function(p) {
  with(p, (anc_base < 1.2) + (hgb_base < 9) + ifelse(plt_base < 75, 2, ifelse(plt_base <= 175, 1, 0)) +
         (crp_base > 30) + ifelse(fer_base > 2000, 2, ifelse(fer_base >= 650, 1, 0)))
}

landmark_features <- function(s, landmark) {
  o <- s$obs[keep == TRUE & day >= 0 & day <= landmark]
  f <- s$patients[, .(patient_id, log_anc_base = log(anc_base), log_plt_base = log(plt_base),
                      hgb_base, log_crp_base = log(crp_base), log_fer_base = log(fer_base),
                      costim_cd28, age, prior_lines, prior_asct)]
  agg <- o[, .(anc_min = log(min(value[var == "anc"])), anc_last = log(value[var == "anc"][which.max(day[var == "anc"])]),
               plt_last = log(value[var == "plt"][which.max(day[var == "plt"])]),
               crp_max = log(max(c(value[var == "crp"], 0.5))), fer_max = log(max(c(value[var == "ferritin"], 1))),
               anc_slope = {
                 a <- var == "anc" & day >= 7
                 if (sum(a) >= 2) unname(stats::coef(stats::lm(log(value[a]) ~ day[a]))[2]) else 0
               }), by = patient_id]
  ev <- s$events[start_day <= landmark, .(toci = as.integer(any(event == "tocilizumab")),
                                          steroid = as.integer(any(event == "steroid"))), by = patient_id]
  f <- merge(f, agg, by = "patient_id", all.x = TRUE)
  f <- merge(f, ev, by = "patient_id", all.x = TRUE)
  f[is.na(toci), toci := 0L]
  f[is.na(steroid), steroid := 0L]
  f
}

# seed: each patient's draws are seeded inside the worker (seed + position),
# because a seed set before mclapply does not reach the forked workers
twin_probabilities <- function(fits, pds, fit, d, cfg, n_draws, seed) {
  st <- fit$pop$struct
  days <- cfg$model$time$sim_days
  par_lapply(names(fits), function(pid) {
    set.seed(seed + match(pid, names(fits)))
    f <- fits[[pid]]
    prow <- d$patients[patient_id == pid]
    draws <- MASS::mvrnorm(n_draws, f$eta, f$cov)
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1, dimnames = list(NULL, names(f$eta)))
    ph <- apply(draws, 1, function(e) {
      e <- stats::setNames(e, names(f$eta))
      m <- marrow_sim(build_par(eta_named_to_indiv(e), st, fit$model), pds[[pid]]$seg, days)
      classify_phenotype(days, m[, "anc_rel"] * prow$anc_base, cfg$phenotype)
    })
    tb <- table(factor(ph, phenotype_levels))
    p <- (as.numeric(tb) + 0.5) / (sum(tb) + 1.5)          # light smoothing avoids 0/1
    data.table::data.table(patient_id = pid, method = "twin", p_quick = p[1], p_second_wave = p[2],
                           p_persistent = p[3])
  }, cfg$runtime$cores)
}

multiclass_brier <- function(P, y) {
  Y <- outer(y, phenotype_levels, `==`) * 1
  rowSums((P - Y)^2)
}

# Predictions of the three methods for the patients `ids` of `te`: the twin
# uses the population `fit`, the comparators are trained on `tr`.
h_predict <- function(tr, te, ids, fit, cfg, s, seed_offset = 0L) {
  te_sub <- te
  te_sub$patients <- te$patients[patient_id %in% ids]

  # twin
  pds <- build_pdats(te_sub, cfg, landmark = s$landmark_day, use_vaf = TRUE)
  fits <- fit_given_population(pds, fit, cfg)
  p_twin <- data.table::rbindlist(twin_probabilities(fits, pds, fit, te, cfg, s$n_draws,
                                                     seed = cfg$runtime$seed + 8L + seed_offset))

  # baseline score
  ytr <- tr$pheno[!is.na(phenotype)]
  btr <- merge(ytr, tr$patients, by = "patient_id")
  btr[, score := hematotox_like_score(btr)]
  bte <- te$patients[match(ids, patient_id)]
  bte[, score := hematotox_like_score(bte)]
  mb <- nnet::multinom(factor(phenotype, phenotype_levels) ~ score, data = btr, trace = FALSE)
  pb <- stats::predict(mb, newdata = bte, type = "probs")
  p_base <- data.table::data.table(patient_id = ids, method = "baseline_score",
                                   p_quick = pb[, 1], p_second_wave = pb[, 2], p_persistent = pb[, 3])

  # elastic-net multinomial
  ftr <- merge(ytr[, .(patient_id, phenotype)], landmark_features(tr, s$landmark_day), by = "patient_id")
  fte <- landmark_features(te_sub, s$landmark_day)[match(ids, patient_id)]
  xcols <- setdiff(names(ftr), c("patient_id", "phenotype"))
  med <- ftr[, lapply(.SD, stats::median, na.rm = TRUE), .SDcols = xcols]
  fill <- function(x) { for (c in xcols) data.table::set(x, which(is.na(x[[c]])), c, med[[c]]); x }
  ftr <- fill(ftr)
  fte <- fill(fte)
  set.seed(cfg$runtime$seed + 9L + seed_offset)
  pm <- tryCatch({
    cvf <- glmnet::cv.glmnet(as.matrix(ftr[, ..xcols]), factor(ftr$phenotype, phenotype_levels),
                             family = "multinomial", alpha = 0.5, nfolds = 5)
    stats::predict(cvf, as.matrix(fte[, ..xcols]), s = "lambda.min", type = "response")[, , 1]
  }, error = function(e) {
    # too few patients in a class to fit: training class frequencies
    fr <- as.numeric(prop.table(table(factor(ftr$phenotype, phenotype_levels))))
    matrix(fr, nrow = length(ids), ncol = 3, byrow = TRUE)
  })
  if (is.null(dim(pm))) pm <- matrix(pm, nrow = 1)
  p_ml <- data.table::data.table(patient_id = ids, method = "ml_elastic_net",
                                 p_quick = pm[, 1], p_second_wave = pm[, 2], p_persistent = pm[, 3])

  rbind(p_twin[match(ids, patient_id)], p_base, p_ml)
}

# Metrics with bootstrap CIs, Brier differences and calibration bins; `pred`
# holds every method's rows in the order of `y`.
h_evaluate <- function(pred, y, cfg, s) {
  # metrics with bootstrap CIs
  metr <- function(P, yy) {
    aucs <- vapply(c("second_wave", "persistent"), function(k) {
      yy_k <- as.integer(yy == k)
      if (length(unique(yy_k)) < 2) return(NA_real_)
      as.numeric(pROC::auc(yy_k, P[, k == phenotype_levels], quiet = TRUE))
    }, numeric(1))
    c(brier = mean(multiclass_brier(P, yy)), auc_second_wave = aucs[[1]], auc_persistent = aucs[[2]])
  }
  pmat <- function(m) as.matrix(pred[method == m, .(p_quick, p_second_wave, p_persistent)])
  methods <- unique(pred$method)
  point <- sapply(methods, function(m) metr(pmat(m), y))
  set.seed(cfg$runtime$seed + 10L)
  boot <- replicate(s$n_boot, {
    i <- sample(length(y), replace = TRUE)
    sapply(methods, function(m) metr(pmat(m)[i, , drop = FALSE], y[i]))
  })
  metrics <- data.table::rbindlist(lapply(methods, function(m) data.table::rbindlist(lapply(rownames(point), function(k)
    data.table::data.table(method = m, metric = k, estimate = round(point[k, m], 3),
                           ci_low = round(stats::quantile(boot[k, m, ], 0.025, na.rm = TRUE), 3),
                           ci_high = round(stats::quantile(boot[k, m, ], 0.975, na.rm = TRUE), 3))))))
  diffs <- data.table::rbindlist(lapply(setdiff(methods, "twin"), function(m) {
    dd <- boot["brier", "twin", ] - boot["brier", m, ]
    data.table::data.table(comparison = paste("twin -", m), metric = "brier",
                           difference = round(point["brier", "twin"] - point["brier", m], 3),
                           ci_low = round(stats::quantile(dd, 0.025, na.rm = TRUE), 3),
                           ci_high = round(stats::quantile(dd, 0.975, na.rm = TRUE), 3))
  }))
  # risk quintiles; one bin when a method gives (nearly) the same probability to everyone
  quintile_bin <- function(p) {
    b <- unique(stats::quantile(p, seq(0, 1, 0.2)))
    if (length(b) < 2) factor(rep(paste0("[", signif(b[1], 3), "]"), length(p)))
    else cut(p, breaks = b, include.lowest = TRUE)
  }
  calib <- data.table::rbindlist(lapply(c("p_second_wave", "p_persistent"), function(col) {
    k <- sub("^p_", "", col)
    pred[, .(method, p = get(col), obs = as.integer(observed == k))][
      , bin := quintile_bin(p), by = method][
      , .(mean_pred = mean(p), obs_rate = mean(obs), n = .N), by = .(method, bin)][, outcome := k]
  }))
  list(metrics = metrics, diffs = diffs, calibration = calib)
}

an_H_holdout <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$H
  if (identical(s$mode, "cv")) return(h_cross_validation(ctx))
  plan <- ctx$holdout_plan
  d <- ctx$data()
  tr <- analysis_set(d, "train")
  te <- analysis_set(d, "test")
  ids <- te$pheno[!is.na(phenotype) & !patient_id %in% plan$excluded_from_holdout, patient_id]
  y <- te$pheno[match(ids, patient_id), phenotype]

  fit <- ctx$results$D4$fits[[plan$chosen_model]]
  ctx$log("  H: twin fits at day ", s$landmark_day, " for ", length(ids), " test patients (", plan$chosen_model, ")")
  pred <- h_predict(tr, te, ids, fit, cfg, s)
  pred[, observed := y[match(patient_id, ids)]]
  ev <- h_evaluate(pred, y, cfg, s)

  list(pred = pred, metrics = ev$metrics, diffs = ev$diffs, calibration = ev$calibration, n_test = length(ids),
       observed = table(factor(y, phenotype_levels)), plan = plan,
       tables = list(predictions = pred, metrics = ev$metrics, brier_differences = ev$diffs,
                     calibration = ev$calibration),
       key_numbers = list(metrics = ev$metrics, diffs = ev$diffs, n = length(ids)))
}

# Stratified k-fold cross-validation over every classifiable patient.
h_cross_validation <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$H
  plan <- ctx$holdout_plan
  d <- ctx$data()
  all <- analysis_set(d, "train")
  ids <- all$pheno[!is.na(phenotype), patient_id]
  y <- all$pheno[match(ids, patient_id), phenotype]
  set.seed(cfg$runtime$seed + 11L)
  fold <- integer(length(ids))
  fold[order(y, stats::runif(length(y)))] <- rep_len(seq_len(s$cv_folds), length(ids))
  start <- ctx$results$D4$fits[[plan$chosen_model]]$pop
  quiet <- function(...) invisible(NULL)

  pred <- data.table::rbindlist(lapply(seq_len(s$cv_folds), function(f) {
    te_ids <- ids[fold == f]
    tr_f <- subset_ids(all, setdiff(all$patients$patient_id, te_ids))
    pds <- build_pdats(tr_f, cfg)
    fit <- fit_population(pds, plan$chosen_model, cfg, init = start, log_fun = quiet)
    ctx$log(sprintf("  H: fold %d/%d - population (%s) on %d patients, predictions for %d",
                    f, s$cv_folds, plan$chosen_model, length(pds), length(te_ids)))
    h_predict(tr_f, subset_ids(all, te_ids), te_ids, fit, cfg, s, seed_offset = 100L * f)[, fold := f]
  }))
  pred <- pred[order(match(method, unique(method)), match(patient_id, ids))]
  pred[, observed := y[match(patient_id, ids)]]
  ev <- h_evaluate(pred, y, cfg, s)
  pred[, predicted := phenotype_levels[max.col(as.matrix(.SD), ties.method = "first")],
       .SDcols = c("p_quick", "p_second_wave", "p_persistent")]
  confusion <- pred[, .N, by = .(method, observed, predicted)][order(method, observed, predicted)]
  accuracy <- pred[, .(accuracy = round(mean(predicted == observed), 3), n = .N), by = method]

  list(pred = pred, metrics = ev$metrics, diffs = ev$diffs, calibration = ev$calibration, n_test = length(ids),
       observed = table(factor(y, phenotype_levels)), plan = plan, mode = "cv", confusion = confusion,
       accuracy = accuracy,
       tables = list(predictions = pred, metrics = ev$metrics, brier_differences = ev$diffs,
                     calibration = ev$calibration, confusion = confusion, accuracy = accuracy),
       key_numbers = list(metrics = ev$metrics, diffs = ev$diffs, n = length(ids), accuracy = accuracy))
}
