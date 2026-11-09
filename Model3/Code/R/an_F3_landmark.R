# =============================================================================
# an_F3_landmark.R
# Purpose : Step F3 - landmark prediction (T4), out of sample through the V1
#           folds: each patient is predicted with the population of the refit of
#           the final model without its fold (no phenotype term, no start from
#           the full-data MAP; comparators trained on the same folds). Without
#           the fold fits the full-data population is used and the rows are
#           flagged oos = FALSE (approximate, Aim4 v1.3). For each patient and landmark L:
#           individual posterior from data up to L, the CAR T-phase
#           driver after L drawn from the ferritin-kernel population
#           conditioned on the patient's ferritin up to L, and the phenotype
#           of [observed ANC up to L] + [model after L]. Compared with
#           CAR-HEMATOTOX (at LD) and a logistic model on counts at L.
# Inputs  : ctx (config, data, results F1).
# Outputs : list(pred, metrics, tables)
# =============================================================================

# CAR-HEMATOTOX (Rejeski 2021) at lymphodepletion; CRP in mg/dL
car_hematotox <- function(plt, anc, hgb, crp_mgdl, ferritin) {
  (ifelse(plt < 75, 2, ifelse(plt <= 175, 1, 0))) + as.numeric(anc < 1.2) + as.numeric(hgb < 9) +
    as.numeric(crp_mgdl >= 3) + (ifelse(ferritin > 2000, 2, ifelse(ferritin >= 650, 1, 0)))
}

auroc <- function(score, y) {
  ok <- !is.na(score) & !is.na(y)
  s <- score[ok]; y <- y[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  r <- rank(s)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

auroc_ci <- function(score, y, n_boot = 1000, seed = 1L) {
  set.seed(seed)
  b <- replicate(n_boot, { i <- sample(length(y), replace = TRUE); auroc(score[i], y[i]) })
  c(stats::quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

# Kernel draws conditioned on the ferritin observed up to L (importance resampling).
kernel_draws <- function(kpop, fer_day, fer_logratio, L, n, sd_obs = 0.3) {
  z <- MASS::mvrnorm(n * 4, kpop$mu, kpop$cov)
  A <- exp(z[, 1]); tp <- exp(z[, 2])
  sel <- fer_day > 0 & fer_day <= L
  w <- if (any(sel)) vapply(seq_along(A), function(j)
    sum(stats::dnorm(fer_logratio[sel], kernel_fun(fer_day[sel], A[j], tp[j], kpop$h), sd_obs, log = TRUE)),
    numeric(1)) else rep(0, length(A))
  w <- exp(w - max(w))
  ix <- sample(length(A), n, replace = TRUE, prob = w)
  cbind(A = A[ix], tp = tp[ix])
}

# driver: observed smooth up to L, kernel draw after L
landmark_driver <- function(d_obs_fer, base, t0, L, A, tp, config) {
  obs <- d_obs_fer[day <= L]
  dr <- if (L >= 3 && nrow(obs) >= 4) driver_from_ferritin(obs$day, obs$ferritin, base, t0, config) else zero_driver(t0, config)
  tt <- dr$t0 + (seq_along(dr$v) - 1) * dr$dt
  k <- kernel_fun(tt, A, tp, config$model$kernel$h)
  dr$v[tt > L] <- k[tt > L]
  dr
}

# Per landmark L: individual posterior from the data up to L (MAP by L-BFGS-B in the search
# box from the prior mean and four prior draws, then one Metropolis chain), CAR T dynamics
# predicted, no G-CSF after L, ferritin after L from the kernel population of the training
# patients conditioned on the patient's ferritin up to L.
forecast_class13 <- function(psi, pd, fit, cfg, drv, obs_upto) {
  spec <- fit$spec
  st <- fit$pop$struct
  ind <- eta_to_indiv(psi, spec)
  ind[["a0"]] <- patient_a0(st, pd$ldh_c)
  b <- pdat_baselines(psi, pd)
  days <- 0:30
  m <- sim_model(build_par(ind, st, spec), make_segments(pd$design, st, 30, cfg), drv, days, spec, cfg)
  if (anyNA(m[, 2])) return(NA_character_)
  phenotype_with_observed(m[, "t"], m[, "anc_rel"] * b[["anc"]], obs_upto, cfg)
}

# starts: extra start vectors (approximate mode only: the full-data MAP)
f3_patient <- function(pid, d, fit, kpop, cfg, starts = list(), seed = 1L) {
  s <- cfg$steps$F3
  prow <- d$patients[patient_id == pid]
  fer <- d$inflam[patient_id == pid & !is.na(ferritin), .(day, ferritin)]
  lr <- log(fer$ferritin / prow$ferritin_ld)
  anc <- d$obs[patient_id == pid & var == "anc" & keep == TRUE]
  tmpl <- par_template(fit$pop$struct, fit$spec)
  out <- lapply(s$landmarks, function(L) {
    set.seed(seed + L)
    kd <- kernel_draws(kpop, fer$day, lr, L, s$n_kernel)
    drv_med <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, stats::median(kd[, 1]), stats::median(kd[, 2]), cfg)
    pd <- build_pdats(d, cfg, landmark = L, ids = pid, drivers = stats::setNames(list(drv_med), pid))[[1]]
    pd$seg <- make_segments(pd$design, fit$pop$struct, max(pd$times), cfg)
    prior <- patient_prior(pd, fit$pop, fit$re)
    nm <- names(prior$mean)
    box <- indiv_box(prior, cfg)
    fn <- function(x) neg_log_post(stats::setNames(x, nm), pd, tmpl, fit$spec, prior, fit$pop$sigma, cfg)
    best <- NULL
    set.seed(seed + 500L + L)
    draws0 <- lapply(1:4, function(j) prior$mean + prior$sd * stats::rnorm(length(nm)))
    for (s0 in c(lapply(starts, function(x) x[nm]), list(prior$mean), draws0)) {
      o <- tryCatch(stats::optim(pmin(pmax(s0, box$lo), box$hi), fn, method = "L-BFGS-B", lower = box$lo,
                                 upper = box$hi, control = list(maxit = 300)), error = function(e) NULL)
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    cv <- solve(pd_floor(stats::optimHess(best$par, fn), 1 / max(prior$sd^2)))
    dimnames(cv) <- list(nm, nm)
    mc <- mcmc_patient13(pd, fit, cfg, stats::setNames(best$par, nm), cv,
                         mc = list(n_chains = 1L, n_iter = 600L, burn = 300L, thin = 1L, scale = 2.38,
                                   rhat_tol = Inf, extend = 0L, max_iter = 600L), seed = seed + 1000L + L)
    ix <- sample(nrow(mc$draws), s$n_draws, replace = TRUE)
    upto <- anc[day <= L & day >= 0]
    obs_upto <- list(day = upto$day, value = upto$value, cens = upto$cens, L = L)
    cl <- vapply(seq_len(s$n_draws), function(j) {
      jk <- (j - 1L) %% nrow(kd) + 1L   # n_draws may exceed n_kernel: reuse the kernel draws cyclically
      drv <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, kd[jk, 1], kd[jk, 2], cfg)
      tryCatch(forecast_class13(mc$draws[ix[j], ], pd, fit, cfg, drv, obs_upto), error = function(e) NA_character_)
    }, character(1))
    n_failed <- sum(is.na(cl))
    cl <- cl[!is.na(cl)]
    data.table::data.table(patient_id = pid, L = L, n_used = length(cl), n_failed = n_failed,
                           p_quick = mean(cl == "quick"), p_second_wave = mean(cl == "second_wave"),
                           p_delayed = mean(cl == "delayed"))
  })
  data.table::rbindlist(out)
}

# Fold fits of the final model (hpc_runs/<final>_fold<k>/fit.rds); NULL when any is missing
f3_fold_fits <- function(cfg, final, runs = file.path(cfg$paths$results, "hpc_runs")) {
  k <- cfg$validation$V1$k
  f <- file.path(runs, paste0(final, "_fold", seq_len(k)), "fit.rds")
  if (!all(file.exists(f))) return(NULL)
  lapply(f, readRDS)
}

# Landmark predictions of every patient: with fold fits (list of k fits) each patient uses the
# fit without its fold and the kernel population of the training patients; without them the
# full-data fit (approximate, oos = FALSE)
landmark_predictions <- function(d, fit_full, fold_fits, folds, cfg) {
  ids <- d$patients$patient_id
  oos <- !is.null(fold_fits)
  psi <- fit_psi13(fit_full)
  data.table::rbindlist(par_lapply(seq_along(ids), function(i) {
    pid <- ids[i]
    k <- folds[patient_id == pid]$fold
    fit <- if (oos) fold_fits[[k]] else fit_full
    train <- if (oos) folds[fold != k]$patient_id else setdiff(ids, pid)
    kpop <- fit_kernel_population(d$drivers[train], cfg)
    # one patient's failure must not discard the others' results: keep it as NA rows
    out <- tryCatch(f3_patient(pid, d, fit, kpop, cfg, starts = if (oos) list() else list(psi[pid, ]),
                               seed = cfg$runtime$seed + i), error = function(e)
      data.table::data.table(patient_id = pid, L = cfg$steps$F3$landmarks, n_used = 0L, n_failed = NA_integer_,
                             p_quick = NA_real_, p_second_wave = NA_real_, p_delayed = NA_real_,
                             error = conditionMessage(e)))
    out[, `:=`(fold = k, oos = oos)]
    out
  }, cfg$runtime$cores), fill = TRUE)
}

an_F3_landmark <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  final <- f1$selected$final
  fit <- f1$fits[[final]]
  folds <- data.table::fread(file.path(cfg$paths$results, "prereg", "V1_folds.csv"))
  fold_fits <- f3_fold_fits(cfg, final)
  if (is.null(fold_fits)) ctx$log("  F3: fold fits of ", final, " missing - APPROXIMATE predictions (full-data population)")
  pred <- landmark_predictions(d, fit, fold_fits, folds, cfg)
  if (!"error" %in% names(pred)) pred[, error := NA_character_]
  pred[d$phenotype, observed := i.primary, on = "patient_id"]
  pts <- d$patients
  pred[pts, gcsf := i.gcsf_after_d14, on = "patient_id"]
  # comparators
  hx <- pts[, .(patient_id, hematotox = car_hematotox(plt_ld, anc_ld, hgb_ld, crp_ld, ferritin_ld))]
  pred[hx, hematotox := i.hematotox, on = "patient_id"]
  lg <- f3_logistic(d, cfg, folds)
  pred[lg, p_nonquick_logit := i.p, on = .(patient_id, L)]
  pred[, y_nonquick := as.numeric(observed != "quick")]
  metrics <- pred[, {
    out <- list()
    for (grp in c("all", "no_gcsf")) {
      sub <- if (grp == "all") .SD else .SD[gcsf == 0]
      ci <- auroc_ci(1 - sub$p_quick, sub$y_nonquick)
      Y <- vapply(cfg$phenotype$classes, function(k) as.numeric(sub$observed == k), numeric(nrow(sub)))
      P <- as.matrix(sub[, .(p_quick, p_second_wave, p_delayed)])
      out[[grp]] <- data.table::data.table(group = grp, n = nrow(sub),
        auroc_model = auroc(1 - sub$p_quick, sub$y_nonquick), auroc_lo = ci[1], auroc_hi = ci[2],
        auroc_hematotox = auroc(sub$hematotox, sub$y_nonquick),
        auroc_logistic = auroc(sub$p_nonquick_logit, sub$y_nonquick),
        brier_3class = mean(rowSums((P - Y)^2), na.rm = TRUE))
    }
    data.table::rbindlist(out)
  }, by = L]
  first_day <- metrics[group == "all" & auroc_lo > auroc_hematotox, min(L)]
  ctx$log("  F3: AUROC by landmark ", paste(sprintf("L%d %.2f", metrics[group == "all"]$L,
                                                    metrics[group == "all"]$auroc_model), collapse = ", "),
          "; CAR-HEMATOTOX ", round(metrics[group == "all"]$auroc_hematotox[1], 2),
          "; failed draws ", sum(pred$n_failed, na.rm = TRUE), "/", sum(pred$n_used + pred$n_failed, na.rm = TRUE),
          "; failed patients ", data.table::uniqueN(pred[!is.na(error)]$patient_id))
  list(pred = pred, metrics = metrics, earliest_day = if (is.finite(first_day)) first_day else NA,
       out_of_sample = !is.null(fold_fits), tables = list(predictions = pred, metrics = metrics))
}

# Ridge logistic model on log ANC, log PLT at L and their slopes, trained on the other V1 folds
# (folds = NULL: leave-one-out).
f3_logistic <- function(d, cfg, folds = NULL) {
  ob <- d$obs[var %in% c("anc", "plt") & keep == TRUE]
  feat <- function(L) {
    data.table::rbindlist(lapply(d$patients$patient_id, function(pid) {
      f <- function(v) {
        o <- ob[patient_id == pid & var == v & day <= L & day >= L - 6][order(day)]
        if (!nrow(o)) return(c(NA, NA))
        lv <- log(pmax(o$value, 0.01))
        c(utils::tail(lv, 1), if (nrow(o) >= 2) (utils::tail(lv, 1) - lv[1]) / max(1, utils::tail(o$day, 1) - o$day[1]) else 0)
      }
      a <- f("anc"); p <- f("plt")
      data.table::data.table(patient_id = pid, L = L, lanc = a[1], sanc = a[2], lplt = p[1], splt = p[2])
    }))
  }
  X <- data.table::rbindlist(lapply(cfg$steps$F3$landmarks, feat))
  X[d$phenotype, y := as.numeric(i.primary != "quick"), on = "patient_id"]
  X[, p := NA_real_]
  if (!is.null(folds)) X[folds, fold := i.fold, on = "patient_id"]
  for (L0 in unique(X$L)) {
    # cv.glmnet assigns its inner folds at random: fixed seed per landmark (reproducible comparator; 2026-11-05)
    set.seed(cfg$runtime$seed + 7919L + L0)
    sub <- X[L == L0]
    for (i in seq_len(nrow(sub))) {
      out <- if (is.null(folds)) -i else which(sub$fold != sub$fold[i])
      tr <- sub[out][stats::complete.cases(sub[out, .(lanc, sanc, lplt, splt)])]
      if (nrow(tr) < 10 || is.na(sub$lanc[i])) next
      m <- tryCatch(glmnet::cv.glmnet(as.matrix(tr[, .(lanc, sanc, lplt, splt)]), tr$y, family = "binomial",
                                      alpha = 0, nfolds = 5), error = function(e) NULL)
      if (is.null(m)) next
      nx <- as.matrix(sub[i, .(lanc, sanc, lplt, splt)])
      nx[is.na(nx)] <- 0
      X[L == L0 & patient_id == sub$patient_id[i], p := as.numeric(stats::predict(m, nx, s = "lambda.min", type = "response"))]
    }
  }
  X[, .(patient_id, L, p)]
}
