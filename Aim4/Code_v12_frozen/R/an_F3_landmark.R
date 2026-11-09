# =============================================================================
# an_F3_landmark.R
# Purpose : Step F3 - leave-one-out landmark prediction (T4). For each patient
#           and landmark L: population refitted without the patient (warm
#           start), individual posterior from data up to L, the CAR T-phase
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

f3_patient <- function(pid, d, fit_loo, kpop, cfg) {
  s <- cfg$steps$F3
  prow <- d$patients[patient_id == pid]
  fer <- d$inflam[patient_id == pid & !is.na(ferritin), .(day, ferritin)]
  lr <- log(fer$ferritin / prow$ferritin_ld)
  anc <- d$obs[patient_id == pid & var == "anc"]
  out <- lapply(s$landmarks, function(L) {
    kd <- kernel_draws(kpop, fer$day, lr, L, s$n_kernel)
    drv_med <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, stats::median(kd[, 1]), stats::median(kd[, 2]), cfg)
    pd <- build_pdats(d, cfg, landmark = L, ids = pid, drivers = stats::setNames(list(drv_med), pid))[[1]]
    pd$seg <- make_segments(pd$design, fit_loo$pop$struct, max(pd$times), cfg)
    mc <- mcmc_patient(pd, fit_loo, cfg, mc = list(n_chains = 1L, n_iter = max(600L, cfg$fit$mcmc$burn + 200L),
                                                     burn = cfg$fit$mcmc$burn, thin = 2L, scale = 2.38))
    ix <- sample(nrow(mc$draws), s$n_draws, replace = TRUE)
    upto <- anc[day <= L & day >= 0]
    obs_upto <- list(day = upto$day, value = upto$value, cens = upto$cens, L = L)
    res <- lapply(seq_len(s$n_draws), function(j) {
      jk <- (j - 1L) %% nrow(kd) + 1L   # n_draws may exceed n_kernel: reuse the kernel draws cyclically
      drv <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, kd[jk, 1], kd[jk, 2], cfg)
      pdj <- pd
      pdj$drv <- drv
      pdj$times <- sort(unique(c(pd$times, 0:30)))
      # a failed counterfactual integration returns NA times, which the observed+model
      # phenotype cannot classify: drop that draw (counted in n_failed)
      tryCatch(counterfactual_draw(mc$draws[ix[j], ], pdj, fit_loo, cfg, 0:30, obs_upto = obs_upto),
               error = function(e) NULL)
    })
    n_failed <- sum(vapply(res, is.null, logical(1)))
    res <- res[!vapply(res, is.null, logical(1))]
    tab <- data.table::as.data.table(do.call(rbind, res))
    tab[, when := when_class(actual, no_a, no_bS, no_abS)]
    data.table::data.table(patient_id = pid, L = L, n_used = nrow(tab), n_failed = n_failed,
                           p_quick = mean(tab$actual == "quick"), p_second_wave = mean(tab$actual == "second_wave"),
                           p_delayed = mean(tab$actual == "delayed"),
                           p_a = mean(tab$when %in% c("a")), p_b = mean(tab$when %in% c("b")),
                           p_ab = mean(tab$when %in% c("ab_joint", "ab_redundant")),
                           p_clone = mean(tab$no_clone != tab$actual & tab$actual != "quick"))
  })
  data.table::rbindlist(out)
}

an_F3_landmark <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  fit <- f1$fits[[f1$selected$final]]
  pdats <- build_pdats(d, cfg)
  ids <- names(pdats)
  set.seed(cfg$runtime$seed)
  inner <- max(1L, floor(cfg$runtime$cores / 2))
  outer <- max(1L, cfg$runtime$cores %/% inner)
  pred <- data.table::rbindlist(par_lapply(seq_along(ids), function(i) {
    pid <- ids[i]
    set.seed(cfg$runtime$seed + i)
    loo <- fit_population(pdats[-i], fit$model, cfg, init = fit$pop, eta_init = fit_eta(fit)[-i, , drop = FALSE],
                          cores = 1L, log_fun = function(...) NULL, max_iter = if (cfg$quick) 1L else 3L)
    kpop <- fit_kernel_population(d$drivers[-i], cfg)
    # one patient's failure must not discard the others' results: keep it as NA rows
    tryCatch(f3_patient(pid, d, loo, kpop, cfg), error = function(e)
      data.table::data.table(patient_id = pid, L = cfg$steps$F3$landmarks, n_used = 0L, n_failed = NA_integer_,
                             p_quick = NA_real_, p_second_wave = NA_real_, p_delayed = NA_real_, p_a = NA_real_,
                             p_b = NA_real_, p_ab = NA_real_, p_clone = NA_real_, error = conditionMessage(e)))
  }, cfg$runtime$cores), fill = TRUE)
  if (!"error" %in% names(pred)) pred[, error := NA_character_]
  pred[d$phenotype, observed := i.primary, on = "patient_id"]
  pts <- d$patients
  pred[pts, gcsf := i.gcsf_after_d14, on = "patient_id"]
  # comparators
  hx <- pts[, .(patient_id, hematotox = car_hematotox(plt_ld, anc_ld, hgb_ld, crp_ld, ferritin_ld))]
  pred[hx, hematotox := i.hematotox, on = "patient_id"]
  lg <- f3_logistic(d, cfg)
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
       tables = list(predictions = pred, metrics = metrics))
}

# Leave-one-out logistic model on log ANC, log PLT at L and their 3-day slopes.
f3_logistic <- function(d, cfg) {
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
  for (L0 in unique(X$L)) {
    sub <- X[L == L0]
    for (i in seq_len(nrow(sub))) {
      tr <- sub[-i][stats::complete.cases(sub[-i, .(lanc, sanc, lplt, splt)])]
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
