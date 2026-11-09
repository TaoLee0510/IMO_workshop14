# =============================================================================
# likelihood.R
# Purpose : Observation model and individual-level objective functions.
#           A "pdat" bundles everything the likelihood needs for one patient;
#           the same constructor serves real data, landmark data (hold-out)
#           and virtual patients (identifiability, virtual cohort).
# Inputs  : standardised data (io_read.R) or simulated observations.
# Outputs : pdat objects; negative log posterior; individual MAP fits with a
#           Laplace approximation of the marginal likelihood.
#
# Observation model (log-normal)
#   ANC  = anc_base * N(t)                 censored below the detection limit
#   PLT  = plt_base * PL(t)
#   CRP  = crp_base + crp_scale * Is(t)    (after tocilizumab: censored upstream)
#   Ferr = fer_base + fer_scale * Is(t)
#   VAF  : log(2*VAF) ~ N(log f0, sd_timepoint); negative -> left-censored
# =============================================================================

obs_vars <- c("anc", "plt", "crp", "ferritin")

# ---- pdat construction ------------------------------------------------------------------------
# obs : data.table(day, var, value, cens) already filtered to usable rows
# base: list(anc, plt, crp, ferritin)
# x   : named list of covariate vectors per random effect (rho, g, f0)
# vaf : NULL or list(status, f_measured, sd, lod)
# vaf_post: NULL or list(day, f_measured, sd)
build_pdat <- function(id, obs, base, seg, x, vaf = NULL, vaf_post = NULL, extra_times = NULL) {
  times <- sort(unique(c(obs$day, extra_times, if (!is.null(vaf_post)) vaf_post$day)))
  y <- lapply(obs_vars, function(v) {
    o <- obs[var == v]
    if (!nrow(o)) return(NULL)
    list(idx = match(o$day, times), logy = log(pmax(o$value, 1e-8)), cens = o$cens)
  })
  names(y) <- obs_vars
  list(id = id, times = times, seg = seg, y = y, base = base, x = x,
       vaf = vaf, vaf_post = vaf_post,
       post_idx = if (!is.null(vaf_post)) match(vaf_post$day, times) else NULL)
}

# Covariate vectors for one patient (names follow config$model$covariates;
# covariates absent from the data, e.g. red-cell indices, are skipped).
patient_covariates <- function(prow, config, cov_f0 = TRUE) {
  cv <- config$model$covariates
  out <- lapply(names(cv), function(k) {
    v <- intersect(cv[[k]], names(prow))
    if (k == "f0" && !cov_f0) v <- character(0)
    stats::setNames(vapply(v, function(c) as.numeric(prow[[c]]), numeric(1)), v)
  })
  stats::setNames(out, names(cv))
}

# pdat for every patient in a standardised data object.
# landmark: if finite, only observations and events up to that day are used
# config$fit$obs_max_day (if set) is the last day of the fitting window
build_pdats <- function(d, config, landmark = Inf, use_vaf = TRUE, cov_f0 = TRUE) {
  t0 <- config$model$time$t0
  tmax <- config$model$time$tmax
  ids <- d$patients$patient_id
  obs_all <- d$obs[keep == TRUE & day >= t0 + 1 & day <= min(landmark, tmax, config$fit$obs_max_day %||% Inf)]
  ev_all <- d$events[start_day <= landmark]
  res <- lapply(ids, function(pid) {
    prow <- d$patients[patient_id == pid]
    o <- obs_all[patient_id == pid]
    seg <- make_segments(t0, tmax, ld_scale = prow$ld_scale, events = ev_all[patient_id == pid])
    base <- list(anc = prow$anc_base, plt = prow$plt_base, crp = prow$crp_base, ferritin = prow$fer_base)
    vaf <- NULL
    if (use_vaf && !is.null(d$vaf) && pid %in% d$vaf$patient_id) {
      v <- d$vaf[patient_id == pid]
      vaf <- list(status = v$status, f_measured = v$f_measured, sd = v$timepoint_sd, lod = v$lod_vaf)
    }
    build_pdat(pid, o, base, seg, patient_covariates(prow, config, cov_f0), vaf = vaf,
               extra_times = min(tmax, config$fit$obs_max_day %||% Inf))
  })
  stats::setNames(res, ids)
}

# ---- parameter template -----------------------------------------------------------------------
# Full parameter vector with structural values and nesting switches applied;
# individual parameters (positions 1-3 = rho, g, f0) are filled per call.
par_template <- function(struct, model) {
  build_par(c(rho = 1, g = 0.1, f0 = 0), struct, model)
}

fill_par <- function(tmpl, eta, model) {
  tmpl[1] <- exp(eta[1])
  tmpl[2] <- exp(eta[2])
  if (model == "M2") tmpl[3] <- stats::plogis(eta[3])
  tmpl
}

# ---- log-likelihood of the data given a simulation ----------------------------------------------
data_loglik <- function(m, pd, sigma, obsm) {
  ll <- 0
  pred <- list(
    anc = m[, "anc_rel"] * pd$base$anc,
    plt = m[, "plt_rel"] * pd$base$plt,
    crp = pd$base$crp + obsm$crp_scale * m[, "is"],
    ferritin = pd$base$ferritin + obsm$fer_scale * m[, "is"])
  for (v in obs_vars) {
    yv <- pd$y[[v]]
    if (is.null(yv)) next
    lp <- log(pmax(pred[[v]][yv$idx], 1e-8))
    s <- sigma[[v]]
    r <- (yv$logy - lp) / s
    if (any(yv$cens)) {
      ll <- ll + sum(stats::pnorm(r[yv$cens], log.p = TRUE))
      r <- r[!yv$cens]
    }
    ll <- ll + sum(-0.5 * r^2 - log(s) - 0.5 * log(2 * pi))
  }
  ll
}

vaf_loglik <- function(eta_f0, pd) {
  v <- pd$vaf
  if (is.null(v)) return(0)
  lf <- stats::plogis(eta_f0, log.p = TRUE)
  if (v$status == "positive")
    return(stats::dnorm(log(v$f_measured), lf, v$sd, log = TRUE))
  stats::pnorm((log(2 * v$lod) - lf) / v$sd, log.p = TRUE)
}

vaf_post_loglik <- function(m, pd) {
  vp <- pd$vaf_post
  if (is.null(vp)) return(0)
  f_pred <- max(m[pd$post_idx, "f_clone"], 1e-8)
  stats::dnorm(log(vp$f_measured), log(f_pred), vp$sd, log = TRUE)
}

# ---- negative log posterior for one patient -----------------------------------------------------
# prior: list(mean, sd) on the eta scale for the model's random effects
neg_log_post <- function(eta, pd, tmpl, model, prior, sigma, config, use_vaf = TRUE) {
  par <- fill_par(tmpl, eta, model)
  m <- marrow_sim(par, pd$seg, pd$times, config$model$solver$rtol,
                  config$model$solver$atol, config$model$solver$hmax)
  if (anyNA(m[, 2])) return(1e10)
  lp <- data_loglik(m, pd, sigma, config$model$obs) +
    sum(stats::dnorm(eta, prior$mean, prior$sd, log = TRUE))
  if (model == "M2") {
    if (use_vaf) lp <- lp + vaf_loglik(eta[3], pd)
    lp <- lp + vaf_post_loglik(m, pd)
  }
  if (!is.finite(lp)) return(1e10)
  -lp
}

# Data-only negative log-likelihood (used for structural updates at fixed eta).
neg_log_lik_data <- function(eta, pd, tmpl, model, sigma, config) {
  par <- fill_par(tmpl, eta, model)
  m <- marrow_sim(par, pd$seg, pd$times, config$model$solver$rtol,
                  config$model$solver$atol, config$model$solver$hmax)
  if (anyNA(m[, 2])) return(1e10)
  -data_loglik(m, pd, sigma, config$model$obs)
}

# Residuals at given eta (for residual-SD updates).
residuals_at <- function(eta, pd, tmpl, model, config) {
  par <- fill_par(tmpl, eta, model)
  m <- marrow_sim(par, pd$seg, pd$times, config$model$solver$rtol,
                  config$model$solver$atol, config$model$solver$hmax)
  pred <- list(anc = m[, "anc_rel"] * pd$base$anc, plt = m[, "plt_rel"] * pd$base$plt,
               crp = pd$base$crp + config$model$obs$crp_scale * m[, "is"],
               ferritin = pd$base$ferritin + config$model$obs$fer_scale * m[, "is"])
  out <- lapply(obs_vars, function(v) {
    yv <- pd$y[[v]]
    if (is.null(yv) || anyNA(m[, 2])) return(numeric(0))
    (yv$logy - log(pmax(pred[[v]][yv$idx], 1e-8)))[!yv$cens]
  })
  stats::setNames(out, obs_vars)
}

# ---- individual MAP + Laplace -----------------------------------------------------------------------
# Positive-definite, well-conditioned version of a numerical Hessian. Eigenvalues
# are bounded to [1e-4, 1e6] (the objective returns 1e10 outside the feasible
# region, which can make finite-difference curvature explode).
make_pd_matrix <- function(H) {
  if (any(!is.finite(H))) return(diag(1, nrow(H)))
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  ev$values <- pmin(pmax(ev$values, 1e-4), 1e6)
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

fit_individual <- function(pd, tmpl, model, prior, sigma, config, eta0, use_vaf = TRUE,
                           extra_starts = NULL) {
  fn <- function(e) neg_log_post(e, pd, tmpl, model, prior, sigma, config, use_vaf)
  starts <- c(list(eta0), extra_starts)
  best <- NULL
  for (s in starts) {
    o <- tryCatch(stats::optim(s, fn, method = "BFGS",
                               control = list(maxit = config$fit$map_maxit, reltol = 1e-8)),
                  error = function(e) NULL)
    if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
  }
  nm <- names(eta0)
  k <- length(eta0)
  if (is.null(best)) {
    return(list(eta = eta0, cov = diag(1, k, dimnames = list(nm, nm)), nll = NA_real_,
                logL = NA_real_, conv = 99L))
  }
  H <- make_pd_matrix(stats::optimHess(best$par, fn))
  logdet <- as.numeric(determinant(H, logarithm = TRUE)$modulus)
  covm <- solve(H)
  dimnames(covm) <- list(nm, nm)
  list(eta = stats::setNames(best$par, nm), cov = covm, nll = best$value,
       logL = -best$value + 0.5 * k * log(2 * pi) - 0.5 * logdet, conv = as.integer(best$convergence))
}

# Prior for one patient from population parameters and covariates.
patient_prior <- function(pd, pop, re) {
  mean <- vapply(re, function(k) {
    b <- pop$beta[[k]]
    x <- pd$x[[k]]
    pop$mu[[k]] + if (length(b)) sum(b * x[names(b)]) else 0
  }, numeric(1))
  list(mean = mean, sd = unname(pop$omega[re]))
}

# Extra starts for the clone fraction (avoid the f0 ~ 0 basin when the data say
# otherwise): a fixed large-clone start, plus the measured clone fraction when a
# positive VAF is available.
clone_starts <- function(eta0, model, pd = NULL) {
  if (model != "M2") return(NULL)
  e <- eta0
  e[3] <- stats::qlogis(0.3)
  out <- list(e)
  v <- if (is.null(pd)) NULL else pd$vaf
  if (!is.null(v) && identical(v$status, "positive")) {
    e2 <- eta0
    e2[3] <- stats::qlogis(min(max(v$f_measured, 0.01), 0.95))
    out <- c(out, list(e2))
  }
  out
}
