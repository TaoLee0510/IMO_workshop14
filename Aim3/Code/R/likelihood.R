# =============================================================================
# likelihood.R
# Purpose : Observation model and individual-level objective functions. A
#           "pdat" = design object (model.R) + the patient's observations; the
#           same constructor serves real data, landmark-truncated data and
#           virtual patients.
# Inputs  : pdes (make_pdes), observation rows (io_read.R build_obs format:
#           day, var, value, ub, keep).
# Outputs : pdat objects; negative log posterior; individual MAP fits with a
#           Laplace approximation of the marginal likelihood.
#
# Observation model (log-normal), lineage l in {anc, plt, amc}
#   y = set_l * B_l(t) * exp(eps),  eps ~ N(0, sigma_l^2)
#   ub = TRUE (no differential, zero on the differential, ANC under G-CSF or
#   steroids): the latent value is <= y  ->  log Phi((log y - log pred) / sigma_l)
#   lb = TRUE (AMC under steroids): latent >= y  ->  log Phi(-(log y - log pred) / sigma_l)
# =============================================================================

obs_vars <- c("anc", "plt", "amc")
rel_col <- c(anc = "anc_rel", plt = "plt_rel", amc = "amc_rel")

build_pdat <- function(pdes, obs, prow = NULL, landmark = Inf) {
  o <- obs[keep == TRUE & day >= pdes$t0 & day <= landmark]
  times <- sort(unique(o$day))
  y <- lapply(obs_vars, function(v) {
    ov <- o[var == v]
    if (!nrow(ov)) return(NULL)
    list(idx = match(ov$day, times), logy = log(ov$value), ub = ov$ub, lb = ov$lb)
  })
  names(y) <- obs_vars
  c(pdes, list(times = times, y = y, prow = prow, n_obs = nrow(o)))
}

# All real patients.
build_pdats <- function(d, struct, config, landmark = Inf) {
  ids <- d$patients$patient_id
  out <- lapply(ids, function(pid) {
    prow <- d$patients[patient_id == pid]
    pdes <- make_pdes(prow, d$inflam[patient_id == pid], struct, config)
    build_pdat(pdes, d$obs[patient_id == pid], prow, landmark)
  })
  stats::setNames(out, ids)
}

# Structural parameters that change the dosing segments or driver.
refresh_pdats <- function(pdats, struct, config) {
  lapply(pdats, function(pd) if (is.null(pd$prow)) pd else refresh_pdes(pd, pd$prow, struct, config))
}

# ---- likelihood ------------------------------------------------------------------------------------
data_loglik <- function(m, pd, sigma) {
  ll <- 0
  for (v in obs_vars) {
    yv <- pd$y[[v]]
    if (is.null(yv)) next
    lp <- log(pmax(m[yv$idx, rel_col[[v]]] * pd$set[[v]], 1e-8))
    s <- sigma[[v]]
    r <- (yv$logy - lp) / s
    if (any(yv$ub)) ll <- ll + sum(stats::pnorm(r[yv$ub], log.p = TRUE))
    if (any(yv$lb)) ll <- ll + sum(stats::pnorm(-r[yv$lb], log.p = TRUE))
    r <- r[!yv$ub & !yv$lb]
    ll <- ll + sum(-0.5 * r^2 - log(s) - 0.5 * log(2 * pi))
  }
  ll
}

# Parameter template: structural values and switches; individual slots filled per call.
par_template <- function(struct, model) build_par(c(th_ap = 1), struct, model)

fill_par <- function(tmpl, eta) {
  tmpl[names(eta)] <- exp(eta)
  tmpl
}

sim_at_obs <- function(par, pd, config) sim_raw(par, pd, pd$times, config)

neg_log_post <- function(eta, pd, tmpl, prior, sigma, config) {
  m <- sim_at_obs(fill_par(tmpl, eta), pd, config)
  if (anyNA(m[, 2])) return(1e10)
  lp <- data_loglik(m, pd, sigma) + sum(stats::dnorm(eta, prior$mean, prior$sd, log = TRUE))
  if (!is.finite(lp)) return(1e10)
  -lp
}

neg_log_lik_data <- function(eta, pd, tmpl, sigma, config) {
  m <- sim_at_obs(fill_par(tmpl, eta), pd, config)
  if (anyNA(m[, 2])) return(1e10)
  -data_loglik(m, pd, sigma)
}

residuals_at <- function(eta, pd, tmpl, config) {
  m <- sim_at_obs(fill_par(tmpl, eta), pd, config)
  out <- lapply(obs_vars, function(v) {
    yv <- pd$y[[v]]
    if (is.null(yv) || anyNA(m[, 2])) return(numeric(0))
    (yv$logy - log(pmax(m[yv$idx, rel_col[[v]]] * pd$set[[v]], 1e-8)))[!yv$ub & !yv$lb]
  })
  stats::setNames(out, obs_vars)
}

# ---- individual MAP + Laplace -------------------------------------------------------------------------
make_pd_matrix <- function(H) {
  if (any(!is.finite(H))) return(diag(1, nrow(H)))
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  ev$values <- pmin(pmax(ev$values, 1e-4), 1e6)
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

fit_individual <- function(pd, tmpl, prior, sigma, config, eta0, extra_starts = NULL) {
  fn <- function(e) neg_log_post(stats::setNames(e, names(eta0)), pd, tmpl, prior, sigma, config)
  best <- NULL
  for (s in c(list(eta0), extra_starts)) {
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

patient_prior <- function(pd, pop, re) {
  mean <- vapply(re, function(k) {
    b <- pop$beta[[k]]
    pop$mu[[k]] + if (length(b)) sum(b * pd$x[[k]][names(b)]) else 0
  }, numeric(1))
  list(mean = mean, sd = unname(pop$omega[re]))
}

# Extra starts so that HSC injury terms are not stuck at the prior basin.
injury_starts <- function(eta0) {
  out <- list()
  for (k in intersect(c("th_as", "th_bs"), names(eta0))) {
    e <- eta0
    e[[k]] <- e[[k]] + 1.5
    out <- c(out, list(e))
  }
  out
}
