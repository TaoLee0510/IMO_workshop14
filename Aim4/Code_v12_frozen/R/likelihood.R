# =============================================================================
# likelihood.R
# Purpose : Observation model and individual-level objective functions. A
#           "pdat" bundles everything the likelihood needs for one patient; the
#           same constructor serves real data, landmark data and virtual
#           patients.
# Inputs  : standardised data (io_read.R) or simulated observations.
# Outputs : pdat objects; negative log posterior; individual MAP fits with a
#           Laplace approximation of the marginal likelihood.
#
# Observation model (log-normal, per lineage l = anc, plt, amc)
#   log C_obs = log(base_l * B_l(t)) + e,  e ~ N(0, sigma_l^2)
#   bound rows (left-censored, or G-CSF upper bounds): Phi((log U - log C)/sigma)
#   optional VAF (simulation designs only): log(2 VAF) ~ N(log f0, sd)
# =============================================================================

lik_vars <- c("anc", "plt", "amc")
rel_col <- c(anc = "anc_rel", plt = "plt_rel", amc = "amc_rel")

# obs : data.table(day, var, value, cens, upper) already filtered to usable rows
build_pdat <- function(id, obs, base, design, drv, x, vaf = NULL, t_end = NULL, landmark = Inf) {
  obs <- obs[day <= landmark]
  times <- sort(unique(c(obs$day, t_end)))
  times <- times[times > design$t0]
  y <- lapply(lik_vars, function(v) {
    o <- obs[var == v & day > design$t0]
    if (!nrow(o)) return(NULL)
    list(idx = match(o$day, times), logy = log(pmax(o$value, 1e-8)), bound = o$cens | o$upper)
  })
  names(y) <- lik_vars
  list(id = id, times = times, design = design, drv = drv, base = base, x = x, y = y, vaf = vaf,
       seg = NULL)
}

# Covariate vectors for one patient (names follow config$model$covariates).
patient_covariates <- function(prow, config) {
  cv <- config$model$covariates
  out <- lapply(names(cv), function(k) {
    v <- intersect(cv[[k]], names(prow))
    stats::setNames(vapply(v, function(c) as.numeric(prow[[c]]), numeric(1)), v)
  })
  stats::setNames(out, names(cv))
}

# pdat for every patient of a data object; landmark: only observations up to that day.
build_pdats <- function(d, config, landmark = Inf, ids = NULL, drivers = NULL) {
  if (is.null(ids)) ids <- d$patients$patient_id
  if (is.null(drivers)) drivers <- d$drivers
  ob <- d$obs[keep == TRUE]
  res <- lapply(ids, function(pid) {
    prow <- d$patients[patient_id == pid]
    o <- ob[patient_id == pid]
    base <- list(anc = prow$anc_base, plt = prow$plt_base, amc = prow$amc_base,
                 alc = if ("alc_base" %in% names(prow)) prow$alc_base else NULL)   # alc: plots only
    vaf <- if (!is.null(d$vaf) && pid %in% d$vaf$patient_id) d$vaf[patient_id == pid] else NULL
    t_end <- if (is.finite(landmark)) max(landmark, 0) else max(o$day, 30)
    build_pdat(pid, o, base, patient_design(prow), drivers[[pid]], patient_covariates(prow, config),
               vaf = vaf, t_end = t_end, landmark = landmark)
  })
  stats::setNames(res, ids)
}

# Attach segment matrices for a given structural vector (they depend on bridging
# and renal-scaling parameters).
prep_segments <- function(pdats, struct, config) {
  lapply(pdats, function(pd) { pd$seg <- make_segments(pd$design, struct, max(pd$times), config); pd })
}

# ---- parameter template -------------------------------------------------------------------
par_template <- function(struct, spec) {
  ind <- if (spec$family == "mf") c(slope = 1, mtt_scale = 1) else c(r_s = 0.1, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0)
  build_par(ind, struct, spec)
}

fill_par <- function(tmpl, eta, spec) {
  ind <- eta_to_indiv(eta, spec)
  tmpl[names(ind)] <- ind
  tmpl
}

simulate_pdat <- function(eta, pd, tmpl, spec, config) {
  sim_model(fill_par(tmpl, eta, spec), pd$seg, pd$drv, pd$times, spec, config)
}

# ---- log-likelihood of the data given a simulation ------------------------------------------
data_loglik <- function(m, pd, sigma) {
  ll <- 0
  for (v in lik_vars) {
    yv <- pd$y[[v]]
    if (is.null(yv)) next
    lp <- log(pmax(m[yv$idx, rel_col[[v]]] * pd$base[[v]], 1e-8))
    s <- sigma[[v]]
    r <- (yv$logy - lp) / s
    if (any(yv$bound)) {
      ll <- ll + sum(stats::pnorm(r[yv$bound], log.p = TRUE))
      r <- r[!yv$bound]
    }
    ll <- ll + sum(-0.5 * r^2 - log(s) - 0.5 * log(2 * pi))
  }
  ll
}

vaf_loglik <- function(eta, pd) {
  v <- pd$vaf
  if (is.null(v) || !"f0" %in% names(eta)) return(0)
  lf <- stats::plogis(eta[["f0"]], log.p = TRUE)
  if (v$status == "positive") return(stats::dnorm(log(v$f_measured), lf, v$sd, log = TRUE))
  stats::pnorm((log(2 * v$lod) - lf) / v$sd, log.p = TRUE)
}

neg_log_post <- function(eta, pd, tmpl, spec, prior, sigma, config) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  if (anyNA(m[, 2])) return(1e10)
  lp <- data_loglik(m, pd, sigma) + sum(stats::dnorm(eta, prior$mean, prior$sd, log = TRUE)) +
    vaf_loglik(eta, pd)
  if (!is.finite(lp)) return(1e10)
  -lp
}

neg_log_lik_data <- function(eta, pd, tmpl, spec, sigma, config) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  if (anyNA(m[, 2])) return(1e10)
  -data_loglik(m, pd, sigma)
}

residuals_at <- function(eta, pd, tmpl, spec, config) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  out <- lapply(lik_vars, function(v) {
    yv <- pd$y[[v]]
    if (is.null(yv) || anyNA(m[, 2])) return(numeric(0))
    (yv$logy - log(pmax(m[yv$idx, rel_col[[v]]] * pd$base[[v]], 1e-8)))[!yv$bound]
  })
  stats::setNames(out, lik_vars)
}

# ---- individual MAP + Laplace -------------------------------------------------------------------
make_pd_matrix <- function(H) {
  if (any(!is.finite(H))) return(diag(1, nrow(H)))
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  ev$values <- pmin(pmax(ev$values, 1e-4), 1e6)
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

fit_individual <- function(pd, tmpl, spec, prior, sigma, config, eta0, extra_starts = NULL) {
  fn <- function(e) neg_log_post(stats::setNames(e, names(eta0)), pd, tmpl, spec, prior, sigma, config)
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

patient_prior <- function(pd, pop, re) {
  mean <- vapply(re, function(k) {
    b <- pop$beta[[k]]
    x <- pd$x[[k]]
    pop$mu[[k]] + if (length(b)) sum(b * x[names(b)]) else 0
  }, numeric(1))
  list(mean = mean, sd = unname(pop$omega[re]))
}

# Extra start for the clone fraction (avoid the f0 ~ 0 basin).
clone_starts <- function(eta0, spec) {
  if (!"f0" %in% names(eta0)) return(NULL)
  e <- eta0
  e[["f0"]] <- stats::qlogis(0.3)
  list(e)
}
