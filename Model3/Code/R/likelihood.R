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
# Observation model (log-normal, per lineage l = anc, plt, amc, alc)
#   log C_obs = log(b_l * B_l(t)) + e,  e ~ N(0, sigma_l^2);  ALC: b_alc * B_L(t) + CT(t)
#   b_l = individual baseline: log b_l = log m_l + d_l, d_l ~ N(0, s_l^2 / n_l + tau_l^2)
#   (m_l = set-point median; the set-point rows are not in the likelihood)
#   bound rows (left-censored): Phi((log U - log C)/sigma)
#   optional VAF (simulation designs only): log(2 VAF) ~ N(log f0, sd)
#   phenotype term (Model1): w log P'(c | psi), c = observed primary-rule class, P' from
#   phenotype_prob() on the model ANC (b_anc * B_N) at the label days; composite likelihood
#   (the day 0-30 ANC enter twice); absent at landmarks (pd$pheno = NULL)
# Individual vector psi = (random effects eta, baseline deviations b_<lineage>).
# =============================================================================

lik_vars <- c("anc", "plt", "amc", "alc")
core_vars <- c("anc", "plt", "amc")          # lineages shared with the Friberg reference MF
rel_col <- c(anc = "anc_rel", plt = "plt_rel", amc = "amc_rel", alc = "alc_rel")

# obs : data.table(day, var, value, cens, upper, setpoint) already filtered to usable rows
# use_upper: v1.2 handling of G-CSF rows as upper bounds (v1.3: FALSE, G-CSF is an input)
# pheno: list(class, days) of the observed phenotype (pheno_label) or NULL
build_pdat <- function(id, obs, base, design, drv, x, vaf = NULL, t_end = NULL, landmark = Inf,
                       base_n = NULL, ldh_c = 0, use_upper = FALSE, pheno = NULL) {
  obs <- obs[day <= landmark]
  if ("setpoint" %in% names(obs)) obs <- obs[setpoint == FALSE]
  if (!is.null(pheno) && (is.na(pheno$class) || !length(pheno$days))) pheno <- NULL
  times <- sort(unique(c(obs$day, t_end, pheno$days)))
  times <- times[times > design$t0]
  if (!is.null(pheno)) pheno$idx <- match(pheno$days, times)
  y <- lapply(lik_vars, function(v) {
    o <- obs[var == v & day > design$t0]
    if (!nrow(o)) return(NULL)
    list(idx = match(o$day, times), logy = log(pmax(o$value, 1e-8)),
         bound = o$cens | (use_upper & o$upper))
  })
  names(y) <- lik_vars
  list(id = id, times = times, design = design, drv = drv, base = base, base_n = base_n, ldh_c = ldh_c,
       x = x, y = y, vaf = vaf, pheno = pheno, seg = NULL)
}

# Observed phenotype of one patient and the days the rule used (all ANC rows, as
# phenotype_all_rules); variant: rule variant of a sensitivity analysis
pheno_label <- function(anc_obs, config, variant = NULL) {
  if (!nrow(anc_obs)) return(NULL)
  phenotype_series(anc_obs$day, anc_obs$value, anc_obs$cens, phenotype_rule(config, variant), return_days = TRUE)
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
# The phenotype term is attached for full-data fits (w > 0) and never at a landmark.
build_pdats <- function(d, config, landmark = Inf, ids = NULL, drivers = NULL) {
  if (is.null(ids)) ids <- d$patients$patient_id
  if (is.null(drivers)) drivers <- d$drivers
  ob <- d$obs[keep == TRUE]
  ph <- config$fit$pheno
  use_ph <- !is.finite(landmark) && !is.null(ph) && ph$w > 0
  bv <- config$model$baseline$vars
  res <- lapply(ids, function(pid) {
    prow <- d$patients[patient_id == pid]
    o <- ob[patient_id == pid]
    base <- lapply(stats::setNames(bv, bv), function(v) prow[[paste0(v, "_base")]])
    base_n <- vapply(bv, function(v) { n <- prow[[paste0("n_base_", v)]]; if (is.null(n) || is.na(n)) 1 else n },
                     numeric(1))
    vaf <- if (!is.null(d$vaf) && pid %in% d$vaf$patient_id) d$vaf[patient_id == pid] else NULL
    t_end <- if (is.finite(landmark)) max(landmark, 0) else max(o$day, 30)
    # landmark forecasts: G-CSF windows only up to the landmark (no future G-CSF)
    des <- patient_design(prow, d$gcsf, until = if (is.finite(landmark)) landmark else Inf)
    lab <- if (use_ph) pheno_label(d$obs[patient_id == pid & var == "anc" & !is.na(value)], config, ph$variant) else NULL
    build_pdat(pid, o, base, des, drivers[[pid]], patient_covariates(prow, config),
               vaf = vaf, t_end = t_end, landmark = landmark, base_n = base_n,
               ldh_c = if ("ldh_c" %in% names(prow)) prow$ldh_c else 0,
               use_upper = isTRUE(config$fit$use_gcsf_upper), pheno = lab)
  })
  stats::setNames(res, ids)
}

# Attach segment matrices for a given structural vector (they depend on bridging
# and renal-scaling parameters).
prep_segments <- function(pdats, struct, config) {
  lapply(pdats, function(pd) { pd$seg <- make_segments(pd$design, struct, max(pd$times), config); pd })
}

# ---- parameter template -------------------------------------------------------------------
# The antigen load a0 = exp(beta_a * ldh_c) is patient-specific: beta_a rides along as an
# attribute and fill_par sets a0 from pd$ldh_c.
par_template <- function(struct, spec) {
  ind <- if (spec$family == "mf") c(slope = 1, mtt_scale = 1) else
    c(r_s = 0.1, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0, rho = 0, a0 = 1)
  out <- build_par(ind, struct, spec)
  attr(out, "beta_a") <- if ("beta_a" %in% names(struct)) struct[["beta_a"]] else 0
  out
}

fill_par <- function(tmpl, eta, spec, pd = NULL) {
  ind <- eta_to_indiv(eta, spec)
  tmpl[names(ind)] <- ind
  if (spec$family != "mf") tmpl[["a0"]] <- exp(attr(tmpl, "beta_a") * (if (is.null(pd$ldh_c)) 0 else pd$ldh_c))
  tmpl
}

simulate_pdat <- function(eta, pd, tmpl, spec, config) {
  sim_model(fill_par(tmpl, eta, spec, pd), pd$seg, pd$drv, pd$times, spec, config)
}

# Natural-scale baselines of a patient from psi (missing deviations = 0)
pdat_baselines <- function(psi, pd) {
  vapply(names(pd$base), function(v) {
    k <- paste0("b_", v)
    pd$base[[v]] * exp(if (k %in% names(psi)) psi[[k]] else 0)
  }, numeric(1))
}

# Model prediction on the observed scale for lineage v at the rows of pd$y[[v]]
predict_obs <- function(m, pd, v, b) {
  yv <- pd$y[[v]]
  pred <- m[yv$idx, rel_col[[v]]] * b[[v]]
  if (v == "alc" && "cart" %in% colnames(m)) pred <- pred + m[yv$idx, "cart"]
  pred
}

# ---- log-likelihood of the data given a simulation ------------------------------------------
# vars: lineages to include (MF has no ALC column and skips it)
data_loglik <- function(m, pd, sigma, b = unlist(pd$base), vars = lik_vars) {
  ll <- 0
  for (v in vars) {
    yv <- pd$y[[v]]
    if (is.null(yv) || !rel_col[[v]] %in% colnames(m) || is.null(sigma[v]) || is.na(sigma[v])) next
    lp <- log(pmax(predict_obs(m, pd, v, b), 1e-8))
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

# w log P'(observed class | psi) from the model ANC at the label days (0 without a label)
pheno_loglik <- function(m, pd, b, config) {
  ph <- pd$pheno
  w <- config$fit$pheno$w
  if (is.null(ph) || is.null(w) || w <= 0) return(0)
  p <- phenotype_prob(pd$times[ph$idx], m[ph$idx, "anc_rel"] * b[["anc"]],
                      phenotype_rule(config, config$fit$pheno$variant), config$fit$pheno$s, config$fit$pheno$eps)
  if (is.null(p)) return(0)
  w * log(p[[ph$class]])
}

vaf_loglik <- function(eta, pd) {
  v <- pd$vaf
  if (is.null(v) || !"f0" %in% names(eta)) return(0)
  lf <- stats::plogis(eta[["f0"]], log.p = TRUE)
  if (v$status == "positive") return(stats::dnorm(log(v$f_measured), lf, v$sd, log = TRUE))
  stats::pnorm((log(2 * v$lod) - lf) / v$sd, log.p = TRUE)
}

# psi: named individual vector (random effects + baseline deviations); prior: list(mean, sd)
# over the same names. vars: lineages in the data term (core_vars for comparisons with MF).
neg_log_post <- function(eta, pd, tmpl, spec, prior, sigma, config, vars = lik_vars) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  if (anyNA(m[, 2])) return(1e10)
  b <- pdat_baselines(eta, pd)
  lp <- data_loglik(m, pd, sigma, b, vars) + pheno_loglik(m, pd, b, config) +
    sum(stats::dnorm(eta, prior$mean[names(eta)], prior$sd[names(eta)], log = TRUE)) + vaf_loglik(eta, pd)
  if (!is.finite(lp)) return(1e10)
  -lp
}

neg_log_lik_data <- function(eta, pd, tmpl, spec, sigma, config) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  if (anyNA(m[, 2])) return(1e10)
  b <- pdat_baselines(eta, pd)
  -data_loglik(m, pd, sigma, b) - pheno_loglik(m, pd, b, config)
}

residuals_at <- function(eta, pd, tmpl, spec, config) {
  m <- simulate_pdat(eta, pd, tmpl, spec, config)
  b <- pdat_baselines(eta, pd)
  out <- lapply(lik_vars, function(v) {
    yv <- pd$y[[v]]
    if (is.null(yv) || anyNA(m[, 2]) || !rel_col[[v]] %in% colnames(m)) return(numeric(0))
    (yv$logy - log(pmax(predict_obs(m, pd, v, b), 1e-8)))[!yv$bound]
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

# Prior of the individual vector: random effects N(mu + beta x, omega^2); baseline
# deviations N(0, s_l^2 / n_l + tau_l^2) (only when pop carries baselines, v1.3)
patient_prior <- function(pd, pop, re) {
  mean <- vapply(re, function(k) {
    b <- pop$beta[[k]]
    x <- pd$x[[k]]
    pop$mu[[k]] + if (length(b)) sum(b * x[names(b)]) else 0
  }, numeric(1))
  sd <- stats::setNames(unname(pop$omega[re]), re)
  if (length(pop$tau)) {
    v <- names(pop$tau)
    sb <- sqrt(pop$base_s[v]^2 / pmax(pd$base_n[v], 1) + pop$tau[v]^2)
    mean <- c(mean, stats::setNames(rep(0, length(v)), paste0("b_", v)))
    sd <- c(sd, stats::setNames(unname(sb), paste0("b_", v)))
  }
  list(mean = mean, sd = sd)
}

# Extra start for the clone fraction (avoid the f0 ~ 0 basin).
clone_starts <- function(eta0, spec) {
  if (!"f0" %in% names(eta0)) return(NULL)
  e <- eta0
  e[["f0"]] <- stats::qlogis(0.3)
  list(e)
}
