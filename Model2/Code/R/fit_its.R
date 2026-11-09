# =============================================================================
# fit_its.R
# Purpose : Population fit by the iterative two-stage method (ITS) with a
#           Laplace approximation of each patient's marginal likelihood
#           (copied from Aim 2 and adapted to model specs, the parent model's
#           random effects and the Friberg reference model).
# Inputs  : list of pdats (likelihood.R), model id, config.
# Outputs : list(model, spec, re, pop, indiv, cov, logLik, n_par, BIC, trace)
#
# Algorithm
#   repeat until convergence or its_max_iter:
#     E  : per-patient MAP of eta given the population prior (parallel), with
#          the Laplace posterior covariance from the Hessian
#     M1 : population means / covariate effects by least squares on eta-hat;
#          variances from residual spread plus posterior variance
#     M2 : residual SDs from the residuals at eta-hat
#     M3 : every `struct_every` iterations, structural parameters by penalised
#          maximum likelihood with eta fixed at eta-hat
#   final E-step -> Laplace marginal log-likelihood -> BIC (n = patients)
# Note    : approximate method (same family as FOCE); adequate for comparing
#           nested models, to be stated when results are reported.
# =============================================================================

par_lapply <- function(X, FUN, cores) {
  out <- if (cores > 1 && .Platform$OS.type == "unix") parallel::mclapply(X, FUN, mc.cores = cores)
         else lapply(X, FUN)
  bad <- vapply(out, function(o) inherits(o, "try-error") || is.null(o), logical(1))
  if (any(bad)) {
    msg <- if (inherits(out[[which(bad)[1]]], "try-error")) as.character(out[[which(bad)[1]]]) else "worker returned NULL"
    stop("parallel worker failed for element ", which(bad)[1], ": ", msg)
  }
  out
}

init_population <- function(spec, config, init = NULL, cov_names = NULL, struct = NULL) {
  re <- model_re_names(spec)
  cv <- if (is.null(cov_names)) config$model$covariates else cov_names
  pop <- list(
    mu = config$model$re_prior$mu[re],
    omega = config$model$re_prior$omega[re],
    beta = lapply(stats::setNames(re, re), function(k) stats::setNames(rep(0, length(cv[[k]])), cv[[k]])),
    struct = if (is.null(struct)) default_struct(spec, config) else struct,
    sigma = config$model$sigma_init)
  if (!is.null(init)) {
    shared <- intersect(re, names(init$mu))
    pop$mu[shared] <- init$mu[shared]
    pop$omega[shared] <- init$omega[shared]
    for (k in shared) {
      b <- intersect(names(pop$beta[[k]]), names(init$beta[[k]]))
      pop$beta[[k]][b] <- init$beta[[k]][b]
    }
    sh <- intersect(names(pop$struct), names(init$struct))
    keep_clone <- if (spec$family == "m2" && spec$clone != "none") c("r_a", "r_b", "delta") else character(0)
    sh <- setdiff(sh, c(keep_clone, if (spec$family == "m2" && spec$q2 != "Mabx") "eta_syn"))
    pop$struct[sh] <- init$struct[sh]
    pop$sigma <- init$sigma
  }
  pop
}

e_step <- function(pdats, spec, pop, config, eta_prev, cores) {
  re <- model_re_names(spec)
  tmpl <- par_template(pop$struct, spec)
  pdats <- prep_segments(pdats, pop$struct, config)
  par_lapply(seq_along(pdats), function(i) {
    pd <- pdats[[i]]
    prior <- patient_prior(pd, pop, re)
    eta0 <- stats::setNames(if (is.null(eta_prev)) prior$mean else eta_prev[i, re], re)
    starts <- clone_starts(eta0, spec)
    if (!is.null(eta_prev)) starts <- c(starts, list(stats::setNames(prior$mean, re)))
    fit_individual(pd, tmpl, spec, prior, pop$sigma, config, eta0, starts)
  }, cores)
}

# Variances are shrunk towards the prior spread with `omega_nu` pseudo-observations
# and capped: a random effect that the data do not inform (e.g. HSC regrowth when
# the stem cells were not hit) otherwise lets the individual estimates drift and
# the variance grow without bound from one iteration to the next.
m_step_re <- function(pdats, fits, pop, re, config) {
  eta <- do.call(rbind, lapply(fits, `[[`, "eta"))
  colnames(eta) <- re
  nu <- config$fit$omega_nu
  for (k in re) {
    b <- names(pop$beta[[k]])
    X <- if (length(b)) do.call(rbind, lapply(pdats, function(pd) pd$x[[k]][b])) else NULL
    fit <- if (is.null(X)) stats::lm(eta[, k] ~ 1) else stats::lm(eta[, k] ~ X)
    cf <- stats::coef(fit)
    cf[is.na(cf)] <- 0
    pop$mu[[k]] <- cf[[1]]
    if (length(b)) pop$beta[[k]][] <- cf[-1]
    post_var <- vapply(fits, function(f) f$cov[k, k], numeric(1))
    ss <- sum(stats::residuals(fit)^2 + post_var)
    w0 <- config$model$re_prior$omega[[k]]
    pop$omega[[k]] <- min(max(sqrt((ss + nu * w0^2) / (length(fits) + nu)), 0.05), config$fit$omega_max)
  }
  pop
}

m_step_sigma <- function(pdats, fits, pop, spec, config, cores) {
  tmpl <- par_template(pop$struct, spec)
  pdats <- prep_segments(pdats, pop$struct, config)
  res <- par_lapply(seq_along(pdats), function(i)
    residuals_at(fits[[i]]$eta, pdats[[i]], tmpl, spec, config), cores)
  for (v in lik_vars) {
    r <- unlist(lapply(res, `[[`, v))
    if (length(r) > 10) pop$sigma[[v]] <- max(sqrt(mean(r^2)), 0.05)
  }
  pop
}

struct_prior_work <- function(est, spec, config) {
  pc <- default_struct(spec, config)[est]
  w <- struct_to_work(pc, config)
  sd <- rep(config$model$struct_prior_sd, length(est))
  names(sd) <- est
  if ("eta_syn" %in% est) sd[["eta_syn"]] <- config$fit$eta_syn_prior[["sd"]]
  if ("e_sr" %in% est) sd[["e_sr"]] <- config$fit$e_sr_prior_sd
  list(w = w, sd = sd)
}

m_step_struct <- function(pdats, fits, pop, spec, config, cores) {
  est <- struct_est_names(spec, config)
  if (!length(est)) return(pop)
  pr <- struct_prior_work(est, spec, config)
  eta <- lapply(fits, `[[`, "eta")
  chunks <- split(seq_along(pdats), rep_len(seq_len(max(1, cores)), length(pdats)))
  obj <- function(w) {
    st <- pop$struct
    st[est] <- work_to_struct(stats::setNames(w, est), config)
    tmpl <- par_template(st, spec)
    pds <- prep_segments(pdats, st, config)
    parts <- par_lapply(chunks, function(ix) sum(vapply(ix, function(i)
      neg_log_lik_data(eta[[i]], pds[[i]], tmpl, spec, pop$sigma, config), numeric(1))), cores)
    sum(unlist(parts)) - sum(stats::dnorm(w, pr$w, pr$sd, log = TRUE))
  }
  w0 <- struct_to_work(pop$struct[est], config)
  o <- stats::optim(w0, obj, method = "L-BFGS-B", lower = pr$w - 4 * pr$sd, upper = pr$w + 4 * pr$sd,
                    control = list(maxit = 60))
  pop$struct[est] <- work_to_struct(stats::setNames(o$par, est), config)
  pop
}

pop_vector <- function(pop, spec, config) {
  est <- struct_est_names(spec, config)
  c(pop$mu, log(pop$omega), unlist(pop$beta), if (length(est)) struct_to_work(pop$struct[est], config))
}

count_pop_params <- function(pop, spec, config) {
  length(pop$mu) + length(pop$omega) + length(unlist(pop$beta)) + length(pop$sigma) +
    length(struct_est_names(spec, config))
}

# init: population list to start from (e.g. a simpler nested fit); eta_init:
# starting eta matrix (rows = patients in pdats order)
fit_population <- function(pdats, model, config, init = NULL, cores = config$runtime$cores,
                           log_fun = message, eta_init = NULL, max_iter = config$fit$its_max_iter,
                           struct = NULL) {
  spec <- model_spec(model)
  re <- model_re_names(spec)
  pop <- init_population(spec, config, init, cov_names = lapply(pdats[[1]]$x, names), struct = struct)
  # starting from a sub-model: an added injury term starts small, i.e. near the sub-model
  if (!is.null(init)) for (k in intersect(setdiff(re, names(init$mu)), c("th_as", "th_bs")))
    pop$mu[[k]] <- pop$mu[[k]] - 2
  eta_prev <- NULL
  if (!is.null(eta_init)) {
    eta_prev <- matrix(NA_real_, length(pdats), length(re), dimnames = list(NULL, re))
    shared <- intersect(re, colnames(eta_init))
    eta_prev[, shared] <- eta_init[, shared]
    for (k in setdiff(re, shared)) eta_prev[, k] <- pop$mu[[k]]
  }
  trace <- list()
  t_start <- Sys.time()
  for (it in seq_len(max_iter)) {
    old <- pop_vector(pop, spec, config)
    fits <- e_step(pdats, spec, pop, config, eta_prev, cores)
    eta_prev <- do.call(rbind, lapply(fits, `[[`, "eta"))
    colnames(eta_prev) <- re
    pop <- m_step_re(pdats, fits, pop, re, config)
    pop <- m_step_sigma(pdats, fits, pop, spec, config, cores)
    if (it %% config$fit$struct_every == 0) pop <- m_step_struct(pdats, fits, pop, spec, config, cores)
    delta <- max(abs(pop_vector(pop, spec, config) - old))
    trace[[it]] <- data.table::data.table(iter = it, max_change = delta,
                                          sum_nll = sum(vapply(fits, `[[`, numeric(1), "nll"), na.rm = TRUE))
    log_fun(sprintf("    %s ITS iter %d: max change %.4f (%.0fs)", model, it, delta,
                    as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    if (it >= 2 && delta < config$fit$its_tol) break
  }
  fits <- e_step(pdats, spec, pop, config, eta_prev, cores)
  logL_i <- vapply(fits, `[[`, numeric(1), "logL")
  n_par <- count_pop_params(pop, spec, config)
  logLik <- sum(logL_i, na.rm = TRUE)
  indiv <- data.table::data.table(patient_id = names(pdats))
  for (k in re) {
    indiv[, (paste0("eta_", k)) := vapply(fits, function(f) f$eta[[k]], numeric(1))]
    indiv[, (paste0("se_", k)) := vapply(fits, function(f) sqrt(f$cov[k, k]), numeric(1))]
  }
  indiv[, `:=`(logL = logL_i, nll = vapply(fits, `[[`, numeric(1), "nll"),
               conv = vapply(fits, `[[`, integer(1), "conv"))]
  list(model = model, spec = spec, re = re, pop = pop, indiv = indiv,
       cov = stats::setNames(lapply(fits, `[[`, "cov"), names(pdats)),
       logLik = logLik, n_par = n_par, n = length(pdats),
       BIC = -2 * logLik + n_par * log(length(pdats)),
       trace = data.table::rbindlist(trace),
       runtime_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

# Fit a nested model starting from the fits of the models it contains (and from
# the default start); keeps the start with the higher marginal likelihood, so a
# richer model is not reported worse than its sub-model because ITS stopped in
# a poorer local optimum.
nested_sources <- function(model) {
  sp <- model_spec(model)
  if (sp$family == "mf") return(character(0))
  q2_src <- switch(sp$q2, M0 = character(0), Ma = "M0", Mb = "M0", Mab = c("Ma", "Mb"), Mabx = "Mab")
  if (sp$clone != "none") c(sp$q2, vapply(q2_src, model_id, character(1), clone = sp$clone)) else q2_src
}

fit_nested <- function(pdats, model, fits, config, cores = config$runtime$cores, log_fun = function(...) NULL) {
  src <- intersect(nested_sources(model), names(fits))
  cands <- list(fit_population(pdats, model, config, cores = cores, log_fun = log_fun))
  if (length(src)) {
    best <- src[which.max(vapply(fits[src], `[[`, numeric(1), "logLik"))]
    cands[[2]] <- fit_population(pdats, model, config, init = fits[[best]]$pop, eta_init = fit_eta(fits[[best]]),
                                 cores = cores, log_fun = log_fun)
    cands[[2]]$start <- best
  }
  cands[[which.max(vapply(cands, `[[`, numeric(1), "logLik"))]]
}

# eta matrix (rows = patients) of a fit
fit_eta <- function(fit) {
  e <- as.matrix(fit$indiv[, paste0("eta_", fit$re), with = FALSE])
  colnames(e) <- fit$re
  rownames(e) <- fit$indiv$patient_id
  e
}

# Individual fits for patients given a fitted population.
fit_given_population <- function(pdats, fit, config, cores = config$runtime$cores) {
  fits <- e_step(pdats, fit$spec, fit$pop, config, NULL, cores)
  stats::setNames(fits, names(pdats))
}
