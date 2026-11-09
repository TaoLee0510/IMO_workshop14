# =============================================================================
# fit_its.R
# Purpose : Population fit of the marrow twin by the iterative two-stage
#           method (ITS) with a Laplace approximation of each patient's
#           marginal likelihood. Dependency-free alternative to nlmixr2.
# Inputs  : list of pdat (likelihood.R), model name, config.
# Outputs : list(model, re, pop, indiv, fits, logLik, n_par, BIC, trace)
#
# Algorithm
#   repeat until convergence or its_max_iter:
#     E  : per-patient MAP of eta given the population prior (parallel),
#          with the Laplace posterior covariance from the Hessian
#     M1 : population means/covariate effects by least squares on eta-hat;
#          variances from residual spread plus posterior variance (ITS)
#     M2 : residual SDs from the residuals at eta-hat
#     M3 : every `struct_every` iterations, structural parameters by
#          penalised maximum likelihood with eta fixed at eta-hat
#   final E-step -> Laplace marginal log-likelihood -> BIC (n = patients)
# Note    : this is an approximate method (same family as FOCE); adequate for
#           comparing nested models, to be stated when results are reported.
# =============================================================================

# Parallel lapply (forked on unix). Worker errors are re-raised with their own
# message instead of surfacing later as an unrelated subscript error.
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

# cov_names: covariate names per random effect as present in the pdats
# (default: the config list)
init_population <- function(model, config, cov_f0 = TRUE, init = NULL, cov_names = NULL) {
  re <- model_re_names(model)
  cv <- if (is.null(cov_names)) config$model$covariates else cov_names
  if (!cov_f0) cv$f0 <- character(0)
  pop <- list(
    mu = config$model$re_prior$mu[re],
    omega = config$model$re_prior$omega[re],
    beta = lapply(stats::setNames(re, re), function(k) stats::setNames(rep(0, length(cv[[k]])), cv[[k]])),
    struct = config$model$struct,
    sigma = config$fit$sigma_init)
  if (!is.null(init)) {
    shared <- intersect(re, names(init$mu))
    pop$mu[shared] <- init$mu[shared]
    pop$omega[shared] <- init$omega[shared]
    for (k in shared) {
      b <- intersect(names(pop$beta[[k]]), names(init$beta[[k]]))
      pop$beta[[k]][b] <- init$beta[[k]][b]
    }
    pop$struct <- init$struct
    pop$sigma <- init$sigma
  }
  pop
}

e_step <- function(pdats, model, pop, config, eta_prev, use_vaf, cores) {
  re <- model_re_names(model)
  tmpl <- par_template(pop$struct, model)
  par_lapply(seq_along(pdats), function(i) {
    pd <- pdats[[i]]
    prior <- patient_prior(pd, pop, re)
    eta0 <- stats::setNames(if (is.null(eta_prev)) prior$mean else eta_prev[i, re], re)
    starts <- clone_starts(eta0, model, if (use_vaf) pd else NULL)
    if (model == "M2" && !is.null(eta_prev)) starts <- c(starts, list(stats::setNames(prior$mean, re)))
    fit_individual(pd, tmpl, model, prior, pop$sigma, config, eta0, use_vaf, starts)
  }, cores)
}

m_step_re <- function(pdats, fits, pop, re) {
  eta <- do.call(rbind, lapply(fits, `[[`, "eta"))
  colnames(eta) <- re
  for (k in re) {
    b <- names(pop$beta[[k]])
    X <- if (length(b)) do.call(rbind, lapply(pdats, function(pd) pd$x[[k]][b])) else NULL
    fit <- if (is.null(X)) stats::lm(eta[, k] ~ 1) else stats::lm(eta[, k] ~ X)
    cf <- stats::coef(fit)
    cf[is.na(cf)] <- 0
    pop$mu[[k]] <- cf[[1]]
    if (length(b)) pop$beta[[k]][] <- cf[-1]
    post_var <- vapply(fits, function(f) f$cov[k, k], numeric(1))
    pop$omega[[k]] <- sqrt(mean(stats::residuals(fit)^2 + post_var))
  }
  pop
}

m_step_sigma <- function(pdats, fits, pop, model, config, cores) {
  tmpl <- par_template(pop$struct, model)
  res <- par_lapply(seq_along(pdats), function(i)
    residuals_at(fits[[i]]$eta, pdats[[i]], tmpl, model, config), cores)
  for (v in obs_vars) {
    r <- unlist(lapply(res, `[[`, v))
    if (length(r) > 10) pop$sigma[[v]] <- max(sqrt(mean(r^2)), 0.05)
  }
  pop
}

m_step_struct <- function(pdats, fits, pop, model, config, cores) {
  est <- config$fit$struct_est[[model]]
  if (!length(est)) return(pop)
  prior_w <- struct_to_work(config$model$struct[est])
  sd_w <- config$model$struct_prior_sd
  eta <- lapply(fits, `[[`, "eta")
  chunks <- split(seq_along(pdats), rep_len(seq_len(max(1, cores)), length(pdats)))
  obj <- function(w) {
    st <- pop$struct
    st[est] <- work_to_struct(stats::setNames(w, est))
    tmpl <- par_template(st, model)
    parts <- par_lapply(chunks, function(ix) sum(vapply(ix, function(i)
      neg_log_lik_data(eta[[i]], pdats[[i]], tmpl, model, pop$sigma, config), numeric(1))), cores)
    sum(unlist(parts)) - sum(stats::dnorm(w, prior_w, sd_w, log = TRUE))
  }
  w0 <- struct_to_work(pop$struct[est])
  o <- stats::optim(w0, obj, method = "L-BFGS-B", lower = prior_w - 4 * sd_w, upper = prior_w + 4 * sd_w,
                    control = list(maxit = 60))
  pop$struct[est] <- work_to_struct(stats::setNames(o$par, est))
  pop
}

pop_vector <- function(pop, model, config) {
  est <- config$fit$struct_est[[model]]
  c(pop$mu, log(pop$omega), unlist(pop$beta), if (length(est)) struct_to_work(pop$struct[est]))
}

count_pop_params <- function(pop, model, config) {
  length(pop$mu) + length(pop$omega) + length(unlist(pop$beta)) + length(pop$sigma) +
    length(config$fit$struct_est[[model]])
}

# eta_init: optional starting eta matrix (rows = patients in pdats order), e.g.
# from a previous fit, so that a refit starts from that fit's individual values.
fit_population <- function(pdats, model, config, cov_f0 = TRUE, use_vaf = TRUE, init = NULL,
                           cores = config$runtime$cores, log_fun = message, eta_init = NULL) {
  re <- model_re_names(model)
  pop <- init_population(model, config, cov_f0, init, cov_names = lapply(pdats[[1]]$x, names))
  if (!cov_f0) pdats <- lapply(pdats, function(pd) { pd$x$f0 <- numeric(0); pd })
  eta_prev <- if (is.null(eta_init)) NULL else eta_init[, re, drop = FALSE]
  trace <- list()
  t_start <- Sys.time()
  for (it in seq_len(config$fit$its_max_iter)) {
    old <- pop_vector(pop, model, config)
    fits <- e_step(pdats, model, pop, config, eta_prev, use_vaf, cores)
    eta_prev <- do.call(rbind, lapply(fits, `[[`, "eta"))
    colnames(eta_prev) <- re
    pop <- m_step_re(pdats, fits, pop, re)
    pop <- m_step_sigma(pdats, fits, pop, model, config, cores)
    if (it %% config$fit$struct_every == 0) pop <- m_step_struct(pdats, fits, pop, model, config, cores)
    delta <- max(abs(pop_vector(pop, model, config) - old))
    trace[[it]] <- data.table::data.table(iter = it, max_change = delta,
                                          sum_nll = sum(vapply(fits, `[[`, numeric(1), "nll"), na.rm = TRUE))
    log_fun(sprintf("    %s ITS iter %d: max change %.4f (%.0fs)", model, it, delta,
                    as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    if (it >= 2 && delta < config$fit$its_tol) break
  }
  fits <- e_step(pdats, model, pop, config, eta_prev, use_vaf, cores)
  logL_i <- vapply(fits, `[[`, numeric(1), "logL")
  n_par <- count_pop_params(pop, model, config)
  logLik <- sum(logL_i, na.rm = TRUE)
  indiv <- data.table::data.table(patient_id = names(pdats))
  for (k in re) {
    indiv[, (paste0("eta_", k)) := vapply(fits, function(f) f$eta[[k]], numeric(1))]
    indiv[, (paste0("se_", k)) := vapply(fits, function(f) sqrt(f$cov[k, k]), numeric(1))]
  }
  indiv[, `:=`(logL = logL_i, nll = vapply(fits, `[[`, numeric(1), "nll"),
               conv = vapply(fits, `[[`, integer(1), "conv"))]
  list(model = model, re = re, cov_f0 = cov_f0, use_vaf = use_vaf, pop = pop, indiv = indiv,
       cov = stats::setNames(lapply(fits, `[[`, "cov"), names(pdats)),
       logLik = logLik, n_par = n_par, n = length(pdats),
       BIC = -2 * logLik + n_par * log(length(pdats)),
       trace = data.table::rbindlist(trace),
       runtime_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

# Profile likelihood over clone traits that trade off with f0. At each grid
# point the clone traits are fixed and the whole M2 population is refitted
# (ITS, started from `fit`), so the surface is a true profile likelihood. Only
# re-optimising the individual parameters at the population of `fit` gives a
# conditional likelihood that favours the traits `fit` was obtained at.
# Returns the profile table, the best traits and the refitted fit at the best
# point (fit_best).
profile_clone_traits <- function(pdats, fit, config, cores = config$runtime$cores, log_fun = message,
                                 keep_eta = FALSE) {
  eta_prev <- as.matrix(fit$indiv[, paste0("eta_", fit$re), with = FALSE])
  colnames(eta_prev) <- fit$re
  if (!fit$cov_f0) pdats <- lapply(pdats, function(pd) { pd$x$f0 <- numeric(0); pd })
  profile_grid(pdats, fit$model, fit$pop, do.call(data.table::CJ, config$fit$profile), config,
               use_vaf = fit$use_vaf, eta_prev = eta_prev, cores = cores, log_fun = log_fun,
               keep_eta = keep_eta, refit = TRUE)
}

# Profile likelihood over a grid of structural parameters (columns of `grid`).
# Every grid point gets the same starting values (eta_prev, or the prior mean
# plus the clone starts when NULL), so the surface does not depend on the order
# in which the grid is visited.
#   refit = FALSE : one E-step at the given population (likelihood conditional
#                   on `pop`; used by D3b, where the population is the truth)
#   refit = TRUE  : full population refit at every grid point (profile
#                   likelihood); the fit at the best point is returned
# keep_eta = TRUE returns the MAP eta matrix of every grid point (rows = patients).
profile_grid <- function(pdats, model, pop, grid, config, use_vaf = TRUE, eta_prev = NULL,
                         cores = config$runtime$cores, log_fun = NULL, keep_eta = FALSE,
                         refit = FALSE) {
  re <- model_re_names(model)
  pars <- names(grid)
  grid <- data.table::copy(grid)[, logLik := NA_real_]
  etas <- vector("list", nrow(grid))
  fit_best <- NULL
  quiet <- function(...) invisible(NULL)
  for (j in seq_len(nrow(grid))) {
    t0 <- Sys.time()
    p <- pop
    for (k in pars) p$struct[[k]] <- grid[[k]][j]
    if (refit) {
      f <- fit_population(pdats, model, config, use_vaf = use_vaf, init = p, cores = cores,
                          log_fun = quiet, eta_init = eta_prev)
      data.table::set(grid, j, "logLik", f$logLik)
      eta_j <- as.matrix(f$indiv[, paste0("eta_", re), with = FALSE])
      if (is.null(fit_best) || f$logLik > fit_best$logLik) fit_best <- f
    } else {
      fits <- e_step(pdats, model, p, config, eta_prev, use_vaf, cores)
      data.table::set(grid, j, "logLik", sum(vapply(fits, `[[`, numeric(1), "logL"), na.rm = TRUE))
      eta_j <- do.call(rbind, lapply(fits, `[[`, "eta"))
    }
    if (keep_eta) etas[[j]] <- `dimnames<-`(eta_j, list(names(pdats), re))
    if (!is.null(log_fun))
      log_fun(sprintf("    profile %s: logLik %.1f (%.0fs)",
                      paste(pars, unlist(grid[j, pars, with = FALSE]), sep = "=", collapse = ", "), grid$logLik[j],
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  grid[, delta_logLik := round(logLik - max(logLik), 2)]
  list(table = grid, best = unlist(grid[which.max(logLik), pars, with = FALSE]),
       eta = if (keep_eta) etas else NULL, fit_best = fit_best)
}

# Individual fits for new patients given a fitted population (hold-out, CH subset).
fit_given_population <- function(pdats, fit, config, use_vaf = TRUE, cores = config$runtime$cores) {
  model <- fit$model
  re <- fit$re
  pop <- fit$pop
  if (!fit$cov_f0) pdats <- lapply(pdats, function(pd) { pd$x$f0 <- numeric(0); pd })
  fits <- e_step(pdats, model, pop, config, NULL, use_vaf, cores)
  stats::setNames(fits, names(pdats))
}

# Natural-scale individual parameters from an eta vector.
eta_named_to_indiv <- function(eta) {
  out <- c(rho = exp(eta[["rho"]]), g = exp(eta[["g"]]))
  if ("f0" %in% names(eta)) out["f0"] <- stats::plogis(eta[["f0"]])
  out
}
