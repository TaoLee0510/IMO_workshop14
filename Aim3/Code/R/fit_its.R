# =============================================================================
# fit_its.R
# Purpose : Population fit by the iterative two-stage method (ITS) with a
#           Laplace approximation of each patient's marginal likelihood.
# Inputs  : list of pdat (likelihood.R), model name, config.
# Outputs : list(model, re, pop, indiv, cov, logLik, n_par, BIC, trace)
#
# Algorithm
#   repeat until convergence or its_max_iter:
#     E  : per-patient MAP of eta given the population prior (parallel), with
#          the Laplace posterior covariance from the Hessian
#     M1 : population means / covariate effects by least squares on eta-hat;
#          SDs from residual spread plus posterior variance
#     M2 : residual SDs from the residuals at eta-hat
#     M3 : every `struct_every` iterations, estimated structural parameters by
#          penalised maximum likelihood with eta fixed at eta-hat
#   final E-step -> Laplace marginal log-likelihood -> BIC (n = patients)
# Note    : approximate method (same family as FOCE); state this when reporting.
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

struct_prior_centre <- function(est, config) {
  c0 <- config$model$struct[est]
  c0[c0 <= 0] <- 1          # eta_syn: estimated on log scale, centred at 1
  log(c0)
}

init_population <- function(model, config, pdats, init = NULL) {
  re <- model_re_names(model)
  cv <- config$model$covariates
  avail <- names(pdats[[1]]$x$r_s)
  pop <- list(
    mu = config$model$re_prior$mu[re],
    omega = config$model$re_prior$omega[re],
    beta = lapply(stats::setNames(re, re), function(k) {
      b <- intersect(cv[[k]], if (k == "r_s") avail else character(0))
      stats::setNames(rep(0, length(b)), b)
    }),
    struct = config$model$struct,
    sigma = config$model$sigma_init)
  if (model == "Mabx" && pop$struct[["eta_syn"]] <= 0) pop$struct[["eta_syn"]] <- 1
  if (!is.null(init)) {
    shared <- intersect(re, names(init$mu))
    pop$mu[shared] <- init$mu[shared]
    pop$omega[shared] <- init$omega[shared]
    for (k in shared) {
      b <- intersect(names(pop$beta[[k]]), names(init$beta[[k]]))
      pop$beta[[k]][b] <- init$beta[[k]][b]
    }
    keep_syn <- pop$struct[["eta_syn"]]
    pop$struct <- init$struct
    if (model == "Mabx" && pop$struct[["eta_syn"]] <= 0) pop$struct[["eta_syn"]] <- keep_syn
    pop$sigma <- init$sigma
  }
  pop
}

e_step <- function(pdats, model, pop, config, eta_prev, cores) {
  re <- model_re_names(model)
  tmpl <- par_template(pop$struct, model)
  par_lapply(seq_along(pdats), function(i) {
    pd <- pdats[[i]]
    prior <- patient_prior(pd, pop, re)
    eta0 <- stats::setNames(if (is.null(eta_prev)) prior$mean else eta_prev[i, re], re)
    starts <- if (is.null(eta_prev)) injury_starts(eta0) else list(stats::setNames(prior$mean, re))
    fit_individual(pd, tmpl, prior, pop$sigma, config, eta0, starts)
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
    # a flat likelihood can give a huge Laplace variance; the posterior
    # variance cannot exceed the prior variance
    post_var <- pmin(vapply(fits, function(f) f$cov[k, k], numeric(1)), pop$omega[[k]]^2)
    pop$omega[[k]] <- max(sqrt(mean(stats::residuals(fit)^2 + post_var)), 0.05)
  }
  pop
}

m_step_sigma <- function(pdats, fits, pop, model, config, cores) {
  tmpl <- par_template(pop$struct, model)
  res <- par_lapply(seq_along(pdats), function(i) residuals_at(fits[[i]]$eta, pdats[[i]], tmpl, config), cores)
  for (v in obs_vars) {
    r <- unlist(lapply(res, `[[`, v))
    if (length(r) > 10) pop$sigma[[v]] <- max(sqrt(mean(r^2)), 0.05)
  }
  pop
}

m_step_struct <- function(pdats, fits, pop, model, config, cores) {
  est <- config$fit$struct_est[[model]]
  if (!length(est)) return(list(pop = pop, pdats = pdats))
  prior_w <- struct_prior_centre(est, config)
  sd_w <- config$model$struct_prior_sd[est]
  eta <- lapply(fits, `[[`, "eta")
  needs_refresh <- any(c("u_br_rgo", "u_br_cy", "beta_renal", "delta_b", "k_dec") %in% est)
  chunks <- split(seq_along(pdats), rep_len(seq_len(max(1, cores)), length(pdats)))
  obj <- function(w) {
    st <- pop$struct
    st[est] <- exp(w)
    tmpl <- par_template(st, model)
    pds <- if (needs_refresh) refresh_pdats(pdats, st, config) else pdats
    parts <- par_lapply(chunks, function(ix) sum(vapply(ix, function(i)
      neg_log_lik_data(eta[[i]], pds[[i]], tmpl, pop$sigma, config), numeric(1))), cores)
    sum(unlist(parts)) - sum(stats::dnorm(w, prior_w, sd_w, log = TRUE))
  }
  w0 <- log(pop$struct[est])
  o <- stats::optim(w0, obj, method = "L-BFGS-B", lower = prior_w - 3 * sd_w, upper = prior_w + 3 * sd_w,
                    control = list(maxit = 40))
  pop$struct[est] <- exp(o$par)
  if (needs_refresh) pdats <- refresh_pdats(pdats, pop$struct, config)
  list(pop = pop, pdats = pdats)
}

pop_vector <- function(pop, model, config) {
  est <- config$fit$struct_est[[model]]
  c(pop$mu, log(pop$omega), unlist(pop$beta), log(pop$struct[est]))
}

count_pop_params <- function(pop, model, config) {
  length(pop$mu) + length(pop$omega) + length(unlist(pop$beta)) + length(pop$sigma) +
    length(config$fit$struct_est[[model]])
}

fit_population <- function(pdats, model, config, init = NULL, cores = config$runtime$cores,
                           log_fun = message, eta_init = NULL) {
  re <- model_re_names(model)
  pop <- init_population(model, config, pdats, init)
  pdats <- refresh_pdats(pdats, pop$struct, config)
  eta_prev <- if (is.null(eta_init)) NULL else eta_init[, re, drop = FALSE]
  trace <- list()
  t_start <- Sys.time()
  for (it in seq_len(config$fit$its_max_iter)) {
    old <- pop_vector(pop, model, config)
    fits <- e_step(pdats, model, pop, config, eta_prev, cores)
    eta_prev <- do.call(rbind, lapply(fits, `[[`, "eta"))
    colnames(eta_prev) <- re
    pop <- m_step_re(pdats, fits, pop, re)
    pop <- m_step_sigma(pdats, fits, pop, model, config, cores)
    if (it %% config$fit$struct_every == 0) {
      ms <- m_step_struct(pdats, fits, pop, model, config, cores)
      pop <- ms$pop
      pdats <- ms$pdats
    }
    delta <- max(abs(pop_vector(pop, model, config) - old))
    trace[[it]] <- data.table::data.table(iter = it, max_change = delta,
                                          sum_nll = sum(vapply(fits, `[[`, numeric(1), "nll"), na.rm = TRUE))
    log_fun(sprintf("    %s ITS iter %d: max change %.4f (%.0fs)", model, it, delta,
                    as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
    if (it >= 2 && delta < config$fit$its_tol) break
  }
  fits <- e_step(pdats, model, pop, config, eta_prev, cores)
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
  list(model = model, re = re, pop = pop, indiv = indiv,
       fits = stats::setNames(fits, names(pdats)),
       logLik = logLik, n_par = n_par, n = length(pdats),
       BIC = -2 * logLik + n_par * log(length(pdats)),
       trace = data.table::rbindlist(trace),
       runtime_s = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
}

# Fit every nested model; Mab starts from Ma/Mb, Mabx from Mab.
fit_all_models <- function(pdats, config, models = config$fit$models, cores = config$runtime$cores,
                           log_fun = message) {
  out <- list()
  for (m in models) {
    init <- switch(m, Mab = out$Ma$pop, Mabx = out$Mab$pop, NULL)
    log_fun("  fitting ", m)
    out[[m]] <- fit_population(pdats, m, config, init = init, cores = cores, log_fun = log_fun)
  }
  out
}

model_table <- function(fits) {
  tb <- data.table::rbindlist(lapply(fits, function(f)
    data.table::data.table(model = f$model, logLik = f$logLik, n_par = f$n_par, BIC = f$BIC,
                           runtime_s = round(f$runtime_s))))
  tb[, delta_BIC := round(BIC - min(BIC), 2)]
  tb[]
}
