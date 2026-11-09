# =============================================================================
# fit_de.R
# Purpose : Population fit by DE-ITS (Aim4 v1.3, plus round 0 for Model1): iterative two-stage estimation in
#           which every optimisation is a global differential-evolution search
#           (DEoptim) with self-detected convergence, plus the state files that
#           let the cluster run it as one Slurm array task per patient (E step)
#           and one job per M step that decides the next round itself.
# Inputs  : data object (read_model1_data / simulate_cohort), model id, config,
#           a run directory.
# Outputs : fit object with the fields of fit_population() (v1.2) plus the
#           importance-sampling marginal likelihoods and DE diagnostics.
#
# Round 0 (Model1; runs without a warm start): every patient fitted on its own with the prior
#   SDs times round0$sd_mult (same search box); mu, beta, omega, tau and sigma are initialised
#   from these MAPs and Hessians by the M-step formulas (M1, M2); the structural parameters keep
#   their prior centres until the M step of round 1.
# Round k
#   E  : per patient, psi = (random effects, baseline deviations) maximising the
#        posterior given population k: DEoptim in chunks of `chunk` generations
#        inside the physiological box, converged when the best value improved by
#        < tol_best (|f| + 1) over `window` generations and the population's
#        objective spread (90% quantile - min) < tol_spread; L-BFGS-B polish;
#        Laplace covariance from the Hessian
#   M1 : mu, beta (least squares on psi_hat), omega and the set-point drift tau
#        (residual + posterior variance, shrinkage towards the prior spread)
#   M2 : residual SDs (RMS of the residuals at psi_hat, incl. ALC)
#   M3 : structural parameters by DEoptim with psi fixed (same convergence rule)
#   convergence: max |change| of the population vector < tol_par and relative
#   change of sum(nll) < tol_obj on n_ok consecutive rounds; components whose
#   change flips sign are damped; max_iter is a safety cap (flag, not a fit).
#   Since 2026-11-04 (deviation): components that keep their direction are
#   over-relaxed; once no component trends over `window` rounds the population
#   is averaged (window mean, then stochastic approximation with step 1/n) and
#   converged when the average is stable (the noise floor of the DE E step)
# Final: E step at the converged population + importance-sampling marginal
#   likelihood (multivariate t proposal at psi_hat), BIC = -2 sum log L + n_par log n.
# Objectives read their data from a global context (.DE) set before the fork
#   cluster is created, so a generation sends only parameter vectors to workers.
# =============================================================================

# ---- population container --------------------------------------------------------------------
init_population13 <- function(spec, config, base_s, cov_names) {
  re <- model_re_names(spec)
  bv <- sub("^b_", "", baseline_names(spec, config))
  list(mu = config$model$re_prior$mu[re], omega = config$model$re_prior$omega[re],
       beta = lapply(stats::setNames(re, re), function(k)
         stats::setNames(rep(0, length(cov_names[[k]])), cov_names[[k]])),
       struct = default_struct(spec, config),
       sigma = config$model$sigma_init[bv],
       tau = stats::setNames(rep(config$model$baseline$tau_init, length(bv)), bv),
       base_s = base_s[bv])
}

pop_vector13 <- function(pop, spec, config) {
  est <- struct_est_names(spec, config)
  b <- unlist(lapply(names(pop$beta), function(k)
    if (length(pop$beta[[k]])) stats::setNames(pop$beta[[k]], paste0("beta_", k, ".", names(pop$beta[[k]])))))
  c(stats::setNames(pop$mu, paste0("mu_", names(pop$mu))),
    stats::setNames(log(pop$omega), paste0("logomega_", names(pop$omega))),
    b,
    stats::setNames(log(pop$tau), paste0("logtau_", names(pop$tau))),
    stats::setNames(log(pop$sigma), paste0("logsigma_", names(pop$sigma))),
    if (length(est)) stats::setNames(struct_to_work(pop$struct[est], config), paste0("struct_", est)))
}

pop_set13 <- function(pop, v, spec, config) {
  for (k in names(pop$mu)) pop$mu[[k]] <- v[[paste0("mu_", k)]]
  for (k in names(pop$omega)) pop$omega[[k]] <- exp(v[[paste0("logomega_", k)]])
  for (k in names(pop$beta)) for (j in names(pop$beta[[k]])) pop$beta[[k]][[j]] <- v[[paste0("beta_", k, ".", j)]]
  for (k in names(pop$tau)) pop$tau[[k]] <- exp(v[[paste0("logtau_", k)]])
  for (k in names(pop$sigma)) pop$sigma[[k]] <- exp(v[[paste0("logsigma_", k)]])
  est <- struct_est_names(spec, config)
  if (length(est)) pop$struct[est] <- work_to_struct(stats::setNames(v[paste0("struct_", est)], est), config)
  pop
}

count_pop_params13 <- function(pop, spec, config) {
  length(pop$mu) + length(pop$omega) + length(unlist(pop$beta)) + length(pop$tau) + length(pop$sigma) +
    length(struct_est_names(spec, config))
}

# ---- search boxes --------------------------------------------------------------------------------
indiv_box <- function(prior, config) {
  nm <- names(prior$mean)
  lim <- config$model$indiv_limits
  bl <- config$model$baseline
  lo <- hi <- stats::setNames(numeric(length(nm)), nm)
  for (k in nm) {
    if (startsWith(k, "b_")) {
      w <- max(bl$bound_min, bl$bound_sd * prior$sd[[k]])
      lo[[k]] <- -w; hi[[k]] <- w
    } else if (k == "f0") {
      lo[[k]] <- stats::qlogis(lim$f0[1]); hi[[k]] <- stats::qlogis(lim$f0[2])
    } else {
      lo[[k]] <- log(lim[[k]][1]); hi[[k]] <- log(lim[[k]][2])
    }
  }
  list(lo = lo, hi = hi)
}

struct_prior_work13 <- function(est, spec, config) {
  pc <- default_struct(spec, config)[est]
  w <- struct_to_work(pc, config)
  sd <- stats::setNames(rep(config$model$struct_prior_sd, length(est)), est)
  ov <- config$model$struct_prior_sd_k
  k <- intersect(names(ov), est)
  sd[k] <- ov[k]
  if ("eta_syn" %in% est) {
    w[["eta_syn"]] <- log(config$fit$eta_syn_prior[["centre"]])
    sd[["eta_syn"]] <- config$fit$eta_syn_prior[["sd"]]
  }
  list(w = w, sd = sd)
}

struct_box <- function(pr, config) {
  lo <- pr$w - 4 * pr$sd
  hi <- pr$w + 4 * pr$sd
  lim <- config$model$struct_limits
  for (k in intersect(names(lim), names(lo))) {
    lw <- struct_to_work(stats::setNames(lim[[k]], c(k, k)), config)
    lo[[k]] <- max(lo[[k]], lw[[1]])
    hi[[k]] <- min(hi[[k]], lw[[2]])
  }
  list(lo = lo, hi = hi)
}

# ---- differential evolution with self-detected convergence -----------------------------------------
# fn: objective of a numeric vector (reads .DE when run on a fork cluster); seeds: start
# vectors placed in the initial population (clipped to the box)
de_run <- function(fn, lo, hi, seeds = list(), config, cl = NULL, np = config$fit$de$np,
                   tol_spread = config$fit$de$tol_spread) {
  s <- config$fit$de
  nm <- names(lo)
  D <- length(lo)
  NP <- max(as.integer(np), 4L * D, 10L)
  pop <- matrix(stats::runif(NP * D), NP, D) * rep(hi - lo, each = NP) + rep(lo, each = NP)
  j <- 0L
  for (sd in seeds) {
    if (is.null(sd) || j >= NP) next
    x <- lo + (hi - lo) / 2
    x[intersect(names(sd), nm)] <- sd[intersect(names(sd), nm)]
    j <- j + 1L
    pop[j, ] <- pmin(pmax(x[nm], lo), hi)
  }
  colnames(pop) <- nm
  eval_pop <- function(P) if (is.null(cl)) apply(P, 1, fn) else as.numeric(parallel::parApply(cl, P, 1, fn))
  best_hist <- numeric(0)
  gen <- 0L
  converged <- FALSE
  spread <- NA_real_
  repeat {
    ctl <- DEoptim::DEoptim.control(NP = NP, itermax = s$chunk, F = s$F, CR = s$CR, strategy = s$strategy,
                                    initialpop = pop, trace = FALSE, cluster = cl)
    r <- suppressWarnings(DEoptim::DEoptim(fn, lo, hi, ctl))
    gen <- gen + as.integer(s$chunk)
    pop <- r$member$pop
    colnames(pop) <- nm
    best_hist <- c(best_hist, r$member$bestvalit)
    f <- utils::tail(best_hist, 1)
    if (length(best_hist) > s$window) {
      prev <- best_hist[length(best_hist) - s$window]
      if (prev - f < s$tol_best * (abs(f) + 1)) {
        vals <- eval_pop(pop)
        spread <- stats::quantile(vals, 0.9, names = FALSE) - min(vals)
        if (spread < tol_spread) { converged <- TRUE; break }
      }
    }
    if (gen >= s$max_gen) break
  }
  list(par = stats::setNames(as.numeric(r$optim$bestmem), nm), value = r$optim$bestval, gen = gen,
       converged = converged, spread = spread, NP = NP)
}

# Fork cluster with a port derived from the process id (several tasks share a node; the
# default random port can collide), retried on another port if it cannot be opened.
de_fork_cluster <- function(cores) {
  for (attempt in 0:9) {
    port <- 11000L + (as.integer(Sys.getpid()) * 7L + attempt * 1009L) %% 40000L
    cl <- tryCatch(parallel::makeForkCluster(cores, port = port), error = function(e) NULL)
    if (!is.null(cl)) return(cl)
    Sys.sleep(1)
  }
  stop("could not create the fork cluster (10 ports tried)")
}

# ---- E step: one patient ------------------------------------------------------------------------------
de_indiv_fn <- function(x) {
  ctx <- get(".DE", envir = globalenv())
  neg_log_post(stats::setNames(x, ctx$nm), ctx$pd, ctx$tmpl, ctx$spec, ctx$prior, ctx$sigma, ctx$config)
}
environment(de_indiv_fn) <- globalenv()

# Sets the global context, creates the fork cluster (cores > 1), fits, and cleans up.
# box: search box (default from the prior; round 0 passes the box of the unwidened prior)
de_fit_individual <- function(pd, tmpl, spec, prior, sigma, config, seeds = list(), cores = 1L,
                              with_is = FALSE, box = NULL) {
  nm <- names(prior$mean)
  if (is.null(box)) box <- indiv_box(prior, config)
  assign(".DE", list(pd = pd, tmpl = tmpl, spec = spec, prior = prior, sigma = sigma, config = config, nm = nm,
                     box = box), envir = globalenv())
  cl <- if (cores > 1) de_fork_cluster(cores) else NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
  t0 <- Sys.time()
  de <- de_run(de_indiv_fn, box$lo, box$hi, c(seeds, list(prior$mean)), config, cl)
  pol <- tryCatch(stats::optim(de$par, de_indiv_fn, method = "L-BFGS-B", lower = box$lo, upper = box$hi,
                               control = list(maxit = config$fit$de$polish_maxit)), error = function(e) NULL)
  par <- de$par
  val <- de$value
  if (!is.null(pol) && is.finite(pol$value) && pol$value < val) {
    par <- stats::setNames(pol$par, nm)
    val <- pol$value
  }
  # Hessian made positive definite with eigenvalues at least the smallest prior precision
  # (at a bound of the box the data Hessian can be indefinite)
  H <- pd_floor(stats::optimHess(par, de_indiv_fn), 1 / max(prior$sd^2))
  k <- length(par)
  logdet <- as.numeric(determinant(H, logarithm = TRUE)$modulus)
  covm <- solve(H)
  dimnames(covm) <- list(nm, nm)
  w <- box$hi - box$lo
  # ALC log-likelihood at psi_hat: the ANC/PLT/AMC-only Laplace value (comparison with MF)
  # is the Laplace value without it
  m_hat <- simulate_pdat(par, pd, tmpl, spec, config)
  ll_alc <- if (anyNA(m_hat[, 2])) NA_real_ else data_loglik(m_hat, pd, sigma, pdat_baselines(par, pd), vars = "alc")
  logL_lap <- -val + 0.5 * k * log(2 * pi) - 0.5 * logdet - box_log_mass(prior, box)
  out <- list(id = pd$id, eta = par, cov = covm, nll = val,
              logL = logL_lap, ll_alc_map = ll_alc, logL_core_laplace = logL_lap - ll_alc,
              conv = if (de$converged) 0L else 1L, gen = de$gen, spread = de$spread,
              at_bound = nm[par <= box$lo + 1e-4 * w | par >= box$hi - 1e-4 * w],
              runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (with_is) out$is <- de_is_logml(par, covm, config, cl)
  out
}

pd_floor <- function(H, floor) {
  if (any(!is.finite(H))) return(diag(floor, nrow(H)))
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE)
  ev$values <- pmin(pmax(ev$values, floor), 1e8)
  ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
}

# The search box is the support of the individual parameters (truncated prior): log of
# the prior mass inside the box (independent normal prior per component)
box_log_mass <- function(prior, box) {
  nm <- names(box$lo)
  m <- prior$mean[nm]
  sd <- prior$sd[nm]
  sum(log(pmax(stats::pnorm((box$hi - m) / sd) - stats::pnorm((box$lo - m) / sd), 1e-300)))
}

# ---- importance-sampling marginal likelihood (uses the .DE context) -----------------------------------
de_is_point <- function(x) {
  ctx <- get(".DE", envir = globalenv())
  if (any(x < ctx$box$lo | x > ctx$box$hi)) return(c(-Inf, -Inf))
  x <- stats::setNames(x, ctx$nm)
  m <- simulate_pdat(x, ctx$pd, ctx$tmpl, ctx$spec, ctx$config)
  if (anyNA(m[, 2])) return(c(-Inf, -Inf))
  b <- pdat_baselines(x, ctx$pd)
  lpr <- sum(stats::dnorm(x, ctx$prior$mean[ctx$nm], ctx$prior$sd[ctx$nm], log = TRUE)) + vaf_loglik(x, ctx$pd) +
    pheno_loglik(m, ctx$pd, b, ctx$config)
  c(data_loglik(m, ctx$pd, ctx$sigma, b) + lpr, data_loglik(m, ctx$pd, ctx$sigma, b, core_vars) + lpr)
}
environment(de_is_point) <- globalenv()

# Two-stage adaptive importance sampling: stage 1 from a multivariate t at psi_hat with the
# scaled Laplace covariance; stage 2 from a t with the weighted mean / covariance of stage 1
# (inflated when stage 1 has too few effective draws); estimate with the deterministic-
# mixture weights of both stages (balance heuristic).
de_is_logml <- function(eta, covm, config, cl = NULL) {
  s <- config$fit$is
  D <- length(eta)
  nu <- s$df
  ctx <- get(".DE", envir = globalenv())
  cap <- max(ctx$prior$sd^2)
  prop <- function(mu, S) {
    S <- solve(pd_floor(solve(pd_floor(S, 1e-8)), 1 / cap))
    L <- chol(S)
    list(mu = mu, L = L, Si = chol2inv(L), ld = sum(log(diag(L))))
  }
  draw <- function(p, n) {
    z <- matrix(stats::rnorm(n * D), n, D) %*% p$L
    sweep(z * sqrt(nu / stats::rchisq(n, nu)), 2, p$mu, `+`)
  }
  logq <- function(p, X) {
    dev <- sweep(X, 2, p$mu)
    q <- rowSums((dev %*% p$Si) * dev)
    lgamma((nu + D) / 2) - lgamma(nu / 2) - D / 2 * log(nu * pi) - p$ld - (nu + D) / 2 * log1p(q / nu)
  }
  evalx <- function(X) if (is.null(cl)) t(apply(X, 1, de_is_point)) else t(parallel::parApply(cl, X, 1, de_is_point))
  lse <- function(a) { mx <- max(a); if (!is.finite(mx)) return(-Inf); mx + log(sum(exp(a - mx))) }
  # centre slightly inside the box when psi_hat lies on a bound
  w <- ctx$box$hi - ctx$box$lo
  p1 <- prop(pmin(pmax(eta, ctx$box$lo + 0.02 * w), ctx$box$hi - 0.02 * w), s$scale * covm)
  X1 <- draw(p1, s$n)
  v1 <- evalx(X1)
  w1 <- v1[, 1] - logq(p1, X1)
  ww <- exp(w1 - max(w1))
  ww <- ww / sum(ww)
  ess1 <- 1 / sum(ww^2)
  p2 <- if (is.finite(ess1) && ess1 >= 5) prop(colSums(X1 * ww), s$scale * stats::cov.wt(X1, wt = ww)$cov) else
    prop(p1$mu, 4 * s$scale * covm)
  n2 <- 2L * s$n
  X2 <- draw(p2, n2)
  v2 <- evalx(X2)
  X <- rbind(X1, X2)
  V <- rbind(v1, v2)
  n1 <- nrow(X1)
  lq1 <- logq(p1, X); lq2 <- logq(p2, X)
  mx <- pmax(lq1, lq2)
  lq <- mx + log((n1 * exp(lq1 - mx) + n2 * exp(lq2 - mx)) / (n1 + n2))
  wf <- V[, 1] - lq
  wc <- V[, 2] - lq
  N <- nrow(X)
  lz <- box_log_mass(get(".DE", envir = globalenv())$prior, get(".DE", envir = globalenv())$box)
  c(logml = lse(wf) - log(N) - lz, logml_core = lse(wc) - log(N) - lz, ess = exp(2 * lse(wf) - lse(2 * wf)),
    ess_stage1 = ess1, n_failed = sum(!is.finite(V[, 1])))
}

# ---- M step ----------------------------------------------------------------------------------------------
m_step_re13 <- function(pdats, fits, pop, re, config) {
  eta <- do.call(rbind, lapply(fits, function(f) f$eta[re]))
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

# set-point drift: tau^2 = shrunk mean of (d_hat^2 + posterior var - s^2 / n)
m_step_tau <- function(pdats, fits, pop, config) {
  if (!length(pop$tau)) return(pop)
  nu <- config$fit$omega_nu
  bl <- config$model$baseline
  for (v in names(pop$tau)) {
    k <- paste0("b_", v)
    d <- vapply(fits, function(f) f$eta[[k]], numeric(1))
    pv <- vapply(fits, function(f) f$cov[k, k], numeric(1))
    n <- vapply(pdats, function(pd) max(pd$base_n[[v]], 1), numeric(1))
    ss <- sum(d^2 + pv - pop$base_s[[v]]^2 / n)
    pop$tau[[v]] <- sqrt(max((ss + nu * bl$tau_init^2) / (length(fits) + nu), bl$tau_min^2))
  }
  pop
}

m_step_sigma13 <- function(pdats, fits, pop, spec, config) {
  tmpl <- par_template(pop$struct, spec)
  pdats <- prep_segments(pdats, pop$struct, config)
  res <- lapply(seq_along(pdats), function(i) residuals_at(fits[[i]]$eta, pdats[[i]], tmpl, spec, config))
  for (v in names(pop$sigma)) {
    r <- unlist(lapply(res, `[[`, v))
    if (length(r) > 10) pop$sigma[[v]] <- max(sqrt(mean(r^2)), 0.05)
  }
  pop
}

de_struct_fn <- function(w) {
  ctx <- get(".DE", envir = globalenv())
  st <- ctx$struct
  st[ctx$est] <- work_to_struct(stats::setNames(w, ctx$est), ctx$config)
  tmpl <- par_template(st, ctx$spec)
  tot <- 0
  for (i in seq_along(ctx$pdats)) {
    pd <- ctx$pdats[[i]]
    pd$seg <- make_segments(pd$design, st, max(pd$times), ctx$config)
    tot <- tot + neg_log_lik_data(ctx$eta[[i]], pd, tmpl, ctx$spec, ctx$sigma, ctx$config)
  }
  tot - sum(stats::dnorm(w, ctx$pr$w, ctx$pr$sd, log = TRUE))
}
environment(de_struct_fn) <- globalenv()

m_step_struct_de <- function(pdats, fits, pop, spec, config, cores = 1L) {
  est <- struct_est_names(spec, config)
  if (!length(est)) return(list(pop = pop, diag = NULL))
  pr <- struct_prior_work13(est, spec, config)
  box <- struct_box(pr, config)
  assign(".DE", list(pdats = pdats, eta = lapply(fits, `[[`, "eta"), struct = pop$struct, est = est, pr = pr,
                     spec = spec, sigma = pop$sigma, config = config), envir = globalenv())
  cl <- if (cores > 1) de_fork_cluster(cores) else NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
  np <- max(config$fit$de$struct_np_min, config$fit$de$struct_np_per_dim * length(est))
  w0 <- struct_to_work(pop$struct[est], config)
  # the objective sums over patients: spread tolerance per patient
  de <- de_run(de_struct_fn, box$lo, box$hi, list(w0, pr$w), config, cl, np = np,
               tol_spread = config$fit$de$tol_spread * length(pdats))
  pol <- tryCatch(stats::optim(de$par, de_struct_fn, method = "L-BFGS-B", lower = box$lo, upper = box$hi,
                               control = list(maxit = 100)), error = function(e) NULL)
  par <- de$par
  val <- de$value
  if (!is.null(pol) && is.finite(pol$value) && pol$value < val) {
    par <- stats::setNames(pol$par, est)
    val <- pol$value
  }
  pop$struct[est] <- work_to_struct(par, config)
  w <- box$hi - box$lo
  list(pop = pop, diag = list(gen = de$gen, converged = de$converged, value = val,
                              at_bound = est[par <= box$lo + 1e-4 * w | par >= box$hi - 1e-4 * w]))
}

# Legal range of every population component (working scale): outer-loop updates are clipped to it.
# xrange: observed range of every covariate (covariate_ranges); a covariate effect may move a
# patient across at most the whole individual box of its random effect over that range.
pop_limits13 <- function(pop, spec, config, xrange = NULL) {
  v <- pop_vector13(pop, spec, config)
  lo <- stats::setNames(rep(-Inf, length(v)), names(v))
  hi <- -lo
  lim <- config$model$indiv_limits
  for (k in names(pop$mu)) {
    r <- if (k == "f0") stats::qlogis(lim$f0) else log(lim[[k]])
    lo[[paste0("mu_", k)]] <- r[1]; hi[[paste0("mu_", k)]] <- r[2]
    for (j in names(pop$beta[[k]])) {
      xr <- xrange[[k]][j]
      if (is.null(xr) || !is.finite(xr) || xr <= 0) next
      lo[[paste0("beta_", k, ".", j)]] <- -diff(r) / xr; hi[[paste0("beta_", k, ".", j)]] <- diff(r) / xr
    }
  }
  lo[startsWith(names(v), "logomega_")] <- log(0.05)
  hi[startsWith(names(v), "logomega_")] <- log(config$fit$omega_max)
  lo[startsWith(names(v), "logtau_")] <- log(config$model$baseline$tau_min)
  lo[startsWith(names(v), "logsigma_")] <- log(0.05)
  est <- struct_est_names(spec, config)
  if (length(est)) {
    box <- struct_box(struct_prior_work13(est, spec, config), config)
    lo[paste0("struct_", est)] <- box$lo
    hi[paste0("struct_", est)] <- box$hi
  }
  list(lo = lo, hi = hi)
}

covariate_ranges <- function(pdats) {
  ks <- names(pdats[[1]]$x)
  stats::setNames(lapply(ks, function(k) {
    X <- do.call(rbind, lapply(pdats, function(pd) pd$x[[k]]))
    if (is.null(X) || !ncol(X)) return(numeric(0))
    apply(X, 2, function(x) diff(range(x)))
  }), ks)
}

# trend of every component over a window of population vectors (rows = rounds): slope per round,
# its SE, and whether the trend is material (|slope| * (W - 1) > tol_drift) and significant (> 2 SE)
window_trend <- function(V, tol_drift) {
  W <- nrow(V)
  t <- seq_len(W) - (W + 1) / 2
  b <- colSums(V * t) / sum(t^2)
  res <- V - matrix(colMeans(V), W, ncol(V), byrow = TRUE) - outer(t, b)
  se <- sqrt(colSums(res^2) / (W - 2) / sum(t^2))
  list(slope = b, se = se, drift = abs(b) * (W - 1) > tol_drift & abs(b) > 2 * se)
}

# Outer-loop bookkeeping: damping of oscillating components, over-relaxation of persistent ones
# (switched off when sum(nll) rose, extrapolation capped), the switch to the averaging phase and
# the convergence decision.
# state : list(round, last_d, lambda, same, last_nll, streak, phase, avg_n, history)
# vwin  : population vectors (rows) of the previous window - 1 rounds (oldest first), or NULL
# xrange: covariate ranges for the legal range of covariate effects (pop_limits13)
de_outer_update <- function(state, pop_old, pop_new, sum_nll, spec, config, vwin = NULL, xrange = NULL) {
  o <- config$fit$outer
  v_old <- pop_vector13(pop_old, spec, config)
  v_new <- pop_vector13(pop_new, spec, config)
  d <- v_new - v_old
  lim <- pop_limits13(pop_new, spec, config, xrange)
  reset <- FALSE
  phase <- if (is.null(state$phase)) 1L else state$phase
  rel_obj <- if (is.null(state$last_nll)) Inf else abs(sum_nll - state$last_nll) / abs(state$last_nll)
  streak0 <- if (is.null(state$streak)) 0 else state$streak
  lam <- state$lambda
  if (is.null(lam)) lam <- stats::setNames(rep(1, length(d)), names(d))
  same <- state$same
  if (is.null(same)) same <- stats::setNames(rep(0, length(d)), names(d))
  avg_n <- state$avg_n
  if (phase == 1L) {
    if (!is.null(state$last_d)) {
      flip <- sign(d) != sign(state$last_d[names(d)]) & abs(d) > o$tol_par
      flip[is.na(flip)] <- FALSE
      # consecutive rounds in which a change (> relax_min) kept its direction
      same <- ifelse(abs(d) > (o$relax_min %||% o$tol_par) & sign(d) == sign(state$last_d[names(d)]),
                     same[names(d)] + 1, 0)
      relax <- !flip & same >= (o$relax_after %||% Inf)
      lam[flip] <- pmax(o$damp_min, pmin(lam[flip], 1) * o$damp)
      lam[relax] <- pmin(o$lam_max, pmax(lam[relax], 1) * o$relax_grow)
      keep <- !flip & !relax
      lam[keep] <- pmin(1, lam[keep] / o$damp)
      # the last over-relaxed step made the fit worse: plain (or damped) steps until directions persist again
      reset <- sum_nll - state$last_nll > (o$reset_rise %||% Inf)
      if (reset) {
        lam <- pmin(lam, 1)
        same[] <- 0
      }
    }
    extra <- (lam - 1) * d
    big <- lam > 1
    extra[big] <- pmax(pmin(extra[big], o$max_extra %||% Inf), -(o$max_extra %||% Inf))
    v_upd <- pmin(pmax(v_old + d + extra, lim$lo), lim$hi)
    step <- v_upd - v_old
    # convergence is judged on the change the M step proposes (undamped)
    max_change <- max(abs(d))
    ok <- max_change < o$tol_par && rel_obj < o$tol_obj
    streak <- if (ok) streak0 + 1 else 0
    status <- if (streak >= o$n_ok) "converged" else if (state$round >= o$max_iter) "cap" else "continue"
    W <- if (is.null(o$window)) 0L else o$window
    if (status == "continue" && W > 2 && !is.null(vwin) && nrow(vwin) == W - 1) {
      V <- rbind(vwin[, names(v_upd), drop = FALSE], v_upd)
      if (!any(window_trend(V, o$tol_drift)$drift)) {
        # stationary: averaging phase from the window mean
        phase <- 2L
        avg_n <- W
        v_upd <- colMeans(V)
        step <- v_upd - v_old
        streak <- 0
      }
    }
    pop_upd <- pop_set13(pop_new, v_upd, spec, config)
  } else {
    # averaging phase (stochastic approximation): the population is the running mean of the
    # window and of every later M-step result; converged when the mean is stable
    step <- pmin(pmax(v_old + d / (avg_n + 1), lim$lo), lim$hi) - v_old
    avg_n <- avg_n + 1L
    max_change <- max(abs(step))
    streak <- if (max_change < o$tol_par) streak0 + 1 else 0
    status <- if (streak >= o$n_ok) "converged" else if (state$round >= o$max_iter) "cap" else "continue"
    pop_upd <- pop_set13(pop_new, v_old + step, spec, config)
  }
  hist <- as.data.frame(data.table::rbindlist(list(state$history, data.frame(
    round = state$round, sum_nll = sum_nll, max_change = max_change, rel_obj = rel_obj,
    n_damped = if (phase == 1L) sum(lam < 1) else 0L, which_max = names(d)[which.max(abs(d))],
    max_step = max(abs(step)), phase = phase, n_relaxed = if (phase == 1L) sum(lam > 1) else 0L,
    reset = as.integer(reset), n_at_limit = sum(v_old + step <= lim$lo | v_old + step >= lim$hi))), fill = TRUE))
  list(pop = pop_upd, status = status,
       state = list(round = state$round + 1L, last_d = d, lambda = lam, same = same, last_nll = sum_nll,
                    streak = streak, phase = phase, avg_n = avg_n, history = hist))
}

# ---- run directory (cluster state) ------------------------------------------------------------------------
run_paths <- function(run) {
  list(dir = run, data = file.path(run, "data.rds"), meta = file.path(run, "meta.rds"),
       status = file.path(run, "status"), fit = file.path(run, "fit.rds"),
       pop = function(k) file.path(run, sprintf("pop_%03d.rds", k)),
       estep = function(k) file.path(run, sprintf("estep_%03d", k)),
       final = file.path(run, "final"))
}

de_write_status <- function(run, ...) writeLines(paste(..., sep = " "), run_paths(run)$status)

# d: data object; config$fit / config$model as used for the fit are stored with the run.
# warm: optional start from another run of the same data, list(pop = population file, maps =
# directory of its per-patient E-step results): shared population components are copied
# (random effects, covariate effects, residual SDs, set-point drift, estimated structural
# parameters), the rest keeps its prior; the MAPs seed round 1 (warm/pNNN.rds)
de_init <- function(run, model, d, config, warm = NULL) {
  rp <- run_paths(run)
  dir.create(run, recursive = TRUE, showWarnings = FALSE)
  spec <- model_spec(model)
  pdats <- build_pdats(d, config)
  cn <- lapply(pdats[[1]]$x, names)
  if (isTRUE(spec$f0_nocov)) cn$f0 <- character(0)
  pop <- init_population13(spec, config, d$base_sd, cn)
  if (!is.null(warm)) {
    src <- readRDS(warm$pop)$pop
    for (k in intersect(names(pop$mu), names(src$mu))) {
      pop$mu[[k]] <- src$mu[[k]]
      pop$omega[[k]] <- src$omega[[k]]
      b <- intersect(names(pop$beta[[k]]), names(src$beta[[k]]))
      pop$beta[[k]][b] <- src$beta[[k]][b]
    }
    pop$sigma[intersect(names(pop$sigma), names(src$sigma))] <- src$sigma[intersect(names(pop$sigma), names(src$sigma))]
    pop$tau[intersect(names(pop$tau), names(src$tau))] <- src$tau[intersect(names(pop$tau), names(src$tau))]
    est <- intersect(struct_est_names(spec, config), names(src$struct))
    pop$struct[est] <- src$struct[est]
    dir.create(file.path(run, "warm"), showWarnings = FALSE)
    # matched by patient id (a fold refit has a subset of the source run's patients)
    src_maps <- lapply(list.files(warm$maps, pattern = "^p[0-9]+\\.rds$", full.names = TRUE), readRDS)
    src_ids <- vapply(src_maps, `[[`, character(1), "id")
    for (i in seq_along(pdats)) {
      j <- match(names(pdats)[i], src_ids)
      if (is.na(j)) next
      saveRDS(list(id = src_ids[j], eta = src_maps[[j]]$eta), file.path(run, "warm", sprintf("p%03d.rds", i)))
    }
  }
  saveRDS(d, rp$data)
  saveRDS(list(model = model, n = length(pdats), ids = names(pdats), version = config$version,
               created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), fit = config$fit, model_cfg = config$model,
               warm = warm),
          rp$meta)
  if (is.null(warm)) {
    # round 0: independent per-patient fits under the widened prior initialise the population
    saveRDS(list(pop = pop, state = list(round = 0L, history = NULL)), rp$pop(0L))
    de_write_status(run, "round0", 0L)
  } else {
    saveRDS(list(pop = pop, state = list(round = 1L, history = NULL)), rp$pop(1L))
    de_write_status(run, "estep", 1L)
  }
  invisible(rp)
}

# Round 0 -> population of round 1 (closed-form M-step parts; structural values unchanged)
de_round0_task <- function(run, config) {
  rp <- run_paths(run)
  config <- de_run_config(run, config)
  meta <- readRDS(rp$meta)
  spec <- model_spec(meta$model)
  es <- de_read_estep(rp$estep(0L), meta$n)
  if (!is.null(es$missing)) {
    de_write_status(run, "missing0", 0L, paste(es$missing, collapse = ","))
    return(invisible("missing"))
  }
  fits <- es$fits
  ps <- readRDS(rp$pop(0L))
  pdats <- build_pdats(readRDS(rp$data), config, ids = meta$ids)
  pop <- m_step_re13(pdats, fits, ps$pop, model_re_names(spec), config)
  pop <- m_step_tau(pdats, fits, pop, config)
  pop <- m_step_sigma13(pdats, fits, pop, spec, config)
  diag <- list(estep = data.table::rbindlist(lapply(fits, function(f) data.table::data.table(
    id = f$id, nll = f$nll, conv = f$conv, gen = f$gen, runtime_s = f$runtime_s,
    at_bound = paste(f$at_bound, collapse = ";")))))
  saveRDS(list(pop = pop, state = list(round = 1L, history = NULL), diag = diag, status = "round0"), rp$pop(1L))
  saveRDS(list(pop = pop, fits = fits, diag = diag), file.path(run, "round0.rds"))
  de_write_status(run, "estep", 1L)
  invisible("estep")
}

# Continue a run that stopped at the safety cap under the current outer-loop settings (the
# settings of a run are otherwise fixed at its start): the cap's final E step and fit are kept
# as final_capNNN / fit_capNNN.rds, and the next round is the one after the cap.
de_continue_task <- function(run, config) {
  rp <- run_paths(run)
  pfiles <- sort(list.files(run, pattern = "^pop_[0-9]+\\.rds$", full.names = TRUE))
  ps <- readRDS(utils::tail(pfiles, 1))
  k <- ps$state$round
  if (!identical(ps$status, "cap")) stop("run did not stop at the cap: ", ps$status)
  meta <- readRDS(rp$meta)
  meta$continued <- list(at_round = k, time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), outer_before = meta$fit$outer)
  meta$fit$outer <- config$fit$outer
  saveRDS(meta, rp$meta)
  if (dir.exists(rp$final)) file.rename(rp$final, file.path(run, sprintf("final_cap%03d", k - 1L)))
  if (file.exists(rp$fit)) file.rename(rp$fit, file.path(run, sprintf("fit_cap%03d.rds", k - 1L)))
  de_write_status(run, "continue", k)
  invisible(k)
}

# Move a running run to the current outer-loop settings from its next M step on (earlier
# settings and the round are recorded in meta$reconfigured; atomic write, tasks read meta).
de_reconfig_task <- function(run, config) {
  rp <- run_paths(run)
  meta <- readRDS(rp$meta)
  pfiles <- sort(list.files(run, pattern = "^pop_[0-9]+\\.rds$"))
  meta$reconfigured <- c(meta$reconfigured, list(list(
    at_round = as.integer(sub("pop_([0-9]+)\\.rds", "\\1", utils::tail(pfiles, 1))),
    time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), outer_before = meta$fit$outer)))
  meta$fit$outer <- config$fit$outer
  tmp <- paste0(rp$meta, ".tmp")
  saveRDS(meta, tmp)
  file.rename(tmp, rp$meta)
  invisible(meta$reconfigured)
}

# run settings override the code defaults (a run keeps the settings it was started with)
de_run_config <- function(run, config) {
  meta <- readRDS(run_paths(run)$meta)
  config$fit <- meta$fit
  config$model <- meta$model_cfg
  config
}

de_estep_task <- function(run, k, i, config, cores = 1L, final = FALSE) {
  rp <- run_paths(run)
  config <- de_run_config(run, config)
  d <- readRDS(rp$data)
  meta <- readRDS(rp$meta)
  spec <- model_spec(meta$model)
  pop <- readRDS(rp$pop(k))$pop
  pd <- prep_segments(build_pdats(d, config, ids = meta$ids[i]), pop$struct, config)[[1]]
  prior <- patient_prior(pd, pop, model_re_names(spec))
  box <- indiv_box(prior, config)
  # round 0: the prior SDs widened, the search box of the unwidened prior
  if (k == 0L) prior$sd <- prior$sd * config$fit$round0$sd_mult
  tmpl <- par_template(pop$struct, spec)
  seeds <- list()
  for (kk in rev(seq_len(k) - 1L)) {
    prev <- file.path(rp$estep(kk), sprintf("p%03d.rds", i))
    if (file.exists(prev)) { seeds <- list(readRDS(prev)$eta); break }
  }
  wf <- file.path(run, "warm", sprintf("p%03d.rds", i))
  if (!length(seeds) && file.exists(wf)) {
    w <- readRDS(wf)$eta
    s <- prior$mean
    s[intersect(names(w), names(s))] <- w[intersect(names(w), names(s))]
    seeds <- list(s)
  }
  set.seed(config$runtime$seed + 1000L * k + i + if (final) 500L else 0L)
  f <- de_fit_individual(pd, tmpl, spec, prior, pop$sigma, config, seeds, cores, with_is = final, box = box)
  out_dir <- if (final) rp$final else rp$estep(k)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(f, file.path(out_dir, sprintf("p%03d.rds", i)))
  invisible(f)
}

de_read_estep <- function(dir, n) {
  files <- file.path(dir, sprintf("p%03d.rds", seq_len(n)))
  miss <- which(!file.exists(files))
  if (length(miss)) return(list(missing = miss))
  list(fits = lapply(files, readRDS))
}

de_mstep_task <- function(run, k, config, cores = 1L) {
  rp <- run_paths(run)
  config <- de_run_config(run, config)
  meta <- readRDS(rp$meta)
  spec <- model_spec(meta$model)
  es <- de_read_estep(rp$estep(k), meta$n)
  if (!is.null(es$missing)) {
    de_write_status(run, "missing", k, paste(es$missing, collapse = ","))
    return(invisible("missing"))
  }
  fits <- es$fits
  ps <- readRDS(rp$pop(k))
  d <- readRDS(rp$data)
  pdats <- build_pdats(d, config, ids = meta$ids)
  re <- model_re_names(spec)
  t0 <- Sys.time()
  pop <- m_step_re13(pdats, fits, ps$pop, re, config)
  pop <- m_step_tau(pdats, fits, pop, config)
  pop <- m_step_sigma13(pdats, fits, pop, spec, config)
  set.seed(config$runtime$seed + 7L * k)
  ms <- m_step_struct_de(prep_segments(pdats, pop$struct, config), fits, pop, spec, config, cores)
  sum_nll <- sum(vapply(fits, `[[`, numeric(1), "nll"))
  W <- if (is.null(config$fit$outer$window)) 0L else config$fit$outer$window
  vwin <- if (W > 2 && k >= W - 1L) do.call(rbind, lapply(seq(k - W + 2L, k), function(j)
    pop_vector13(readRDS(rp$pop(j))$pop, spec, config))) else NULL
  up <- de_outer_update(ps$state, ps$pop, ms$pop, sum_nll, spec, config, vwin, covariate_ranges(pdats))
  diag <- list(estep = data.table::rbindlist(lapply(fits, function(f) data.table::data.table(
    id = f$id, nll = f$nll, conv = f$conv, gen = f$gen, runtime_s = f$runtime_s,
    at_bound = paste(f$at_bound, collapse = ";")))), struct = ms$diag,
    mstep_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(list(pop = up$pop, state = up$state, diag = diag, status = up$status), rp$pop(k + 1L))
  utils::write.csv(up$state$history, file.path(run, "trace.csv"), row.names = FALSE)
  de_write_status(run, up$status, k + 1L)
  invisible(up$status)
}

de_assemble_task <- function(run, config) {
  rp <- run_paths(run)
  config <- de_run_config(run, config)
  meta <- readRDS(rp$meta)
  spec <- model_spec(meta$model)
  es <- de_read_estep(rp$final, meta$n)
  if (!is.null(es$missing)) {
    de_write_status(run, "missing_final", paste(es$missing, collapse = ","))
    return(invisible("missing"))
  }
  fits <- es$fits
  pfiles <- sort(list.files(run, pattern = "^pop_[0-9]+\\.rds$", full.names = TRUE))
  ps <- readRDS(utils::tail(pfiles, 1))
  pop <- ps$pop
  re <- model_re_names(spec)
  ids <- meta$ids
  nm <- names(fits[[1]]$eta)
  indiv <- data.table::data.table(patient_id = ids)
  for (k in nm) {
    col <- if (k %in% re) paste0("eta_", k) else k
    indiv[, (col) := vapply(fits, function(f) f$eta[[k]], numeric(1))]
    indiv[, (paste0("se_", k)) := vapply(fits, function(f) sqrt(f$cov[k, k]), numeric(1))]
  }
  ess_min <- if (is.null(config$fit$is$ess_min)) 30 else config$fit$is$ess_min
  is_ok <- vapply(fits, function(f) is.finite(f$is[["logml"]]) && is.finite(f$is[["ess"]]) &&
                    f$is[["ess"]] >= ess_min, logical(1))
  indiv[, `:=`(logL_laplace = vapply(fits, `[[`, numeric(1), "logL"),
               logL = vapply(seq_along(fits), function(i) if (is_ok[i]) fits[[i]]$is[["logml"]] else fits[[i]]$logL,
                             numeric(1)),
               logL_core = vapply(seq_along(fits), function(i) if (is_ok[i]) fits[[i]]$is[["logml_core"]] else
                 if (is.null(fits[[i]]$logL_core_laplace)) NA_real_ else fits[[i]]$logL_core_laplace, numeric(1)),
               logL_source = ifelse(is_ok, "IS", "Laplace"),
               is_ess = vapply(fits, function(f) f$is[["ess"]], numeric(1)),
               nll = vapply(fits, `[[`, numeric(1), "nll"),
               conv = vapply(fits, `[[`, integer(1), "conv"),
               gen = vapply(fits, `[[`, integer(1), "gen"),
               at_bound = vapply(fits, function(f) paste(f$at_bound, collapse = ";"), character(1)))]
  n_par <- count_pop_params13(pop, spec, config)
  logLik <- sum(indiv$logL)
  logLik_core <- sum(indiv$logL_core, na.rm = TRUE)
  n <- length(ids)
  fit <- list(model = meta$model, spec = spec, re = re, pop = pop, indiv = indiv,
              cov = stats::setNames(lapply(fits, `[[`, "cov"), ids),
              logLik = logLik, logLik_laplace = sum(indiv$logL_laplace), logLik_core = logLik_core,
              n_par = n_par, n = n, BIC = -2 * logLik + n_par * log(n), BIC_core = -2 * logLik_core + n_par * log(n),
              trace = ps$state$history, rounds = ps$state$round - 1L,
              converged = identical(ps$status, "converged"),
              converged_by = if (!identical(ps$status, "converged")) NA_character_ else
                if (identical(ps$state$phase, 2L)) "averaging" else "change",
              warm = meta$warm, continued = meta$continued,
              is_fallback = sum(!is_ok), version = meta$version, method = "DE-ITS")
  saveRDS(fit, rp$fit)
  de_write_status(run, "done")
  invisible(fit)
}

# ---- local driver (tests and small runs; same task functions as the cluster) ------------------------------
de_fit_local <- function(run, model, d, config, cores = 1L, log_fun = message, warm = NULL) {
  de_init(run, model, d, config, warm)
  n <- readRDS(run_paths(run)$meta)$n
  if (is.null(warm)) {
    for (i in seq_len(n)) de_estep_task(run, 0L, i, config, cores)
    de_round0_task(run, config)
    log_fun(sprintf("    %s round 0: independent fits done", model))
  }
  k <- 1L
  repeat {
    for (i in seq_len(n)) de_estep_task(run, k, i, config, cores)
    st <- de_mstep_task(run, k, config, cores)
    tr <- utils::read.csv(file.path(run, "trace.csv"))
    log_fun(sprintf("    %s round %d: sum nll %.2f, max change %.4f (%s)", model, k,
                    utils::tail(tr$sum_nll, 1), utils::tail(tr$max_change, 1), st))
    k <- k + 1L
    if (st != "continue") break
  }
  for (i in seq_len(n)) de_estep_task(run, k, i, config, cores, final = TRUE)
  de_assemble_task(run, config)
}
