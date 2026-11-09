# =============================================================================
# fit_mcmc.R
# Purpose : Individual posteriors by adaptive Metropolis given the population
#           parameters (Laplace posteriors are too narrow and push attribution
#           probabilities to 0 or 1 - Aim 2 D3c).
# Inputs  : pdat, population fit (fit_its.R), config.
# Outputs : list(draws = matrix (rows = draws, cols = random effects), accept,
#           rhat)
# =============================================================================

mcmc_patient <- function(pd, fit, config, eta_map = NULL, cov_map = NULL, mc = config$fit$mcmc,
                         seed = 1L) {
  spec <- fit$spec
  re <- fit$re
  pop <- fit$pop
  tmpl <- par_template(pop$struct, spec)
  if (is.null(pd$seg)) pd$seg <- make_segments(pd$design, pop$struct, max(pd$times), config)
  prior <- patient_prior(pd, pop, re)
  lp <- function(e) -neg_log_post(stats::setNames(e, re), pd, tmpl, spec, prior, pop$sigma, config)
  if (is.null(eta_map)) {
    f <- fit_individual(pd, tmpl, spec, prior, pop$sigma, config, stats::setNames(prior$mean, re),
                        clone_starts(stats::setNames(prior$mean, re), spec))
    eta_map <- f$eta
    cov_map <- f$cov
  }
  k <- length(re)
  set.seed(seed)
  base_cov <- (mc$scale^2 / k) * (cov_map + diag(1e-6, k))
  chains <- lapply(seq_len(mc$n_chains), function(ch) {
    x <- as.numeric(eta_map) + if (ch > 1) as.numeric(MASS::mvrnorm(1, rep(0, k), cov_map)) else 0
    lx <- lp(x)
    if (!is.finite(lx) || lx < -1e9) { x <- as.numeric(eta_map); lx <- lp(x) }
    prop <- base_cov
    out <- matrix(NA_real_, mc$n_iter, k)
    acc <- 0
    for (it in seq_len(mc$n_iter)) {
      y <- x + as.numeric(MASS::mvrnorm(1, rep(0, k), prop))
      ly <- lp(y)
      if (log(stats::runif(1)) < ly - lx) { x <- y; lx <- ly; acc <- acc + 1 }
      out[it, ] <- x
      # adapt during burn-in (Haario et al. 2001)
      if (it <= mc$burn && it >= 100 && it %% 50 == 0) {
        cv <- stats::cov(out[max(1, it - 300):it, , drop = FALSE])
        if (all(is.finite(cv))) prop <- (mc$scale^2 / k) * (cv + diag(1e-6, k))
      }
    }
    keep <- seq(mc$burn + 1, mc$n_iter, by = mc$thin)
    list(draws = out[keep, , drop = FALSE], accept = acc / mc$n_iter)
  })
  draws <- do.call(rbind, lapply(chains, `[[`, "draws"))
  colnames(draws) <- re
  rhat <- if (mc$n_chains > 1) vapply(seq_len(k), function(j)
    gelman_rhat(lapply(chains, function(c) c$draws[, j])), numeric(1)) else rep(NA_real_, k)
  list(draws = draws, accept = vapply(chains, `[[`, numeric(1), "accept"),
       rhat = stats::setNames(rhat, re), eta_map = eta_map)
}

gelman_rhat <- function(chs) {
  n <- min(lengths(chs))
  chs <- lapply(chs, function(x) x[seq_len(n)])
  m <- length(chs)
  means <- vapply(chs, mean, numeric(1))
  W <- mean(vapply(chs, stats::var, numeric(1)))
  B <- n * stats::var(means)
  v <- (n - 1) / n * W + B / n
  if (W <= 0) return(NA_real_)
  sqrt(v / W)
}

# MCMC for every patient of a fit (parallel over patients).
mcmc_all <- function(pdats, fit, config, cores = config$runtime$cores, seed = 1L) {
  pdats <- prep_segments(pdats, fit$pop$struct, config)
  ids <- names(pdats)
  out <- par_lapply(seq_along(ids), function(i) {
    eta_map <- stats::setNames(as.numeric(fit$indiv[patient_id == ids[i], paste0("eta_", fit$re), with = FALSE]), fit$re)
    cv <- fit$cov[[ids[i]]]
    mcmc_patient(pdats[[i]], fit, config, eta_map, cv, seed = seed + i)
  }, cores)
  stats::setNames(out, ids)
}
