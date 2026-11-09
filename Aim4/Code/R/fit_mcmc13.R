# =============================================================================
# fit_mcmc13.R
# Purpose : v1.3 individual posteriors by adaptive Metropolis given the
#           population of a DE-ITS fit: the individual vector psi (random
#           effects + baseline deviations) is confined to the search box (the
#           support of the truncated prior), and chains are extended until the
#           largest R-hat is below the tolerance (safety cap max_iter).
# Inputs  : pdat, fit (fit_de.R), psi_hat and its Laplace covariance, config.
# Outputs : list(draws, accept, rhat, n_iter, converged)
# =============================================================================

mcmc_patient13 <- function(pd, fit, config, psi_map, cov_map, mc = config$fit$mcmc, seed = 1L) {
  spec <- fit$spec
  pop <- fit$pop
  tmpl <- par_template(pop$struct, spec)
  if (is.null(pd$seg)) pd$seg <- make_segments(pd$design, pop$struct, max(pd$times), config)
  prior <- patient_prior(pd, pop, fit$re)
  nm <- names(prior$mean)
  box <- indiv_box(prior, config)
  psi_map <- pmin(pmax(psi_map[nm], box$lo), box$hi)
  lp <- function(e) {
    if (any(e < box$lo | e > box$hi)) return(-Inf)
    -neg_log_post(stats::setNames(e, nm), pd, tmpl, spec, prior, pop$sigma, config)
  }
  k <- length(nm)
  S0 <- solve(pd_floor(solve(pd_floor(cov_map[nm, nm], 1e-8)), 1 / max(prior$sd^2)))
  set.seed(seed)
  step <- function(st, n, adapt) {
    out <- matrix(NA_real_, n, k)
    for (it in seq_len(n)) {
      y <- st$x + as.numeric(MASS::mvrnorm(1, rep(0, k), st$prop))
      ly <- lp(y)
      if (log(stats::runif(1)) < ly - st$lx) { st$x <- y; st$lx <- ly; st$acc <- st$acc + 1 }
      out[it, ] <- st$x
      st$n <- st$n + 1
      if (adapt && st$n >= 100 && st$n %% 50 == 0) {
        hist <- rbind(st$hist, out[seq_len(it), , drop = FALSE])
        cv <- stats::cov(utils::tail(hist, 300))
        if (all(is.finite(cv))) st$prop <- (mc$scale^2 / k) * (cv + diag(1e-6, k))
      }
    }
    st$hist <- rbind(st$hist, out)
    st
  }
  chains <- lapply(seq_len(mc$n_chains), function(ch) {
    x <- psi_map
    if (ch > 1) x <- pmin(pmax(psi_map + as.numeric(MASS::mvrnorm(1, rep(0, k), S0)), box$lo), box$hi)
    lx <- lp(x)
    if (!is.finite(lx)) { x <- psi_map; lx <- lp(x) }
    st <- list(x = x, lx = lx, prop = (mc$scale^2 / k) * (S0 + diag(1e-6, k)), acc = 0, n = 0, hist = NULL)
    st <- step(st, mc$burn, adapt = TRUE)
    step(st, mc$n_iter - mc$burn, adapt = FALSE)
  })
  kept <- function(st) {
    h <- st$hist[-seq_len(mc$burn), , drop = FALSE]
    h[seq(1, nrow(h), by = mc$thin), , drop = FALSE]
  }
  rhat_all <- function() if (mc$n_chains > 1) vapply(seq_len(k), function(j)
    gelman_rhat(lapply(chains, function(c) kept(c)[, j])), numeric(1)) else rep(NA_real_, k)
  rh <- rhat_all()
  total <- mc$n_iter
  ext <- if (is.null(mc$extend)) 0L else mc$extend
  while (mc$n_chains > 1 && any(rh > mc$rhat_tol, na.rm = TRUE) && ext > 0 && total + ext <= mc$max_iter) {
    chains <- lapply(chains, function(st) step(st, ext, adapt = FALSE))
    total <- total + ext
    rh <- rhat_all()
  }
  draws <- do.call(rbind, lapply(chains, kept))
  colnames(draws) <- nm
  list(draws = draws, accept = vapply(chains, function(c) c$acc / c$n, numeric(1)),
       rhat = stats::setNames(rh, nm), n_iter = total,
       converged = mc$n_chains < 2 || all(rh <= mc$rhat_tol, na.rm = TRUE))
}

# MCMC for every patient of a v1.3 fit (parallel over patients)
mcmc_all13 <- function(pdats, fit, config, cores = config$runtime$cores, seed = 1L) {
  pdats <- prep_segments(pdats, fit$pop$struct, config)
  psi <- fit_psi13(fit)
  ids <- names(pdats)
  out <- par_lapply(seq_along(ids), function(i) {
    mcmc_patient13(pdats[[i]], fit, config, psi[ids[i], ], fit$cov[[ids[i]]], seed = seed + i)
  }, cores)
  stats::setNames(out, ids)
}
