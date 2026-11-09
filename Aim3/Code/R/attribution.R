# =============================================================================
# attribution.R
# Purpose : Individual posteriors by random-walk Metropolis (given the fitted
#           population), and counterfactual attribution of each patient's
#           recovery phenotype to chemotherapy (a) and/or the CAR T phase (b).
# Inputs  : pdat, population fit (fit_its.R), config.
# Outputs : posterior draws of eta; per-patient P(mechanism class) and the
#           continuous read-outs f_a, f_b.
#
# Attribution rule (Aim3_plan.md 4.7); "improves" = phenotype at least one
# level better (delayed -> second wave -> quick)
#   none          actual phenotype is quick
#   a             removing a improves, removing b does not
#   b             removing b improves, removing a does not
#   ab_joint      removing either improves
#   ab_redundant  neither alone improves, removing both does
#   other         removing both does not improve
# =============================================================================

mech_classes <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
pheno_rank <- c(quick = 1, second_wave = 2, delayed = 3)

# Random-walk Metropolis started at the MAP, proposal from the Laplace
# covariance, adapted once from the burn-in sample.
mcmc_individual <- function(pd, tmpl, prior, sigma, config, fit_i, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ctl <- config$fit$mcmc
  re <- names(fit_i$eta)
  k <- length(re)
  lp <- function(e) -neg_log_post(stats::setNames(e, re), pd, tmpl, prior, sigma, config)
  cur <- fit_i$eta
  lcur <- lp(cur)
  # proposal: Laplace covariance with eigenvalues capped at the prior variance
  ev <- eigen((fit_i$cov + t(fit_i$cov)) / 2, symmetric = TRUE)
  ev$values <- pmin(pmax(ev$values, 1e-6), max(prior$sd)^2)
  base <- ev$vectors %*% diag(ev$values, k) %*% t(ev$vectors)
  scale <- ctl$scale^2 / k
  L <- t(chol(base * scale + diag(1e-10, k)))
  keep <- matrix(NA_real_, 0, k, dimnames = list(NULL, re))
  burn_store <- matrix(NA_real_, ctl$burn, k)
  acc <- 0
  acc_win <- 0
  for (it in seq_len(ctl$n_iter)) {
    prop <- cur + as.vector(L %*% stats::rnorm(k))
    lprop <- lp(prop)
    if (log(stats::runif(1)) < lprop - lcur) {
      cur <- prop
      lcur <- lprop
      acc <- acc + 1
      acc_win <- acc_win + 1
    }
    if (it <= ctl$burn) {
      burn_store[it, ] <- cur
      # step-size adaptation during burn-in (target acceptance 0.2-0.4)
      if (it %% 50 == 0) {
        r <- acc_win / 50
        if (r < 0.15) scale <- scale * 0.5 else if (r > 0.45) scale <- scale * 1.6
        acc_win <- 0
        L <- t(chol(base * scale + diag(1e-10, k)))
      }
      if (it == ctl$burn) {
        S <- stats::cov(burn_store[(ctl$burn %/% 2):ctl$burn, , drop = FALSE])
        if (all(is.finite(S)) && min(eigen(S, only.values = TRUE)$values) > 1e-6)
          L <- t(chol(S * ctl$scale^2 / k + diag(1e-10, k)))
        acc <- 0
      }
    } else if ((it - ctl$burn) %% ctl$thin == 0) {
      keep <- rbind(keep, cur)
    }
  }
  attr(keep, "accept") <- acc / (ctl$n_iter - ctl$burn)
  keep
}

# Phenotypes of the four counterfactual worlds for one natural-scale parameter vector.
counterfactual_phenotypes <- function(par, pd, config) {
  worlds <- list(actual = par,
                 no_a = replace(par, "th_as", 0),
                 no_b = replace(par, "th_bs", 0),
                 no_ab = replace(par, c("th_as", "th_bs"), 0),
                 no_bp = replace(par, c("q_b", "eps"), 0))
  vapply(worlds, function(p) sim_phenotype(p, pd, config), character(1))
}

mech_class <- function(ph) {
  if (anyNA(ph)) return(NA_character_)
  r <- stats::setNames(pheno_rank[ph], names(ph))
  if (r[["actual"]] == 1) return("none")
  imp_a <- r[["no_a"]] < r[["actual"]]
  imp_b <- r[["no_b"]] < r[["actual"]]
  if (imp_a && imp_b) return("ab_joint")
  if (imp_a) return("a")
  if (imp_b) return("b")
  if (r[["no_ab"]] < r[["actual"]]) "ab_redundant" else "other"
}

# Attribution over posterior draws (rows of eta) for one patient.
attribute_patient <- function(draws, pd, tmpl, config, n_use = config$fit$attrib_draws) {
  idx <- if (nrow(draws) > n_use) round(seq(1, nrow(draws), length.out = n_use)) else seq_len(nrow(draws))
  res <- lapply(idx, function(j) {
    par <- fill_par(tmpl, draws[j, ])
    ph <- counterfactual_phenotypes(par, pd, config)
    if (anyNA(ph)) return(c(pheno = NA, mech = NA, f_a = NA, f_b = NA, bp = NA))
    kf <- kill_fractions(par, pd, config)
    c(pheno = ph[["actual"]], mech = mech_class(ph), f_a = kf[["f_a"]], f_b = kf[["f_b"]],
      bp = as.numeric(pheno_rank[[ph[["no_bp"]]]] < pheno_rank[[ph[["actual"]]]]))
  })
  m <- do.call(rbind, res)
  mech <- factor(m[, "mech"], levels = mech_classes)
  ph <- factor(m[, "pheno"], levels = names(pheno_rank))
  fa <- as.numeric(m[, "f_a"])
  fb <- as.numeric(m[, "f_b"])
  out <- data.table::data.table(patient_id = pd$id, n_draws = length(idx))
  for (k in mech_classes) out[, (paste0("p_", k)) := mean(mech == k, na.rm = TRUE)]
  for (k in names(pheno_rank)) out[, (paste0("pp_", k)) := mean(ph == k, na.rm = TRUE)]
  # removing the CAR T-phase progenitor-level effects (suppression q_b and
  # accelerated release eps) improves the phenotype
  out[, p_bp := mean(as.numeric(m[, "bp"]), na.rm = TRUE)]
  out[, `:=`(f_a = stats::median(fa, na.rm = TRUE), f_a_lo = stats::quantile(fa, 0.1, na.rm = TRUE),
             f_a_hi = stats::quantile(fa, 0.9, na.rm = TRUE), f_b = stats::median(fb, na.rm = TRUE),
             f_b_lo = stats::quantile(fb, 0.1, na.rm = TRUE), f_b_hi = stats::quantile(fb, 0.9, na.rm = TRUE))]
  out[, map_class := mech_classes[which.max(unlist(.SD))], .SDcols = paste0("p_", mech_classes)]
  out[]
}

# MCMC + attribution for every patient of a population fit.
attribute_cohort <- function(pdats, fit, config, cores = config$runtime$cores, seed = 1L) {
  re <- fit$re
  tmpl <- par_template(fit$pop$struct, fit$model)
  pdats <- refresh_pdats(pdats, fit$pop$struct, config)
  rows <- par_lapply(seq_along(pdats), function(i) {
    pd <- pdats[[i]]
    prior <- patient_prior(pd, fit$pop, re)
    draws <- mcmc_individual(pd, tmpl, prior, fit$pop$sigma, config, fit$fits[[i]], seed = seed + i)
    a <- attribute_patient(draws, pd, tmpl, config)
    a[, accept := attr(draws, "accept")]
    a
  }, cores)
  data.table::rbindlist(rows)
}

# True mechanism class of a virtual patient (known parameters).
true_mech_class <- function(par, pd, config) mech_class(counterfactual_phenotypes(par, pd, config))
