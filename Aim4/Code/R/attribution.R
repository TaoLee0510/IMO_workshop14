# =============================================================================
# attribution.R
# Purpose : Counterfactual attribution of each patient's recovery phenotype
#           from posterior draws, with two readouts:
#             when (Q2): a-driven / b-driven / a+b joint / a+b redundant /
#                        none (quick) / other
#             who  (Q1): clone / reserve / both / other (non-quick patients)
#           plus continuous readouts f_a, f_b (HSC loss attributable to a and
#           b) and the clone fraction on day 30.
# Inputs  : draws (mcmc), pdat, fit, config.
# Outputs : per-patient probabilities and per-draw table.
#
# Counterfactuals (each simulated per draw and classified by the same rule)
#   actual | no_a: th_aS = 0 | no_bS: th_bS = 0 | no_abS: both (HSC level) |
#   no_b: D = 0 and th_bS = 0 (whole CAR T phase, incl. the reversible
#   suppression formerly H2) | no_ab: th_aS = 0 and no CAR T phase |
#   no_clone: f0 = 0 | reserve: r_S = population median
# when       = HSC-level readout (primary; matches the Q2 model comparison)
# when_total = CAR T phase removed entirely (shows the reversible part, ex-H2)
# "Improves" = phenotype rank decreases (delayed > second wave > quick;
# persistent > delayed when the late class is assessed).
# =============================================================================

when_class <- function(act, no_a, no_b, no_ab) {
  ra <- pheno_rank(act)
  ia <- pheno_rank(no_a) < ra
  ib <- pheno_rank(no_b) < ra
  iab <- pheno_rank(no_ab) < ra
  ifelse(act == "quick", "none",
         ifelse(ia & !ib, "a",
                ifelse(ib & !ia, "b",
                       ifelse(ia & ib, "ab_joint",
                              ifelse(iab, "ab_redundant", "other")))))
}

who_class <- function(act, no_clone, reserve) {
  ra <- pheno_rank(act)
  ic <- pheno_rank(no_clone) < ra
  ir <- pheno_rank(reserve) < ra
  ifelse(act == "quick", NA_character_,
         ifelse(ic & ir, "both", ifelse(ic, "clone", ifelse(ir, "reserve", "other"))))
}

trapz <- function(x, y) sum(diff(x) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# Simulate the counterfactual set for one draw.
# days: daily grid (0..end); with_persistent: classify the late class too
counterfactual_draw <- function(eta, pd, fit, config, days, with_persistent = FALSE,
                                obs_upto = NULL) {
  spec <- fit$spec
  st <- fit$pop$struct
  ind <- eta_to_indiv(eta, spec)
  all_days <- sort(unique(c(days, pd$times[pd$times > pd$design$t0])))
  all_days <- all_days[all_days > pd$design$t0]
  seg <- make_segments(pd$design, st, max(all_days), config)
  zero <- zero_driver(pd$drv$t0, config)
  zero$v <- zero$v[seq_along(pd$drv$v)]
  run <- function(i2, drv) {
    m <- sim_model(build_par(i2, st, spec), seg, drv, all_days, spec, config)
    m[all_days %in% days, , drop = FALSE]
  }
  cls <- function(m) {
    anc <- m[, "anc_rel"] * pd$base$anc
    if (!is.null(obs_upto)) return(phenotype_with_observed(m[, "t"], anc, obs_upto, config, with_persistent))
    phenotype_sim(m[, "t"], anc, config, with_persistent = with_persistent)
  }
  i_a <- ind; i_a[["th_as"]] <- 0
  i_bs <- ind; i_bs[["th_bs"]] <- 0
  i_c <- ind; i_c[["f0"]] <- 0
  i_r <- ind
  if ("r_s" %in% names(fit$pop$mu)) i_r[["r_s"]] <- exp(fit$pop$mu[["r_s"]])
  m_act <- run(ind, pd$drv)
  if (anyNA(m_act[, 2])) return(NULL)
  i_abs <- i_a; i_abs[["th_bs"]] <- 0
  m_noa <- run(i_a, pd$drv)
  m_nobs <- run(i_bs, pd$drv)
  m_noabs <- run(i_abs, pd$drv)
  m_nob <- run(i_bs, zero)
  m_noab <- run(i_abs, zero)
  m_c <- if (ind[["f0"]] > 0) run(i_c, pd$drv) else m_act
  m_r <- run(i_r, pd$drv)
  # HSC loss attributable to each phase (Aim 3 continuous readouts)
  tt <- m_act[, "t"]
  f_a <- 1 - exp(-ind[["th_as"]] * trapz(tt, m_act[, "C"]))
  f_b <- 1 - exp(-ind[["th_bs"]] * trapz(tt[tt >= 0], m_act[tt >= 0, "I"]))
  i30 <- which.min(abs(tt - 30))
  out <- c(cls(m_act), cls(m_noa), cls(m_nobs), cls(m_noabs), cls(m_nob), cls(m_noab), cls(m_c), cls(m_r),
           unname(f_a), unname(f_b), unname(m_act[i30, "f_clone"]), unname(m_act[i30, "S"]))
  names(out) <- c("actual", "no_a", "no_bS", "no_abS", "no_b", "no_ab", "no_clone", "reserve",
                  "f_a", "f_b", "f30", "S30")
  out
}

# Phenotype of [observed ANC up to the landmark] + [model after it] (landmark forecasts).
# obs_upto: list(day, value, cens, L)
phenotype_with_observed <- function(day, anc, obs_upto, config, with_persistent = FALSE) {
  keep <- day > obs_upto$L
  d2 <- c(obs_upto$day, day[keep])
  a2 <- c(obs_upto$value, anc[keep])
  c2 <- c(obs_upto$cens, rep(FALSE, sum(keep)))
  phenotype_series(d2, a2, c2, phenotype_rule(config),
                   if (with_persistent) config$phenotype$persistent else NULL)
}

attribute_patient <- function(draws, pd, fit, config, n_draws = config$fit$attrib_draws,
                              days = 0:30, with_persistent = FALSE, seed = 1L) {
  set.seed(seed)
  ix <- if (nrow(draws) > n_draws) sample(nrow(draws), n_draws) else seq_len(nrow(draws))
  res <- lapply(ix, function(j) counterfactual_draw(stats::setNames(draws[j, ], colnames(draws)), pd, fit,
                                                    config, days, with_persistent))
  res <- res[!vapply(res, is.null, logical(1))]
  tab <- data.table::as.data.table(do.call(rbind, res))
  for (v in c("f_a", "f_b", "f30", "S30")) tab[, (v) := as.numeric(get(v))]
  tab[, when := when_class(actual, no_a, no_bS, no_abS)]          # HSC level (primary)
  tab[, when_total := when_class(actual, no_a, no_b, no_ab)]      # whole CAR T phase incl. reversible suppression
  tab[, who := who_class(actual, no_clone, reserve)]
  tab[, patient_id := pd$id]
  tab
}

summarise_attribution <- function(tab, config) {
  cls <- if (any(tab$actual == "persistent")) config$phenotype$classes4 else config$phenotype$classes
  p_ph <- vapply(cls, function(k) mean(tab$actual == k), numeric(1))
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  p_when <- vapply(wc, function(k) mean(tab$when == k), numeric(1))
  p_when_tot <- vapply(wc, function(k) mean(tab$when_total == k), numeric(1))
  hc <- c("clone", "reserve", "both", "other")
  nq <- tab[actual != "quick"]
  p_who <- vapply(hc, function(k) if (nrow(nq)) mean(nq$who == k) else NA_real_, numeric(1))
  q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  data.table::data.table(
    patient_id = tab$patient_id[1], n_draws = nrow(tab),
    t(stats::setNames(p_ph, paste0("p_", cls))),
    t(stats::setNames(p_when, paste0("p_when_", wc))),
    t(stats::setNames(p_when_tot, paste0("p_whentot_", wc))),
    t(stats::setNames(p_who, paste0("p_who_", hc))),
    f_a_lo = q(tab$f_a)[1], f_a = q(tab$f_a)[2], f_a_hi = q(tab$f_a)[3],
    f_b_lo = q(tab$f_b)[1], f_b = q(tab$f_b)[2], f_b_hi = q(tab$f_b)[3],
    f30 = stats::median(tab$f30), S30 = stats::median(tab$S30))
}
