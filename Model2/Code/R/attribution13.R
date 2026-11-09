# =============================================================================
# attribution13.R
# Purpose : v1.3 counterfactual attribution per posterior draw. As attribution.R
#           (v1.2), with the CAR T phase (b) split into its two HSC routes:
#           the ferritin-driven inflammation I (th_bS) and the CAR T - target
#           signal Z (th_zS); individual baselines scale the ANC that is
#           classified; exposures are integrated on a 0.25-day grid from the
#           start of the simulation (incl. lymphodepletion).
# Inputs  : psi draw, pdat, fit (fit_de.R), config.
# Outputs : per-draw classes and per-patient probabilities.
#
# Counterfactuals (each simulated per draw and classified by the same rule)
#   actual | no_a: th_aS = 0 | no_bI: th_bS = 0 | no_bZ: th_zS = 0 | no_bS: both |
#   no_abS: th_aS = th_bS = th_zS = 0 | no_b: no CAR T phase at all (no ferritin driver,
#   no infusion; th_bS = th_zS = 0) | no_ab: no_b and th_aS = 0 | no_clone: f0 = 0 |
#   reserve: r_S = population median
# when       = when_class(actual, no_a, no_bS, no_abS)      (HSC level, b = both routes)
# when_total = when_class(actual, no_a, no_b, no_ab)        (whole CAR T phase)
# b_route    = for b-involved draws: ferritin (only no_bI improves), cart_target (only
#              no_bZ improves), both (each alone improves), joint (only removing both)
# =============================================================================

b_route_class <- function(act, when, no_bi, no_bz) {
  ra <- pheno_rank(act)
  ii <- pheno_rank(no_bi) < ra
  iz <- pheno_rank(no_bz) < ra
  ifelse(!when %in% c("b", "ab_joint"), NA_character_,
         ifelse(ii & iz, "both", ifelse(ii, "ferritin", ifelse(iz, "cart_target", "joint"))))
}

counterfactual_draw13 <- function(psi, pd, fit, config, days = 0:30, obs_upto = NULL, drv = pd$drv) {
  spec <- fit$spec
  st <- fit$pop$struct
  ind <- eta_to_indiv(psi, spec)
  ind[["a0"]] <- patient_a0(st, pd$ldh_c)
  b <- pdat_baselines(psi, pd)
  t_end <- max(days)
  fine <- seq(ceiling(pd$design$t0) + 1, t_end, by = 0.25)
  zero <- zero_driver(drv$t0, config)
  zero$v <- zero$v[seq_along(drv$v)]
  run <- function(i2, s2 = st, dr = drv, no_car = FALSE, grid = days) {
    des <- pd$design
    des$no_car <- no_car
    seg <- make_segments(des, s2, t_end, config)
    sim_model(build_par(i2, s2, spec), seg, dr, grid, spec, config)
  }
  cls <- function(m) {
    anc <- m[, "anc_rel"] * b[["anc"]]
    if (!is.null(obs_upto)) return(phenotype_with_observed(m[, "t"], anc, obs_upto, config))
    phenotype_sim(m[, "t"], anc, config)
  }
  st0 <- st
  st0[["th_zs"]] <- 0
  i_a <- ind; i_a[["th_as"]] <- 0
  i_bi <- ind; i_bi[["th_bs"]] <- 0
  i_abs <- i_a; i_abs[["th_bs"]] <- 0
  i_c <- ind; i_c[["f0"]] <- 0
  i_r <- ind
  if ("r_s" %in% names(fit$pop$mu)) i_r[["r_s"]] <- exp(fit$pop$mu[["r_s"]])
  m_fine <- run(ind, grid = fine)
  if (anyNA(m_fine[, 2])) return(NULL)
  m_act <- m_fine[m_fine[, "t"] %in% days, , drop = FALSE]
  tt <- m_fine[, "t"]
  pos <- tt >= 0
  f_a <- 1 - exp(-ind[["th_as"]] * trapz(tt, m_fine[, "C"]))
  f_bi <- 1 - exp(-ind[["th_bs"]] * trapz(tt[pos], m_fine[pos, "I"]))
  f_bz <- 1 - exp(-(if (isTRUE(spec$z_hsc)) st[["th_zs"]] else 0) * trapz(tt[pos], m_fine[pos, "z"]))
  i30 <- which.min(abs(tt - 30))
  out <- c(actual = cls(m_act), no_a = cls(run(i_a)), no_bI = cls(run(i_bi)), no_bZ = cls(run(ind, st0)),
           no_bS = cls(run(i_bi, st0)), no_abS = cls(run(i_abs, st0)),
           no_b = cls(run(i_bi, st0, zero, TRUE)), no_ab = cls(run(i_abs, st0, zero, TRUE)),
           no_clone = if (ind[["f0"]] > 0) cls(run(i_c)) else cls(m_act), reserve = cls(run(i_r)),
           f_a = f_a, f_bi = f_bi, f_bz = f_bz, f30 = unname(m_fine[i30, "f_clone"]), S30 = unname(m_fine[i30, "S"]),
           peak_cart = max(m_fine[pos, "cart"]))
  out
}

attribute_patient13 <- function(draws, pd, fit, config, n_draws = config$fit$attrib_draws, days = 0:30, seed = 1L) {
  set.seed(seed)
  ix <- if (nrow(draws) > n_draws) sample(nrow(draws), n_draws) else seq_len(nrow(draws))
  res <- lapply(ix, function(j) tryCatch(counterfactual_draw13(stats::setNames(draws[j, ], colnames(draws)), pd, fit,
                                                               config, days), error = function(e) NULL))
  n_failed <- sum(vapply(res, is.null, logical(1)))
  res <- res[!vapply(res, is.null, logical(1))]
  tab <- data.table::as.data.table(do.call(rbind, res))
  for (v in c("f_a", "f_bi", "f_bz", "f30", "S30", "peak_cart")) tab[, (v) := as.numeric(get(v))]
  tab[, when := when_class(actual, no_a, no_bS, no_abS)]
  tab[, when_total := when_class(actual, no_a, no_b, no_ab)]
  tab[, b_route := b_route_class(actual, when, no_bI, no_bZ)]
  tab[, who := who_class(actual, no_clone, reserve)]
  tab[, patient_id := pd$id]
  attr(tab, "n_failed") <- n_failed
  tab
}

summarise_attribution13 <- function(tab, config) {
  cls <- config$phenotype$classes
  p_ph <- vapply(cls, function(k) mean(tab$actual == k), numeric(1))
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  p_when <- vapply(wc, function(k) mean(tab$when == k), numeric(1))
  p_when_tot <- vapply(wc, function(k) mean(tab$when_total == k), numeric(1))
  rc <- c("ferritin", "cart_target", "both", "joint")
  br <- tab[!is.na(b_route)]
  p_route <- vapply(rc, function(k) if (nrow(br)) mean(br$b_route == k) else NA_real_, numeric(1))
  hc <- c("clone", "reserve", "both", "other")
  nq <- tab[actual != "quick"]
  p_who <- vapply(hc, function(k) if (nrow(nq)) mean(nq$who == k) else NA_real_, numeric(1))
  q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  data.table::data.table(
    patient_id = tab$patient_id[1], n_draws = nrow(tab), n_failed = attr(tab, "n_failed") %||% 0L,
    t(stats::setNames(p_ph, paste0("p_", cls))),
    t(stats::setNames(p_when, paste0("p_when_", wc))),
    t(stats::setNames(p_when_tot, paste0("p_whentot_", wc))),
    p_b_involved = mean(!is.na(tab$b_route)),
    t(stats::setNames(p_route, paste0("p_route_", rc))),
    t(stats::setNames(p_who, paste0("p_who_", hc))),
    f_a_lo = q(tab$f_a)[1], f_a = q(tab$f_a)[2], f_a_hi = q(tab$f_a)[3],
    f_bi_lo = q(tab$f_bi)[1], f_bi = q(tab$f_bi)[2], f_bi_hi = q(tab$f_bi)[3],
    f_bz_lo = q(tab$f_bz)[1], f_bz = q(tab$f_bz)[2], f_bz_hi = q(tab$f_bz)[3],
    # both CAR T-phase routes together (for the f_a / f_b phenotype map)
    f_b_lo = q(1 - (1 - tab$f_bi) * (1 - tab$f_bz))[1], f_b = q(1 - (1 - tab$f_bi) * (1 - tab$f_bz))[2],
    f_b_hi = q(1 - (1 - tab$f_bi) * (1 - tab$f_bz))[3],
    f30 = stats::median(tab$f30), S30 = stats::median(tab$S30), peak_cart = stats::median(tab$peak_cart))
}
