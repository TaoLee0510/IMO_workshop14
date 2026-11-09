# =============================================================================
# model_checks.R
# Purpose : Checks of Model3 used by step K1 and its unit tests: the
#           steady-state Jacobian (stability), early nadir / recovery summaries,
#           the comparison with the Aim4 v1.3 reference engine, and the
#           calibrations of the th_aP prior centre, of the lymphocyte constants
#           (th_aL, c_L) and of the CAR T - target module.
# Inputs  : config, structural vectors.
# Outputs : parameter vectors, eigenvalues, summary tables.
# =============================================================================

# Calibration targets
#   th_aP prior centre: typical-design ANC nadir (M0, no CAR T phase) = Aim4 v1.3 at its
#          prior centre (aim4_v13_th_ap_centre)
#   th_aL, c_L: typical-design model ALC on days -2..0 = median per-patient ALC / base on
#          days -2..0 in the real extract (35 patients; Aim4 K1), and ALC on day 30 = Aim4 v1.3
#   CAR T (Aim4 v1.3 K1; the module does not involve the haematopoietic layers): typical
#          design (a0 = 1, k_t prior centre) peaks on day 8 at 0.05 x 10^9/L; z50 = typical
#          peak of Z; th_zS centre: typical HSC loss through Z = through the ferritin-driven I
#          at the th_bS prior centre, with the typical real ferritin kernel (A 0.88, tp 7.6 d)
# Model3: th_aM such that the typical-design AMC on days -1..0 = 0.078 x baseline (median per-patient value
#          of the real extract, 2026-11-04)
calib_targets <- list(alc_rel_ld = 0.042, alc_day = 30, amc_rel_ld = 0.078,
                      cart_peak_day = 8, cart_peak = 0.05, kernel = c(A = 0.88, tp = 7.6))

typical_m0_run <- function(config, struct = config$model$struct, th_ap = exp(config$model$re_prior$mu[["th_ap"]]),
                           days = seq(-30, 90, by = 0.25)) {
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  seg <- make_segments(des, struct, max(days), config)
  m3_sim(build_par(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = th_ap, th_bs = 0), struct, model_spec("M0")),
         seg, z$v, z$t0, z$dt, days)
}

typical_aim4_run <- function(config, days = seq(-30, 90, by = 0.25)) {
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  seg <- make_segments(des, aim4_v13_struct, max(days), config)
  aim4_sim(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = aim4_v13_th_ap_centre, th_bs = 0), aim4_v13_struct,
           model_spec("M0"), seg, z, days, config)
}

# Eigenvalues of the Jacobian at the pre-treatment steady state (f0 = 0; the neutral S_X
# direction of the shared niche and the antigen load before any CAR T cells are removed).
# Central differences.
steady_eigen <- function(struct, r_s = 0.08, h = 1e-6) {
  par <- build_par(c(r_s = r_s, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0), struct, model_spec("M0"))
  nm <- as.character(m3_state_names())
  y0 <- stats::setNames(rep(1, length(nm)), nm)
  y0[c("SX", "C", "I1", "I", "CT", "Z1", "Z", "GC")] <- 0
  f <- function(y) m3_rhs(0, y, par, 0, 0, 0, c(0, 0), -1, 1)
  # (unused neutrophil compartments beyond n_tr_n are constant: zero eigenvalues, left out)
  unused <- paste0("TN", seq(struct[["n_tr_n"]] + 1, length.out = 8 - struct[["n_tr_n"]]))
  keep <- setdiff(seq_along(nm), which(nm %in% c("SX", "AG", unused)))
  J <- vapply(keep, function(j) {
    e <- numeric(length(nm)); e[j] <- h
    (f(y0 + e) - f(y0 - e)) / (2 * h)
  }, numeric(length(nm)))[keep, , drop = FALSE]
  eigen(J, only.values = TRUE)$values
}

# Stability grid: feedback exponent gamma (up to the upper limit of its search box) x scaling
# of the MPP self-renewal rate (capped at 0.95 k_p, since sr_p < k_p). With sr_p = 0.21 the
# steady state loses stability at gamma = 2.09 (Model3 chain, mtt_n 5; 2.37 with the Model2 chain);
# with sr_p/k_p = 0.95 already at gamma ~ 1.
stability_grid <- function(config, gammas = c(0.25, 0.5, 1, 1.5, 2, config$model$struct_limits$gamma[2]),
                           scales = c(0.5, 1, 1.2)) {
  st <- config$model$struct
  g <- data.table::CJ(gamma = unique(gammas), sr_scale = scales)
  g[, c("max_re", "slowest_rate", "period_days") := {
    s <- st
    s[["gamma"]] <- gamma
    s[["sr_p"]] <- min(sr_scale * st[["sr_p"]], 0.95 * st[["k_p"]])
    ev <- steady_eigen(s)
    osc <- ev[abs(Im(ev)) > 1e-8]
    list(max(Re(ev)), -max(Re(ev)), if (length(osc)) 2 * pi / abs(Im(osc[which.max(Re(osc))])) else NA_real_)
  }, by = .(gamma, sr_scale)]
  g
}

# Early nadir and recovery per lineage from a trajectory matrix / table with
# columns t (or day) and *_rel; recovery = first day after the nadir >= thr.
nadir_summary <- function(m, cols = c(ANC = "anc_rel", PLT = "plt_rel", AMC = "amc_rel", ALC = "alc_rel"),
                          from = -5, thr = 0.8) {
  m <- as.data.frame(m)
  tt <- if ("t" %in% names(m)) m$t else m$day
  data.table::rbindlist(lapply(names(cols), function(k) {
    sel <- tt >= from
    x <- m[[cols[[k]]]][sel]; d <- tt[sel]
    i <- which.min(x)
    r <- which(seq_along(x) > i & x >= thr)
    data.table::data.table(lineage = k, nadir_day = d[i], nadir_rel = x[i],
                           recovery_day = if (length(r)) d[r[1]] else NA_real_)
  }))
}

# Early nadir and recovery of the typical patient: Aim4 v1.3 vs Model3 at the calibrated
# prior centres (M0, no CAR T phase).
aim4_comparison <- function(config) {
  ap <- exp(config$model$re_prior$mu[["th_ap"]])
  rbind(nadir_summary(typical_aim4_run(config))[, engine := sprintf("Aim4 v1.3 (th_aP = %.2f)", aim4_v13_th_ap_centre)],
        nadir_summary(typical_m0_run(config))[, engine := sprintf("Model3 (th_aP = %.2f)", ap)])
}

# ---- calibrations (rerun by K1; config.R holds the rounded results) -------------------------
# th_aP prior centre such that the typical-design ANC nadir (M0, no CAR T phase) equals Aim4 v1.3
calibrate_th_ap <- function(config) {
  days <- seq(-30, 90, by = 0.25)
  target <- min(typical_aim4_run(config, days)[days >= -5, "anc_rel"])
  f <- function(lx) log(min(typical_m0_run(config, th_ap = exp(lx), days = days)[days >= -5, "anc_rel"])) - log(target)
  exp(stats::uniroot(f, c(log(0.2), log(30)), tol = 1e-8)$root)
}

# th_aL and c_L such that the typical-design ALC on days -2..0 equals the observed median and
# the ALC on day 30 equals Aim4 v1.3
calibrate_lymph <- function(config) {
  days <- seq(-30, 90, by = 0.25)
  a4 <- typical_aim4_run(config, days)
  target <- c(ld = calib_targets$alc_rel_ld, d30 = unname(a4[which(days == calib_targets$alc_day), "alc_rel"]))
  st <- config$model$struct
  obj <- function(w) {
    s2 <- st
    s2[c("th_al", "c_l")] <- exp(w)
    m <- typical_m0_run(config, s2, days = days)
    (log(stats::median(m[days >= -2 & days <= 0, "alc_rel"])) - log(target[["ld"]]))^2 +
      (log(m[days == calib_targets$alc_day, "alc_rel"]) - log(target[["d30"]]))^2
  }
  o <- stats::optim(log(st[c("th_al", "c_l")]), obj, control = list(reltol = 1e-14, maxit = 2000))
  c(th_al = exp(o$par[[1]]), c_l = exp(o$par[[2]]), alc_d30_target = target[["d30"]], loss = o$value)
}

# Model3: th_aM such that the typical-design AMC on days -1..0 equals the observed median
calibrate_th_am <- function(config) {
  days <- seq(-30, 5, by = 0.25)
  st <- config$model$struct
  f <- function(lx) {
    s2 <- st
    s2[["th_am"]] <- exp(lx)
    m <- typical_m0_run(config, s2, days = days)
    log(stats::median(m[days >= -1 & days <= 0, "amc_rel"])) - log(calib_targets$amc_rel_ld)
  }
  exp(stats::uniroot(f, c(log(1e-3), log(20)), tol = 1e-8)$root)
}

# ---- CAR T - target interaction (Aim4 v1.3 calibrations, rechecked on the Model3 engine) -----
# Typical design with the CAR T infusion (no chemotherapy injury of the stem cells), quarter-
# day grid; returns the trajectory matrix.
typical_cart_run <- function(config, struct = config$model$struct, rho = exp(config$model$re_prior$mu[["rho"]]),
                             drv = NULL, spec = model_spec("Mb"), th_bs = exp(config$model$re_prior$mu[["th_bs"]]),
                             days = seq(-30, 60, by = 0.25)) {
  des <- typical_design(config)
  if (is.null(drv)) drv <- driver_from_kernel(des$t0, calib_targets$kernel[["A"]], calib_targets$kernel[["tp"]], config)
  seg <- make_segments(des, struct, max(days), config)
  ind <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = exp(config$model$re_prior$mu[["th_ap"]]), th_bs = th_bs,
           rho = rho, a0 = 1)
  m3_sim(build_par(ind, struct, spec), seg, drv$v, drv$t0, drv$dt, days)
}

# rho prior centre and kappa: typical CAR T peak on day cart_peak_day at cart_peak
calibrate_cart <- function(config) {
  st <- config$model$struct
  obj <- function(w) {
    s2 <- st
    s2[["kappa"]] <- exp(w[2])
    m <- typical_cart_run(config, s2, rho = exp(w[1]))
    i <- which.max(m[, "cart"])
    (m[i, "t"] - calib_targets$cart_peak_day)^2 + (10 * log(m[i, "cart"] / calib_targets$cart_peak))^2
  }
  o <- stats::optim(c(log(0.7), log(20)), obj, control = list(reltol = 1e-12, maxit = 2000))
  c(rho = exp(o$par[1]), kappa = exp(o$par[2]), loss = o$value)
}

# z50 and the th_zS prior centre (given the calibrated rho and kappa)
calibrate_z <- function(config) {
  st <- config$model$struct
  m <- typical_cart_run(config, st)
  tt <- m[, "t"]
  trap <- function(y) sum(diff(tt) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
  sel <- tt >= 0
  auc_i <- trap(ifelse(sel, m[, "I"], 0))
  auc_z <- trap(ifelse(sel, m[, "z"], 0))
  th_bs0 <- exp(config$model$re_prior$mu[["th_bs"]])
  c(z50 = max(m[, "z"]), th_zs = th_bs0 * auc_i / auc_z, auc_i = auc_i, auc_z = auc_z,
    hsc_loss_route = 1 - exp(-th_bs0 * auc_i))
}

# ---- reachable phenotypes ----------------------------------------------------------------------
# Phenotypes that the CAR T phase alone (b) and one chemotherapy hit alone (a) can give over a
# grid of driver strength / HSC injury (x emax for the CAR T phase), at the prior centres otherwise (typical design,
# ANC baseline 3.5). Used to set the emax prior centre and checked by the unit tests.
phenotype_screen <- function(config, emax = c(0.3, 0.5, 0.7, 0.9)) {
  des <- typical_design(config)
  ap <- exp(config$model$re_prior$mu[["th_ap"]])
  days <- seq(-30, 90)
  b <- data.table::CJ(emax = emax, A = c(1, 1.5, 2.5, 3.5), th_bs = c(0, 0.02, 0.06))
  b[, phenotype := {
    s2 <- config$model$struct
    s2[["emax"]] <- emax
    drv <- driver_from_kernel(des$t0, A, 9, config)
    m <- sim_daily(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = th_bs, rho = exp(config$model$re_prior$mu[["rho"]]),
                     a0 = 1), s2, model_spec("Mb"), des, drv, config, days = days)
    phenotype_sim(m$day, m$anc_rel * 3.5, config)
  }, by = .(emax, A, th_bs)]
  # chemotherapy alone (no CAR T phase: emax has no effect here)
  a <- data.table::CJ(r_s = c(0.02, 0.04, 0.08), th_as = c(0.4, 1, 1.5, 2.5))
  des_nc <- des
  des_nc$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  a[, phenotype := {
    m <- sim_daily(c(r_s = r_s, f0 = 0, th_as = th_as, th_ap = ap, th_bs = 0), config$model$struct, model_spec("Ma"),
                   des_nc, z, config, days = days)
    phenotype_sim(m$day, m$anc_rel * 3.5, config)
  }, by = .(r_s, th_as)]
  list(b_only = b, a_only = a)
}
