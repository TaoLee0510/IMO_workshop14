# =============================================================================
# model_checks.R
# Purpose : Checks of the v1.2 parent model used by step K1 and its unit tests:
#           the frozen v1.1 engine (regression reference), the steady-state
#           Jacobian (stability), early nadir / recovery summaries, and the
#           calibration targets of th_aP and th_aL.
# Inputs  : config, structural vectors.
# Outputs : parameter vectors, eigenvalues, summary tables.
# =============================================================================

# Calibration targets (2026-11-03, written into config.R)
#   th_aP prior centre: typical-design ANC nadir = v1.1 engine at th_aP = 1
#   th_aL: typical-design model ALC on days -2..0 = median per-patient ALC / base
#          on days -2..0 in the real extract (35 patients)
#   v1.3 CAR T: typical design (a0 = 1, k_t prior centre) peaks on day 8 at 0.05 x 10^9/L
#          (the observed median ALC excess on days 8-11 over the lymphocyte recovery);
#          z50 = typical peak of Z; th_zS centre: typical HSC loss through Z = through the
#          ferritin-driven I at the th_bS prior centre, with the typical real ferritin kernel
#          (kernel population means of the 35 drivers: A 0.88, tp 7.6 d)
calib_targets <- list(v11_th_ap = 1, alc_rel_ld = 0.042,
                      cart_peak_day = 8, cart_peak = 0.05, kernel = c(A = 0.88, tp = 7.6))

load_v11_reference <- function(config) {
  if (exists("v11_a4_sim", mode = "function")) return(invisible(TRUE))
  dir.create(config$paths$rcpp_cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(file.path(config$paths$code, "tests", "ref", "aim4_model_v11.cpp"),
                  cacheDir = config$paths$rcpp_cache, env = globalenv())
  invisible(TRUE)
}

# v1.1 defaults (frozen copy of the v1.1 config), individual values on top.
v11_par <- function(indiv = c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0), overrides = NULL) {
  p <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0,
         alpha = 1, k_m = 0.3, phi = 0.7, gmax = 3, gamma_m = 0.5, r_a = 1, r_b = 1, delta = 1,
         ktr_n = 4 / 5, ktr_pl = 4 / 6, ktr_mo = 4 / 3.5, c_n = 2.3, c_pl = 0.12, c_mo = 0.5,
         gam_n = 0.2, gam_pl = 0.2, gam_mo = 0.2, k_i = 0.2, emax = 0.9, ic50 = 0.5, hill = 2,
         e_pl = 0.7, e_mo = 0.7, eta_syn = 0, k_e = 0.8, gcsf_fold = 2.5, tpo_fold = 1.8)
  p[names(indiv)] <- indiv
  if (!is.null(overrides)) p[names(overrides)] <- overrides
  p
}

# v1.2 structure that reduces to a v1.1 parameter vector: CMP an instantaneous
# relay, GMP / MEP the v1.1 committed progenitors (no self-renewal), feedback on
# the MPP only, no production feedback (gam_l = 0 in both), N and Mo with the same
# transit rate (they share the GMP). No CAR T-phase driver.
v12_as_v11 <- function(p11, struct) {
  s <- struct
  s[c("k_mpp", "sr_mpp", "gmax", "gamma_m", "alpha")] <-
    c(p11[["k_m"]], p11[["phi"]] * p11[["k_m"]], p11[["gmax"]], p11[["gamma_m"]], p11[["alpha"]])
  s[c("k_cmp", "sr_cmp", "k_gmp", "sr_gmp", "k_mep", "sr_mep")] <- c(1000, 0, p11[["ktr_n"]], 0, p11[["ktr_pl"]], 0)
  s[c("fb_cmp", "fb_gmp", "fb_mep", "fb_mpp")] <- c(0, 0, 0, 1)
  s[c("mtt_n", "mtt_mo", "mtt_pl")] <- c(3 / p11[["ktr_n"]], 3 / p11[["ktr_mo"]], 3 / p11[["ktr_pl"]])
  s[c("c_n", "c_mo", "c_pl")] <- p11[c("c_n", "c_mo", "c_pl")]
  s[c("gam_n", "gam_mo", "gam_pl")] <- 0
  s
}

# Eigenvalues of the Jacobian at the pre-treatment steady state (f0 = 0; the
# neutral S_X direction of the shared niche is removed). Central differences.
steady_eigen <- function(struct, r_s = 0.08, h = 1e-6) {
  par <- build_par(c(r_s = r_s, f0 = 0, th_as = 0, th_ap = 1, th_bs = 0), struct, model_spec("M0"))
  nm <- as.character(a4_state_names())
  y0 <- stats::setNames(rep(1, length(nm)), nm)
  y0[c("SX", "C", "I1", "I", "CT", "Z1", "Z", "GC")] <- 0
  f <- function(y) a4_rhs(0, y, par, 0, 0, 0, c(0, 0), -1, 1)
  # neutral directions: S_X (shared niche) and the antigen load before any CAR T cells
  keep <- setdiff(seq_along(nm), which(nm %in% c("SX", "AG")))
  J <- vapply(keep, function(j) {
    e <- numeric(length(nm)); e[j] <- h
    (f(y0 + e) - f(y0 - e)) / (2 * h)
  }, numeric(length(nm)))[keep, , drop = FALSE]
  eigen(J, only.values = TRUE)$values
}

# Stability grid: feedback gain gamma_M x scaling of every self-renewal rate
# (capped at 0.95 k, since sr < k).
stability_grid <- function(config, gammas = c(0.25, 0.5, 1, 1.5), scales = c(0.5, 1, 1.5)) {
  st <- config$model$struct
  prog <- c("mpp", "cmp", "gmp", "mep", "clp", "ly")
  g <- data.table::CJ(gamma_m = gammas, sr_scale = scales)
  g[, c("max_re", "slowest_rate", "period_days") := {
    s <- st
    s[["gamma_m"]] <- gamma_m
    for (x in prog) s[[paste0("sr_", x)]] <- min(sr_scale * st[[paste0("sr_", x)]], 0.95 * st[[paste0("k_", x)]])
    ev <- steady_eigen(s)
    lead <- ev[which.max(Re(ev))]
    osc <- ev[abs(Im(ev)) > 1e-8]
    list(max(Re(ev)), -max(Re(ev)), if (length(osc)) 2 * pi / abs(Im(osc[which.max(Re(osc))])) else NA_real_)
  }, by = .(gamma_m, sr_scale)]
  g
}

# Early nadir and recovery per lineage from a trajectory matrix / table with
# columns t (or day) and *_rel; recovery = first day after the nadir >= thr.
nadir_summary <- function(m, cols = c(ANC = "anc_rel", PLT = "plt_rel", AMC = "amc_rel"), from = -5, thr = 0.8) {
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

# Early nadir and recovery of the typical patient: v1.1 at th_aP = 1 vs v1.2 at the
# calibrated prior centre (M0, no CAR T phase).
v11_comparison <- function(config) {
  load_v11_reference(config)
  st <- config$model$struct
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  days <- seq(-30, 90, by = 0.25)
  seg <- make_segments(des, st, max(days), config)
  m11 <- v11_a4_sim(v11_par(c(th_ap = calib_targets$v11_th_ap)), seg[, 1:6, drop = FALSE], z$v, z$t0, z$dt, days)
  ap <- exp(config$model$re_prior$mu[["th_ap"]])
  m12 <- a4_sim(build_par(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0), st, model_spec("M0")),
                seg, z$v, z$t0, z$dt, days)
  rbind(nadir_summary(m11)[, engine := "v1.1 (th_aP = 1)"],
        nadir_summary(m12, c(ANC = "anc_rel", PLT = "plt_rel", AMC = "amc_rel", ALC = "alc_rel"))[
          , engine := sprintf("v1.2 (th_aP = %.2f)", ap)])
}

# ---- calibrations (rerun by K1; config.R holds the rounded results) -------------------------
# th_aP prior centre such that the typical-design ANC nadir (M0, no CAR T phase) equals
# v1.1 at th_aP = 1
calibrate_th_ap <- function(config) {
  load_v11_reference(config)
  st <- config$model$struct
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  days <- seq(-30, 90, by = 0.25)
  seg <- make_segments(des, st, max(days), config)
  target <- min(v11_a4_sim(v11_par(c(th_ap = calib_targets$v11_th_ap)), seg[, 1:6, drop = FALSE], z$v, z$t0, z$dt,
                           days)[days >= -5, "anc_rel"])
  f <- function(lx) {
    m <- a4_sim(build_par(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = exp(lx), th_bs = 0), st, model_spec("M0")),
                seg, z$v, z$t0, z$dt, days)
    log(min(m[days >= -5, "anc_rel"])) - log(target)
  }
  exp(stats::uniroot(f, c(log(0.2), log(20)), tol = 1e-8)$root)
}

# th_aL such that the typical-design model ALC on days -2..0 equals the observed median
calibrate_th_al <- function(config) {
  st <- config$model$struct
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  days <- seq(-30, 5, by = 0.25)
  seg <- make_segments(des, st, max(days), config)
  ap <- exp(config$model$re_prior$mu[["th_ap"]])
  f <- function(lx) {
    s2 <- st
    s2[["th_al"]] <- exp(lx)
    m <- a4_sim(build_par(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0), s2, model_spec("M0")),
                seg, z$v, z$t0, z$dt, days)
    log(stats::median(m[days >= -2 & days <= 0, "alc_rel"])) - log(calib_targets$alc_rel_ld)
  }
  exp(stats::uniroot(f, c(log(0.05), log(10)), tol = 1e-8)$root)
}

# Screen of the CMP/GMP/MEP self-renewal rate (one shared value): for each value,
# th_aP recalibrated, the typical ANC nadir / recovery, and the phenotypes a single
# chemotherapy hit can give over a grid of reserve x HSC injury (gamma_M 0.5 and 1).
progenitor_screen <- function(config, values = c(0.5, 0.75, 1)) {
  des <- typical_design(config)
  des$no_car <- TRUE
  z <- zero_driver(des$t0, config)
  days <- seq(-30, 90)
  data.table::rbindlist(lapply(values, function(sp) {
    cfg <- config
    cfg$model$struct[c("sr_cmp", "sr_gmp", "sr_mep")] <- sp
    ap <- calibrate_th_ap(cfg)
    cfg$model$re_prior$mu[["th_ap"]] <- log(ap)
    st <- cfg$model$struct
    vc <- v11_comparison(cfg)[grepl("v1.2", engine) & lineage == "ANC"]
    g <- data.table::CJ(r_s = c(0.01, 0.02, 0.04, 0.08), th_as = c(0.4, 1, 1.5, 2, 2.5, 3.5), gamma_m = c(0.5, 1))
    ph <- g[, .(ph = {
      s2 <- st
      s2[["gamma_m"]] <- gamma_m
      m <- sim_daily(c(r_s = r_s, f0 = 0, th_as = th_as, th_ap = ap, th_bs = 0), s2, model_spec("Ma"), des, z, cfg,
                     days = days)
      phenotype_sim(m$day, m$anc_rel * 3.5, cfg)
    }), by = names(g)]$ph
    data.table::data.table(sr_progenitor = sp, th_ap_centre = ap, anc_nadir_day = vc$nadir_day,
                           anc_recovery_day = vc$recovery_day, n_quick = sum(ph == "quick"),
                           n_second_wave = sum(ph == "second_wave"), n_delayed = sum(ph == "delayed"))
  }))
}

# ---- v1.3: CAR T - target interaction --------------------------------------------------------
# Typical design with the CAR T infusion (no chemotherapy injury of the stem cells), daily
# grid; returns the trajectory matrix.
typical_cart_run <- function(config, struct = config$model$struct, rho = exp(config$model$re_prior$mu[["rho"]]),
                             drv = NULL, spec = model_spec("Mb"), th_bs = exp(config$model$re_prior$mu[["th_bs"]]),
                             days = seq(-30, 60, by = 0.25)) {
  des <- typical_design(config)
  if (is.null(drv)) drv <- driver_from_kernel(des$t0, calib_targets$kernel[["A"]], calib_targets$kernel[["tp"]], config)
  seg <- make_segments(des, struct, max(days), config)
  ind <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = exp(config$model$re_prior$mu[["th_ap"]]), th_bs = th_bs,
           rho = rho, a0 = 1)
  a4_sim(build_par(ind, struct, spec), seg, drv$v, drv$t0, drv$dt, days)
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
