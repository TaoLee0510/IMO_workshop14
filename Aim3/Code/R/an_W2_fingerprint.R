# =============================================================================
# an_W2_fingerprint.R
# Purpose : Step W2 - what the blood counts can show about when the stem cell
#           pool was hit. Uses no patient counts: a design patient with the
#           cohort's median ferritin driver and median set points.
#   1 trajectories  same total HSC loss delivered in a, in b, or split
#   2 phase map     phenotype over (f_a, f_b), two ANC set points, eta 0 / 3
#   3 pseudo wave   can a chemotherapy-only hit produce a second wave? scanned
#                   over multipotent-progenitor regimes (phi * gmax < or >= 1)
#   4 timing        day of each lineage's largest HSC footprint as the driver
#                   peak is shifted: b follows the driver, a does not
# Inputs  : ctx (W0 derived data, W1 design driver).
# Outputs : tables for plotting and the summary.
# =============================================================================

shift_driver <- function(drv, days) {
  k <- round(days / drv$dt)
  v <- drv$v
  n <- length(v)
  if (k > 0) v <- c(rep(0, k), v[seq_len(n - k)])
  if (k < 0) v <- c(v[(1 - k):n], rep(v[n], -k))
  drv$v <- v
  drv
}

w2_phase_map <- function(pd, st, config, set_name) {
  g <- config$steps$W2$f_grid
  mu <- exp(config$model$design$mu)
  grid <- data.table::CJ(f_a = g, f_b = g, eta = config$steps$W2$eta_levels)
  grid[, pheno := mapply(function(fa, fb, eta) {
    s <- replace(st, "eta_syn", eta)
    ind <- c(th_ap = mu[["th_ap"]], th_as = theta_for_loss(fa, pd, s, "a"),
             th_bs = theta_for_loss(fb, pd, s, "b"), r_s = mu[["r_s"]])
    sim_phenotype(build_par(ind, s, "Mabx"), pd, config)
  }, f_a, f_b, eta)]
  grid[, set := set_name]
  grid
}

# Second-wave depth of a trajectory: after the first value >= 1.5 that follows
# the early nadir, the lowest later ANC (days <= 30).
echo_depth <- function(tr) {
  x <- tr[day >= 0 & day <= 30]
  nd <- x[day <= 14][which.min(anc), day]
  r <- x[day > nd & anc >= 1.5, day]
  if (!length(r)) return(c(rec_day = NA_real_, min_after = NA_real_))
  rec <- min(r)
  c(rec_day = rec, min_after = x[day > rec, if (.N) min(anc) else NA_real_])
}

w2_pseudo_wave <- function(pd, config) {
  mu <- exp(config$model$design$mu)
  f_seq <- 1 - 10^seq(0, -3, length.out = 31)
  data.table::rbindlist(lapply(names(config$steps$W2$regimes), function(rg) {
    st <- config$model$struct
    st[names(config$steps$W2$regimes[[rg]])] <- config$steps$W2$regimes[[rg]]
    data.table::rbindlist(lapply(f_seq, function(fa) {
      ind <- c(th_ap = mu[["th_ap"]], th_as = theta_for_loss(fa, pd, st, "a"), th_bs = 0, r_s = mu[["r_s"]])
      par <- build_par(ind, st, "Ma")
      tr <- sim_daily(par, pd, config, days = 0:45)
      e <- echo_depth(tr)
      data.table::data.table(regime = rg, phi_gmax = st[["phi"]] * st[["gmax"]], f_a = fa,
                             pheno = classify_sim(tr$day, tr$anc, config),
                             rec_day = e[["rec_day"]], min_after = e[["min_after"]])
    }))
  }))
}

w2_timing <- function(pd, st, config) {
  mu <- exp(config$model$design$mu)
  base_ind <- c(th_ap = mu[["th_ap"]], th_as = 0, th_bs = 0, r_s = mu[["r_s"]])
  days <- seq(0, 50, by = 0.25)
  data.table::rbindlist(lapply(config$steps$W2$peak_shifts, function(s) {
    pds <- pd
    pds$drv <- shift_driver(pd$drv, s)
    ref <- sim_daily(build_par(base_ind, st, "Mab"), pds, config, days = days)
    cases <- list(a = replace(base_ind, "th_as", theta_for_loss(0.97, pds, st, "a")),
                  b = replace(base_ind, "th_bs", theta_for_loss(0.97, pds, st, "b")))
    data.table::rbindlist(lapply(names(cases), function(k) {
      tr <- sim_daily(build_par(cases[[k]], st, "Mab"), pds, config, days = days)
      data.table::rbindlist(lapply(c("amc", "anc", "plt"), function(v) {
        r <- tr[[v]] / ref[[v]]
        data.table::data.table(phase = k, shift = s, lineage = v, footprint_day = days[which.min(r)],
                               footprint_depth = min(r))
      }))
    }))
  }))
}

an_W2_fingerprint <- function(ctx) {
  cfg <- ctx$config
  w1 <- ctx$results$W1
  st <- cfg$model$struct
  pd_med <- design_pdes(st, cfg, w1$drv, set = w1$set_med)
  pd_low <- design_pdes(st, cfg, w1$drv, set = w1$set_low)

  traj <- typical_cases(pd_med, st, cfg, loss = cfg$steps$W2$loss_demo)
  phase <- rbind(w2_phase_map(pd_med, st, cfg, "median"), w2_phase_map(pd_low, st, cfg, "lower quartile"))
  pseudo <- w2_pseudo_wave(pd_low, cfg)
  timing <- w2_timing(pd_med, st, cfg)

  slopes <- timing[, .(slope = stats::coef(stats::lm(footprint_day ~ shift))[[2]]), by = .(phase, lineage)]
  thr <- pseudo[pheno == "second_wave", .(min_f_a = if (.N) min(f_a) else NA_real_), by = regime]
  thr <- merge(pseudo[, .(phi_gmax = phi_gmax[1], max_dip = suppressWarnings(min(min_after, na.rm = TRUE))), by = regime],
               thr, by = "regime", all.x = TRUE)
  ctx$log("  W2: pseudo second wave from chemotherapy alone: ",
          paste(thr$regime, ifelse(is.na(thr$min_f_a), "none", sprintf("f_a >= %.3f", thr$min_f_a)), collapse = "; "))
  ctx$log("  W2: footprint slope vs driver shift (ANC): a ", round(slopes[phase == "a" & lineage == "anc", slope], 2),
          ", b ", round(slopes[phase == "b" & lineage == "anc", slope], 2))
  list(tables = list(trajectories = traj, phase_map = phase, pseudo_wave = pseudo,
                     pseudo_threshold = thr, timing = timing, timing_slopes = slopes))
}
