# =============================================================================
# sim_cohort.R
# Purpose : Synthetic cohorts with known truth, in the same format as
#           read_aim4_data(). Two designs:
#             real35 - the 35 real patients' sampling days, ferritin curves,
#                      eGFR, bridging and set points (their counts are NOT used)
#             n350   - the planned design: 350 patients, daily counts to day 14,
#                      every 2 days to day 90, ferritin from the kernel model
#           Used by K3 (identifiability), TestData and the virtual cohorts.
# Inputs  : design rows, truth population, mechanism scenario, config.
# Outputs : list(patients, obs, inflam, drivers, truth, vaf)
# =============================================================================

# ---- designs ---------------------------------------------------------------------------------
# One row per design patient: patient table row, driver, sampling days per variable.
design_real35 <- function(d) {
  lapply(stats::setNames(d$patients$patient_id, d$patients$patient_id), function(pid) {
    o <- d$obs[patient_id == pid]
    list(prow = d$patients[patient_id == pid], drv = d$drivers[[pid]],
         days = lapply(stats::setNames(lik_vars, lik_vars), function(v) sort(unique(o[var == v]$day))),
         inflam = d$inflam[patient_id == pid])
  })
}

# Population of ferritin kernels fitted to the real smoothed drivers.
fit_kernel_population <- function(drivers, config) {
  h <- config$model$kernel$h
  pars <- t(vapply(drivers, function(dr) {
    tt <- dr$t0 + (seq_along(dr$v) - 1) * dr$dt
    sel <- tt > 0 & tt <= 35
    y <- dr$v[sel]
    if (max(y) < 0.05) return(c(A = 0.05, tp = 7))
    obj <- function(w) sum((y - kernel_fun(tt[sel], exp(w[1]), exp(w[2]), h))^2)
    o <- stats::optim(c(log(max(y)), log(tt[sel][which.max(y)])), obj)
    c(A = exp(o$par[1]), tp = exp(o$par[2]))
  }, numeric(2)))
  lp <- log(pmax(pars, 1e-3))
  list(pars = pars, mu = colMeans(lp), cov = stats::cov(lp) + diag(1e-4, 2), h = h)
}

design_n350 <- function(n, d, kpop, config, rng_seed) {
  set.seed(rng_seed)
  base_rows <- d$patients[sample(nrow(d$patients), n, replace = TRUE)]
  days_cbc <- c(-30, -20, -10, -7:14, seq(16, 90, by = 2))
  kd <- MASS::mvrnorm(n, kpop$mu, kpop$cov)
  out <- lapply(seq_len(n), function(i) {
    prow <- data.table::copy(base_rows[i])
    prow[, `:=`(patient_id = sprintf("V%03d", i), bridging = "none", t0 = config$model$t_min)]
    A <- exp(kd[i, 1]); tp <- exp(kd[i, 2])
    drv <- driver_from_kernel(prow$t0, A, tp, config)
    fd <- c(-5, 0:30)
    fer <- prow$ferritin_ld * exp(kernel_fun(fd, A, tp, kpop$h) + stats::rnorm(length(fd), 0, 0.15))
    list(prow = prow, drv = drv, days = list(anc = days_cbc, plt = days_cbc, amc = days_cbc),
         inflam = data.table::data.table(patient_id = prow$patient_id, day = fd, ferritin = fer, crp = NA_real_))
  })
  stats::setNames(out, vapply(out, function(x) x$prow$patient_id, character(1)))
}

# ---- truth -------------------------------------------------------------------------------------
# mech: a_only | b_only | mix | synergy; p_clone: share of clone carriers (ddr type)
draw_truth <- function(ids, mech, p_clone, pop, config, shift = c(r_s = 0, inj = 0)) {
  n <- length(ids)
  z <- function(k) pop$mu[[k]] + pop$omega[[k]] * stats::rnorm(n)
  r_s <- exp(z("r_s") + shift[["r_s"]])
  th_ap <- exp(z("th_ap"))
  th_as <- exp(z("th_as") + shift[["inj"]])
  th_bs <- exp(z("th_bs") + shift[["inj"]])
  cat_ <- switch(mech,
                 a_only = rep("a", n), b_only = rep("b", n), synergy = rep("ab", n),
                 mix = sample(c("a", "b", "ab"), n, replace = TRUE, prob = config$steps$K3$mix))
  th_bs[cat_ == "a"] <- 0
  th_as[cat_ == "b"] <- 0
  carrier <- stats::runif(n) < p_clone
  f0 <- ifelse(carrier, stats::plogis(z("f0")), 0)
  data.table::data.table(patient_id = ids, true_cat = cat_, carrier = carrier,
                         r_s = r_s, f0 = f0, th_as = th_as, th_ap = th_ap, th_bs = th_bs)
}

truth_struct <- function(mech, config) {
  st <- config$model$struct
  st[names(config$model$clone_types$ddr)] <- config$model$clone_types$ddr
  st[["eta_syn"]] <- if (mech == "synergy") config$model$design$eta_syn else 0
  st
}
truth_spec <- function(mech) model_spec(if (mech == "synergy") "MabxC_ddr" else "MabC_ddr")

# ---- one virtual patient -------------------------------------------------------------------------
# Latent daily trajectory, G-CSF given when the latent ANC falls below the threshold after
# day `from` (probability p), observations with noise, censoring and transfusion exclusions.
simulate_virtual <- function(dz, tr, struct, spec, config, gcsf_allowed = TRUE) {
  ds <- config$model$design
  prow <- dz$prow
  design <- patient_design(prow)
  base <- list(anc = prow$anc_base, plt = prow$plt_base, amc = prow$amc_base)
  ind <- c(r_s = tr$r_s, f0 = tr$f0, th_as = tr$th_as, th_ap = tr$th_ap, th_bs = tr$th_bs)
  all_days <- sort(unique(c(unlist(dz$days), seq(0, config$model$sim_end))))
  all_days <- all_days[all_days > design$t0]
  lat <- sim_daily(ind, struct, spec, design, dz$drv, config, days = all_days, base = base)
  gcsf <- 0L
  obs_traj <- lat
  if (gcsf_allowed) {
    low <- lat[day >= ds$gcsf$from & day <= 30 & anc < ds$gcsf$anc_thr]
    if (nrow(low) && stats::runif(1) < ds$gcsf$p) {
      st2 <- struct
      st2[["gcsf_fold"]] <- stats::runif(1, ds$gcsf$fold[1], ds$gcsf$fold[2])
      design$gcsf <- cbind(min(low$day), min(low$day) + ds$gcsf$days)
      obs_traj <- sim_daily(ind, st2, spec, design, dz$drv, config, days = all_days, base = base)
      gcsf <- 1L
    }
  }
  if (anyNA(lat$anc_rel)) return(NULL)
  ob <- data.table::rbindlist(lapply(lik_vars, function(v) {
    dd <- dz$days[[v]]
    dd <- dd[dd > design$t0]
    val <- obs_traj[match(dd, day)][[v]] * exp(stats::rnorm(length(dd), 0, ds$sigma[[v]]))
    data.table::data.table(patient_id = prow$patient_id, day = dd, var = v, value = val,
                           cens = FALSE, upper = FALSE, keep = TRUE, reason = "")
  }))
  ob[var == "anc" & value < ds$cens$anc, `:=`(value = ds$cens$anc_bound, cens = TRUE, reason = "no_diff_low_wbc")]
  ob[var == "amc" & value < ds$cens$amc, `:=`(value = ds$cens$amc_bound, cens = TRUE, reason = "manual_diff_zero")]
  if (gcsf == 1L) ob[var == "anc" & day > config$data$gcsf_from_day, `:=`(upper = TRUE, reason = "gcsf_after_d14")]
  txd <- ob[var == "plt" & value < ds$plt_tx_below, day]
  for (t in txd) ob[var == "plt" & day > t & day <= t + config$data$plt_tx$drop_days,
                    `:=`(keep = FALSE, reason = "after_plt_tx")]
  lat_ph <- phenotype_sim(lat$day, lat$anc, config, with_persistent = max(lat$day) >= 75)
  list(obs = ob, gcsf = gcsf, latent = lat, latent_phenotype = lat_ph)
}

# ---- a whole cohort -------------------------------------------------------------------------------
# dz_list: design rows (design_real35 / design_n350); ids map onto designs in order
simulate_cohort <- function(dz_list, mech, p_clone, pop, config, seed, shift = c(r_s = 0, inj = 0),
                            ids = NULL, vaf_n = 0) {
  set.seed(seed)
  if (is.null(ids)) ids <- names(dz_list)
  truth <- draw_truth(ids, mech, p_clone, pop, config, shift)
  st <- truth_struct(mech, config)
  spec <- truth_spec(mech)
  sims <- lapply(seq_along(ids), function(i) simulate_virtual(dz_list[[i]], truth[i], st, spec, config))
  ok <- !vapply(sims, is.null, logical(1))
  truth <- truth[ok]; sims <- sims[ok]; dz_list <- dz_list[ok]
  pts <- data.table::rbindlist(lapply(seq_along(sims), function(i) {
    r <- data.table::copy(dz_list[[i]]$prow)
    r[, `:=`(patient_id = truth$patient_id[i], gcsf_after_d14 = sims[[i]]$gcsf)]
    r
  }), fill = TRUE)
  obs <- data.table::rbindlist(lapply(seq_along(sims), function(i) sims[[i]]$obs[, patient_id := truth$patient_id[i]]))
  inflam <- data.table::rbindlist(lapply(seq_along(sims), function(i)
    data.table::copy(dz_list[[i]]$inflam)[, patient_id := truth$patient_id[i]]), fill = TRUE)
  drivers <- stats::setNames(lapply(dz_list, `[[`, "drv"), truth$patient_id)
  truth[, `:=`(latent_phenotype = vapply(sims, `[[`, character(1), "latent_phenotype"), mech = mech)]
  vaf <- NULL
  if (vaf_n > 0) {
    sel <- if (is.finite(vaf_n)) sample(nrow(truth), min(vaf_n, nrow(truth))) else seq_len(nrow(truth))
    sd <- config$steps$K3$n350$vaf_sd
    vaf <- truth[sel, .(patient_id, status = ifelse(f0 > 0.01, "positive", "negative"),
                        f_measured = pmax(f0 * exp(stats::rnorm(.N, 0, sd)), 0.01), sd = sd, lod = 0.01)]
  }
  d <- list(patients = pts, obs = obs, inflam = inflam, drivers = drivers, truth = truth, vaf = vaf)
  d$phenotype <- phenotype_all_rules(obs, config)
  truth[d$phenotype, obs_phenotype := i.primary, on = "patient_id"]
  d$truth <- truth
  d
}

# Shift the population so that the observed phenotype mix matches the target
# (two knobs: reserve r_s and injury scale; damped fixed-point updates).
calibrate_design <- function(dz_list, mech, p_clone, pop, config, seed = 1L, log_fun = NULL) {
  s <- config$steps$K3
  tgt <- config$model$design$target
  shift <- c(r_s = 0, inj = 0)
  hist <- list()
  for (it in seq_len(s$calib_iter)) {
    ix <- rep_len(seq_along(dz_list), s$calib_n)
    ids <- sprintf("C%05d", seq_along(ix))
    co <- simulate_cohort(dz_list[ix], mech, p_clone, pop, config, seed + it, shift, ids = ids)
    p <- prop.table(table(factor(co$phenotype$primary, levels = names(tgt))))
    hist[[it]] <- data.table::data.table(iter = it, r_s_shift = shift[["r_s"]], inj_shift = shift[["inj"]],
                                         quick = p[["quick"]], second_wave = p[["second_wave"]],
                                         delayed = p[["delayed"]])
    if (!is.null(log_fun)) log_fun(sprintf("    calib %s it %d: %.2f/%.2f/%.2f (shift r_s %.2f, inj %.2f)",
                                           mech, it, p[[1]], p[[2]], p[[3]], shift[["r_s"]], shift[["inj"]]))
    shift[["r_s"]] <- shift[["r_s"]] - 2.0 * (tgt[["delayed"]] - p[["delayed"]])
    shift[["inj"]] <- shift[["inj"]] + 2.5 * ((1 - tgt[["quick"]]) - (1 - p[["quick"]]))
  }
  list(shift = shift, history = data.table::rbindlist(hist))
}

design_population <- function(config) {
  list(mu = config$model$design$mu, omega = config$model$design$omega)
}
