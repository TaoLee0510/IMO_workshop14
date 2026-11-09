# =============================================================================
# an_G_export.R
# Purpose : Step G for Model3 (public name MarrowTwin) - precompute everything GreenLight
#           (Model3/GreenLight) displays, with the R engine:
#             - per patient (35 real, de-identified, local use only; 16 Model3
#               digital patients with known truth) and per as-of day (D-5..D+30
#               daily, weekly to D+90): individual posterior from the data up to
#               that day without the phenotype term (MAP + Laplace), posterior
#               bands of ANC/PLT/AMC/ALC, P(quick / second wave / delayed),
#               predicted ANC recovery day and nadir, the when / who readouts,
#               a plan library ranked by an illustrative utility, and what the
#               in-browser what-if needs (MAP, Cholesky factor, driver);
#             - the public index (digital roster, prediction evolution, population,
#               parameter table from config, K1 results, results status);
#             - reference trajectories for the JavaScript port of Model3 (digital
#               patients' designs; the Model2 cases plus cases that exercise the Model3 terms:
#               other neutrophil chains n_tr_n / nk_n incl. the Model2 chain 3 / 2, P1 th_cp > 0
#               (engine term, fixed 0 in Model3), P2 d_zpl != 0).
#           Population: the final model's fit (Results/real/hpc_runs/_chain/final)
#           when it exists, otherwise the prior centres ("population: prior").
# Access control (Model1 plan section 13, carried over by Model3 plan 6.5): real patients are never written into
#           GreenLight. Each real patient gets a clinical payload (all fields) and a
#           research payload minimised at generation time (13.5: HMAC pseudonym,
#           age band, prior-lines band, ASCT, bridging class, binned LD labs, CRS
#           grade; no exact age, regimen, drugs, design, G-CSF windows, ferritin,
#           plans or what-if inputs), written to Results/real/greenlight_stage/;
#           GreenLight/tools/seal.js encrypts them. The pseudonym key lives in
#           Results/real/greenlight_admin/ only. Digital patients stay public.
# Inputs  : config (real), the real data, optional fit.
# Outputs : GreenLight: data/model3_index.js, data/patients/DP*.js, tools/ref_model3.json;
#           Results/real/greenlight_stage: clinical/<id>.json, research/<id>.json,
#           cohort_research.json
# =============================================================================

g_schema <- "glm3-1.0"
g_public_name <- "MarrowTwin"
g_grid_days <- c(-10:30, seq(32, 90, by = 2))
g_asof_days <- c(-5:30, seq(37, 86, by = 7))

# ---- population ---------------------------------------------------------------------------
g_find_final_fit <- function(cfg) {
  runs <- file.path(cfg$paths$results, "hpc_runs")
  fin <- file.path(runs, "_chain", "final")
  if (!file.exists(fin)) return(NULL)
  f <- file.path(runs, readLines(fin)[1], "fit.rds")
  if (file.exists(f)) readRDS(f) else NULL
}

# fit-like object (spec, pop, re) of the final fit or of the prior centres
g_fit_object <- function(d, cfg, fit = NULL, model = "MabC_ddr") {
  if (!is.null(fit))
    return(list(spec = fit$spec, pop = fit$pop, re = fit$re, model = fit$model, source = "fit",
                rounds = fit$rounds, converged = fit$converged))
  spec <- model_spec(model)
  cn <- lapply(build_pdats(d, cfg, ids = d$patients$patient_id[1])[[1]]$x, names)
  list(spec = spec, pop = init_population13(spec, cfg, d$base_sd, cn), re = spec$re, model = model, source = "prior")
}

# ---- digital patients ------------------------------------------------------------------------
# 16 Model3 digital patients (MabC_ddr structure, prior-centre population): mechanism,
# status and data horizon (in-treatment patients), individual truth, ferritin kernel.
g_dp_table <- function() {
  data.table::data.table(
    patient_id = sprintf("DP%02d", 1:16),
    status = rep(c("active", "completed"), each = 8),
    horizon = c(-3, 1, 5, 8, 12, 16, 21, 27, rep(NA, 8)),
    mechanism = c("typical", "car_t_phase_second_wave", "low_reserve", "clone_ddr", "chemo_hsc_injury",
                  "car_t_target_route", "typical_strong_reserve", "inflammation_mixed",
                  "typical", "car_t_phase_second_wave", "low_reserve_delayed", "clone_ddr_chemo",
                  "chemo_and_car_t", "car_t_target_route", "strong_inflammation", "quick_young"),
    r_s   = c(0.08, 0.08, 0.02, 0.06, 0.05, 0.08, 0.15, 0.06, 0.09, 0.07, 0.012, 0.05, 0.04, 0.08, 0.06, 0.12),
    th_as = c(0.10, 0.08, 0.25, 0.30, 1.20, 0.10, 0.08, 0.15, 0.12, 0.10, 0.80, 1.00, 0.60, 0.10, 0.12, 0.05),
    th_ap = c(2.7, 2.7, 3.2, 2.7, 3.0, 2.7, 2.2, 2.9, 2.6, 2.8, 3.4, 2.9, 3.0, 2.7, 2.8, 2.0),
    th_bs = c(0.01, 0.07, 0.02, 0.02, 0.01, 0.02, 0.01, 0.05, 0.015, 0.08, 0.02, 0.02, 0.06, 0.02, 0.09, 0.005),
    f0    = c(0, 0, 0, 0.35, 0, 0, 0, 0.1, 0, 0, 0, 0.45, 0, 0, 0, 0),
    rho   = c(0.81, 0.85, 0.8, 0.8, 0.8, 1.15, 0.75, 0.9, 0.8, 0.9, 0.8, 0.8, 0.85, 1.2, 0.9, 0.7),
    A     = c(0.7, 2.6, 0.9, 0.9, 0.8, 1.2, 0.5, 1.8, 0.8, 2.8, 1.0, 0.9, 2.2, 1.3, 3.4, 0.4),
    tp    = c(7, 9, 8, 7, 8, 8, 7, 9, 7, 9, 8, 8, 9, 8, 10, 6),
    age   = c(62, 58, 74, 71, 66, 55, 45, 68, 60, 63, 77, 72, 69, 52, 64, 38),
    female = c(0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0),
    prior_lines = c(3, 2, 6, 4, 5, 2, 2, 3, 3, 2, 7, 5, 4, 2, 3, 2),
    prior_asct = c(1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0),
    bridging = c("none", "none", "cy", "none", "rgo", "none", "none", "cy", "none", "none", "rgo", "cy", "rgo", "none", "none", "none"),
    anc_base = c(3.6, 3.1, 2.2, 2.9, 3.4, 4.2, 4.8, 3.0, 3.9, 2.8, 1.9, 2.6, 3.0, 4.4, 3.3, 5.2),
    plt_base = c(190, 175, 120, 160, 185, 210, 240, 150, 200, 170, 95, 140, 165, 220, 180, 260),
    amc_base = c(0.62, 0.58, 0.45, 0.55, 0.6, 0.7, 0.75, 0.6, 0.65, 0.55, 0.4, 0.5, 0.58, 0.72, 0.6, 0.8),
    alc_base = c(0.7, 0.6, 0.4, 0.5, 0.65, 0.9, 1.1, 0.6, 0.75, 0.55, 0.35, 0.45, 0.6, 0.85, 0.7, 1.2),
    ferritin_ld = c(320, 650, 480, 400, 300, 520, 180, 900, 280, 780, 600, 450, 700, 500, 1400, 150),
    ldh_ld = c(240, 380, 260, 230, 300, 420, 190, 520, 220, 450, 280, 260, 400, 480, 600, 180),
    crp_ld = c(0.8, 2.5, 1.5, 1.0, 0.9, 2.0, 0.4, 4.0, 0.7, 3.2, 1.8, 1.2, 2.6, 2.2, 5.0, 0.3),
    hgb_ld = c(10.5, 9.8, 8.6, 9.9, 10.2, 11.0, 12.1, 9.4, 10.8, 9.5, 8.4, 9.1, 9.6, 11.3, 9.9, 12.8),
    egfr = c(88, 95, 62, 70, 85, 102, 110, 76, 92, 90, 55, 68, 80, 100, 84, 118),
    rdw_base = c(49, 52, 58, 56, 50, 47, 45, 53, 48, 51, 60, 57, 54, 46, 52, 44),
    mcv_base = c(91, 93, 97, 99, 92, 89, 88, 94, 90, 92, 100, 98, 95, 89, 93, 87))
}

# Data object of the digital patients (same format as read_model3_data) with the truth.
g_digital_patients <- function(d, cfg, seed = 20261104L) {
  set.seed(seed)
  tab <- g_dp_table()
  rp <- d$patients
  sp <- model_spec("MabC_ddr")
  st <- default_struct(sp, cfg)
  ds <- cfg$model$design
  rule <- cfg$data$gcsf_rule
  bv <- cfg$model$baseline$vars
  days_obs <- c(-28, -21, -14, -5:30, seq(33, 90, by = 3))
  out <- lapply(seq_len(nrow(tab)), function(i) {
    r <- tab[i]
    pid <- r$patient_id
    prow <- data.table::data.table(
      patient_id = pid, age = r$age, female = r$female, prior_lines = r$prior_lines, prior_asct = r$prior_asct,
      bridging = r$bridging, bridging_regimen = c(none = "none", cy = "oral cyclophosphamide", rgo = "R-Gem-Ox")[[r$bridging]],
      apheresis_day = -27, egfr = r$egfr, ferritin_ld = r$ferritin_ld, crp_ld = r$crp_ld, ldh_ld = r$ldh_ld,
      hgb_ld = r$hgb_ld, anc_ld = r$anc_base * 0.9, plt_ld = r$plt_base * 0.9, rdw_base = r$rdw_base, mcv_base = r$mcv_base,
      gcsf_after_d14 = 0L, tocilizumab = as.integer(r$A > 2), steroids = as.integer(r$A > 3),
      max_crs = if (r$A > 3) 2L else if (r$A > 1.5) 1L else 0L, last_lab_day = 90, t0 = cfg$model$t_min)
    prow[, `:=`(age_c = (age - 60) / 10, prior_lines_c = prior_lines - stats::median(rp$prior_lines),
                log_rdw_c = log(rdw_base) - mean(log(rp$rdw_base)), mcv_c = (mcv_base - mean(rp$mcv_base)) / 10,
                ldh_c = log(ldh_ld) - stats::median(log(rp$ldh_ld), na.rm = TRUE))]
    fd <- c(-5, seq(0, 30, by = 2), 35, 42)
    lr <- kernel_fun(fd, r$A, r$tp, cfg$model$kernel$h) + stats::rnorm(length(fd), 0, 0.15)
    inflam <- data.table::data.table(patient_id = pid, day = fd, ferritin = r$ferritin_ld * exp(lr), crp = NA_real_)
    drv <- driver_from_ferritin(inflam$day, inflam$ferritin, r$ferritin_ld, prow$t0, cfg)
    ind <- c(r_s = r$r_s, f0 = r$f0, th_as = r$th_as, th_ap = r$th_ap, th_bs = r$th_bs, rho = r$rho,
             a0 = patient_a0(st, prow$ldh_c))
    b_true <- c(anc = r$anc_base, plt = r$plt_base, amc = r$amc_base, alc = r$alc_base)
    des <- patient_design(prow)
    sim <- function(des) sim_daily(ind, st, sp, des, drv, cfg, days = seq(-30, 90), base = as.list(b_true))
    lat <- sim(des)
    win <- NULL
    for (course in 1:3) {
      from <- if (is.null(win)) rule$from else max(win[, 2])
      low <- lat[day >= from & day <= rule$until & anc < rule$thr, day]
      if (!length(low)) break
      win <- rbind(win, c(min(low), min(low) + rule$days))
      des$gcsf <- win
      lat <- sim(des)
      prow[, gcsf_after_d14 := 1L]
    }
    ob <- data.table::rbindlist(lapply(bv, function(v) {
      sp_rows <- days_obs < -5
      val <- ifelse(sp_rows, b_true[[v]] * exp(stats::rnorm(length(days_obs), 0, d$base_sd[[v]])),
                    lat[match(days_obs, day)][[v]] * exp(stats::rnorm(length(days_obs), 0, ds$sigma[[v]])))
      data.table::data.table(patient_id = pid, day = days_obs, var = v, value = val, cens = FALSE, upper = FALSE,
                             keep = TRUE, reason = "", setpoint = sp_rows)
    }))
    ob[var == "anc" & value < ds$cens$anc, `:=`(value = ds$cens$anc_bound, cens = TRUE, reason = "no_diff_low_wbc")]
    ob[var %in% c("amc", "alc") & value < ds$cens$amc, `:=`(value = ds$cens$amc_bound, cens = TRUE, reason = "manual_diff_zero")]
    for (v in bv) {
      x <- ob[var == v & setpoint & !cens]$value
      prow[, (paste0(v, "_base")) := stats::median(x)]
      prow[, (paste0("n_base_", v)) := length(x)]
    }
    tt <- lat$day
    pos <- tt >= 0
    truth <- data.table::data.table(patient_id = pid, mechanism = r$mechanism, r_s = r$r_s, f0 = r$f0, th_as = r$th_as,
                                    th_ap = r$th_ap, th_bs = r$th_bs, rho = r$rho, A = r$A, tp = r$tp,
                                    b_anc = r$anc_base, b_plt = r$plt_base, b_amc = r$amc_base, b_alc = r$alc_base,
                                    f_bi = 1 - exp(-r$th_bs * trapz(tt[pos], lat$I[pos])),
                                    f_a = 1 - exp(-r$th_as * trapz(tt, lat$C)),
                                    latent_phenotype = phenotype_sim(lat$day, lat$anc, cfg))
    list(prow = prow, obs = ob, inflam = inflam, drv = drv, truth = truth, latent = lat[, .(day, anc, plt, amc, alc)])
  })
  pts <- data.table::rbindlist(lapply(out, `[[`, "prow"), fill = TRUE)
  data.table::setkey(pts, patient_id)
  obs <- data.table::rbindlist(lapply(out, `[[`, "obs"))
  dd <- list(patients = pts, obs = obs, inflam = data.table::rbindlist(lapply(out, `[[`, "inflam")),
             drivers = stats::setNames(lapply(out, `[[`, "drv"), pts$patient_id),
             truth = data.table::rbindlist(lapply(out, `[[`, "truth")),
             latent = stats::setNames(lapply(out, `[[`, "latent"), pts$patient_id),
             status = tab[, .(patient_id, status, horizon)], base_sd = d$base_sd)
  dd$gcsf <- impute_gcsf_windows(obs, pts, rule)
  dd$phenotype <- phenotype_all_rules(obs, cfg)
  dd
}

# ---- plans -------------------------------------------------------------------------------------
# Illustrative plan library: decisions from the day after the as-of day L.
g_plans <- list(
  standard      = list(burden = 0.00),
  gcsf_early    = list(burden = 0.04),
  gcsf_reactive = list(burden = 0.02),
  anti_inflam   = list(burden = 0.05, scale = 0.5),
  tpo_ra        = list(burden = 0.04),
  cd34_boost    = list(burden = 0.08, add = 0.3),
  combo         = list(burden = 0.07, scale = 0.5))
g_utility <- c(delayed = 1.0, second_wave = 0.6, plt20 = 0.4)

# Simulate one draw under plan `plan` from day L: returns the daily matrix on 0..90
g_plan_sim <- function(plan, ind, st, spec, design, drv, L, b, cfg, days = 0:90) {
  dec <- max(floor(L) + 1, 1)
  des <- design
  dr <- drv
  if (plan %in% c("anti_inflam", "combo")) dr <- scale_driver(drv, max(L, 0), g_plans[[plan]]$scale)
  if (plan == "gcsf_early") des$gcsf <- rbind(des$gcsf, c(dec, dec + 5))
  if (plan == "tpo_ra") des$tpo <- rbind(des$tpo, c(max(dec, 7), 60))
  if (plan == "cd34_boost") des$cd34 <- data.frame(day = max(dec, 7), add = g_plans$cd34_boost$add)
  run <- function(des) sim_model(build_par(ind, st, spec), make_segments(des, st, max(days), cfg), dr, days, spec, cfg)
  m <- run(des)
  if (plan %in% c("gcsf_reactive", "combo")) {
    rule <- cfg$data$gcsf_rule
    for (course in 1:2) {
      if (anyNA(m[, 2])) break
      from <- max(rule$from, dec, if (!is.null(des$gcsf) && nrow(des$gcsf)) max(des$gcsf[, 2]) else -Inf)
      low <- m[m[, "t"] >= from & m[, "t"] <= 60 & m[, "anc_rel"] * b[["anc"]] < rule$thr, "t"]
      if (!length(low)) break
      des$gcsf <- rbind(des$gcsf, c(min(low), min(low) + rule$days))
      m <- run(des)
    }
  }
  m
}

# ---- one patient, every as-of day -----------------------------------------------------------------
g_q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE, names = FALSE))
g_sig <- function(x, k = 4) signif(x, k)

# recovery day and nadir of [observed ANC up to L] + [model ANC after L] (days 0..90)
g_recovery <- function(day, anc, obs_upto, rec_thr = 1.5) {
  keep <- day > obs_upto$L
  dd <- c(obs_upto$day, day[keep])
  aa <- c(ifelse(obs_upto$cens, 0.25, obs_upto$value), anc[keep])
  o <- order(dd); dd <- dd[o]; aa <- aa[o]
  nw <- which(dd >= 0 & dd <= 14)
  if (!length(nw)) return(c(rec = NA, nadir_day = NA, nadir = NA))
  nad <- nw[which.min(aa[nw])]
  r <- which(seq_along(aa) > nad & aa >= rec_thr)
  c(rec = if (length(r)) dd[r[1]] else Inf, nadir_day = dd[nad], nadir = aa[nad])
}

g_patient <- function(pid, d, fobj, kpop, cfg, days = g_asof_days, n_draws = 40L, n_cf = 16L, seed = 1L) {
  spec <- fobj$spec
  st <- fobj$pop$struct
  prow <- d$patients[patient_id == pid]
  fer <- d$inflam[patient_id == pid & !is.na(ferritin), .(day, ferritin)]
  lr <- log(fer$ferritin / prow$ferritin_ld)
  obs_all <- d$obs[patient_id == pid & keep == TRUE]
  anc <- obs_all[var == "anc"]
  plt <- obs_all[var == "plt"]
  tmpl <- par_template(st, spec)
  prev <- NULL
  res <- vector("list", length(days))
  for (li in seq_along(days)) {
    L <- days[li]
    set.seed(seed + 101L * li)
    kd <- kernel_draws(kpop, fer$day, lr, L, n_draws)
    drv_med <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, stats::median(kd[, 1]), stats::median(kd[, 2]), cfg)
    pd <- build_pdats(d, cfg, landmark = L, ids = pid, drivers = stats::setNames(list(drv_med), pid))[[1]]
    pd$seg <- make_segments(pd$design, st, max(pd$times), cfg)
    prior <- patient_prior(pd, fobj$pop, fobj$re)
    nm <- names(prior$mean)
    box <- indiv_box(prior, cfg)
    fn <- function(x) neg_log_post(stats::setNames(x, nm), pd, tmpl, spec, prior, fobj$pop$sigma, cfg)
    starts <- if (is.null(prev)) list(prior$mean) else list(prev[nm])
    if (L %in% c(0, 7, 14, 21, 30)) starts <- c(starts, list(prior$mean))
    best <- NULL
    for (s0 in starts) {
      o <- tryCatch(stats::optim(pmin(pmax(s0, box$lo), box$hi), fn, method = "L-BFGS-B", lower = box$lo,
                                 upper = box$hi, control = list(maxit = 300)), error = function(e) NULL)
      if (!is.null(o) && (is.null(best) || o$value < best$value)) best <- o
    }
    map <- stats::setNames(best$par, nm)
    prev <- map
    cv <- solve(pd_floor(stats::optimHess(map, fn), 1 / max(prior$sd^2)))
    cv <- (cv + t(cv)) / 2
    ch <- tryCatch(t(chol(cv)), error = function(e) diag(sqrt(pmax(diag(cv), 1e-10))))
    z <- matrix(stats::rnorm(n_draws * length(nm)), n_draws)
    draws <- sweep(z %*% t(ch), 2, map, `+`)
    draws <- sweep(sweep(draws, 2, box$lo, pmax), 2, box$hi, pmin)
    colnames(draws) <- nm
    draws[1, ] <- map
    upto <- anc[day <= L & day >= 0]
    obs_upto <- list(day = upto$day, value = upto$value, cens = upto$cens, L = L)
    plt_obs20 <- any(plt[day >= 14 & day <= min(30, L)]$value < 20)
    gd <- seq(-10, 90)
    traj <- array(NA_real_, c(n_draws, length(gd), 4), dimnames = list(NULL, NULL, c("anc", "plt", "amc", "alc")))
    cls <- character(n_draws)
    rec <- matrix(NA_real_, n_draws, 3)
    for (j in seq_len(n_draws)) {
      psi <- draws[j, ]
      ind <- eta_to_indiv(psi, spec)
      ind[["a0"]] <- patient_a0(st, pd$ldh_c)
      b <- pdat_baselines(psi, pd)
      drv_j <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, kd[j, 1], kd[j, 2], cfg)
      m <- tryCatch(sim_model(build_par(ind, st, spec), make_segments(pd$design, st, 90, cfg), drv_j, gd, spec, cfg),
                    error = function(e) NULL)
      if (is.null(m) || anyNA(m[, 2])) next
      traj[j, , "anc"] <- m[, "anc_rel"] * b[["anc"]]
      traj[j, , "plt"] <- m[, "plt_rel"] * b[["plt"]]
      traj[j, , "amc"] <- m[, "amc_rel"] * b[["amc"]]
      traj[j, , "alc"] <- m[, "alc_rel"] * b[["alc"]] + m[, "cart"]
      sel <- gd >= 0 & gd <= 30
      cls[j] <- phenotype_with_observed(gd[sel], traj[j, sel, "anc"], obs_upto, cfg)
      rec[j, ] <- g_recovery(gd[gd >= 0], traj[j, gd >= 0, "anc"], obs_upto)
    }
    ok <- nzchar(cls)
    gi <- match(g_grid_days, gd)
    band <- lapply(c(anc = "anc", plt = "plt", amc = "amc", alc = "alc"), function(v) {
      x <- traj[ok, gi, v, drop = FALSE][, , 1]
      list(g_sig(apply(x, 2, g_q, 0.1), 3), g_sig(apply(x, 2, g_q, 0.5), 3), g_sig(apply(x, 2, g_q, 0.9), 3))
    })
    p_cls <- vapply(cfg$phenotype$classes, function(k) mean(cls[ok] == k), numeric(1))
    recd <- rec[ok, 1]
    # counterfactual readouts and plans on the first n_cf draws (each with its own driver)
    cf <- list()
    plans <- lapply(names(g_plans), function(k) list(cls = character(0), plt20 = logical(0)))
    names(plans) <- names(g_plans)
    for (j in which(ok)[seq_len(min(n_cf, sum(ok)))]) {
      psi <- draws[j, ]
      drv_j <- landmark_driver(fer, prow$ferritin_ld, prow$t0, L, kd[j, 1], kd[j, 2], cfg)
      x <- tryCatch(counterfactual_draw13(psi, pd, fobj, cfg, days = 0:30, obs_upto = obs_upto, drv = drv_j),
                    error = function(e) NULL)
      if (!is.null(x)) cf[[length(cf) + 1]] <- x
      ind <- eta_to_indiv(psi, spec)
      ind[["a0"]] <- patient_a0(st, pd$ldh_c)
      b <- pdat_baselines(psi, pd)
      for (k in names(g_plans)) {
        m <- tryCatch(g_plan_sim(k, ind, st, spec, pd$design, drv_j, L, b, cfg), error = function(e) NULL)
        if (is.null(m) || anyNA(m[, 2])) next
        a <- m[, "anc_rel"] * b[["anc"]]
        sel <- m[, "t"] <= 30
        plans[[k]]$cls <- c(plans[[k]]$cls, phenotype_with_observed(m[sel, "t"], a[sel], obs_upto, cfg))
        pm <- m[, "plt_rel"] * b[["plt"]]
        plans[[k]]$plt20 <- c(plans[[k]]$plt20, plt_obs20 || any(pm[m[, "t"] >= max(14, L + 1) & m[, "t"] <= 30] < 20))
      }
    }
    cft <- if (length(cf)) data.table::as.data.table(do.call(rbind, cf)) else NULL
    when <- who <- NULL
    fl <- NULL
    if (!is.null(cft)) {
      w <- when_class(cft$actual, cft$no_a, cft$no_bS, cft$no_abS)
      when <- vapply(c("none", "a", "b", "ab_joint", "ab_redundant", "other"), function(k) mean(w == k), numeric(1))
      h <- who_class(cft$actual, cft$no_clone, cft$reserve)
      nq <- !is.na(h)
      who <- vapply(c("reserve", "clone", "both", "other"), function(k) if (any(nq)) mean(h[nq] == k) else 0, numeric(1))
      fl <- c(a = stats::median(as.numeric(cft$f_a)), bi = stats::median(as.numeric(cft$f_bi)),
              bz = stats::median(as.numeric(cft$f_bz)), f30 = stats::median(as.numeric(cft$f30)))
    }
    pl <- lapply(names(g_plans), function(k) {
      cl <- plans[[k]]$cls
      p <- vapply(cfg$phenotype$classes, function(c) if (length(cl)) mean(cl == c) else NA_real_, numeric(1))
      p20 <- if (length(cl)) mean(plans[[k]]$plt20) else NA_real_
      score <- -(g_utility[["delayed"]] * p[["delayed"]] + g_utility[["second_wave"]] * p[["second_wave"]] +
                   g_utility[["plt20"]] * p20) - g_plans[[k]]$burden
      list(id = k, p = g_sig(unname(p), 3), plt20 = g_sig(p20, 3), score = g_sig(score, 3), n = length(cl))
    })
    ord <- order(-vapply(pl, function(x) if (is.na(x$score)) -Inf else x$score, numeric(1)))
    tg <- drv_med$t0 + (seq_along(drv_med$v) - 1) * drv_med$dt
    dsel <- abs(tg - round(tg)) < 1e-9
    res[[li]] <- list(
      L = L, n = sum(ok),
      map = as.list(g_sig(map, 6)), chol = g_sig(ch[lower.tri(ch, diag = TRUE)], 6), names = nm,
      lo = unname(g_sig(box$lo, 6)), hi = unname(g_sig(box$hi, 6)),
      kern = g_sig(c(stats::median(kd[, 1]), stats::median(kd[, 2])), 4),
      drv = list(t0 = tg[dsel][1], dt = 1, v = g_sig(drv_med$v[dsel], 4)),
      p = g_sig(unname(p_cls), 3),
      rec = g_sig(c(stats::median(recd), g_q(recd, 0.1), g_q(recd, 0.9), mean(is.finite(recd) & recd <= 30)), 3),
      nadir = g_sig(c(stats::median(rec[ok, 2]), stats::median(rec[ok, 3])), 3),
      band = band, when = if (is.null(when)) NULL else as.list(g_sig(when, 3)),
      who = if (is.null(who)) NULL else as.list(g_sig(who, 3)), f = if (is.null(fl)) NULL else as.list(g_sig(fl, 3)),
      plans = pl[ord], top3 = vapply(pl[ord][1:3], `[[`, character(1), "id"))
  }
  res
}

# ---- patient record (observations, design, summary) ---------------------------------------------
g_patient_record <- function(pid, d, cfg, type, status, horizon = NA, truth = NULL) {
  p <- d$patients[patient_id == pid]
  o <- d$obs[patient_id == pid]
  obs <- lapply(c(anc = "anc", plt = "plt", amc = "amc", alc = "alc"), function(v) {
    x <- o[var == v][order(day)]
    lapply(seq_len(nrow(x)), function(i) list(x$day[i], g_sig(x$value[i], 4), as.integer(x$cens[i]),
                                               as.integer(x$keep[i]), as.integer(isTRUE(x$setpoint[i]))))
  })
  fer <- d$inflam[patient_id == pid & !is.na(ferritin)][order(day)]
  gw <- if (!is.null(d$gcsf) && nrow(d$gcsf)) d$gcsf[patient_id == pid] else NULL
  ph <- d$phenotype[patient_id == pid]
  hx <- car_hematotox(p$plt_ld, p$anc_ld, p$hgb_ld, p$crp_ld, p$ferritin_ld)
  list(id = pid, type = type, status = status, horizon = if (is.na(horizon)) NULL else horizon,
       summary = list(age = p$age, female = p$female, prior_lines = p$prior_lines, prior_asct = p$prior_asct,
                      bridging = p$bridging, bridging_regimen = p$bridging_regimen, apheresis_day = p$apheresis_day,
                      egfr = g_sig(p$egfr, 3), ferritin_ld = p$ferritin_ld, crp_ld = p$crp_ld, ldh_ld = p$ldh_ld,
                      hgb_ld = p$hgb_ld, plt_ld = p$plt_ld, anc_ld = p$anc_ld, max_crs = p$max_crs,
                      tocilizumab = p$tocilizumab, steroids = p$steroids, gcsf_after_d14 = p$gcsf_after_d14,
                      hematotox = if (is.na(hx)) NULL else hx),
       base = list(anc = g_sig(p$anc_base, 4), plt = g_sig(p$plt_base, 4), amc = g_sig(p$amc_base, 4),
                   alc = g_sig(p$alc_base, 4)),
       design = list(t0 = p$t0, egfr = g_sig(p$egfr, 6), bridging = p$bridging, apheresis_day = p$apheresis_day,
                     ldh_c = g_sig(p$ldh_c, 6),
                     gcsf = if (is.null(gw) || !nrow(gw)) list() else lapply(seq_len(nrow(gw)), function(i) c(gw$start[i], gw$end[i]))),
       cov = list(age_c = p$age_c, prior_lines_c = p$prior_lines_c, prior_asct = p$prior_asct,
                  log_rdw_c = g_sig(p$log_rdw_c, 4), mcv_c = g_sig(p$mcv_c, 4)),
       obs = obs, ferritin = lapply(seq_len(nrow(fer)), function(i) c(fer$day[i], g_sig(fer$ferritin[i], 4))),
       label = if (nrow(ph)) as.list(ph[, -1]) else NULL,
       truth = if (is.null(truth)) NULL else as.list(truth))
}

# ---- index pieces --------------------------------------------------------------------------------
g_param_table <- function(cfg, fobj) {
  st <- fobj$pop$struct
  big <- model_spec("MabxC_ddr")
  est <- struct_est_names(big, cfg)
  pr <- struct_prior_work13(est, big, cfg)
  box <- struct_box(pr, cfg)
  rows <- lapply(names(cfg$model$struct), function(k) {
    e <- k %in% est
    list(id = k, value = g_sig(unname(st[[k]]), 4), prior = g_sig(unname(cfg$model$struct[[k]]), 4),
         kind = if (e) "estimated" else "fixed",
         scale = if (k %in% cfg$model$unit_interval) "logit" else if (k %in% cfg$model$real_line) "identity" else "log",
         prior_sd = if (e) g_sig(unname(pr$sd[[k]]), 3) else NULL,
         box = if (e) g_sig(unname(work_to_struct(stats::setNames(c(box$lo[[k]], box$hi[[k]]), c(k, k)), cfg)), 3) else NULL)
  })
  re <- lapply(names(cfg$model$re_prior$mu)[1:6], function(k) {
    mu <- cfg$model$re_prior$mu[[k]]
    list(id = k, kind = "individual", link = if (k == "f0") "logit" else "log",
         centre = g_sig(if (k == "f0") stats::plogis(mu) else exp(mu), 4),
         omega = cfg$model$re_prior$omega[[k]], box = cfg$model$indiv_limits[[k]],
         covariates = cfg$model$covariates[[k]],
         value = if (fobj$source == "fit" && k %in% names(fobj$pop$mu))
           g_sig(if (k == "f0") stats::plogis(fobj$pop$mu[[k]]) else exp(fobj$pop$mu[[k]]), 4) else NULL)
  })
  list(struct = rows, re = re)
}

g_k1 <- function(cfg) {
  tdir <- file.path(cfg$paths$results, "tables")
  rd <- function(f) { p <- file.path(tdir, f); if (file.exists(p)) data.table::fread(p) else NULL }
  # public name in what GreenLight shows (MarrowTwin; the predecessor engine is not named)
  pub <- function(x) {
    x <- gsub("frozen Model1 v1.0 engine", "frozen first-version engine", x, fixed = TRUE)
    x <- gsub("frozen Model2 v1.0 engine", "frozen previous-version engine", x, fixed = TRUE)
    x <- gsub("the Model2 values", "the previous-version values", x, fixed = TRUE)
    x <- gsub("\\bModel[12]\\b", "an earlier version", x, perl = TRUE)
    gsub("\\bModel3\\b", g_public_name, x, perl = TRUE)
  }
  cmp <- rd("K1_aim4_comparison.csv")
  if (!is.null(cmp) && "engine" %in% names(cmp)) cmp[, engine := pub(engine)]
  ut <- rd("K1_unit_tests.csv")
  if (!is.null(ut) && "test" %in% names(ut)) ut[, test := pub(test)]
  list(unit_tests = ut, calibration = rd("K1_calibration.csv"),
       aim4_comparison = cmp, stability = rd("K1_stability.csv"),
       cart_calibration = rd("K1_cart_calibration.csv"))
}

# step S tables (synthetic cohorts S1-S4: summaries only; S5: cohort-level counts), for the accuracy tab
g_s_results <- function(cfg) {
  tabs <- c("S1_struct", "S1_indiv", "S1_loss", "S1_convergence", "S2_struct", "S2_indiv", "S3_auroc",
            "S4_indiv", "S4_loss", "S4_convergence", "S5_vpc")
  out <- lapply(stats::setNames(tabs, tabs), function(n) {
    f <- file.path(cfg$paths$results, "tables", paste0("S_", n, ".csv"))
    if (!file.exists(f)) return(NULL)
    x <- data.table::fread(f)
    # (ratios of a parameter whose truth is 0, e.g. d_zpl, are not finite: written as null)
    lapply(seq_len(nrow(x)), function(i) lapply(as.list(x[i]), function(v)
      if (is.double(v)) { if (is.finite(v)) g_sig(v, 3) else NA } else v))
  })
  out <- Filter(Negate(is.null), out)
  if (length(out)) out else NULL
}

# models with an assembled fit (copied back from the cluster), for the researcher's accuracy tab
g_results_status <- function(cfg) {
  runs <- file.path(cfg$paths$results, "hpc_runs")
  ids <- if (dir.exists(runs)) list.dirs(runs, recursive = FALSE, full.names = FALSE) else character(0)
  # group: only candidates and drop-one runs share the likelihood and data of the final model (ranked by
  # BIC); sensitivity runs change the phenotype term, fold refits the data, MF has no ALC; synthetic
  # cohorts (S1-S4) have other data altogether
  run_group <- function(id) if (id == "MF") "reference" else if (grepl("^(a4)?syn[0-9]+", id)) "synthetic" else
    if (grepl("_fold[0-9]+$", id)) "fold" else
    if (grepl("_(s|w)[0-9.]+$", id)) "sensitivity" else if (grepl("_no(N|P1|P2)$", id)) "drop_one" else "candidate"
  fits <- lapply(ids, function(m) {
    f <- file.path(runs, m, "fit.rds")
    if (!file.exists(f)) return(NULL)
    x <- readRDS(f)
    # (a provisional fold refit, built from a population for a preliminary V1, has no BIC)
    sig <- function(v) if (is.numeric(v)) g_sig(v, 6) else NULL
    list(id = m, model = x$model, BIC = sig(x$BIC), BIC_core = sig(x$BIC_core), n_par = x$n_par,
         rounds = x$rounds, converged = isTRUE(x$converged), group = run_group(m), provisional = isTRUE(x$provisional))
  })
  ff <- file.path(runs, "_chain", "final")
  final <- if (file.exists(ff)) readLines(ff)[1] else NULL
  # V1 (F3): cohort-level metrics only; provisional while a fold refit of the final model is provisional
  v1 <- {
    fm <- file.path(cfg$paths$results, "tables", "F3_metrics.csv")
    fp <- file.path(cfg$paths$results, "tables", "F3_predictions.csv")
    if (file.exists(fm) && file.exists(fp)) {
      m <- data.table::fread(fm)
      obs <- unique(data.table::fread(fp)[, .(patient_id, observed)])
      pk <- obs[, .N, by = observed]$N / nrow(obs)
      fold_dirs <- if (is.null(final)) character(0) else file.path(runs, paste0(final, "_fold", seq_len(cfg$validation$V1$k)))
      list(metrics = lapply(seq_len(nrow(m)), function(i) lapply(as.list(m[i]), function(x) if (is.numeric(x)) g_sig(x, 3) else x)),
           n = nrow(obs), n_event = sum(obs$observed != "quick"), brier_ref = g_sig(1 - sum(pk^2), 3),
           provisional = any(file.exists(file.path(fold_dirs, "PROVISIONAL.txt"))))
    } else NULL
  }
  list(fits = Filter(Negate(is.null), fits),
       final = final, v1 = v1, s = g_s_results(cfg),
       prereg = file.exists(file.path(cfg$paths$results, "prereg", "PR.sha256")),
       cohorts = {
         cdir <- file.path(cfg$paths$results, "recovery", "cohorts")
         f <- file.path(cdir, "cohorts_overview.csv")
         # (no overview table: phenotype counts from the cohort files themselves)
         if (file.exists(f)) data.table::fread(f) else {
           cf <- list.files(cdir, pattern = "^(a4)?syn[0-9]+\\.rds$", full.names = TRUE)
           if (!length(cf)) NULL else data.table::rbindlist(lapply(cf, function(p) {
             co <- readRDS(p)
             ph <- table(factor(co$phenotype$primary, levels = cfg$phenotype$classes))
             data.table::data.table(cohort = sub("\\.rds$", "", basename(p)), generator = co$generator,
                                    quick = ph[["quick"]], second_wave = ph[["second_wave"]], delayed = ph[["delayed"]])
           }))
         }
       })
}

g_data_description <- function(d, cfg) {
  o <- d$obs[keep == TRUE & setpoint == FALSE & day <= cfg$data$obs_max_day]
  per_var <- o[, .(n = .N, censored = sum(cens), patients = data.table::uniqueN(patient_id)), by = var]
  excl <- d$obs[keep == FALSE, .N, by = reason]
  cov <- d$obs[var == "anc" & day >= -30 & day <= 90, .(days = list(sort(unique(day)))), by = patient_id]
  ph <- table(factor(d$phenotype$primary, levels = cfg$phenotype$classes))
  list(n_patients = nrow(d$patients), per_var = per_var, excluded = excl,
       gcsf_windows = if (is.null(d$gcsf)) 0L else nrow(d$gcsf),
       gcsf_patients = if (is.null(d$gcsf)) 0L else data.table::uniqueN(d$gcsf$patient_id),
       phenotype = as.list(stats::setNames(as.integer(ph), names(ph))),
       bridging = as.list(table(d$patients$bridging)),
       coverage = stats::setNames(lapply(cov$days, as.integer), cov$patient_id),
       rules = cfg$data[c("wbc_no_diff", "manual_diff_cells", "gcsf_rule", "plt_tx", "setpoint_window", "obs_max_day")],
       phenotype_rule = cfg$phenotype$primary)
}

# ridge-logistic example module (counts at the landmark -> P(non-quick)), fitted in sample on the
# real cohort; for the ML-module interface demonstration only (not a validated predictor)
g_logistic_module <- function(d, cfg, landmarks = c(0, 7, 10, 14, 21), lambda = 1) {
  y <- as.integer(d$phenotype[match(d$patients$patient_id, patient_id)]$primary != "quick")
  lapply(landmarks, function(L) {
    X <- t(vapply(d$patients$patient_id, function(pid) {
      vapply(c("anc", "plt"), function(v) {
        x <- d$obs[patient_id == pid & var == v & keep == TRUE & day <= L & day >= L - 3]
        if (!nrow(x)) x <- d$obs[patient_id == pid & var == v & keep == TRUE & day <= L]
        log(max(utils::tail(x[order(day)]$value, 1), 0.01))
      }, numeric(1))
    }, numeric(2)))
    mu <- colMeans(X); sdv <- apply(X, 2, stats::sd)
    Z <- sweep(sweep(X, 2, mu), 2, sdv, `/`)
    nll <- function(b) { eta <- b[1] + Z %*% b[-1]; -sum(y * eta - log1p(exp(eta))) + lambda * sum(b[-1]^2) / 2 }
    b <- stats::optim(c(0, 0, 0), nll, method = "BFGS")$par
    list(L = L, intercept = g_sig(b[1], 5), coef = g_sig(b[-1], 5), centre = g_sig(unname(mu), 5),
         scale = g_sig(unname(sdv), 5), vars = c("log_anc", "log_plt"))
  })
}

# ---- reference trajectories for the JavaScript port ------------------------------------------------
g_ref_trajectories <- function(d, cfg, fobj, ids = d$patients$patient_id[c(1, 2, 4)]) {   # d: digital patients (public designs)
  st <- fobj$pop$struct
  cases <- list(
    list(name = "M0 typical", model = "M0", ind = c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = 2.7, th_bs = 0, rho = 0.81)),
    list(name = "MabC_ddr clone + G-CSF + CD34 + TPO", model = "MabC_ddr",
         ind = c(r_s = 0.04, f0 = 0.3, th_as = 0.4, th_ap = 3, th_bs = 0.05, rho = 0.9), extra = TRUE),
    list(name = "Mabx_Zm synergy", model = "Mabx_Zm", ind = c(r_s = 0.06, f0 = 0, th_as = 0.2, th_ap = 2.5, th_bs = 0.04, rho = 1.1)),
    list(name = "Mb_Zh strong CAR T", model = "Mb_Zh", ind = c(r_s = 0.1, f0 = 0, th_as = 0, th_ap = 2.2, th_bs = 0.08, rho = 1.3)),
    list(name = "MabC_ddr strong chemo kill", model = "MabC_ddr", ind = c(r_s = 0.05, f0 = 0.2, th_as = 0.3, th_ap = 3.2, th_bs = 0.03, rho = 0.9),
         struct = c(kappa_c = 2.5, th_am = 1.2)),
    list(name = "Mab nested (kappa_c = th_am = 0)", model = "Mab", ind = c(r_s = 0.08, f0 = 0, th_as = 0.2, th_ap = 2.7, th_bs = 0.04, rho = 0.9),
         struct = c(kappa_c = 0, th_am = 0)),
    # Model3 terms: neutrophil chain variants (N), platelet consumption (P1, engine only), platelet e_Z (P2)
    list(name = "MabC_ddr N 8/6 + P1 + P2 + G-CSF", model = "MabC_ddr", ind = c(r_s = 0.05, f0 = 0.25, th_as = 0.3, th_ap = 2.9, th_bs = 0.05, rho = 1),
         struct = c(n_tr_n = 8, nk_n = 6, th_cp = 0.3, d_zpl = 1.5), extra = TRUE),
    list(name = "Mab_Zm Model2 chain 3/2, d_zpl -1", model = "Mab_Zm", ind = c(r_s = 0.07, f0 = 0, th_as = 0.15, th_ap = 2.6, th_bs = 0.03, rho = 0.9),
         struct = c(n_tr_n = 3, nk_n = 2, d_zpl = -1)),
    list(name = "Mb N 5/1, MTT_N 12, gamma 0.9, P1", model = "Mb", ind = c(r_s = 0.09, f0 = 0, th_as = 0, th_ap = 2.4, th_bs = 0.06, rho = 1.1),
         struct = c(n_tr_n = 5, nk_n = 1, mtt_n = 12, gamma = 0.9, th_cp = 0.15, d_zpl = 2.5)),
    list(name = "Mab nested Model2 (N 3/2, th_cp = d_zpl = 0)", model = "Mab", ind = c(r_s = 0.08, f0 = 0, th_as = 0.2, th_ap = 2.7, th_bs = 0.04, rho = 0.9),
         struct = c(n_tr_n = 3, nk_n = 2, th_cp = 0, d_zpl = 0)))
  lapply(seq_along(cases), function(i) {
    cs <- cases[[i]]
    pid <- ids[(i - 1) %% length(ids) + 1]
    p <- d$patients[patient_id == pid]
    sp <- model_spec(cs$model)
    s2 <- default_struct(sp, cfg)
    s2[names(st)[names(st) %in% names(s2)]] <- st[names(st)[names(st) %in% names(s2)]]
    if (sp$clone != "none") s2[names(cfg$model$clone_types[[sp$clone]])] <- cfg$model$clone_types[[sp$clone]]
    if (sp$q2 == "Mabx") s2[["eta_syn"]] <- 1.5
    if (!is.null(cs$struct)) s2[names(cs$struct)] <- cs$struct
    des <- patient_design(p, d$gcsf)
    if (isTRUE(cs$extra)) { des$gcsf <- rbind(des$gcsf, c(18, 23)); des$tpo <- rbind(des$tpo, c(10, 40)); des$cd34 <- data.frame(day = 12, add = 0.3) }
    ind <- c(cs$ind, a0 = patient_a0(s2, p$ldh_c))
    drv <- d$drivers[[pid]]
    days <- seq(ceiling(des$t0), 90)
    seg <- make_segments(des, s2, max(days), cfg)
    par <- build_par(ind, s2, sp)
    m <- m3_sim(par, seg, drv$v, drv$t0, drv$dt, days, cfg$model$solver$rtol, cfg$model$solver$atol, cfg$model$solver$hmax)
    list(name = cs$name, model = cs$model, par = as.list(par), struct = as.list(s2),
         design = list(t0 = des$t0, egfr = des$egfr, bridging = des$bridging, apheresis_day = des$apheresis_day,
                       gcsf = if (is.null(des$gcsf)) list() else lapply(seq_len(nrow(des$gcsf)), function(k) des$gcsf[k, ]),
                       tpo = if (is.null(des$tpo)) list() else lapply(seq_len(nrow(des$tpo)), function(k) des$tpo[k, ]),
                       cd34 = if (is.null(des$cd34)) list() else lapply(seq_len(nrow(des$cd34)), function(k) as.list(des$cd34[k, ]))),
         drv = list(t0 = drv$t0, dt = drv$dt, v = drv$v), days = days,
         out = lapply(as.data.frame(m[, c("anc_rel", "plt_rel", "amc_rel", "alc_rel", "S", "f_clone", "M", "cart", "z", "I")]), identity))
  })
}

# ---- access control: staging, pseudonyms, minimisation (plan 13.5 / 13.6) ------------------------
g_access_dirs <- function(cfg) list(stage = file.path(cfg$paths$results, "greenlight_stage"),
                                    admin = file.path(cfg$paths$results, "greenlight_admin"))

# steward-held HMAC key (32 random bytes, hex); created once, never copied into GreenLight
g_pseudonym_key <- function(cfg) {
  f <- file.path(g_access_dirs(cfg)$admin, "pseudonym_key.hex")
  if (!file.exists(f)) {
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines(paste(as.character(readBin("/dev/urandom", "raw", 32L)), collapse = ""), f)
    Sys.chmod(f, "600")
  }
  readLines(f, warn = FALSE)[1]
}
g_pseudonym <- function(id, key_hex) {
  key <- as.raw(strtoi(substring(key_hex, seq(1, 63, 2), seq(2, 64, 2)), 16L))
  paste0("R-", toupper(substr(digest::hmac(key, id, algo = "sha256", serialize = FALSE), 1, 6)))
}

# cohort tertiles of the LD labs (research bins); cut points are cohort summaries
g_research_bins <- function(patients) {
  labs <- c("ferritin_ld", "ldh_ld", "crp_ld", "hgb_ld", "plt_ld", "anc_ld")
  stats::setNames(lapply(labs, function(k) unname(stats::quantile(patients[[k]], c(1 / 3, 2 / 3), na.rm = TRUE))), labs)
}
g_band <- function(x, cuts, labels) if (is.null(x) || is.na(x)) NA_character_ else labels[findInterval(x, cuts) + 1]
g_research_summary <- function(s, bins) {
  age <- if (is.na(s$age)) NA_character_ else if (s$age >= 90) "90+" else sprintf("%d-%d", 10 * (s$age %/% 10), 10 * (s$age %/% 10) + 9)
  list(sex = if (isTRUE(s$female == 1)) "F" else "M", age_band = age,
       lines_band = g_band(s$prior_lines, c(2.5, 4.5), c("1-2", "3-4", ">=5")),
       asct = isTRUE(s$prior_asct == 1), bridging_class = if (identical(s$bridging, "none")) "none" else "chemo",
       crs = s$max_crs,
       labs = c(list(egfr = g_band(s$egfr, c(60, 90), c("<60", "60-89", ">=90"))),
                stats::setNames(lapply(names(bins), function(k) g_band(s[[k]], bins[[k]], c("T1", "T2", "T3"))),
                                sub("_ld$", "", names(bins)))))
}
# research payload: generated from the clinical record, exact fields never copied
g_research_record <- function(rec, pseudo, bins) {
  keep <- c("L", "n", "map", "chol", "names", "lo", "hi", "p", "rec", "nadir", "band", "when", "who", "f")
  list(kind = "research", id = pseudo, type = "real", status = "completed", summary = g_research_summary(rec$summary, bins),
       base = rec$base, obs = rec$obs, label = rec$label, asof = lapply(rec$asof, function(a) a[intersect(keep, names(a))]))
}
# small cells (1-4) of research aggregates are shown as "<5"
g_mask_small <- function(n) ifelse(n >= 1 & n < 5, "<5", as.character(n))
g_cohort_research <- function(d, cfg, pseudo, bins) {
  o <- d$obs[keep == TRUE & setpoint == FALSE & day <= cfg$data$obs_max_day]
  per_var <- o[, .(n = .N, censored = sum(cens), patients = data.table::uniqueN(patient_id)), by = var]
  excl <- d$obs[keep == FALSE, .N, by = reason]
  cls <- cfg$phenotype$classes
  ph <- factor(d$phenotype[match(d$patients$patient_id, patient_id)]$primary, levels = cls)
  bc <- factor(ifelse(d$patients$bridging == "none", "none", "chemo"), levels = c("chemo", "none"))
  lb <- factor(vapply(d$patients$prior_lines, function(x) g_band(x, c(2.5, 4.5), c("1-2", "3-4", ">=5")), ""), levels = c("1-2", "3-4", ">=5"))
  xt <- function(f) { tb <- table(ph, f); lapply(stats::setNames(colnames(tb), colnames(tb)), function(j) as.list(stats::setNames(g_mask_small(tb[, j]), rownames(tb)))) }
  cov <- d$obs[var == "anc" & day >= -30 & day <= 90, .(days = list(sort(unique(day)))), by = patient_id]
  rows <- lapply(order(pseudo[cov$patient_id]), function(i) list(id = unname(pseudo[cov$patient_id[i]]),
                                                                 label = as.character(ph[match(cov$patient_id[i], d$patients$patient_id)]),
                                                                 days = as.integer(cov$days[[i]])))
  list(kind = "cohort", n_patients = nrow(d$patients), per_var = per_var,
       excluded = excl[, .(reason, N = g_mask_small(N))], gcsf_windows = g_mask_small(if (is.null(d$gcsf)) 0L else nrow(d$gcsf)),
       phenotype = as.list(stats::setNames(g_mask_small(as.integer(table(ph))), cls)),
       bridging_class = as.list(stats::setNames(g_mask_small(as.integer(table(bc))), levels(bc))),
       crosstab = list(bridging_class = xt(bc), lines_band = xt(lb)),
       coverage = rows,
       bins = list(egfr = c(60, 90), tertiles = stats::setNames(lapply(bins, function(x) signif(x, 3)), sub("_ld$", "", names(bins)))),
       rules = cfg$data[c("wbc_no_diff", "manual_diff_cells", "gcsf_rule", "plt_tx", "setpoint_window", "obs_max_day")],
       phenotype_rule = cfg$phenotype$primary)
}
g_write_json <- function(obj, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  jsonlite::write_json(obj, path, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  Sys.chmod(path, "600")
}

g_write_js <- function(obj, path, prefix, suffix = ";") {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  js <- jsonlite::toJSON(obj, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  writeLines(paste0(prefix, js, suffix), path, useBytes = TRUE)
}

# ---- main ------------------------------------------------------------------------------------------
# ids: optional subset. Real patients -> Results/real/greenlight_stage (clinical + research payloads);
# run GreenLight/tools/seal.js afterwards to encrypt them for the users (plan 13.6).
g_export <- function(cfg, cores = cfg$runtime$cores, gl_dir = cfg$paths$greenlight, ids = NULL, fit = g_find_final_fit(cfg),
                     log_fun = message, n_draws = 40L, n_cf = 16L) {
  d <- read_model3_data(cfg)
  fobj <- g_fit_object(d, cfg, fit)
  dp <- g_digital_patients(d, cfg)
  kpop <- fit_kernel_population(d$drivers, cfg)
  dirs <- g_access_dirs(cfg)
  pkey <- g_pseudonym_key(cfg)
  pseudo <- vapply(stats::setNames(d$patients$patient_id, d$patients$patient_id), g_pseudonym, "", key_hex = pkey)
  if (anyDuplicated(pseudo)) stop("pseudonym collision: lengthen the pseudonym")
  bins <- g_research_bins(d$patients)
  saveRDS(dp, file.path(gl_dir, "tools", "digital_patients.rds"))
  # one data object for both groups (patient ids are distinct)
  both <- list(patients = rbind(d$patients, dp$patients, fill = TRUE), obs = rbind(d$obs, dp$obs, fill = TRUE),
               inflam = rbind(d$inflam, dp$inflam, fill = TRUE), drivers = c(d$drivers, dp$drivers),
               gcsf = rbind(d$gcsf, dp$gcsf, fill = TRUE), phenotype = rbind(d$phenotype, dp$phenotype, fill = TRUE),
               base_sd = d$base_sd)
  data.table::setkey(both$patients, patient_id)
  all_ids <- c(d$patients$patient_id, dp$patients$patient_id)
  todo <- if (is.null(ids)) all_ids else intersect(ids, all_ids)
  t0 <- Sys.time()
  log_fun("GreenLight export: ", length(todo), " patients, population ", fobj$source, " (", fobj$model, "), ", cores, " cores")
  done <- par_lapply(seq_along(todo), function(i) {
    pid <- todo[i]
    is_dp <- startsWith(pid, "DP")
    st <- if (is_dp) dp$status[patient_id == pid] else NULL
    rec <- g_patient_record(pid, both, cfg, if (is_dp) "digital" else "real",
                            if (is_dp) st$status else "completed", if (is_dp) st$horizon else NA,
                            if (is_dp) dp$truth[patient_id == pid] else NULL)
    t1 <- Sys.time()
    rec$asof <- tryCatch(g_patient(pid, both, fobj, kpop, cfg, n_draws = n_draws, n_cf = n_cf,
                                   seed = cfg$runtime$seed + match(pid, all_ids)),
                         error = function(e) { rec$error <<- conditionMessage(e); list() })
    rec$runtime_s <- round(as.numeric(difftime(Sys.time(), t1, units = "secs")))
    if (is_dp) {
      g_write_js(rec, file.path(gl_dir, "data", "patients", paste0(pid, ".js")), "GLM3.addPatient(", ");")
    } else {
      rec$kind <- "clinical"
      g_write_json(rec, file.path(dirs$stage, "clinical", paste0(pid, ".json")))
      g_write_json(g_research_record(rec, pseudo[[pid]], bins), file.path(dirs$stage, "research", paste0(pid, ".json")))
    }
    list(id = pid, runtime = rec$runtime_s, error = rec$error)
  }, cores)
  log_fun("  patients done in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  # no real-patient plaintext inside GreenLight (earlier exports wrote P*.js there)
  unlink(list.files(file.path(gl_dir, "data", "patients"), pattern = "^P[0-9]+\\.js$", full.names = TRUE))
  g_write_json(g_cohort_research(d, cfg, pseudo, bins), file.path(dirs$stage, "cohort_research.json"))
  g_write_index(cfg, gl_dir, d, dp, fobj)
  ref <- g_ref_trajectories(dp, cfg, fobj, ids = dp$patients$patient_id[c(2, 5, 13)])
  jsonlite::write_json(ref, file.path(gl_dir, "tools", "ref_model3.json"), auto_unbox = TRUE, digits = NA)
  log_fun("  real patients staged in ", dirs$stage, "; run: node GreenLight/tools/seal.js")
  invisible(data.table::rbindlist(lapply(done, function(x) data.table::data.table(id = x$id, runtime_s = x$runtime,
                                                                               error = x$error %||% NA_character_))))
}

# public index: digital roster + aggregate model information (no real-patient data)
g_write_index <- function(cfg, gl_dir, d, dp, fobj) {
  pdir <- file.path(gl_dir, "data", "patients")
  files <- list.files(pdir, pattern = "^DP[0-9]+\\.js$", full.names = TRUE)
  roster <- lapply(files, function(f) {
    x <- jsonlite::fromJSON(sub(");\\s*$", "", sub("^GLM3.addPatient\\(", "", paste(readLines(f, warn = FALSE), collapse = ""))),
                            simplifyVector = FALSE)
    ev <- lapply(x$asof, function(a) c(a$L, unlist(a$p)))
    # gcsf: the G-CSF windows, so the roster can list who received G-CSF (for a patient in treatment: before its data horizon)
    list(id = x$id, type = x$type, status = x$status, horizon = x$horizon, summary = x$summary,
         label = x$label$primary, mechanism = x$truth$mechanism, evolution = ev, gcsf = x$design$gcsf %||% list(),
         n_obs = sum(vapply(x$obs, length, integer(1))), error = x$error)
  })
  st <- fobj$pop$struct
  idx <- list(
    schema = g_schema, generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), model_name = g_public_name,
    model_version = cfg$version, model_label = sprintf("%s (%s)", g_public_name, cfg$version),
    population = list(source = fobj$source, model = fobj$model, rounds = fobj$rounds,
                      struct = as.list(g_sig(st, 6)), mu = as.list(g_sig(fobj$pop$mu, 6)),
                      omega = as.list(g_sig(fobj$pop$omega, 6)),
                      beta = lapply(fobj$pop$beta, function(b) as.list(g_sig(b, 6))),
                      sigma = as.list(g_sig(fobj$pop$sigma, 4)), tau = as.list(g_sig(fobj$pop$tau, 4)),
                      base_s = as.list(g_sig(fobj$pop$base_s, 4)), re = fobj$re,
                      clone_types = cfg$model$clone_types),
    asof_days = g_asof_days, grid_days = g_grid_days,
    plans = lapply(names(g_plans), function(k) c(list(id = k), g_plans[[k]])), utility = as.list(g_utility),
    phenotype = list(rule = cfg$phenotype$primary, term = cfg$fit$pheno[c("w", "s", "eps")]),
    params = g_param_table(cfg, fobj), k1 = g_k1(cfg), results = g_results_status(cfg),
    logistic = g_logistic_module(d, cfg), real_cohort = list(n = nrow(d$patients), sealed = TRUE),
    roster = roster)
  g_write_js(idx, file.path(gl_dir, "data", "model3_index.js"), "window.GLM3_INDEX = ")
  invisible(idx)
}

# Step G (run_all.R): uses the F1 final fit when present, otherwise the prior population
an_G_export <- function(ctx) {
  cfg <- ctx$config
  fit <- if (!is.null(ctx$results$F1)) ctx$results$F1$fits[[ctx$results$F1$selected$final]] else g_find_final_fit(cfg)
  tab <- g_export(cfg, fit = fit, log_fun = function(...) ctx$log(...))
  list(tables = list(export = tab))
}
