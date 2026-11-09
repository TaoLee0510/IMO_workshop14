# =============================================================================
# an_D3_identifiability.R
# Purpose : D3 - can the data tell the mechanisms apart? Virtual patients are
#           drawn from the design population, observed with the real sampling
#           schedule and noise, and refitted under three information
#           scenarios:
#             counts              blood counts + CRP/ferritin only
#             counts_pre_vaf      + pre-treatment VAF (what we have)
#             counts_pre_post_vaf + one VAF at day +30 (what Aim 2 asks for)
#           Individual fits use the true design population as prior, which
#           isolates identifiability from population estimation. The design
#           population is first calibrated to the target phenotype mix under
#           the clinical process (calibrate_design_pop). Virtual
#           patients receive the same supportive care as the test data
#           (CRS treatment, G-CSF; sim_clinical.R) and the same sampling
#           (sparse follow-up after day 30 for some, random missed draws);
#           persistently low patients may get TPO-RA or a CD34+ boost, and
#           relapse, death or a boost censors later data. The fit knows the
#           interventions, as it would with real data.
# Inputs  : ctx$config (design values; no patient data).
# Outputs : list(patients, tables, key_numbers)
# =============================================================================

# Sampling schedule of one virtual patient. clin = NULL gives the full schedule.
d3_schedule <- function(clin = NULL) {
  labs <- c(-5:14, seq(16, 90, by = 2))
  inflam <- c(-5:14, 17, 21, 24, 28)
  if (!is.null(clin)) {
    sm <- clin$sampling
    if (stats::runif(1) < sm$sparse_frac)
      labs <- c(labs[labs < sm$sparse_after], seq(sm$sparse_after, 90, by = 7))
    labs <- labs[stats::runif(length(labs)) >= sm$missing | labs == -5]
    inflam <- inflam[stats::runif(length(inflam)) >= sm$inflam_missing | inflam == -5]
  }
  list(labs = labs, inflam = inflam)
}

# Simulate one virtual patient: clinical process, then noisy observations.
# CRP after tocilizumab is dropped, as io_read.R does for real data.
d3_virtual_patient <- function(id, indiv, struct, cfg, sched) {
  b <- cfg$model$design$base
  b$anc <- baseline_with_clone(b$anc, indiv[["f0"]], struct)
  sig <- cfg$fit$sigma_init
  clin <- cfg$model$design$clinical
  s <- simulate_truth_patient(c(indiv[c("rho", "g", "f0")], anc_base = b$anc), struct, cfg, clin, ld_scale = 1)
  days <- cfg$model$time$sim_days
  m <- s$m
  ev <- s$events
  if (identical(classify_phenotype(days, m[, "anc_rel"] * b$anc, cfg$phenotype), "persistent")) {
    ev <- add_persistent_support(ev, clin$events)
    if (any(ev$event == "tpo_ra"))
      m <- marrow_sim(build_par(indiv, struct, "M2"),
                      make_segments(cfg$model$time$t0, max(days), ld_scale = 1, events = ev), days)
  }
  # censoring: CD34+ boost, relapse, death (the clinical course is truncated, not the truth)
  cens <- min(c(ev[event == "cd34_boost", start_day],
                if (stats::runif(1) < clin$events$p_relapse) sample(30:90, 1),
                if (stats::runif(1) < clin$events$p_death) sample(20:90, 1), Inf))
  sched$labs <- sched$labs[sched$labs <= cens]
  sched$inflam <- sched$inflam[sched$inflam <= cens]
  seg <- make_segments(cfg$model$time$t0, cfg$model$time$tmax, ld_scale = 1, events = ev)
  toci <- ev[event == "tocilizumab", start_day]
  if (length(toci)) sched$inflam_crp <- sched$inflam[sched$inflam <= min(toci)]
  else sched$inflam_crp <- sched$inflam
  sub <- function(dd) m[match(dd, days), , drop = FALSE]
  ml <- sub(sched$labs)
  mi <- sub(sched$inflam)
  mc <- sub(sched$inflam_crp)
  anc <- ml[, "anc_rel"] * b$anc * exp(stats::rnorm(nrow(ml), 0, sig[["anc"]]))
  obs <- rbind(
    data.table::data.table(day = sched$labs, var = "anc", value = pmax(anc, cfg$data$lod$anc),
                           cens = anc <= cfg$data$lod$anc),
    data.table::data.table(day = sched$labs, var = "plt", cens = FALSE,
                           value = ml[, "plt_rel"] * b$plt * exp(stats::rnorm(nrow(ml), 0, sig[["plt"]]))),
    data.table::data.table(day = sched$inflam_crp, var = "crp", cens = FALSE,
                           value = (b$crp + cfg$model$obs$crp_scale * mc[, "is"]) *
                             exp(stats::rnorm(nrow(mc), 0, sig[["crp"]]))),
    data.table::data.table(day = sched$inflam, var = "ferritin", cens = FALSE,
                           value = (b$ferritin + cfg$model$obs$fer_scale * mi[, "is"]) *
                             exp(stats::rnorm(nrow(mi), 0, sig[["ferritin"]]))))
  lod <- cfg$data$lod$vaf
  nz <- cfg$model$vaf_noise_sd
  f_pre <- indiv[["f0"]] * exp(stats::rnorm(1, 0, nz))
  f_post_true <- m[match(cfg$s3_post_day, days), "f_clone"]
  f_post <- f_post_true * exp(stats::rnorm(1, 0, nz))
  list(id = id, obs = obs, seg = seg, base = b, events = ev, censor_day = cens,
       vaf_pre = if (f_pre / 2 >= lod) list(status = "positive", f_measured = min(f_pre, 1),
                                            sd = cfg$model$vaf_prior_sd$pre_LD, lod = lod)
                 else list(status = "negative", f_measured = NA_real_,
                           sd = cfg$model$vaf_prior_sd$negative, lod = lod),
       vaf_post = if (f_post / 2 >= lod && cens >= cfg$s3_post_day) list(day = cfg$s3_post_day, f_measured = min(f_post, 1),
                                              sd = cfg$model$vaf_post_sd) else NULL)
}

# Design population calibrated so that, under the design clone traits and the
# clinical process, the phenotype mix matches config$model$design$target. Same
# update rule as the TestData generator (mu_g moves the share with a second
# dip or worse, mu_rho the persistent share). Reseeds the RNG.
calibrate_design_pop <- function(cfg, st) {
  des <- cfg$model$design
  cal <- des$calib
  tgt <- des$target
  pop <- des$pop
  days <- cfg$model$time$sim_days
  set.seed(cfg$runtime$seed + 30L)
  n <- cal$n
  z <- matrix(stats::rnorm(3 * n), n, 3)
  ch <- stats::runif(n) < pop$p_ch
  trace <- list()
  for (it in seq_len(cal$max_iter)) {
    rho <- exp(pop$mu[["rho"]] + pop$omega[["rho"]] * z[, 1])
    g <- exp(pop$mu[["g"]] + pop$omega[["g"]] * z[, 2])
    f0 <- ifelse(ch, stats::plogis(pop$mu[["f0"]] + pop$omega[["f0"]] * z[, 3]), pop$f0_non_ch)
    ph <- vapply(seq_len(n), function(i) {
      a0 <- baseline_with_clone(des$base$anc, f0[i], st)
      s <- simulate_truth_patient(c(rho = rho[i], g = g[i], f0 = f0[i], anc_base = a0), st, cfg,
                                  des$clinical, ld_scale = 1)
      classify_phenotype(days, s$m[, "anc_rel"] * a0, cfg$phenotype)
    }, character(1))
    pr <- prop.table(table(factor(ph, phenotype_levels)))
    trace[[it]] <- data.table::data.table(iter = it, mu_rho = pop$mu[["rho"]], mu_g = pop$mu[["g"]],
                                          quick = pr[["quick"]], second_wave = pr[["second_wave"]],
                                          persistent = pr[["persistent"]])
    err_hit <- (pr[["second_wave"]] + pr[["persistent"]]) - (tgt[["second_wave"]] + tgt[["persistent"]])
    err_pers <- pr[["persistent"]] - tgt[["persistent"]]
    if ((abs(err_hit) < cal$tol && abs(err_pers) < cal$tol) || it == cal$max_iter) break
    pop$mu[["g"]] <- pop$mu[["g"]] - 1.5 * err_hit
    pop$mu[["rho"]] <- pop$mu[["rho"]] + 2.5 * err_pers
  }
  list(pop = pop, trace = data.table::rbindlist(trace))
}

# Virtual cohort drawn from a design population (parameters only).
d3_truth_cohort <- function(n, cfg, st, pop, prefix = "V") {
  ch <- stats::runif(n) < pop$p_ch
  data.table::data.table(
    id = sprintf("%s%04d", prefix, seq_len(n)),
    rho = exp(pop$mu[["rho"]] + pop$omega[["rho"]] * stats::rnorm(n)),
    g = exp(pop$mu[["g"]] + pop$omega[["g"]] * stats::rnorm(n)),
    f0 = ifelse(ch, stats::plogis(pop$mu[["f0"]] + pop$omega[["f0"]] * stats::rnorm(n)), pop$f0_non_ch),
    ch = ch)
}

# True phenotype and mechanism under the true clone traits and each patient's
# actual supportive care (events of the simulated clinical course).
d3_add_truth_mechanism <- function(truth, vp, st, cfg) {
  tm <- lapply(seq_len(nrow(truth)), function(i)
    infer_mechanism(c(rho = truth$rho[i], g = truth$g[i], f0 = truth$f0[i]), st, "M2",
                    vp[[i]]$base$anc, 1, cfg, vp[[i]]$events))
  truth[, `:=`(phenotype_true = vapply(tm, `[[`, character(1), "phenotype"),
               mechanism_true = vapply(tm, `[[`, character(1), "mechanism"))]
  truth
}

# Prior used for the individual fits in D3/D3b: the design population, with f0
# centred on the carrier / non-carrier mixture mean.
d3_design_pop <- function(cfg, st, dp) {
  mu_f0 <- stats::qlogis(dp$p_ch * stats::plogis(dp$mu[["f0"]]) + (1 - dp$p_ch) * dp$f0_non_ch)
  list(mu = c(rho = dp$mu[["rho"]], g = dp$mu[["g"]], f0 = mu_f0),
       omega = c(rho = dp$omega[["rho"]], g = dp$omega[["g"]], f0 = 2.0),
       beta = list(rho = numeric(0), g = numeric(0), f0 = numeric(0)),
       sigma = cfg$fit$sigma_init, struct = st)
}

# Observations of one virtual patient as a pdat; vaf / vaf_post as supplied.
d3_pdat <- function(v, cfg, vaf = NULL, vaf_post = NULL) {
  build_pdat(v$id, v$obs, v$base, v$seg, list(rho = numeric(0), g = numeric(0), f0 = numeric(0)),
             vaf = vaf, vaf_post = vaf_post, extra_times = cfg$model$time$tmax)
}

an_D3_identifiability <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$D3
  cfg$s3_post_day <- s$post_vaf_day
  st <- design_struct(cfg)
  n <- s$n_virtual
  cal <- calibrate_design_pop(cfg, st)
  set.seed(cfg$runtime$seed + 3L)
  truth <- d3_truth_cohort(n, cfg, st, cal$pop)

  vp <- lapply(seq_len(n), function(i)
    d3_virtual_patient(truth$id[i], c(rho = truth$rho[i], g = truth$g[i], f0 = truth$f0[i]), st, cfg,
                       d3_schedule(cfg$model$design$clinical)))
  truth <- d3_add_truth_mechanism(truth, vp, st, cfg)

  pop <- d3_design_pop(cfg, st, cal$pop)

  fit_one <- function(v, model, scenario) {
    pd <- d3_pdat(v, cfg, vaf = if (scenario != "counts") v$vaf_pre else NULL,
                  vaf_post = if (scenario == "counts_pre_post_vaf") v$vaf_post else NULL)
    re <- model_re_names(model)
    prior <- patient_prior(pd, pop, re)
    eta0 <- stats::setNames(prior$mean, re)
    f <- fit_individual(pd, par_template(st, model), model, prior, pop$sigma, cfg, eta0,
                        use_vaf = scenario != "counts",
                        extra_starts = clone_starts(eta0, model, if (scenario != "counts") pd else NULL))
    indiv <- eta_named_to_indiv(f$eta)
    mech <- infer_mechanism(indiv, st, model, v$base$anc, 1, cfg, v$events)
    data.table::data.table(id = v$id, model = model, scenario = scenario, logL = f$logL,
                           f0_est = if (model == "M2") indiv[["f0"]] else NA_real_,
                           f0_sd_logit = if (model == "M2") sqrt(f$cov["f0", "f0"]) else NA_real_,
                           rho_est = indiv[["rho"]], g_est = indiv[["g"]],
                           phenotype_pred = mech$phenotype, mechanism_pred = mech$mechanism)
  }
  jobs <- rbind(data.table::CJ(i = seq_len(n), model = "M2", scenario = s$scenarios),
                data.table::CJ(i = seq_len(n), model = c("M0", "M1"), scenario = "counts"))
  res <- par_lapply(seq_len(nrow(jobs)), function(j)
    fit_one(vp[[jobs$i[j]]], jobs$model[j], jobs$scenario[j]), cfg$runtime$cores)
  fits <- merge(data.table::rbindlist(res), truth, by = "id")

  m2 <- fits[model == "M2"]
  acc <- m2[!is.na(mechanism_true), .(n = .N, mechanism_accuracy = round(mean(mechanism_pred == mechanism_true, na.rm = TRUE), 3)),
            by = scenario]
  acc_pers <- m2[mechanism_true %in% c("H1", "H3"),
                 .(n = .N, h1_vs_h3_accuracy = round(mean(mechanism_pred == mechanism_true, na.rm = TRUE), 3),
                   h3_sensitivity = round(mean(mechanism_pred[mechanism_true == "H3"] == "H3"), 3),
                   h1_specificity = round(mean(mechanism_pred[mechanism_true == "H1"] == "H1"), 3)),
                 by = scenario]
  f0_rec <- m2[ch == TRUE, .(n_ch = .N, cor_logit_f0 = round(stats::cor(stats::qlogis(f0_est), stats::qlogis(f0)), 3),
                             median_post_sd_logit = round(stats::median(f0_sd_logit), 3)), by = scenario]
  confusion <- m2[!is.na(mechanism_true), .N, by = .(scenario, mechanism_true, mechanism_pred)]
  sel <- fits[scenario == "counts", .SD[which.max(logL)], by = id][, .(id, selected = model)]
  sel <- merge(sel, truth, by = "id")
  sel[, needed := data.table::fcase(mechanism_true == "H3", "M2", mechanism_true %in% c("H1", "H2"), "M1",
                                    default = "M0")]
  model_sel <- sel[!is.na(mechanism_true), .N, by = .(mechanism_true, selected)]

  list(patients = fits, truth = truth,
       tables = list(design_calibration = cal$trace, mechanism_accuracy = acc, persistent_h1_vs_h3 = acc_pers, f0_recovery = f0_rec,
                     confusion = confusion, model_selection_counts = model_sel,
                     truth_phenotypes = truth[, .N, by = .(phenotype_true, mechanism_true)]),
       key_numbers = list(acc = acc, acc_pers = acc_pers, f0_rec = f0_rec))
}
