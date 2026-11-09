# =============================================================================
# an_D6_virtual_cohort.R
# Purpose : D6 - virtual cohort from the fitted M2 population (training split).
#           (1) How much persistent cytopenia is attributable to clones?
#               Counterfactual: the same virtual patients with f0 = 0.
#           Each virtual patient inherits the CRS treatment (tocilizumab,
#           steroids) of the training patient whose covariates it reuses;
#           response-driven support (G-CSF, transfusions) is not simulated.
#           (2) Paired sampling design: power to detect clonal expansion
#               (pre-treatment vs day +30 VAF) in persistent patients, and the
#               sequencing cost of a stratified day +30 sampling rule.
# Inputs  : ctx$data() (train covariates/baselines), ctx$results$D4.
# Outputs : list(cohort, tables)
# =============================================================================

d6_draw_cohort <- function(fit, tr, n, cfg) {
  rows <- tr$patients[sample(.N, n, replace = TRUE)]
  pop <- fit$pop
  eta <- sapply(fit$re, function(k) {
    b <- pop$beta[[k]]
    mean_i <- pop$mu[[k]] + if (length(b)) as.matrix(rows[, names(b), with = FALSE]) %*% b else 0
    as.numeric(mean_i) + pop$omega[[k]] * stats::rnorm(n)
  })
  data.table::data.table(vid = seq_len(n), source_id = rows$patient_id, anc_base = rows$anc_base,
                         ld_scale = rows$ld_scale,
                         eta_rho = eta[, "rho"], eta_g = eta[, "g"], eta_f0 = eta[, "f0"])
}

d6_simulate <- function(coh, shift, fit, cfg, events) {
  st <- fit$pop$struct
  crs_ev <- events[event %in% c("tocilizumab", "steroid")]
  ev_by_id <- split(crs_ev, by = "patient_id")
  post_day <- cfg$steps$D6$post_vaf_day
  out <- par_lapply(seq_len(nrow(coh)), function(i) {
    r <- coh[i]
    indiv <- c(rho = exp(r$eta_rho), g = exp(r$eta_g), f0 = stats::plogis(r$eta_f0 + shift))
    ev <- ev_by_id[[r$source_id]]
    a <- simulate_anc_nat(indiv, st, "M2", r$anc_base, r$ld_scale, cfg, ev)
    cf <- indiv
    cf["f0"] <- 0
    b <- simulate_anc_nat(cf, st, "M2", r$anc_base, r$ld_scale, cfg, ev)
    k <- which(a$day == post_day)
    c(ph = classify_phenotype(a$day, a$anc, cfg$phenotype),
      ph_cf = classify_phenotype(b$day, b$anc, cfg$phenotype),
      f0 = indiv[["f0"]], f_post = a$f_clone[k], anc_post = a$anc[k])
  }, cfg$runtime$cores)
  out <- do.call(rbind, out)
  data.table::data.table(vid = coh$vid, phenotype = out[, "ph"], phenotype_cf = out[, "ph_cf"],
                         f0 = as.numeric(out[, "f0"]), f_post = as.numeric(out[, "f_post"]),
                         anc_post = as.numeric(out[, "anc_post"]))
}

d6_power <- function(sim, cfg, scenario) {
  s <- cfg$steps$D6
  lod <- cfg$data$lod$vaf
  nz <- cfg$model$vaf_noise_sd
  n <- nrow(sim)
  pre <- sim$f0 / 2 * exp(stats::rnorm(n, 0, nz))
  post <- pmax(sim$f_post / 2 * exp(stats::rnorm(n, 0, nz)), lod / 2)
  sampled <- sim$anc_post < 1.0 | stats::runif(n) < s$quick_sample_frac
  elig <- sim$phenotype %in% "persistent" & pre >= lod
  diffs <- log(post) - log(pmax(pre, lod / 2))
  eff <- data.table::data.table(scenario = scenario, n_eligible_in_cohort = sum(elig),
                                mean_log_change = round(mean(diffs[elig]), 3),
                                sd_log_change = round(stats::sd(diffs[elig]), 3))
  eff[, cohen_d := round(mean_log_change / sd_log_change, 2)]
  pw <- data.table::rbindlist(lapply(s$n_enrol, function(N) {
    res <- replicate(s$n_rep, {
      i <- sample(n, N, replace = TRUE)
      e <- i[elig[i] & sampled[i]]
      ok <- length(e) >= 3 && stats::sd(diffs[e]) > 0
      c(sig = if (ok) { tt <- stats::t.test(diffs[e]); tt$p.value < 0.05 && mean(diffs[e]) > 0 } else FALSE,
        n_pairs = length(e), assays_strat = N + sum(sampled[i]))
    })
    data.table::data.table(scenario = scenario, n_enrolled = N, power = round(mean(res["sig", ]), 3),
                           median_pairs = stats::median(res["n_pairs", ]),
                           assays_all = 2 * N, assays_stratified = round(mean(res["assays_strat", ])))
  }))
  pw[, assay_saving_pct := round(100 * (1 - assays_stratified / assays_all), 1)]
  list(effect = eff, power = pw)
}

an_D6_virtual_cohort <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$D6
  fit <- ctx$results$D4$fits$M2
  tr <- analysis_set(ctx$data(), "train")
  set.seed(cfg$runtime$seed + 6L)
  coh <- d6_draw_cohort(fit, tr, s$n_virtual, cfg)
  thr <- stats::qlogis(2 * cfg$data$lod$vaf)     # clone detectable at VAF >= LOD
  p_fit <- mean(coh$eta_f0 >= thr)
  scen <- c(as_fitted = 0, stats::setNames(
    vapply(s$ch_prevalence, function(p) thr - stats::quantile(coh$eta_f0, 1 - p, names = FALSE), numeric(1)),
    paste0("prevalence_", round(100 * s$ch_prevalence))))

  sims <- list()
  af <- list()
  power <- list()
  effects <- list()
  for (nm in names(scen)) {
    ctx$log("  D6: scenario ", nm)
    x <- d6_simulate(coh, scen[[nm]], fit, cfg, tr$events)
    sims[[nm]] <- x
    p_pers <- mean(x$phenotype %in% "persistent")
    p_pers_cf <- mean(x$phenotype_cf %in% "persistent")
    af[[nm]] <- data.table::data.table(
      scenario = nm, ch_prevalence = round(mean(x$f0 >= 2 * cfg$data$lod$vaf), 3),
      p_quick = round(mean(x$phenotype %in% "quick"), 3),
      p_second_wave = round(mean(x$phenotype %in% "second_wave"), 3),
      p_persistent = round(p_pers, 3), p_persistent_without_clones = round(p_pers_cf, 3),
      attributable_fraction = round(if (p_pers > 0) (p_pers - p_pers_cf) / p_pers else NA_real_, 3))
    pw <- d6_power(x, cfg, nm)
    power[[nm]] <- pw$power
    effects[[nm]] <- pw$effect
  }
  af <- data.table::rbindlist(af)
  power <- data.table::rbindlist(power)
  effects <- data.table::rbindlist(effects)
  list(cohort = coh, sims = sims, fitted_prevalence = p_fit,
       tables = list(attributable_fraction = af, power = power, effect_size = effects),
       key_numbers = list(af = af, power = power, effects = effects))
}
