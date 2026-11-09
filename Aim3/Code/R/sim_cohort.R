# =============================================================================
# sim_cohort.R
# Purpose : Virtual cohorts at the real sampling design (W3). Each virtual
#           patient inherits one real patient's sampling days, censoring and
#           G-CSF pattern, transfusion exclusions, ferritin driver, eGFR,
#           bridging and set points; only the counts are simulated. Mechanism
#           scenarios: a only, b only, mixture, synergy.
# Inputs  : converted data (io_read.R), design population (config$model$design).
# Outputs : list(obs, truth) in the build_obs format; obs replaces the real counts.
# =============================================================================

# Injury type per virtual patient under a scenario.
scenario_types <- function(scenario, n, config) {
  st <- config$steps$W3
  switch(scenario,
         a_only = rep("a", n), b_only = rep("b", n), synergy = rep("ab", n),
         mix = sample(c("a", "b", "ab"), n, replace = TRUE, prob = st$mix),
         stop("unknown scenario: ", scenario))
}

# Natural-scale individual parameters for given types and a shift of the
# injury means (used to calibrate the phenotype mix).
draw_indiv <- function(types, shift, config) {
  ds <- config$model$design
  n <- length(types)
  e <- sapply(re_names_all, function(k) stats::rnorm(n, ds$mu[[k]], ds$omega[[k]]))
  if (n == 1) e <- matrix(e, nrow = 1, dimnames = list(NULL, re_names_all))
  e[, "th_as"] <- e[, "th_as"] + shift
  e[, "th_bs"] <- e[, "th_bs"] + shift
  x <- exp(e)
  x[types == "b", "th_as"] <- 0
  x[types == "a", "th_bs"] <- 0
  x
}

scenario_struct <- function(scenario, struct, config) {
  if (scenario == "synergy") struct[["eta_syn"]] <- config$steps$W3$eta_syn
  struct
}

# Shift of the injury means so that the share of non-quick phenotypes matches
# the target (bisection over virtual patients built on the real designs).
calibrate_shift <- function(pdes_list, scenario, config, n = config$steps$W3$calib_n) {
  st <- scenario_struct(scenario, config$model$struct, config)
  target <- 1 - config$model$design$target[["quick"]]
  frac <- function(shift) {
    set.seed(config$runtime$seed)
    types <- scenario_types(scenario, n, config)
    x <- draw_indiv(types, shift, config)
    ids <- sample(seq_along(pdes_list), n, replace = TRUE)
    ph <- vapply(seq_len(n), function(j)
      sim_phenotype(build_par(x[j, ], st, "Mabx"), pdes_list[[ids[j]]], config), character(1))
    mean(ph != "quick", na.rm = TRUE)
  }
  lo <- -3
  hi <- 4
  for (i in 1:9) {
    mid <- (lo + hi) / 2
    if (frac(mid) < target) lo <- mid else hi <- mid
  }
  (lo + hi) / 2
}

# One virtual cohort on the real designs.
simulate_cohort <- function(d, pdes_list, scenario, shift, config, seed) {
  set.seed(seed)
  st <- scenario_struct(scenario, config$model$struct, config)
  ds <- config$model$design
  ids <- d$patients$patient_id
  types <- scenario_types(scenario, length(ids), config)
  x <- draw_indiv(types, shift, config)
  obs_out <- list()
  truth <- list()
  for (j in seq_along(ids)) {
    pid <- ids[j]
    pd <- pdes_list[[pid]]
    o <- data.table::copy(d$obs[patient_id == pid & keep == TRUE])
    par <- build_par(x[j, ], st, "Mabx")
    m <- sim_raw(par, pd, sort(unique(o$day)), config)
    mi <- match(o$day, m[, "t"])
    pred <- m[cbind(mi, match(rel_col[o$var], colnames(m)))] * pd$set[o$var]
    sig <- ds$sigma[o$var]
    eps <- stats::rnorm(nrow(o), 0, sig)
    val <- pred * exp(eps)
    cen <- o$reason %in% c("no_differential", "zero_on_differential")
    val[cen] <- pred[cen] * exp(abs(eps[cen])) * 3
    gc <- o$reason == "gcsf_upper_bound"
    val[gc] <- val[gc] * stats::runif(sum(gc), ds$gcsf_fold[1], ds$gcsf_fold[2])
    sa <- o$reason == "steroid_upper_bound"
    val[sa] <- val[sa] * stats::runif(sum(sa), ds$steroid_anc_fold[1], ds$steroid_anc_fold[2])
    sm <- o$reason == "steroid_lower_bound"
    val[sm] <- val[sm] * stats::runif(sum(sm), ds$steroid_amc_fold[1], ds$steroid_amc_fold[2])
    o[, value := pmax(val, 1e-4)]
    obs_out[[pid]] <- o
    kf <- kill_fractions(par, pd, config)
    truth[[pid]] <- data.table::data.table(
      patient_id = pid, type = types[j], th_ap = x[j, "th_ap"], th_as = x[j, "th_as"],
      th_bs = x[j, "th_bs"], r_s = x[j, "r_s"],
      pheno = sim_phenotype(par, pd, config), mech = true_mech_class(par, pd, config),
      f_a = kf[["f_a"]], f_b = kf[["f_b"]])
  }
  list(obs = data.table::rbindlist(obs_out), truth = data.table::rbindlist(truth),
       scenario = scenario, shift = shift, seed = seed)
}

# pdats for a virtual cohort (real designs, simulated observations).
cohort_pdats <- function(d, sim, config) {
  d2 <- d
  d2$obs <- sim$obs
  build_pdats(d2, config$model$struct, config)
}

# Which model is "correct" for a scenario (strict) and acceptable alternatives.
scenario_truth_model <- function(scenario) {
  switch(scenario, a_only = "Ma", b_only = "Mb", mix = "Mab", synergy = "Mabx")
}
