# =============================================================================
# sim_cohort.R
# Purpose : Simulate one synthetic cohort with known ground truth, in memory.
#           Patients are simulated from the marrow twin itself (model M2),
#           then passed through realistic clinical processes: CRS treatment,
#           G-CSF, transfusions, CD34 boosts, relapse, death, sparse
#           follow-up, missingness, measurement noise and detection limits.
#           Used by testdata/make_testdata.R (which writes TestData/) and by
#           step D3c (repeated cohorts at the TestData scale).
# Inputs  : config; scenario; seed; optionally a pre-computed calibration.
# Outputs : list(raw, post_muts, truth, spec, cal) - `raw` has the same
#           tables as read_raw() returns for files on disk.
#
# Scenarios
#   main   : data-generating model = analysis model (checks code correctness)
#   stress : deliberately misspecified - two clone types (DDR vs DTA) with
#            different traits, sparser follow-up, more missingness and noise,
#            and half of the marrow-involvement cases left unflagged
# =============================================================================

# ---- truth specification ------------------------------------------------------
truth_spec <- function(scenario = "main") {
  spec <- list(
    scenario = scenario,
    n = 350L, n_sequenced = 15L, n_post = 3L, train_frac = 0.7,
    # structural truth: config defaults, overridden here
    struct_override = c(r_ld = 0.3, r_inf = 0.2, delta = 0.03, q = 0.12),
    clone_types = NULL,
    pop = list(mu_rho = 0, sd_rho = 0.7, mu_g = -1.2, sd_g = 0.6,
               mu_f = -0.9, sd_f = 0.9,
               beta_rho = c(prior_lines = -0.10, prior_asct = -0.25,
                            prior_bendamustine = -0.20, age = -0.008),
               beta_g = c(costim_cd28 = 0.35),
               ch_logit = c(int = -0.3, age = 0.04, lines = 0.25, asct = 0.6),
               beta_f_lines = 0.15,
               f0_non_ch = 0.001),
    target = c(quick = 0.50, second_wave = 0.40, persistent = 0.10),
    calib = list(n = 2000L, max_iter = 8L, tol = 0.02),
    noise = c(anc = 0.20, plt = 0.15, hgb = 0.04, crp = 0.25, ferritin = 0.25, vaf = 0.15),
    sampling = list(sparse_frac = 0.30, sparse_after = 30, missing = 0.08,
                    inflam_days = c(-5:14, 17, 21, 24, 28), inflam_missing = 0.05),
    events = list(p_relapse = 0.15, p_death = 0.05, p_marrow_inv = 0.05,
                  p_iec_hs = 0.02, p_viral = 0.05, p_gcsf = 0.6, p_gcsf_second = 0.5,
                  p_toci_g2 = 0.8, p_toci_g1 = 0.2, p_steroid_g3 = 0.7,
                  p_steroid_g2 = 0.15, p_anakinra_g3 = 0.2,
                  p_cd34_persistent = 0.3, p_tpo_persistent = 0.3,
                  marrow_inv_unflagged = 0),
    crs_cut = c(0.3, 1.2, 2.4, 3.6),      # Is peak thresholds for CRS grade 1..4
    # red-cell indices (assumed effects): log RDW and MCV per unit clone
    # fraction, per prior line, per log marrow reserve; between-patient and
    # day-to-day noise; relative RDW rise from day 7
    red_cell = list(rdw = c(base = 14.0, f0 = 0.15, lines = 0.02, log_rho = -0.04, sd = 0.08,
                            day_sd = 0.03, post = 0.06),
                    mcv = c(base = 91, f0 = 6, lines = 1.0, sd = 5, day_sd = 1.2))
  )
  if (scenario == "stress") {
    spec$clone_types <- list(
      DDR = c(r_ld = 0.3, r_inf = 0.8, delta = 0.03),
      DTA = c(r_ld = 0.9, r_inf = 0.2, delta = 0.9))
    spec$noise <- spec$noise + c(anc = 0.05, plt = 0.05, hgb = 0, crp = 0.05, ferritin = 0.05, vaf = 0.05)
    spec$sampling$sparse_frac <- 0.5
    spec$sampling$sparse_after <- 21
    spec$sampling$missing <- 0.15
    spec$events$marrow_inv_unflagged <- 0.5
  }
  spec
}

ch_hotspots <- list(DNMT3A = "p.R882H", TET2 = "p.Q1034*", ASXL1 = "p.G646Wfs*12",
                    PPM1D = "p.R458*", TP53 = "p.R248Q")

# ---- covariates and individual parameters ----------------------------------------
draw_covariates <- function(n) {
  data.table::data.table(
    age = round(pmin(pmax(stats::rnorm(n, 62, 11), 25), 85)),
    sex = ifelse(stats::rbinom(n, 1, 0.38) == 1, "F", "M"),
    product = sample(c("A", "B", "C"), n, TRUE, prob = c(0.55, 0.35, 0.10)),
    prior_lines = pmin(2L + stats::rpois(n, 0.9), 7L),
    prior_asct = stats::rbinom(n, 1, 0.25),
    prior_bendamustine = stats::rbinom(n, 1, 0.20),
    months_since_last_chemo = round(stats::rlnorm(n, log(3), 0.6), 1),
    bridging = sample(c("none", "chemo", "radiation"), n, TRUE, prob = c(0.4, 0.4, 0.2)),
    ld_dose_reduced = stats::rbinom(n, 1, 0.10)
  )[, `:=`(costim = ifelse(product == "A", "CD28", "4-1BB"),
           ld_scale = ifelse(ld_dose_reduced == 1, 0.75, 1.0))]
}

draw_individual <- function(cov, spec, mu_rho, mu_g, z) {
  p <- spec$pop
  lin_rho <- mu_rho + p$beta_rho["prior_lines"] * (cov$prior_lines - 2) +
    p$beta_rho["prior_asct"] * cov$prior_asct +
    p$beta_rho["prior_bendamustine"] * cov$prior_bendamustine +
    p$beta_rho["age"] * (cov$age - 62)
  lin_ch <- p$ch_logit["int"] + p$ch_logit["age"] * (cov$age - 62) +
    p$ch_logit["lines"] * (cov$prior_lines - 2) + p$ch_logit["asct"] * cov$prior_asct
  ch <- as.integer(z$u_ch < stats::plogis(lin_ch))
  f0 <- ifelse(ch == 1, stats::plogis(p$mu_f + p$beta_f_lines * (cov$prior_lines - 2) + p$sd_f * z$f),
               p$f0_non_ch)
  clone_class <- ifelse(ch == 1, ifelse(z$u_type < 0.5, "DDR", "DTA"), "none")
  data.table::data.table(
    rho = unname(exp(lin_rho + p$sd_rho * z$rho)),
    g = unname(exp(mu_g + p$beta_g["costim_cd28"] * (cov$costim == "CD28") + p$sd_g * z$g)),
    f0 = unname(f0), ch = ch, clone_class = clone_class)
}

draw_z <- function(n) {
  list(rho = stats::rnorm(n), g = stats::rnorm(n), f = stats::rnorm(n),
       u_ch = stats::runif(n), u_type = stats::runif(n))
}

# Pre-treatment ANC: low reserve and large low-output clones lower the baseline.
draw_anc_base <- function(rho, f0, struct) {
  exp(stats::rnorm(length(rho), log(3.2) + 0.4 * log(rho), 0.30)) *
    ((1 - f0) + struct[["delta"]] * f0)^(1 / (1 + struct[["gamma_n"]]))
}

# Structural truth for one patient (clone traits may depend on clone class).
patient_struct <- function(struct, spec, clone_class) {
  if (!is.null(spec$clone_types) && clone_class %in% names(spec$clone_types)) {
    tr <- spec$clone_types[[clone_class]]
    struct[names(tr)] <- tr
  }
  struct
}

# Clinical process (CRS treatment, G-CSF): crs_from_is() and
# simulate_truth_patient() live in R/sim_clinical.R (shared with D3/D3b).

phenotype_truth <- function(m, anc_base, rule, days) {
  classify_phenotype(days, m[, "anc_rel"] * anc_base, rule)
}

# ---- population calibration -------------------------------------------------------
calibrate_population <- function(cfg, spec, struct) {
  n <- spec$calib$n
  cov <- draw_covariates(n)
  if (isTRUE(spec$all_cd28)) cov[, `:=`(product = "A", costim = "CD28")]
  z <- draw_z(n)
  days <- cfg$model$time$sim_days
  mu_rho <- spec$pop$mu_rho
  mu_g <- spec$pop$mu_g
  tgt <- spec$target
  for (it in seq_len(spec$calib$max_iter)) {
    ind <- draw_individual(cov, spec, mu_rho, mu_g, z)
    anc0 <- draw_anc_base(ind$rho, ind$f0, struct)
    # calibrate on trajectories that include the clinical process (CRS treatment, G-CSF)
    ph <- vapply(seq_len(n), function(i) {
      st <- patient_struct(struct, spec, ind$clone_class[i])
      s <- simulate_truth_patient(c(rho = ind$rho[i], g = ind$g[i], f0 = ind$f0[i], anc_base = anc0[i]),
                                  st, cfg, spec, cov$ld_scale[i])
      classify_phenotype(days, s$m[, "anc_rel"] * anc0[i], cfg$phenotype)
    }, character(1))
    pr <- prop.table(table(factor(ph, phenotype_levels)))
    msg <- sprintf("  calibration %d: mu_rho %.3f mu_g %.3f -> quick %.3f wave %.3f persistent %.3f",
                   it, mu_rho, mu_g, pr[["quick"]], pr[["second_wave"]], pr[["persistent"]])
    message(msg)
    err_hit <- (pr[["second_wave"]] + pr[["persistent"]]) - (tgt[["second_wave"]] + tgt[["persistent"]])
    err_pers <- pr[["persistent"]] - tgt[["persistent"]]
    if (abs(err_hit) < spec$calib$tol && abs(err_pers) < spec$calib$tol) break
    mu_g <- mu_g - 1.5 * err_hit
    mu_rho <- mu_rho + 2.5 * err_pers
  }
  list(mu_rho = mu_rho, mu_g = mu_g, proportions = as.numeric(pr))
}

# ---- observation process -------------------------------------------------------------
# spec$sampling$lab_days (optional) replaces the default schedule; a late_frac
# share of patients then get one more draw on a day from late_days (D3d).
lab_schedule <- function(sparse, spec, death_day) {
  sm <- spec$sampling
  if (!is.null(sm$lab_days)) {
    d <- sm$lab_days[stats::runif(length(sm$lab_days)) >= sm$missing | sm$lab_days == -5]
    if (stats::runif(1) < sm$late_frac) d <- c(d, sample(sm$late_days, 1))
  } else {
    d <- c(-5:14, seq(16, 90, by = 2))
    if (sparse) d <- c(d[d < sm$sparse_after], seq(sm$sparse_after, 90, by = 7))
    d <- d[stats::runif(length(d)) >= sm$missing | d == -5]
  }
  if (!is.na(death_day)) d <- d[d < death_day]
  sort(unique(d))
}

# Non-model perturbations applied to the observed (not ODE) trajectory.
perturb_counts <- function(day, anc, plt, pt) {
  f <- rep(1, length(day))
  if (pt$marrow_involvement == 1) f <- f * ifelse(day >= 30, pmax(0.15, 1 - (day - 30) / 40), 1)
  if (pt$iec_hs == 1) f <- f * ifelse(day >= 14 & day <= 40, 0.3, 1)
  if (pt$viral_infection == 1) f <- f * ifelse(day >= pt$viral_day & day <= pt$viral_day + 10, 0.6, 1)
  if (!is.na(pt$cd34_day)) {
    rec <- ifelse(day > pt$cd34_day, pmin(1, (day - pt$cd34_day) / 14), 0)
    anc <- anc * (1 - rec) + pt$anc_base * 0.8 * rec
    plt <- plt * (1 - rec) + pt$plt_base * 0.7 * rec
  }
  list(anc = anc * f, plt = plt * f)
}

observe_labs <- function(pid, truth_daily, pt, spec, days_obs) {
  td <- truth_daily[match(days_obs, truth_daily$day)]
  pc <- perturb_counts(days_obs, td$anc, td$plt, pt)
  nz <- spec$noise
  anc <- pc$anc * exp(stats::rnorm(length(days_obs), 0, nz[["anc"]]))
  anc[anc < 0.01] <- 0
  plt_true <- pc$plt * exp(stats::rnorm(length(days_obs), 0, nz[["plt"]]))

  # platelet transfusions: triggered by an observed count < 10, transient increment
  plt <- plt_true
  tx_days <- numeric(0)
  for (i in seq_along(days_obs)) {
    inc <- sum(25 * exp(-0.35 * (days_obs[i] - tx_days[tx_days < days_obs[i]])))
    plt[i] <- plt_true[i] + inc
    if (plt[i] < 10 && days_obs[i] >= 0) tx_days <- c(tx_days, days_obs[i])
  }
  # haemoglobin: slow follower of marrow output, RBC transfusions below 7.5
  hgb_true <- pt$hgb_base * (0.78 + 0.22 * pmin(1, td$f_rel_smooth))
  hgb <- hgb_true
  rbc_days <- numeric(0)
  for (i in seq_along(days_obs)) {
    inc <- sum(1.2 * exp(-0.03 * (days_obs[i] - rbc_days[rbc_days < days_obs[i]])))
    hgb[i] <- hgb_true[i] * exp(stats::rnorm(1, 0, nz[["hgb"]])) + inc
    if (hgb[i] < 7.5 && days_obs[i] >= 0) rbc_days <- c(rbc_days, days_obs[i])
  }
  # lymphocytes: depleted by LD, CAR T expansion bump, slow partial recovery
  alc <- ifelse(days_obs < -3, pt$alc_base,
         ifelse(days_obs < 0, pt$alc_base * 0.3, 0.1 + 0.6 * pt$g_true * exp(-((days_obs - 9) / 4)^2) +
                  0.5 * pt$alc_base * pmin(1, pmax(0, days_obs) / 90)))
  alc <- alc * exp(stats::rnorm(length(days_obs), 0, 0.25))
  wbc <- anc + alc + 0.3 * exp(stats::rnorm(length(days_obs), 0, 0.3))

  labs <- data.table::data.table(patient_id = pid, day = days_obs, anc = round(anc, 2),
                                 plt = round(plt), wbc = round(wbc, 2), hgb = round(hgb, 1),
                                 alc = round(alc, 2))
  ev <- rbind(
    if (length(tx_days)) data.table::data.table(patient_id = pid, event = "plt_tx",
                                                start_day = tx_days, end_day = tx_days, dose = 1),
    if (length(rbc_days)) data.table::data.table(patient_id = pid, event = "rbc_tx",
                                                 start_day = rbc_days, end_day = rbc_days, dose = 1))
  list(labs = labs, events = ev)
}

observe_inflam <- function(pid, truth_daily, pt, spec, toci_day, cfg) {
  d <- spec$sampling$inflam_days
  d <- d[stats::runif(length(d)) >= spec$sampling$inflam_missing | d == -5]
  if (!is.na(pt$death_day)) d <- d[d < pt$death_day]
  td <- truth_daily[match(d, truth_daily$day)]
  crp_infl <- cfg$model$obs$crp_scale * td$is
  # tocilizumab blocks hepatic CRP production for ~3 weeks: observation artefact
  if (!is.na(toci_day)) crp_infl <- crp_infl * ifelse(d >= toci_day & d <= toci_day + 21, 0.15, 1)
  crp <- (pt$crp_base + crp_infl) * exp(stats::rnorm(length(d), 0, spec$noise[["crp"]]))
  fer <- (pt$fer_base + cfg$model$obs$fer_scale * td$is) * exp(stats::rnorm(length(d), 0, spec$noise[["ferritin"]]))
  data.table::data.table(patient_id = pid, day = d, crp = round(crp, 1), ferritin = round(fer))
}

# ---- mutation data -----------------------------------------------------------------------
draw_gene <- function(clone_class) {
  if (clone_class == "DDR") sample(c("PPM1D", "TP53"), 1, prob = c(0.6, 0.4))
  else sample(c("DNMT3A", "TET2", "ASXL1"), 1, prob = c(0.55, 0.30, 0.15))
}

observe_mutations_pre <- function(seq_ids, pts, spec, lod) {
  tested <- list()
  muts <- list()
  for (pid in seq_ids) {
    pt <- pts[patient_id == pid]
    tp <- sample(c("apheresis", "post_bridging", "pre_LD"), 1, prob = c(0.5, 0.25, 0.25))
    sday <- switch(tp, apheresis = -30 + sample(-5:5, 1), post_bridging = -10 + sample(-2:2, 1), pre_LD = -6)
    tested[[pid]] <- data.table::data.table(patient_id = pid, sample_day = sday, sample_timepoint = tp,
                                            sample_type = "PB", lod_vaf = lod,
                                            panel_version = "Myeloid-45 (synthetic)")
    f_at <- pt$f0_true
    # bridging chemotherapy between apheresis and LD selects DDR-type clones
    if (tp == "apheresis" && pt$bridging == "chemo") f_at <- f_at * 0.7
    vaf <- f_at / 2 * exp(stats::rnorm(1, 0, spec$noise[["vaf"]]))
    if (pt$ch_true == 1 && vaf >= lod) {
      gene <- pt$gene_true
      muts[[length(muts) + 1]] <- data.table::data.table(
        patient_id = pid, gene = gene, protein_change = ch_hotspots[[gene]],
        vaf = round(vaf, 4), depth = sample(1500:3000, 1))
      if (stats::runif(1) < 0.3) {
        g2 <- sample(setdiff(names(ch_hotspots), gene), 1)
        muts[[length(muts) + 1]] <- data.table::data.table(
          patient_id = pid, gene = g2, protein_change = ch_hotspots[[g2]],
          vaf = round(stats::runif(1, lod, 0.03), 4), depth = sample(1500:3000, 1))
      }
    }
  }
  list(tested = data.table::rbindlist(tested),
       muts = if (length(muts)) data.table::rbindlist(muts) else
         data.table::data.table(patient_id = character(0), gene = character(0),
                                protein_change = character(0), vaf = numeric(0), depth = integer(0)))
}

# ---- red-cell indices -------------------------------------------------------------------
# Optional CBC columns. A large clone raises RDW and MCV modestly; prior lines
# and low marrow reserve raise RDW too, so both are confounded proxies of the
# clone fraction. Effect sizes are assumptions for testing (spec$red_cell).
add_red_cell_indices <- function(labs, pts, spec) {
  rc <- spec$red_cell
  rdw0 <- exp(log(rc$rdw[["base"]]) + rc$rdw[["f0"]] * pts$f0_true + rc$rdw[["lines"]] * (pts$prior_lines - 2) +
                rc$rdw[["log_rho"]] * log(pts$rho_true) + stats::rnorm(nrow(pts), 0, rc$rdw[["sd"]]))
  mcv0 <- rc$mcv[["base"]] + rc$mcv[["f0"]] * pts$f0_true + rc$mcv[["lines"]] * (pts$prior_lines - 2) +
    stats::rnorm(nrow(pts), 0, rc$mcv[["sd"]])
  i <- match(labs$patient_id, pts$patient_id)
  out <- data.table::copy(labs)
  # RDW widens after lymphodepletion while the marrow regenerates
  out[, mcv := round(mcv0[i] + stats::rnorm(.N, 0, rc$mcv[["day_sd"]]))]
  out[, rdw := round(rdw0[i] * (1 + rc$rdw[["post"]] * (day >= 7)) * exp(stats::rnorm(.N, 0, rc$rdw[["day_sd"]])), 1)]
  out[]
}

# ---- one cohort ------------------------------------------------------------------------------
# cal: list(mu_rho, mu_g) from calibrate_population(); computed here when NULL
simulate_cohort <- function(config, scenario = "main", seed = NULL, cal = NULL, spec = truth_spec(scenario)) {
  if (is.null(seed)) seed <- config$runtime$seed + if (scenario == "stress") 1L else 0L
  set.seed(seed)
  load_marrow_model(config)

  struct <- config$model$struct
  struct[names(spec$struct_override)] <- spec$struct_override
  days <- config$model$time$sim_days
  rule <- config$phenotype
  if (is.null(cal)) {
    message("Calibrating population (", scenario, ") ...")
    cal <- calibrate_population(config, spec, struct)
  }

  n <- spec$n
  ids <- sprintf("P%03d", seq_len(n))
  cov <- draw_covariates(n)
  # every patient received the same CD28 product (real-data extract)
  if (isTRUE(spec$all_cd28)) cov[, `:=`(product = "A", costim = "CD28")]
  ind <- draw_individual(cov, spec, cal$mu_rho, cal$mu_g, draw_z(n))
  ev_p <- spec$events

  pts <- data.table::data.table(patient_id = ids, cov, rho_true = ind$rho, g_true = ind$g,
                                f0_true = ind$f0, ch_true = ind$ch, clone_class = ind$clone_class)
  pts[, gene_true := vapply(clone_class, function(cc) if (cc == "none") NA_character_ else draw_gene(cc), character(1))]
  if (is.null(spec$clone_types)) pts[ch_true == 1, clone_class := ifelse(gene_true %in% c("PPM1D", "TP53"), "DDR", "DTA")]
  # baseline values: low reserve and large clones lower the pre-treatment counts
  pts[, `:=`(
    anc_base = round(draw_anc_base(rho_true, f0_true, struct), 2),
    plt_base = round(exp(stats::rnorm(.N, log(185) + 0.5 * log(rho_true), 0.30))),
    hgb_base = round(pmax(stats::rnorm(.N, 10.6, 1.3), 7.8), 1),
    crp_base = round(exp(stats::rnorm(.N, log(8), 0.7)), 1),
    fer_base = round(exp(stats::rnorm(.N, log(600), 0.6))),
    alc_base = round(exp(stats::rnorm(.N, log(0.8), 0.5)), 2),
    marrow_involvement = stats::rbinom(.N, 1, ev_p$p_marrow_inv),
    iec_hs = stats::rbinom(.N, 1, ev_p$p_iec_hs),
    viral_infection = stats::rbinom(.N, 1, ev_p$p_viral),
    sparse = stats::runif(.N) < spec$sampling$sparse_frac)]
  pts[, viral_day := ifelse(viral_infection == 1, sample(20:60, .N, TRUE), NA_integer_)]
  pts[, relapse_day := ifelse(stats::runif(.N) < ev_p$p_relapse, sample(30:90, .N, TRUE), NA_integer_)]
  pts[, subsequent_tx_day := ifelse(is.na(relapse_day), NA_integer_, relapse_day + sample(5:14, .N, TRUE))]
  pts[, death_day := ifelse(stats::runif(.N) < ev_p$p_death, sample(20:90, .N, TRUE), NA_integer_)]
  pts[, last_followup_day := ifelse(is.na(death_day), 180L, death_day)]

  # ---- simulate each patient's truth with clinical process ----
  message("Simulating ", n, " patients ...")
  truth_list <- vector("list", n)
  events_list <- vector("list", n)
  crs_list <- vector("list", n)
  for (i in seq_len(n)) {
    pt <- pts[i]
    st <- patient_struct(struct, spec, pt$clone_class)
    indiv <- c(rho = pt$rho_true, g = pt$g_true, f0 = pt$f0_true, anc_base = pt$anc_base)
    s <- simulate_truth_patient(indiv[c("rho", "g", "f0", "anc_base")], st, config, spec, pt$ld_scale)
    # TPO-RA and CD34 boost for truly persistent patients
    ph <- phenotype_truth(s$m, pt$anc_base, rule, days)
    ev <- s$events
    if (identical(ph, "persistent")) {
      ev <- add_persistent_support(ev, ev_p)
      if (any(ev$event == "tpo_ra")) {
        s$m <- marrow_sim(build_par(indiv[c("rho", "g", "f0")], st, "M2"),
                          make_segments(config$model$time$t0, max(days), pt$ld_scale, ev), days)
      }
    }
    truth_list[[i]] <- data.table::data.table(patient_id = pt$patient_id, day = days,
                                              anc = s$m[, "anc_rel"] * pt$anc_base,
                                              plt = s$m[, "plt_rel"] * pt$plt_base,
                                              is = s$m[, "is"], im = s$m[, "im"],
                                              f_clone = s$m[, "f_clone"], f_rel = s$m[, "f_rel"])
    if (nrow(ev)) {
      ev <- data.table::copy(ev)[, patient_id := pt$patient_id]
      events_list[[i]] <- ev[, .(patient_id, event, start_day, end_day, dose)]
    }
    crs_list[[i]] <- data.table::data.table(patient_id = pt$patient_id, crs_grade = s$crs$grade,
                                            crs_onset_day = s$crs$onset, crs_end_day = s$crs$end)
  }
  truth <- data.table::rbindlist(truth_list)
  events <- data.table::rbindlist(events_list, fill = TRUE)
  pts <- merge(pts, data.table::rbindlist(crs_list), by = "patient_id")
  pts[, cd34_day := events[event == "cd34_boost"][match(pts$patient_id, patient_id), start_day]]
  truth[, f_rel_smooth := data.table::frollmean(f_rel, 10, align = "right", fill = 1), by = patient_id]

  # ---- truth phenotype and mechanism ----
  ph_true <- truth[, .(phenotype_true = classify_phenotype(day, anc, rule)), by = patient_id]
  mech <- vapply(seq_len(n), function(i) {
    pt <- pts[i]
    ph <- ph_true$phenotype_true[i]
    if (is.na(ph)) return(NA_character_)
    if (ph == "quick") return("none")
    if (ph == "second_wave") return("H2")
    st <- patient_struct(struct, spec, pt$clone_class)
    m0 <- marrow_sim(build_par(c(rho = pt$rho_true, g = pt$g_true, f0 = 0), st, "M2"),
                     make_segments(config$model$time$t0, max(days), pt$ld_scale), days)
    if (identical(classify_phenotype(days, m0[, "anc_rel"] * pt$anc_base, rule), "persistent")) "H1" else "H3"
  }, character(1))
  ph_true[, mechanism_true := mech]

  # ---- observed data ----
  message("Applying observation process ...")
  toci_day <- events[event == "tocilizumab"][, .(toci_day = min(start_day)), by = patient_id]
  pts <- merge(pts, toci_day, by = "patient_id", all.x = TRUE)
  labs_l <- vector("list", n); inflam_l <- vector("list", n); tx_l <- vector("list", n)
  for (i in seq_len(n)) {
    pt <- pts[i]
    td <- truth[patient_id == pt$patient_id]
    ol <- observe_labs(pt$patient_id, td, pt, spec, lab_schedule(pt$sparse, spec, pt$death_day))
    labs_l[[i]] <- ol$labs
    tx_l[[i]] <- ol$events
    inflam_l[[i]] <- observe_inflam(pt$patient_id, td, pt, spec, pt$toci_day, config)
  }
  labs <- data.table::rbindlist(labs_l)
  inflam <- data.table::rbindlist(inflam_l)
  events <- data.table::rbindlist(list(events, data.table::rbindlist(tx_l)), fill = TRUE)
  data.table::setorder(events, patient_id, start_day, event)

  # ---- split (post-sample patients forced into train) ----
  # post-treatment samples: patients with persistent cytopenia and a clone
  cand <- ph_true[phenotype_true == "persistent", patient_id]
  cand_clone <- intersect(cand, pts[ch_true == 1 & f0_true > 0.05, patient_id])
  post_ids <- head(c(sample(cand_clone), sample(setdiff(cand, cand_clone))), spec$n_post)
  n_train <- round(spec$train_frac * n)
  train_ids <- c(post_ids, sample(setdiff(ids, post_ids), n_train - length(post_ids)))
  split <- data.table::data.table(patient_id = ids, split = ifelse(ids %in% train_ids, "train", "test"),
                                  locked_at = "2026-11-02 19:00", seed = seed)

  # ---- mutation data ----
  seq_ids <- c(post_ids, sample(setdiff(ids, post_ids), spec$n_sequenced - length(post_ids)))
  mp <- observe_mutations_pre(sort(seq_ids), pts, spec, config$data$lod$vaf)
  post_planned <- data.table::data.table(patient_id = post_ids,
                                         sample_day = sample(60:85, length(post_ids), TRUE),
                                         sample_type = "BM", indication = "persistent cytopenia work-up")
  post_muts <- data.table::rbindlist(lapply(seq_len(nrow(post_planned)), function(k) {
    pid <- post_planned$patient_id[k]
    sd <- post_planned$sample_day[k]
    pt <- pts[patient_id == pid]
    f_t <- truth[patient_id == pid & day == sd, f_clone]
    vaf <- f_t / 2 * exp(stats::rnorm(1, 0, spec$noise[["vaf"]]))
    if (pt$ch_true == 1 && vaf >= config$data$lod$vaf)
      data.table::data.table(patient_id = pid, sample_day = sd, gene = pt$gene_true,
                             protein_change = ch_hotspots[[pt$gene_true]], vaf = round(vaf, 4),
                             depth = sample(1500:3000, 1))
  }))
  if (!nrow(post_muts))
    post_muts <- data.table::data.table(patient_id = character(0), sample_day = integer(0), gene = character(0),
                                        protein_change = character(0), vaf = numeric(0), depth = integer(0))

  # ---- clinical phenotype label (rule on observed data, ~3% disagreement) ----
  obs_ph <- labs[, .(ph = classify_phenotype(day, anc, rule)), by = patient_id]
  flip <- stats::runif(nrow(obs_ph)) < 0.03 & !is.na(obs_ph$ph)
  obs_ph[flip, ph := vapply(ph, function(x) sample(setdiff(phenotype_levels, x), 1), character(1))]

  # marrow involvement flags (stress: some left unflagged)
  pts[, marrow_involvement_flag := ifelse(marrow_involvement == 1 &
                                            stats::runif(.N) < spec$events$marrow_inv_unflagged, 0L,
                                          marrow_involvement)]

  patients <- pts[, .(patient_id, age, sex, product, costim, prior_lines, prior_asct, prior_bendamustine,
                      months_since_last_chemo, bridging,
                      ld_flu_total = ifelse(ld_dose_reduced == 1, 67.5, 90),
                      ld_cy_total = ifelse(ld_dose_reduced == 1, 1125, 1500), ld_dose_reduced,
                      crs_grade, crs_onset_day, crs_end_day,
                      marrow_involvement = marrow_involvement_flag, iec_hs, viral_infection,
                      relapse_day, subsequent_tx_day, death_day, last_followup_day)]
  patients[, phenotype_clinical := obs_ph$ph[match(patient_id, obs_ph$patient_id)]]

  # ---- red-cell indices (drawn last, so the other tables do not depend on them) ----
  labs <- add_red_cell_indices(labs, pts, spec)

  truth_pt <- merge(pts[, .(patient_id, rho = rho_true, g = g_true, f0 = f0_true, ch = ch_true,
                            clone_class, gene = gene_true, anc_base, plt_base, crp_base,
                            fer_base, marrow_involvement_true = marrow_involvement,
                            iec_hs, crs_grade)], ph_true, by = "patient_id")
  list(raw = list(patients = patients, labs = labs, inflam = inflam, events = events, split = split,
                  mut_tested = mp$tested, mut_pre = mp$muts, post_planned = post_planned),
       post_muts = post_muts,
       truth = list(patient = truth_pt, traj = truth,
                    model = list(scenario = scenario, seed = seed, struct = as.list(struct),
                                 clone_types = spec$clone_types, pop = spec$pop,
                                 calibrated = list(mu_rho = cal$mu_rho, mu_g = cal$mu_g))),
       spec = spec, cal = cal)
}
