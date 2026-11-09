# =============================================================================
# an_K1_model.R
# Purpose : Step K1 - unit tests of the parent model (v1.2 lineage tree) and
#           MF; typical trajectories of every nested model (prior-centre
#           parameters, typical CAR T-phase driver) with the mean time and
#           number of divisions since MPP entry of the blood cells; the
#           pre-treatment ages / generations per compartment; steady-state
#           stability; the calibrations (v1.1 comparison, progenitor screen).
# Inputs  : ctx (config).
# Outputs : list(tables, tests, traj, steady, stability, v11, screen, all_passed)
# =============================================================================

an_K1_model <- function(ctx) {
  cfg <- ctx$config
  assign("aim4_test_config", cfg, envir = globalenv())
  load_v11_reference(cfg)
  rep <- testthat::test_file(file.path(cfg$paths$code, "tests", "test_model.R"),
                             reporter = testthat::ListReporter$new(), stop_on_failure = FALSE)
  tests <- data.table::rbindlist(lapply(rep, function(r) {
    ok <- all(vapply(r$results, function(x) inherits(x, "expectation_success"), logical(1)))
    data.table::data.table(test = r$test, passed = ok)
  }))
  ctx$log("  K1: ", sum(tests$passed), "/", nrow(tests), " unit tests passed")

  des <- typical_design(cfg)
  drv <- driver_from_kernel(des$t0, 1.2, 8, cfg)
  ind <- c(r_s = 0.08, f0 = 0.25, th_as = 0.3, th_ap = 1.5 * exp(cfg$model$re_prior$mu[["th_ap"]]), th_bs = 0.04)
  ids <- c("M0", "Ma", "Mb", "Mab", "MabC_ddr", "MabC_dta", "MF")
  traj <- data.table::rbindlist(lapply(ids, function(id) {
    sp <- model_spec(id)
    i2 <- if (sp$family == "mf") c(slope = 1.5, mtt_scale = 1) else {
      x <- ind
      if (sp$clone == "none") x[["f0"]] <- 0
      if (!"th_as" %in% sp$re) x[["th_as"]] <- 0
      if (!"th_bs" %in% sp$re) x[["th_bs"]] <- 0
      x
    }
    m <- sim_daily(i2, default_struct(sp, cfg), sp, des, drv, cfg, days = seq(-10, 90),
                   base = list(anc = 3.5, plt = 200, amc = 0.6, alc = 1.5), track = sp$family == "a4")
    m[, model := id]
    m[, phenotype := phenotype_sim(day, anc, cfg, with_persistent = TRUE)]
    m
  }), fill = TRUE)

  bk <- cfg$model$bookkeeping
  steady <- steady_age_generation(cfg$model$struct)
  steady[, note := data.table::fcase(
    compartment == "PLT", sprintf("generation of the megakaryocyte; + %d endomitotic rounds (%d-%d, ploidy %d-%dN), ~%g platelets per megakaryocyte",
                                  bk$endomitosis_rounds, bk$endomitosis_range[1], bk$endomitosis_range[2],
                                  bk$ploidy_range[1], bk$ploidy_range[2], bk$platelets_per_megakaryocyte),
    default = "")]
  stability <- stability_grid(cfg)
  v11 <- v11_comparison(cfg)
  screen <- progenitor_screen(cfg)
  ctx$log("  K1: steady state stable over the grid: ", all(stability$max_re < 0),
          "; ANC recovery day v1.1 / v1.2: ", paste(v11[lineage == "ANC", recovery_day], collapse = " / "))

  list(tables = list(unit_tests = tests,
                     typical = traj[, .(phenotype = phenotype[1], anc_min = min(anc), anc_d21 = anc[day == 21],
                                        anc_d30 = anc[day == 30], plt_d30 = plt[day == 30],
                                        alc_d30 = alc[day == 30],
                                        S_d30 = S[day == 30]), by = model],
                     steady_age_generation = steady, stability = stability, v11_comparison = v11,
                     progenitor_screen = screen),
       tests = tests, traj = traj, steady = steady, stability = stability, v11 = v11, screen = screen,
       all_passed = all(tests$passed))
}
