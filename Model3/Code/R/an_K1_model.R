# =============================================================================
# an_K1_model.R
# Purpose : Step K1 - unit tests of Model3 and MF; typical trajectories of every
#           nested model (prior-centre parameters, typical CAR T-phase driver);
#           steady-state stability; the comparison with Aim4 v1.3; the
#           calibrations (th_aP centre, lymphocytes, CAR T module) and the
#           reachable phenotypes.
# Inputs  : ctx (config).
# Outputs : list(tables, tests, traj, stability, aim4, screen, all_passed)
# =============================================================================

an_K1_model <- function(ctx) {
  cfg <- ctx$config
  assign("model3_test_config", cfg, envir = globalenv())
  load_aim4_reference(cfg)
  rep <- testthat::test_file(file.path(cfg$paths$code, "tests", "test_model.R"),
                             reporter = testthat::ListReporter$new(), stop_on_failure = FALSE)
  tests <- data.table::rbindlist(lapply(rep, function(r) {
    ok <- all(vapply(r$results, function(x) inherits(x, "expectation_success"), logical(1)))
    data.table::data.table(test = r$test, passed = ok)
  }))
  ctx$log("  K1: ", sum(tests$passed), "/", nrow(tests), " unit tests passed")

  des <- typical_design(cfg)
  drv <- driver_from_kernel(des$t0, 1.2, 8, cfg)
  ind <- c(r_s = 0.08, f0 = 0.25, th_as = 0.3, th_ap = 1.5 * exp(cfg$model$re_prior$mu[["th_ap"]]), th_bs = 0.04,
           rho = exp(cfg$model$re_prior$mu[["rho"]]), a0 = 1)
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
                   base = list(anc = 3.5, plt = 200, amc = 0.6, alc = 1.5))
    m[, model := id]
    m[, phenotype := phenotype_sim(day, anc, cfg, with_persistent = TRUE)]
    m
  }), fill = TRUE)

  stability <- stability_grid(cfg)
  aim4 <- aim4_comparison(cfg)
  th_ap <- calibrate_th_ap(cfg)
  lymph <- calibrate_lymph(cfg)
  th_am <- calibrate_th_am(cfg)
  cart <- data.table::as.data.table(as.list(c(calibrate_cart(cfg), calibrate_z(cfg))))
  screen <- phenotype_screen(cfg)
  ctx$log("  K1: steady state stable for gamma up to its box limit at sr_p <= its value: ",
          all(stability[sr_scale <= 1]$max_re < 0), " (with sr_p x 1.2 from gamma ",
          suppressWarnings(min(stability[sr_scale > 1 & max_re >= 0]$gamma)), ")",
          "; ANC recovery day Aim4 / Model3: ", paste(aim4[lineage == "ANC", recovery_day], collapse = " / "),
          "; th_aP matching the Aim4 nadir ", round(th_ap, 3), " (not used: centre kept at 2.71)", "; th_aL ", round(lymph[["th_al"]], 3), ", c_L ", round(lymph[["c_l"]], 3), "; th_aM ", round(th_am, 3))

  list(tables = list(unit_tests = tests,
                     typical = traj[, .(phenotype = phenotype[1], anc_min = min(anc), anc_d21 = anc[day == 21],
                                        anc_d30 = anc[day == 30], plt_d30 = plt[day == 30],
                                        alc_d30 = alc[day == 30],
                                        S_d30 = S[day == 30]), by = model],
                     stability = stability, aim4_comparison = aim4,
                     calibration = data.table::data.table(th_ap_aim4_nadir = th_ap, th_am = th_am, th_al = lymph[["th_al"]],
                                                          c_l = lymph[["c_l"]], alc_d30_target = lymph[["alc_d30_target"]]),
                     cart_calibration = cart, screen_b_only = screen$b_only, screen_a_only = screen$a_only),
       tests = tests, traj = traj, stability = stability, aim4 = aim4, screen = screen, cart = cart,
       all_passed = all(tests$passed))
}
