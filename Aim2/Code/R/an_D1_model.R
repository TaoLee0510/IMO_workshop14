# =============================================================================
# an_D1_model.R
# Purpose : D1 - model correctness. Runs the unit tests and simulates the
#           archetypal trajectories of the three nested models.
# Inputs  : ctx$config (design values; no patient data).
# Outputs : list(tables = tests/archetypes, traj, all_passed)
# =============================================================================

an_D1_model <- function(ctx) {
  cfg <- ctx$config
  assign("aim2_test_config", cfg, envir = globalenv())
  res <- testthat::test_file(file.path(cfg$paths$code, "tests", "test_model.R"),
                             reporter = testthat::SilentReporter, stop_on_failure = FALSE)
  tests <- data.table::as.data.table(as.data.frame(res))[, .(test, n_expectations = nb, failed,
                                                             error, skipped, real_time = round(real, 2))]
  tests[, passed := failed == 0 & !error & !skipped]

  st <- design_struct(cfg)
  b <- cfg$model$design$base
  archetypes <- list(
    quick       = c(rho = 1.0, g = 0.15, f0 = 0),
    second_wave = c(rho = 1.0, g = 0.70, f0 = 0),
    persistent_H1 = c(rho = 0.2, g = 0.80, f0 = 0),
    persistent_H3 = c(rho = 1.0, g = 1.0, f0 = 0.6))
  traj <- data.table::rbindlist(lapply(names(archetypes), function(a) {
    data.table::rbindlist(lapply(cfg$fit$models, function(mdl) {
      a0 <- baseline_with_clone(b$anc, archetypes[[a]][["f0"]], st)
      s <- sim_daily(build_par(archetypes[[a]], st, mdl), cfg,
                     base = list(anc = a0, plt = b$plt, crp = b$crp, ferritin = b$ferritin))
      s[, `:=`(archetype = a, model = mdl,
               phenotype = classify_phenotype(day, anc, cfg$phenotype))]
    }))
  }))
  arch_tab <- traj[, .(phenotype = phenotype[1], anc_d14 = round(anc[day == 14], 2),
                       anc_min_21_60 = round(min(anc[day >= 21 & day <= 60]), 2),
                       anc_d90 = round(anc[day == 90], 2), f_clone_d90 = round(f_clone[day == 90], 2)),
                   by = .(archetype, model)]

  list(tables = list(unit_tests = tests, archetypes = arch_tab),
       traj = traj, all_passed = all(tests$passed),
       key_numbers = list(tests_passed = sum(tests$passed), tests_total = nrow(tests)))
}
