# =============================================================================
# test_model.R
# Purpose : Unit tests for the marrow twin (run by step D1 via testthat).
# Inputs  : functions sourced by run_all.R; config in `aim2_test_config`.
# Outputs : testthat results.
# =============================================================================

cfg <- get("aim2_test_config", envir = globalenv())
st <- design_struct(cfg)
days <- cfg$model$time$sim_days
seg_ld <- make_segments(cfg$model$time$t0, max(days), ld_scale = 1)
sim <- function(indiv, model, seg = seg_ld, struct = st) {
  marrow_sim(build_par(indiv, struct, model), seg, days)
}

test_that("M2 with f0 = 0 is identical to M1", {
  a <- sim(c(rho = 0.8, g = 0.6, f0 = 0), "M2")
  b <- sim(c(rho = 0.8, g = 0.6), "M1")
  expect_lt(max(abs(a[, "anc_rel"] - b[, "anc_rel"])), 1e-10)
})

test_that("M2 with a vanishing clone stays within 1% of M1", {
  a <- sim(c(rho = 0.8, g = 0.6, f0 = 1e-6), "M2")
  b <- sim(c(rho = 0.8, g = 0.6), "M1")
  expect_lt(max(abs(a[, "anc_rel"] / b[, "anc_rel"] - 1)), 0.01)
})

test_that("M1 without marrow effects of inflammation equals M0", {
  s1 <- st
  s1[c("q", "emax")] <- 0
  a <- sim(c(rho = 0.8, g = 0.6), "M1", struct = s1)
  b <- sim(c(rho = 0.8, g = 0.6), "M0")
  expect_lt(max(abs(a[, "anc_rel"] - b[, "anc_rel"])), 1e-10)
})

test_that("without treatment the system stays at its pre-treatment steady state", {
  seg0 <- make_segments(cfg$model$time$t0, max(days), ld_scale = 0)
  m <- sim(c(rho = 0.7, g = 0, f0 = 0.3), "M2", seg = seg0)
  expect_lt(max(abs(m[, "anc_rel"] - 1)), 1e-6)
  expect_lt(max(abs(m[, "plt_rel"] - 1)), 1e-6)
  expect_lt(max(abs(m[, "f_clone"] - 0.3)), 1e-6)
})

test_that("a neutral clone conserves the stem pool and its own fraction", {
  sn <- st
  sn[c("r_ld", "r_inf", "delta")] <- 1
  a <- sim(c(rho = 0.8, g = 0.6, f0 = 0.25), "M2", struct = sn)
  b <- sim(c(rho = 0.8, g = 0.6, f0 = 0), "M2", struct = sn)
  expect_lt(max(abs(a[, "x_tot"] - b[, "x_tot"])), 1e-6)
  expect_lt(max(abs(a[, "f_clone"] - 0.25)), 1e-6)
})

test_that("C++ integrator agrees with deSolve::lsoda", {
  ev <- data.table::data.table(event = c("gcsf", "tocilizumab", "steroid"),
                               start_day = c(20, 5, 6), end_day = c(24, 5, 9))
  seg <- make_segments(cfg$model$time$t0, max(days), ld_scale = 1, events = ev)
  par <- build_par(c(rho = 0.6, g = 0.8, f0 = 0.3), st, "M2")
  cpp <- marrow_sim(par, seg, days, rtol = 1e-9, atol = 1e-12, hmax = 0.1)
  rhs <- function(t, y, p) {
    k <- max(which(seg[, "start"] <= t + 1e-12))
    list(marrow_rhs(t, y, p, seg[k, input_columns]))
  }
  times <- c(cfg$model$time$t0, sort(unique(c(days, seg[, "end"]))))
  ref <- deSolve::lsoda(marrow_init(par), times, rhs, par, rtol = 1e-10, atol = 1e-12,
                        tcrit = max(times), hmax = 0.1)
  ref_anc <- ref[match(days, ref[, 1]), 1 + which(marrow_state_names() == "N")]
  expect_lt(max(abs(cpp[, "anc_rel"] / ref_anc - 1)), 1e-4)
})

test_that("counts fall after lymphodepletion and recover with little inflammation", {
  m <- sim(c(rho = 1, g = 0.05, f0 = 0), "M2")
  expect_lt(min(m[days %in% 0:10, "anc_rel"]), 0.5)
  expect_gt(m[days == 28, "anc_rel"], 0.8)
})

test_that("phenotype rule labels canonical trajectories", {
  rule <- cfg$phenotype
  d <- 0:90
  quick <- ifelse(d < 7, 0.3, 3)
  wave <- ifelse(d < 7, 0.3, ifelse(d >= 25 & d <= 35, 0.4, 3))
  pers <- ifelse(d < 7, 0.3, 0.2)
  expect_equal(classify_phenotype(d, quick, rule), "quick")
  expect_equal(classify_phenotype(d, wave, rule), "second_wave")
  expect_equal(classify_phenotype(d, pers, rule), "persistent")
  short <- if (identical(rule$scheme, "early")) 0:15 else 0:30
  expect_true(is.na(classify_phenotype(short, rep(3, length(short)), rule)))
  # early scheme: high counts before the nadir are not a recovery
  if (identical(rule$scheme, "early"))
    expect_equal(classify_phenotype(d, ifelse(d < 2, 3, 0.3), rule), "persistent")
})

test_that("segments place events on the correct intervals", {
  ev <- data.table::data.table(event = c("gcsf", "tocilizumab"), start_day = c(20, 5), end_day = c(24, 5))
  seg <- make_segments(-6, 90, ld_scale = 0.75, events = ev)
  expect_equal(unname(seg[seg[, "start"] == -5, "ld"]), 0.75)
  expect_equal(unname(seg[seg[, "start"] == 20, "gcsf"]), 1)
  expect_equal(unname(seg[seg[, "start"] >= 25, "gcsf"]), rep(0, sum(seg[, "start"] >= 25)))
  expect_true(all(seg[seg[, "start"] >= 5, "toci"] == 1))
})
