# =============================================================================
# test_model.R
# Purpose : Unit tests for the parent model (v1.2 lineage tree) and the
#           Friberg reference model (run by step K1 via testthat).
# Inputs  : functions sourced by run_all.R; config in `aim4_test_config`.
# Outputs : testthat results.
# =============================================================================

cfg <- get("aim4_test_config", envir = globalenv())
st <- cfg$model$struct
des <- typical_design(cfg)
days <- seq(-30, 90)
drv <- driver_from_kernel(des$t0, 1.5, 8, cfg)
ap <- exp(cfg$model$re_prior$mu[["th_ap"]])
ind <- c(r_s = 0.08, f0 = 0, th_as = 0.2, th_ap = 1.5 * ap, th_bs = 0.03)
sim <- function(indiv, id, struct = st, design = des, driver = drv, track = FALSE) {
  sp <- model_spec(id)
  s2 <- struct
  if (sp$clone != "none") s2[names(cfg$model$clone_types[[sp$clone]])] <- cfg$model$clone_types[[sp$clone]]
  sim_daily(indiv, s2, sp, design, driver, cfg, days = days, track = track)
}
phen <- function(m) phenotype_sim(m$day, m$anc_rel * 3.5, cfg)

test_that("the parameter vector matches the C++ layout and self-renewal must stay below the exit rate", {
  expect_identical(names(build_par(ind, st, model_spec("Mab"))), as.character(a4_param_names()))
  expect_error(build_par(ind, replace(st, "sr_gmp", st[["k_gmp"]]), model_spec("Mab")), "self-renewal")
})

test_that("without treatment the system stays at its pre-treatment steady state (with a clone, with tracking)", {
  d0 <- des
  d0$egfr <- 1e12
  m <- sim(c(r_s = 0.08, f0 = 0.3, th_as = 0.2, th_ap = 1, th_bs = 0.05), "MabC_ddr", design = d0,
           driver = zero_driver(des$t0, cfg), track = TRUE)
  expect_lt(max(abs(as.matrix(m[, .(anc_rel, plt_rel, amc_rel, alc_rel, S, M, CMP, GMP, MEP, CLP, L)]) - 1)), 1e-5)
  expect_lt(max(abs(m$f_clone - 0.3)), 1e-6)
  # mean time and number of divisions since MPP entry equal the closed forms
  ss <- steady_age_generation(st)
  ref <- ss[match(c("ANC", "PLT", "AMC", "ALC"), compartment)]
  expect_lt(max(abs(as.matrix(m[, .(age_anc, age_plt, age_amc, age_alc)]) - rep(ref$age_days, each = nrow(m)))), 1e-4)
  expect_lt(max(abs(as.matrix(m[, .(gen_anc, gen_plt, gen_amc, gen_alc)]) - rep(ref$generations, each = nrow(m)))), 1e-4)
})

test_that("age and generation means stay finite and positive under treatment", {
  m <- sim(ind, "Mab", track = TRUE)
  x <- as.matrix(m[, .(age_anc, age_plt, age_amc, age_alc, gen_anc, gen_plt, gen_amc, gen_alc)])
  expect_true(all(is.finite(x)) && all(x > 0))
})

test_that("the clone model with f0 = 0 is identical to the model without a clone", {
  a <- sim(ind, "MabC_ddr")
  b <- sim(ind, "Mab")
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-10)
})

test_that("Mab without HSC injury equals M0", {
  i0 <- ind
  i0[c("th_as", "th_bs")] <- 0
  a <- sim(i0, "Mab")
  b <- sim(ind, "M0")
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-10)
})

test_that("Mabx with eta = 0 equals Mab", {
  s0 <- st
  s0[["eta_syn"]] <- 0
  a <- sim(ind, "Mabx", struct = s0)
  b <- sim(ind, "Mab")
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-10)
})

test_that("the lymphoid branch feeds nothing upstream: ANC, PLT, AMC ignore every lymphoid parameter", {
  s2 <- st
  s2[c("k_clp", "sr_clp", "k_ly", "sr_ly", "th_al", "c_l", "e_l", "fb_clp", "fb_ly")] <-
    c(0.4, 0.1, 0.2, 0.15, 3, 0.5, 0.2, 0, 0)
  seg <- make_segments(des, st, 90, cfg)
  run <- function(s) a4_sim(build_par(ind, s, model_spec("Mab")), seg, drv$v, drv$t0, drv$dt, days, 1e-10, 1e-12, 0.1)
  a <- run(st)
  b <- run(s2)
  expect_lt(max(abs(a[, c("anc_rel", "plt_rel", "amc_rel", "S", "M")] - b[, c("anc_rel", "plt_rel", "amc_rel", "S", "M")])), 1e-7)
  expect_gt(max(abs(a[, "alc_rel"] - b[, "alc_rel"])), 0.01)
})

test_that("the C++ integrator agrees with deSolve", {
  par <- build_par(ind, st, model_spec("Mab"))
  nm <- as.character(a4_state_names())
  y0 <- stats::setNames(rep(1, length(nm)), nm)
  y0[c("SX", "C", "I1", "I")] <- 0
  seg <- make_segments(des, st, 60, cfg)
  rhs <- function(t, y, p) {
    j <- max(which(seg[, 1] <= t))
    list(a4_rhs(t, y, par, seg[j, "dose"], 0, 0, drv$v, drv$t0, drv$dt))
  }
  tt <- seq(des$t0, 60, by = 0.5)
  ref <- deSolve::lsoda(y0, tt, rhs, NULL, rtol = 1e-9, atol = 1e-11, hmax = 0.02)
  cpp <- a4_sim(par, seg, drv$v, drv$t0, drv$dt, tt, 1e-9, 1e-11, 0.1)
  expect_lt(max(abs(ref[, "BN"] - cpp[, "anc_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BP"] - cpp[, "plt_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BM"] - cpp[, "amc_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BL"] - cpp[, "alc_rel"])), 1e-4)
})

test_that("regression: with v1.1 settings the v1.2 tree reproduces the v1.1 engine (< 1%)", {
  # CMP an instantaneous relay, GMP / MEP without self-renewal, feedback on the MPP
  # only, gam_l = 0, equal N / Mo transit (they share the GMP), no CAR T phase
  load_v11_reference(cfg)
  p11 <- v11_par(c(r_s = 0.08, f0 = 0, th_as = 0.3, th_ap = 1.5, th_bs = 0),
                 c(gam_n = 0, gam_pl = 0, gam_mo = 0, ktr_mo = 4 / 5))
  s12 <- v12_as_v11(p11, st)
  z <- zero_driver(des$t0, cfg)
  tt <- seq(-30, 90, by = 0.5)
  seg <- make_segments(des, st, 90, cfg)
  a <- v11_a4_sim(p11, seg, z$v, z$t0, z$dt, tt)
  b <- a4_sim(build_par(c(r_s = 0.08, f0 = 0, th_as = 0.3, th_ap = 1.5, th_bs = 0), s12, model_spec("Ma")),
              seg, z$v, z$t0, z$dt, tt)
  for (v in c("anc_rel", "plt_rel", "amc_rel", "S", "M")) expect_lt(max(abs(a[, v] / b[, v] - 1)), 0.01)
})

test_that("the pre-treatment steady state is stable over the gamma_M x self-renewal grid", {
  g <- stability_grid(cfg)
  expect_true(all(g$max_re < 0))
})

test_that("the calibrated th_aP centre and th_aL in config.R match the calibrations (2%)", {
  expect_lt(abs(calibrate_th_ap(cfg) / ap - 1), 0.02)
  expect_lt(abs(calibrate_th_al(cfg) / st[["th_al"]] - 1), 0.02)
})

test_that("one chemotherapy hit can give a pseudo second wave, and a deeper one delayed recovery", {
  z <- zero_driver(des$t0, cfg)
  expect_identical(phen(sim(c(r_s = 0.03, f0 = 0, th_as = 1.5, th_ap = ap, th_bs = 0), "Ma", driver = z)), "second_wave")
  expect_identical(phen(sim(c(r_s = 0.02, f0 = 0, th_as = 2.5, th_ap = ap, th_bs = 0), "Ma", driver = z)), "delayed")
})

test_that("the CAR T phase alone can give a second wave through suppressed self-renewal (e_sr), not without it", {
  d2 <- driver_from_kernel(des$t0, 2.5, 9, cfg)
  i2 <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0)
  expect_identical(phen(sim(i2, "Mb", driver = d2)), "second_wave")
  expect_identical(phen(sim(i2, "Mb", struct = replace(st, "e_sr", 0), driver = d2)), "quick")
})

test_that("MF agrees with a deSolve implementation of the original Friberg model", {
  sp <- model_spec("MF")
  ms <- cfg$model$mf_struct
  par <- build_par(c(slope = 1.5, mtt_scale = 1), ms, sp)
  seg <- make_segments(des, ms, 60, cfg)
  k <- 4 / ms[["mtt_n"]]
  fr <- function(t, y, p) {
    j <- max(which(seg[, 1] <= t))
    e <- 1.5 * y[1]
    list(c(seg[j, "dose"] - ms[["k_e"]] * y[1],
           k * y[2] * (1 - e) * y[6]^-ms[["gam_mf"]] - k * y[2],
           k * (y[2] - y[3]), k * (y[3] - y[4]), k * (y[4] - y[5]), k * (y[5] - y[6])))
  }
  tt <- seq(des$t0, 60, by = 0.5)
  ref <- deSolve::lsoda(c(0, 1, 1, 1, 1, 1), tt, fr, NULL, rtol = 1e-9, atol = 1e-11, hmax = 0.05)
  cpp <- mf_sim(par, seg, tt, 1e-9, 1e-11, 0.1)
  expect_lt(max(abs(ref[, 7] - cpp[, "anc_rel"])), 1e-4)
})
