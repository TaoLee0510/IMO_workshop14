# =============================================================================
# test_model.R
# Purpose : Unit tests for the parent model (v1.3: v1.2 lineage tree + CAR T -
#           target module + G-CSF input), the likelihood, DE-ITS pieces and the
#           Friberg reference model (run by step K1 via testthat).
# Inputs  : functions sourced by run_all.R; config in `aim4_test_config`.
# Outputs : testthat results.
# =============================================================================

cfg <- get("aim4_test_config", envir = globalenv())
st <- cfg$model$struct
des <- typical_design(cfg)
des_nc <- des
des_nc$no_car <- TRUE          # no CAR T infusion: the v1.2 dynamics
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
  d0 <- des_nc
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

test_that("Mab without HSC injury (incl. the CAR T - target route) equals M0", {
  i0 <- ind
  i0[c("th_as", "th_bs")] <- 0
  a <- sim(i0, "Mab", struct = replace(st, "th_zs", 0))
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
  y0[c("SX", "C", "I1", "I", "CT", "Z1", "Z", "GC")] <- 0
  seg <- make_segments(des_nc, st, 60, cfg)
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
  seg <- make_segments(des_nc, st, 90, cfg)
  a <- v11_a4_sim(p11, seg[, 1:6, drop = FALSE], z$v, z$t0, z$dt, tt)
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
  expect_identical(phen(sim(c(r_s = 0.03, f0 = 0, th_as = 1.5, th_ap = ap, th_bs = 0), "Ma", driver = z, design = des_nc)),
                   "second_wave")
  expect_identical(phen(sim(c(r_s = 0.02, f0 = 0, th_as = 2.5, th_ap = ap, th_bs = 0), "Ma", driver = z, design = des_nc)),
                   "delayed")
})

test_that("the CAR T phase alone can give a second wave through suppressed self-renewal (e_sr), not without it", {
  d2 <- driver_from_kernel(des$t0, 2.5, 9, cfg)
  i2 <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0)
  expect_identical(phen(sim(i2, "Mb", driver = d2, design = des_nc)), "second_wave")
  expect_identical(phen(sim(i2, "Mb", struct = replace(st, "e_sr", 0), driver = d2, design = des_nc)), "quick")
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


# ---- v1.3 --------------------------------------------------------------------------------------
test_that("v1.3 without the CAR T infusion and G-CSF reproduces the frozen v1.2 engine (< 1%)", {
  e12 <- new.env()
  Rcpp::sourceCpp(file.path(cfg$paths$code, "tests", "ref", "aim4_model_v12.cpp"),
                  cacheDir = cfg$paths$rcpp_cache, env = e12)
  for (id in c("Mab", "MabC_ddr", "Mabx")) {
    sp <- model_spec(id)
    s2 <- st
    if (sp$clone != "none") s2[names(cfg$model$clone_types[[sp$clone]])] <- cfg$model$clone_types[[sp$clone]]
    if (sp$q2 == "Mabx") s2[["eta_syn"]] <- 2
    p13 <- build_par(c(ind, f0 = 0.2), s2, sp)
    seg <- make_segments(des_nc, s2, 90, cfg)
    a <- a4_sim(p13, seg, drv$v, drv$t0, drv$dt, days, track = TRUE)
    b <- e12$a4_sim(p13[as.character(e12$a4_param_names())], seg[, 1:6, drop = FALSE], drv$v, drv$t0, drv$dt, days,
                    track = TRUE)
    cols <- setdiff(intersect(colnames(a), colnames(b)), "t")
    expect_lt(max(abs(a[, cols] / b[, cols] - 1), na.rm = TRUE), 0.01)
  }
})

test_that("CAR T compartments stay empty without the infusion; the typical design peaks on day 8 at 0.05", {
  m0 <- sim(ind, "Mab", design = des_nc)
  expect_equal(max(abs(m0$cart)), 0)
  expect_equal(max(abs(m0$z)), 0)
  m <- typical_cart_run(cfg)
  i <- which.max(m[, "cart"])
  expect_equal(unname(m[i, "t"]), calib_targets$cart_peak_day, tolerance = 0.5)
  expect_lt(abs(m[i, "cart"] / calib_targets$cart_peak - 1), 0.05)
})

test_that("the CAR T - target signal acts on the marrow only through th_zS (HSC) and e_Z (MPP)", {
  i2 <- c(ind, rho = 0.9)
  a <- sim(i2, "Mab", struct = replace(replace(st, "th_zs", 0), "e_z", 0))
  b <- sim(i2, "Mab", design = des_nc)
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-6)
  c1 <- sim(i2, "Mab")
  expect_gt(max(abs(c1$S - b$S)), 1e-3)                          # HSC route on
  # route variants switch the two C++ parameters (S still moves under _Zm through the faster
  # HSC regrowth while the MPP is depleted, which is not a loss)
  expect_equal(build_par(i2, st, model_spec("Mab_Zm"))[["th_zs"]], 0)
  expect_equal(build_par(i2, st, model_spec("Mab_Zh"))[["e_z"]], 0)
  expect_equal(build_par(i2, st, model_spec("Ma"))[["th_zs"]], 0)
  d1 <- sim(i2, "Mab_Zm")
  expect_gt(max(abs(d1$M - b$M)), 1e-3)
})

test_that("G-CSF releases neutrophils within a day and switches off completely", {
  g <- des_nc
  g$gcsf <- cbind(40, 41)
  i0 <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0)
  m <- sim(i0, "M0", design = g, track = TRUE)
  m0 <- sim(i0, "M0", design = des_nc)
  expect_gt(m[day == 41]$anc_rel / m[day == 40]$anc_rel, 1.5)
  expect_true(all(is.finite(as.matrix(m[, .(age_anc, gen_anc)]))))
  s0 <- replace(replace(st, "gcsf_fold", 1), "gcsf_rel", 0)
  # extra segment breakpoints only (integration path): differences at solver tolerance
  expect_lt(max(abs(sim(i0, "M0", struct = s0, design = g)$anc_rel - sim(i0, "M0", struct = s0, design = des_nc)$anc_rel)), 1e-6)
})

test_that("set-point rows are only in the baseline prior; ALC = b_alc * lymphocytes + CAR T", {
  dd <- read_aim4_data(cfg)
  pd <- build_pdats(dd, cfg, ids = dd$patients$patient_id[1])[[1]]
  o <- dd$obs[patient_id == pd$id & keep == TRUE & setpoint == FALSE & day > pd$design$t0]
  expect_equal(sum(vapply(pd$y, function(y) length(y$idx), integer(1))), nrow(o))
  pd$seg <- make_segments(pd$design, st, max(pd$times), cfg)
  sp <- model_spec("Mab")
  tm <- par_template(st, sp)
  psi <- c(r_s = log(0.08), th_ap = log(ap), th_as = log(0.1), th_bs = log(0.02), rho = log(0.9), b_alc = 0.2)
  mm <- simulate_pdat(psi, pd, tm, sp, cfg)
  b <- pdat_baselines(psi, pd)
  expect_equal(predict_obs(mm, pd, "alc", b), mm[pd$y$alc$idx, "alc_rel"] * pd$base$alc * exp(0.2) + mm[pd$y$alc$idx, "cart"])
})

test_that("DE with the chunked convergence rule finds the minimum of a quadratic and stops by itself", {
  c2 <- cfg
  c2$fit$de <- utils::modifyList(cfg$fit$de, list(np = 30L, chunk = 20L, window = 40L, max_gen = 4000L))
  set.seed(1)
  r <- de_run(function(x) sum((x - c(1, -2))^2), c(a = -5, b = -5), c(a = 5, b = 5), list(), c2)
  expect_true(r$converged)
  expect_lt(max(abs(r$par - c(1, -2))), 1e-2)
  expect_lt(r$gen, 4000)
})

test_that("the calibrated CAR T constants in config.R match the calibrations (2%)", {
  cc <- calibrate_cart(cfg)
  expect_lt(abs(cc[["kappa"]] / st[["kappa"]] - 1), 0.02)
  expect_lt(abs(cc[["rho"]] / exp(cfg$model$re_prior$mu[["rho"]]) - 1), 0.02)
  zz <- calibrate_z(cfg)
  expect_lt(abs(zz[["z50"]] / st[["z50"]] - 1), 0.02)
  expect_lt(abs(zz[["th_zs"]] / st[["th_zs"]] - 1), 0.02)
})
