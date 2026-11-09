# =============================================================================
# test_model.R
# Purpose : Unit tests for Model3 (Model1 + chemotherapy kill of the chain
#           precursors and circulating monocytes; Friberg-like haematopoietic layer + CAR T -
#           target module + G-CSF input), the phenotype term, the likelihood,
#           DE-ITS pieces (incl. round 0) and the Friberg reference model (run by
#           step K1 via testthat).
# Inputs  : functions sourced by run_all.R; config in `model3_test_config`.
# Outputs : testthat results.
# =============================================================================

cfg <- get("model3_test_config", envir = globalenv())
st <- cfg$model$struct
des <- typical_design(cfg)
des_nc <- des
des_nc$no_car <- TRUE          # no CAR T infusion
days <- seq(-30, 90)
drv <- driver_from_kernel(des$t0, 1.5, 8, cfg)
ap <- exp(cfg$model$re_prior$mu[["th_ap"]])
ind <- c(r_s = 0.08, f0 = 0, th_as = 0.2, th_ap = 1.5 * ap, th_bs = 0.03)
sim <- function(indiv, id, struct = st, design = des, driver = drv) {
  sp <- model_spec(id)
  s2 <- struct
  if (sp$clone != "none") s2[names(cfg$model$clone_types[[sp$clone]])] <- cfg$model$clone_types[[sp$clone]]
  sim_daily(indiv, s2, sp, design, driver, cfg, days = days)
}
phen <- function(m) phenotype_sim(m$day, m$anc_rel * 3.5, cfg)
# Model2 values of the Model3 changes (N: 3 compartments, 2 killed; P1 off; P2: e_z,Pl = e_z)
m2_values <- function(s) replace(s, c("n_tr_n", "nk_n", "th_cp", "d_zpl"), c(3, 2, 0, 0))

test_that("the parameter vector matches the C++ layout and self-renewal must stay below the exit rate", {
  expect_identical(names(build_par(ind, st, model_spec("Mab"))), as.character(m3_param_names()))
  expect_error(build_par(ind, replace(st, "sr_p", st[["k_p"]]), model_spec("Mab")), "self-renewal")
  expect_length(st, 43L)
})

test_that("without treatment the system stays at its pre-treatment steady state (with a clone)", {
  d0 <- des_nc
  d0$egfr <- 1e12
  m <- sim(c(r_s = 0.08, f0 = 0.3, th_as = 0.2, th_ap = 1, th_bs = 0.05), "MabC_ddr", design = d0,
           driver = zero_driver(des$t0, cfg))
  expect_lt(max(abs(as.matrix(m[, .(anc_rel, plt_rel, amc_rel, alc_rel, S, PN, PMO, PPL, PL)]) - 1)), 1e-5)
  expect_lt(max(abs(m$f_clone - 0.3)), 1e-6)
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
  a <- sim(ind, "Mabx", struct = replace(st, "eta_syn", 0))
  b <- sim(ind, "Mab")
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-10)
})

test_that("the lymphoid component feeds nothing else: ANC, PLT, AMC ignore c_L and th_aL", {
  s2 <- st
  s2[c("th_al", "c_l")] <- c(3, 0.5)
  seg <- make_segments(des, st, 90, cfg)
  run <- function(s) m3_sim(build_par(ind, s, model_spec("Mab")), seg, drv$v, drv$t0, drv$dt, days, 1e-10, 1e-12, 0.1)
  a <- run(st)
  b <- run(s2)
  cols <- c("anc_rel", "plt_rel", "amc_rel", "S", "M", "PL")
  expect_lt(max(abs(a[, cols] - b[, cols])), 1e-7)
  expect_gt(max(abs(a[, "alc_rel"] - b[, "alc_rel"])), 0.01)
})

test_that("the C++ integrator agrees with deSolve", {
  par <- build_par(ind, st, model_spec("Mab"))
  nm <- as.character(m3_state_names())
  y0 <- stats::setNames(rep(1, length(nm)), nm)
  y0[c("SX", "C", "I1", "I", "CT", "Z1", "Z", "GC")] <- 0
  seg <- make_segments(des_nc, st, 60, cfg)
  rhs <- function(t, y, p) {
    j <- max(which(seg[, 1] <= t))
    list(m3_rhs(t, y, par, seg[j, "dose"], 0, 0, drv$v, drv$t0, drv$dt))
  }
  tt <- seq(des$t0, 60, by = 0.5)
  ref <- deSolve::lsoda(y0, tt, rhs, NULL, rtol = 1e-9, atol = 1e-11, hmax = 0.02)
  cpp <- m3_sim(par, seg, drv$v, drv$t0, drv$dt, tt, 1e-9, 1e-11, 0.1)
  expect_lt(max(abs(ref[, "BN"] - cpp[, "anc_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BP"] - cpp[, "plt_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BM"] - cpp[, "amc_rel"])), 1e-4)
  expect_lt(max(abs(ref[, "BL"] - cpp[, "alc_rel"])), 1e-4)
})

test_that("Friberg limit: with sr_P -> k_P an MPP component is Friberg's proliferating compartment", {
  # no stem-cell input; feedback min(gmax, B^-gamma) on proliferation at rate k_P; chemotherapy
  # as a kill term; chain of three transit compartments and blood
  s2 <- replace(m2_values(st), "sr_p", st[["k_p"]] * (1 - 1e-12))   # three compartments as Friberg
  s2[c("kappa_c", "th_am")] <- 0                      # Model2 kill terms off: Friberg's chain has no kill
  par <- build_par(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0), s2, model_spec("M0"))
  seg <- make_segments(des_nc, s2, 60, cfg)
  z <- zero_driver(des$t0, cfg)
  k <- s2[["k_p"]]; kt <- 3 / s2[["mtt_n"]]
  fr <- function(t, y, p) {
    j <- max(which(seg[, 1] <= t))
    g <- min(s2[["gmax"]], y[6]^-s2[["gamma"]])
    list(c(seg[j, "dose"] - s2[["k_e"]] * y[1],
           k * y[2] * (g - 1) - ap * y[1] * y[2],
           kt * (y[2] - y[3]), kt * (y[3] - y[4]), kt * (y[4] - y[5]), s2[["c_n"]] * (y[5] - y[6])))
  }
  tt <- seq(des$t0, 60, by = 0.5)
  ref <- deSolve::lsoda(c(0, 1, 1, 1, 1, 1), tt, fr, NULL, rtol = 1e-9, atol = 1e-11, hmax = 0.05)
  cpp <- m3_sim(par, seg, z$v, z$t0, z$dt, tt, 1e-9, 1e-11, 0.1)
  expect_lt(max(abs(ref[, 7] - cpp[, "anc_rel"])), 1e-4)
})

test_that("the pre-treatment steady state is stable for gamma up to its box limit (sr_p at or below its value)", {
  # sr_p is fixed; the rows with sr_p x 1.2 (in the K1 table) lose stability from gamma ~ 1.5
  g <- stability_grid(cfg)
  expect_true(all(g[sr_scale <= 1]$max_re < 0))
})

test_that("the calibrated th_aM, th_aL and c_L in config.R match the calibrations (2%)", {
  # th_aP keeps the Model1 centre (2.71); with kappa_c = 1 the typical ANC nadir is deeper than in Model1
  expect_lt(abs(calibrate_th_am(cfg) / st[["th_am"]] - 1), 0.02)
  ly <- calibrate_lymph(cfg)
  expect_lt(abs(ly[["th_al"]] / st[["th_al"]] - 1), 0.02)
  expect_lt(abs(ly[["c_l"]] / st[["c_l"]] - 1), 0.02)
})

test_that("reachable phenotypes: the CAR T phase can give a second wave (through emax), a chemotherapy hit alone gives quick or delayed", {
  d2 <- driver_from_kernel(des$t0, 2.5, 9, cfg)
  i2 <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0.06, rho = exp(cfg$model$re_prior$mu[["rho"]]), a0 = 1)
  expect_identical(phen(sim(i2, "Mb", driver = d2)), "second_wave")
  expect_identical(phen(sim(i2, "Mb", struct = replace(replace(st, "emax", 0), "e_z", 0), driver = d2)), "quick")
  z <- zero_driver(des$t0, cfg)
  ph <- vapply(c(0.4, 1, 1.5, 2, 2.5, 3), function(a)
    phen(sim(c(r_s = 0.03, f0 = 0, th_as = a, th_ap = ap, th_bs = 0), "Ma", driver = z, design = des_nc)), character(1))
  expect_false("second_wave" %in% ph)
  expect_identical(ph[length(ph)], "delayed")
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
  # extra breakpoint at the infusion (integration path): differences at solver tolerance
  expect_lt(max(abs(a$anc_rel - b$anc_rel)), 1e-5)
  c1 <- sim(i2, "Mab")
  expect_gt(max(abs(c1$S - b$S)), 1e-3)                          # HSC route on
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
  m <- sim(i0, "M0", design = g)
  expect_gt(m[day == 41]$anc_rel / m[day == 40]$anc_rel, 1.5)
  s0 <- replace(replace(st, "gcsf_fold", 1), "gcsf_rel", 0)
  # extra segment breakpoints only (integration path): differences at solver tolerance
  expect_lt(max(abs(sim(i0, "M0", struct = s0, design = g)$anc_rel - sim(i0, "M0", struct = s0, design = des_nc)$anc_rel)), 1e-6)
})

test_that("set-point rows are only in the baseline prior; ALC = b_alc * lymphocytes + CAR T", {
  dd <- read_model3_data(cfg)
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

# ---- phenotype term ------------------------------------------------------------------------------
test_that("phenotype_prob: probabilities sum to one, s -> 0 is the rule on the model series, continuous at the thresholds", {
  r <- phenotype_rule(cfg)
  i2 <- c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0.06, rho = 0.811, a0 = 1)
  runs <- list(sim(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = 0), "M0", design = des_nc,
                   driver = zero_driver(des$t0, cfg)),
               sim(i2, "Mb", driver = driver_from_kernel(des$t0, 2.5, 9, cfg)),
               sim(c(r_s = 0.03, f0 = 0, th_as = 2.5, th_ap = ap, th_bs = 0), "Ma", design = des_nc,
                   driver = zero_driver(des$t0, cfg)))
  for (m in runs) {
    w <- m[day >= 0 & day <= 30]
    p <- phenotype_prob(w$day, w$anc_rel * 3.5, r, cfg$fit$pheno$s, cfg$fit$pheno$eps)
    expect_equal(sum(p), 1, tolerance = 1e-12)
    p0 <- phenotype_prob(w$day, w$anc_rel * 3.5, r, 1e-6)
    expect_identical(names(which.max(p0)), phen(m))
  }
  dd <- 0:30
  step <- function(v) ifelse(dd < 10, 0.3, v)
  a <- phenotype_prob(dd, step(1.5 * exp(-1e-6)), r, 0.3)
  b <- phenotype_prob(dd, step(1.5 * exp(1e-6)), r, 0.3)
  expect_lt(max(abs(a - b)), 1e-5)
})

test_that("the phenotype term enters the individual objective with weight w and not at landmarks", {
  dd <- read_model3_data(cfg)
  pid <- dd$phenotype[primary == "second_wave"]$patient_id[1]
  sp <- model_spec("Mab")
  tm <- par_template(st, sp)
  psi <- c(r_s = log(0.08), th_ap = log(ap), th_as = log(0.1), th_bs = log(0.02), rho = log(0.9))
  pd <- prep_segments(build_pdats(dd, cfg, ids = pid), st, cfg)[[1]]
  expect_identical(pd$pheno$class, "second_wave")
  pop <- list(mu = cfg$model$re_prior$mu[sp$re], omega = cfg$model$re_prior$omega[sp$re],
              beta = lapply(stats::setNames(sp$re, sp$re), function(k) numeric(0)))
  prior <- patient_prior(pd, pop, sp$re)
  sg <- cfg$model$sigma_init
  c0 <- cfg
  c0$fit$pheno$w <- 0
  f1 <- neg_log_post(psi, pd, tm, sp, prior, sg, cfg)
  f0 <- neg_log_post(psi, pd, tm, sp, prior, sg, c0)
  mm <- simulate_pdat(psi, pd, tm, sp, cfg)
  expect_equal(f0 - f1, pheno_loglik(mm, pd, pdat_baselines(psi, pd), cfg), tolerance = 1e-10)
  expect_lt(pheno_loglik(mm, pd, pdat_baselines(psi, pd), cfg), 0)
  expect_null(build_pdats(dd, cfg, landmark = 14, ids = pid)[[1]]$pheno)
  expect_null(build_pdats(dd, c0, ids = pid)[[1]]$pheno)
})

# ---- estimation ----------------------------------------------------------------------------------
test_that("DE with the chunked convergence rule finds the minimum of a quadratic and stops by itself", {
  c2 <- cfg
  c2$fit$de <- utils::modifyList(cfg$fit$de, list(np = 30L, chunk = 20L, window = 40L, max_gen = 4000L))
  set.seed(1)
  r <- de_run(function(x) sum((x - c(1, -2))^2), c(a = -5, b = -5), c(a = 5, b = 5), list(), c2)
  expect_true(r$converged)
  expect_lt(max(abs(r$par - c(1, -2))), 1e-2)
  expect_lt(r$gen, 4000)
})

test_that("round 0: widened prior, unchanged box, and the population of round 1 from the M-step formulas", {
  dd <- read_model3_data(cfg)
  ids <- dd$patients$patient_id[1:3]
  d3 <- list(patients = dd$patients[patient_id %in% ids], obs = dd$obs[patient_id %in% ids], inflam = dd$inflam,
             drivers = dd$drivers[ids], phenotype = dd$phenotype[patient_id %in% ids], gcsf = dd$gcsf,
             base_sd = dd$base_sd)
  c2 <- cfg
  c2$fit$de <- utils::modifyList(cfg$fit$de, list(np = 16L, chunk = 5L, window = 5L, max_gen = 10L, polish_maxit = 5L))
  run <- file.path(tempdir(), "m1_round0_test")
  unlink(run, recursive = TRUE)
  de_init(run, "M0", d3, c2)
  expect_identical(readLines(file.path(run, "status")), "round0 0")
  for (i in 1:3) de_estep_task(run, 0L, i, c2)
  expect_identical(de_round0_task(run, c2), "estep")
  expect_identical(readLines(file.path(run, "status")), "estep 1")
  r0 <- readRDS(file.path(run, "round0.rds"))
  p1 <- readRDS(file.path(run, "pop_001.rds"))$pop
  pdats <- build_pdats(d3, c2)
  ref <- m_step_re13(pdats, r0$fits, readRDS(file.path(run, "pop_000.rds"))$pop, model_spec("M0")$re, c2)
  expect_equal(p1$mu, ref$mu)
  expect_equal(p1$struct, readRDS(file.path(run, "pop_000.rds"))$pop$struct)
  # round 1 is seeded by the round-0 MAPs
  f1 <- de_estep_task(run, 1L, 1L, c2)
  expect_true(is.finite(f1$nll))
  unlink(run, recursive = TRUE)
})

test_that("the calibrated CAR T constants in config.R match the calibrations (2%)", {
  cc <- calibrate_cart(cfg)
  expect_lt(abs(cc[["kappa"]] / st[["kappa"]] - 1), 0.02)
  expect_lt(abs(cc[["rho"]] / exp(cfg$model$re_prior$mu[["rho"]]) - 1), 0.02)
  zz <- calibrate_z(cfg)
  expect_lt(abs(zz[["z50"]] / st[["z50"]] - 1), 0.02)
  expect_lt(abs(zz[["th_zs"]] / st[["th_zs"]] - 1), 0.02)
})

# ---- Model3 --------------------------------------------------------------------------------------

# (the unused neutrophil compartments are constant but enter the step-size control: differences
# of ~1e-9 at rtol 1e-9 come from the adaptive step size)
test_that("Model3 with the Model2 values of N and P2 and kappa_c = th_aM = 0 reproduces the frozen Model1 v1.0 engine", {
  e1 <- new.env()
  Rcpp::sourceCpp(file.path(cfg$paths$code, "tests", "ref", "model1_v10.cpp"), cacheDir = cfg$paths$rcpp_cache, env = e1)
  s0 <- m2_values(st)
  s0[c("kappa_c", "th_am")] <- 0
  for (id in c("Mab", "MabC_ddr", "Mabx_Zm")) {
    sp <- model_spec(id)
    s2 <- s0
    if (sp$clone != "none") s2[names(cfg$model$clone_types[[sp$clone]])] <- cfg$model$clone_types[[sp$clone]]
    if (sp$q2 == "Mabx") s2[["eta_syn"]] <- 2
    p2 <- build_par(c(ind, f0 = 0.2, rho = 0.9), s2, sp)
    seg <- make_segments(des, s2, 90, cfg)
    a <- m3_sim(p2, seg, drv$v, drv$t0, drv$dt, days, 1e-9, 1e-11, 0.1)
    b <- e1$m1_sim(p2[as.character(e1$m1_param_names())], seg, drv$v, drv$t0, drv$dt, days, 1e-9, 1e-11, 0.1)
    expect_lt(max(abs(a[, colnames(b)] - b), na.rm = TRUE), 1e-6)
  }
})

test_that("Model3 with the Model2 values of N and P2 reproduces the frozen Model2 v1.0 engine", {
  e2 <- new.env()
  Rcpp::sourceCpp(file.path(cfg$paths$code, "tests", "ref", "model2_v10.cpp"), cacheDir = cfg$paths$rcpp_cache, env = e2)
  for (id in c("Mab", "MabC_ddr", "Mabx_Zm", "Mb")) {
    sp <- model_spec(id)
    s2 <- m2_values(default_struct(sp, cfg))
    if (sp$q2 == "Mabx") s2[["eta_syn"]] <- 2
    p3 <- build_par(c(ind, f0 = 0.2, rho = 0.9), s2, sp)
    seg <- make_segments(des, s2, 90, cfg)
    a <- m3_sim(p3, seg, drv$v, drv$t0, drv$dt, days, 1e-9, 1e-11, 0.1)
    b <- e2$m2_sim(p3[as.character(e2$m2_param_names())], seg, drv$v, drv$t0, drv$dt, days, 1e-9, 1e-11, 0.1)
    expect_lt(max(abs(a - b), na.rm = TRUE), 1e-6)
  }
  # the drop-one ids carry exactly these values
  expect_equal(default_struct(model_spec("Mb_noN"), cfg)[c("n_tr_n", "nk_n")], c(n_tr_n = 3, nk_n = 2))
  expect_equal(default_struct(model_spec("Mb_noP2"), cfg)[["d_zpl"]], 0)
  expect_true("d_zpl" %in% struct_est_names(model_spec("Mb"), cfg))
  expect_false("d_zpl" %in% struct_est_names(model_spec("Mb_noP2"), cfg))
  expect_false("d_zpl" %in% struct_est_names(model_spec("Mb_Zh"), cfg))
  expect_identical(model_spec("MabC_ddr_noN")$drop, "N")
})

test_that("N: the neutrophil chain changes only the neutrophil lineage; its storage pool delays the fall of ANC", {
  seg <- make_segments(des_nc, st, 90, cfg)
  z <- zero_driver(des$t0, cfg)
  run <- function(s) m3_sim(build_par(ind, s, model_spec("Mab")), seg, z$v, z$t0, z$dt, days, 1e-10, 1e-12, 0.1)
  a <- run(st)
  b <- run(m2_values(st))
  # (PN senses B_N, so only the neutrophil lineage and the mean M change)
  cols <- c("plt_rel", "amc_rel", "alc_rel", "S", "PMO", "PPL", "PL")
  expect_lt(max(abs(a[, cols] - b[, cols])), 1e-6)
  # with the same MTT, ANC holds longer during lymphodepletion (days -3..0) with the storage pool
  i <- days >= -3 & days <= 0
  expect_gt(mean(a[i, "anc_rel"]), mean(b[i, "anc_rel"]))
  # killing more compartments of the same chain deepens the nadir
  c3 <- run(replace(st, "nk_n", 2))
  expect_lt(min(a[days >= 0, "anc_rel"]), min(c3[days >= 0, "anc_rel"]))
})

test_that("P2: d_zPl acts only on the platelet MPP component and only through the CAR T - target signal", {
  sp <- model_spec("Mab")
  seg <- make_segments(des, st, 90, cfg)
  run <- function(s, d = drv) m3_sim(build_par(c(ind, rho = 0.9), s, sp), seg, d$v, d$t0, d$dt, days, 1e-10, 1e-12, 0.1)
  a <- run(st)
  b <- run(replace(st, "d_zpl", 2))
  cols <- c("anc_rel", "amc_rel", "alc_rel", "S", "PN", "PMO", "PL")
  expect_lt(max(abs(a[, cols] - b[, cols])), 1e-6)
  expect_lt(min(b[days > 0, "PPL"]), min(a[days > 0, "PPL"]))
  # without the MPP route e_z,Pl is 0 whatever d_zpl
  expect_equal(build_par(c(ind, rho = 0.9), replace(st, "d_zpl", 2), model_spec("Mab_Zh"))[["e_z_pl"]], 0)
})

test_that("th_aM acts only on the monocyte lineage; kappa_c deepens the nadir of every lineage", {
  seg <- make_segments(des_nc, st, 90, cfg)
  z <- zero_driver(des$t0, cfg)
  run <- function(s) m3_sim(build_par(ind, s, model_spec("Mab")), seg, z$v, z$t0, z$dt, days, 1e-10, 1e-12, 0.1)
  a <- run(st)
  b <- run(replace(st, "th_am", 0))
  # (the monocyte MPP component senses B_Mo, so only the monocyte lineage and the mean M change)
  # (differences of ~1e-8 come from the adaptive step size, which also sees the monocyte states)
  cols <- c("anc_rel", "plt_rel", "alc_rel", "S", "PN", "PPL", "PL")
  expect_lt(max(abs(a[, cols] - b[, cols])), 1e-6)
  expect_lt(min(a[, "amc_rel"]), min(b[, "amc_rel"]))
  c0 <- run(replace(st, "kappa_c", 0))
  for (v in c("anc_rel", "plt_rel", "amc_rel")) expect_lt(min(a[days >= 0, v]), min(c0[days >= 0, v]))
})
