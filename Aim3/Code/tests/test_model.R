# =============================================================================
# test_model.R
# Purpose : Unit tests of the marrow model, the nesting switches, the
#           integrator (against deSolve), the injury read-outs and the
#           phenotype rule. Called by step W1; every test returns a row.
# Inputs  : config, a design patient (pdes).
# Outputs : data.table(test, passed, detail)
# =============================================================================

run_unit_tests <- function(config, pd) {
  st <- config$model$struct
  rows <- list()
  add <- function(test, passed, detail) rows[[length(rows) + 1]] <<-
    data.table::data.table(test = test, passed = isTRUE(passed), detail = detail)
  days <- seq(-10, 60, by = 1)

  # 1. steady state without any input
  pd0 <- pd
  pd0$seg[, "dose"] <- 0
  pd0$drv$v[] <- 0
  m <- sim_raw(build_par(c(th_ap = 4, th_as = 0.2, th_bs = 0.1, r_s = 0.05), st, "Mab"), pd0, days, config)
  dev <- max(abs(m[, c("anc_rel", "plt_rel", "amc_rel", "S", "M")] - 1))
  add("steady state without input", dev < 1e-8, sprintf("max deviation %.2e", dev))

  # 2. nesting: Mab with no HSC injury equals M0
  a <- sim_raw(build_par(c(th_ap = 4, th_as = 0, th_bs = 0, r_s = 0.05), st, "Mab"), pd, days, config)
  b <- sim_raw(build_par(c(th_ap = 4, th_as = 0.2, th_bs = 0.1, r_s = 0.05), st, "M0"), pd, days, config)
  dev <- max(abs(a[, -1] - b[, -1]))
  add("M0 = Mab without HSC injury", dev < 1e-10, sprintf("max deviation %.2e", dev))

  # 3. switches
  p_ma <- build_par(c(th_as = 0.2, th_bs = 0.1), replace(st, "eta_syn", 2), "Ma")
  p_mb <- build_par(c(th_as = 0.2, th_bs = 0.1), replace(st, "eta_syn", 2), "Mb")
  p_ab <- build_par(c(th_as = 0.2, th_bs = 0.1), replace(st, "eta_syn", 2), "Mab")
  p_abx <- build_par(c(th_as = 0.2, th_bs = 0.1), replace(st, "eta_syn", 2), "Mabx")
  ok <- p_ma[["th_bs"]] == 0 && p_ma[["th_as"]] > 0 && p_mb[["th_as"]] == 0 && p_mb[["th_bs"]] > 0 &&
    p_ab[["eta_syn"]] == 0 && p_abx[["eta_syn"]] == 2
  add("model switches", ok, "Ma: th_bS = 0; Mb: th_aS = 0; Mab: eta = 0; Mabx keeps eta")

  # 4. integrator against deSolve (LSODA) for a patient with both injuries
  par <- build_par(c(th_ap = 4, th_as = 0.2, th_bs = 0.08, r_s = 0.05), replace(st, "eta_syn", 1), "Mabx")
  seg <- pd$seg
  f <- function(t, y, p) {
    i <- max(which(seg[, "start"] <= t + 1e-12))
    list(hsc_rhs(t, y, par, seg[i, "dose"], pd$drv$v, pd$drv$t0, pd$drv$dt))
  }
  y0 <- rep(1, length(hsc_state_names()))
  y0[3] <- 0
  tt <- sort(unique(c(seq(pd$t0, 60, by = 0.25), seg[, "start"], seg[, "end"])))
  tt <- tt[tt <= 60]
  ds <- deSolve::lsoda(y0, tt, f, NULL, rtol = 1e-8, atol = 1e-10, hmax = 0.1, tcrit = 60)
  ours <- sim_raw(par, pd, days, config)
  ref <- ds[match(days, ds[, 1]), 1 + c(8, 13, 18)]
  rel <- max(abs(ours[, c("anc_rel", "plt_rel", "amc_rel")] - ref) / pmax(abs(ref), 1e-3))
  add("RK45 agrees with deSolve", rel < 2e-3, sprintf("max relative difference %.2e", rel))

  # 5. monotonicity of the injury terms
  anc_at <- function(ind, day) sim_daily(build_par(ind, st, "Mab"), pd, config, days = day)$anc
  a1 <- anc_at(c(th_ap = 4, th_as = 0.1, th_bs = 0, r_s = 0.05), 14)
  a2 <- anc_at(c(th_ap = 4, th_as = 0.3, th_bs = 0, r_s = 0.05), 14)
  b1 <- anc_at(c(th_ap = 4, th_as = 0, th_bs = 0.05, r_s = 0.05), 25)
  b2 <- anc_at(c(th_ap = 4, th_as = 0, th_bs = 0.2, r_s = 0.05), 25)
  add("more injury -> lower ANC", a2 < a1 && b2 < b1,
      sprintf("a: %.2f -> %.2f (day 14); b: %.2f -> %.2f (day 25)", a1, a2, b1, b2))

  # 6. kill fraction read-out: f_a = 1 - exp(-th_aS * integral C)
  p6 <- build_par(c(th_ap = 4, th_as = 0.2, th_bs = 0, r_s = 0.05), st, "Ma")
  kf <- kill_fractions(p6, pd, config)
  dose_int <- sum((seg[, "end"] - seg[, "start"]) * seg[, "dose"])
  expect <- 1 - exp(-0.2 * dose_int / st[["k_e"]])
  add("f_a matches the exposure integral", abs(kf[["f_a"]] - expect) < 0.01,
      sprintf("f_a %.3f vs analytic %.3f", kf[["f_a"]], expect))

  # 7. synergy raises the CAR T-phase stem cell loss when chemotherapy has hit first
  p7 <- build_par(c(th_ap = 4, th_as = 0.2, th_bs = 0.05, r_s = 0.05), replace(st, "eta_syn", 3), "Mabx")
  p7b <- replace(p7, "eta_syn", 0)
  add("synergy increases f_b", kill_fractions(p7, pd, config)[["f_b"]] > kill_fractions(p7b, pd, config)[["f_b"]],
      sprintf("f_b %.3f (eta 3) vs %.3f (eta 0)", kill_fractions(p7, pd, config)[["f_b"]],
              kill_fractions(p7b, pd, config)[["f_b"]]))

  # 8. phenotype rule on constructed series
  r <- pheno_rule(config)
  d0 <- 0:30
  quick <- c(rep(0.1, 10), rep(2, 21))
  wave <- c(rep(0.1, 10), rep(2, 8), rep(0.6, 5), rep(2, 8))
  late <- c(rep(0.1, 25), rep(1.2, 6))
  spike <- c(rep(0.1, 20), 0.5, 3.0, 0.6, rep(0.8, 8))
  ok <- classify_anc(d0, quick, r) == "quick" && classify_anc(d0, wave, r) == "second_wave" &&
    classify_anc(d0, late, r) == "delayed" && classify_anc(d0, spike, r) == "delayed"
  add("phenotype rule on constructed series", ok, "quick / second wave / delayed / G-CSF spike ignored")

  data.table::rbindlist(rows)
}
