# =============================================================================
# an_K2_phasemap.R
# Purpose : Step K2 - fingerprints and phenotype maps (no patient data).
#   1. (f_a, f_b) map: phenotype over HSC loss from chemotherapy x CAR T phase,
#      for eta = 0 and eta > 0 (synergy)
#   2. reserve x clone fraction map (Aim 2's question), DDR- and DTA-type clones
#   3. the link: the same a/b injuries with no clone, a DDR clone and a DTA clone
#   4. pseudo second wave: phenotype vs HSC loss from a single chemotherapy hit
#   5. fingerprints: the same 70% HSC loss placed in a vs b (ANC, PLT, AMC)
# Inputs  : ctx (config).
# Outputs : list(tables, maps ...)
# =============================================================================

k2_base <- list(anc = 3.5, plt = 200, amc = 0.6)

k2_run <- function(ind, spec, struct, design, drv, cfg) {
  m <- sim_daily(ind, struct, spec, design, drv, cfg, days = seq(-10, 90), base = k2_base)
  tt <- m$day
  list(ph = phenotype_sim(m$day, m$anc, cfg, with_persistent = TRUE),
       f_a = 1 - exp(-ind[["th_as"]] * trapz(tt, m$C)),
       f_b = 1 - exp(-ind[["th_bs"]] * trapz(tt[tt >= 0], m$I[tt >= 0])),
       S30 = m$S[m$day == 30])
}

an_K2_phasemap <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$K2
  des <- typical_design(cfg)
  drv <- driver_from_kernel(des$t0, 1.2, 8, cfg)
  n <- s$n_grid
  lg <- function(r) exp(seq(log(r[1]), log(r[2]), length.out = n))
  ap <- 1.5 * exp(cfg$model$re_prior$mu[["th_ap"]])   # progenitor hit: 1.5 x the prior centre

  # 1. (f_a, f_b) map
  g1 <- data.table::CJ(th_as = lg(s$th_as_range), th_bs = lg(s$th_bs_range), eta = s$eta_levels)
  map_ab <- g1[, {
    sp <- model_spec(if (eta > 0) "Mabx" else "Mab")
    st <- cfg$model$struct
    st[["eta_syn"]] <- eta
    r <- k2_run(c(r_s = 0.08, f0 = 0, th_as = th_as, th_ap = ap, th_bs = th_bs), sp, st, des, drv, cfg)
    .(phenotype = r$ph, f_a = r$f_a, f_b = r$f_b, S30 = r$S30)
  }, by = .(th_as, th_bs, eta)]

  # 2. reserve x clone fraction
  g2 <- data.table::CJ(r_s = lg(s$r_s_range), f0 = seq(s$f0_range[1], s$f0_range[2], length.out = n),
                       clone = c("ddr", "dta"))
  map_rf <- g2[, {
    sp <- model_spec(model_id("Mab", clone))
    r <- k2_run(c(r_s = r_s, f0 = f0, th_as = 0.3, th_ap = ap, th_bs = 0.03), sp,
                default_struct(sp, cfg), des, drv, cfg)
    .(phenotype = r$ph, S30 = r$S30)
  }, by = .(r_s, f0, clone)]

  # 3. the link: the a/b map under no clone / DDR clone / DTA clone (f0 = 0.3)
  g3 <- data.table::CJ(th_as = lg(s$th_as_range), th_bs = lg(s$th_bs_range), clone = c("none", "ddr", "dta"))
  map_link <- g3[, {
    sp <- model_spec(model_id("Mab", clone))
    f0 <- if (clone == "none") 0 else 0.3
    r <- k2_run(c(r_s = 0.08, f0 = f0, th_as = th_as, th_ap = ap, th_bs = th_bs), sp,
                default_struct(sp, cfg), des, drv, cfg)
    .(phenotype = r$ph, f_a = r$f_a, f_b = r$f_b)
  }, by = .(th_as, th_bs, clone)]

  # 4. pseudo second wave: single chemotherapy hit, no CAR T phase, three self-renewal feedback gains
  g4 <- data.table::CJ(th_as = lg(c(0.01, 3)), gamma_m = c(0.5, 1, 1.5))
  pseudo <- g4[, {
    st <- cfg$model$struct
    st[["gamma_m"]] <- gamma_m
    r <- k2_run(c(r_s = 0.05, f0 = 0, th_as = th_as, th_ap = ap, th_bs = 0), model_spec("Ma"), st, des,
                zero_driver(des$t0, cfg), cfg)
    .(phenotype = r$ph, hsc_loss = r$f_a)
  }, by = .(th_as, gamma_m)]
  thr <- pseudo[phenotype == "second_wave", .(threshold_hsc_loss = min(hsc_loss)), by = gamma_m]

  # 5. fingerprints: 70% HSC loss from a vs from b
  fp <- local({
    st <- cfg$model$struct
    ia <- uniroot(function(x) k2_run(c(r_s = 0.08, f0 = 0, th_as = x, th_ap = ap, th_bs = 0), model_spec("Ma"),
                                     st, des, drv, cfg)$f_a - 0.7, c(0.01, 5))$root
    ib <- uniroot(function(x) k2_run(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = x), model_spec("Mb"),
                                     st, des, drv, cfg)$f_b - 0.7, c(1e-4, 2))$root
    a <- sim_daily(c(r_s = 0.08, f0 = 0, th_as = ia, th_ap = ap, th_bs = 0), st, model_spec("Ma"), des, drv, cfg,
                   days = seq(-10, 60), base = k2_base)[, route := "a: chemotherapy"]
    b <- sim_daily(c(r_s = 0.08, f0 = 0, th_as = 0, th_ap = ap, th_bs = ib), st, model_spec("Mb"), des, drv, cfg,
                   days = seq(-10, 60), base = k2_base)[, route := "b: CAR T phase"]
    rbind(a, b)
  })
  ctx$log("  K2: pseudo second-wave threshold (HSC loss) by self-renewal feedback gain: ",
          paste(sprintf("gamma_M %.1f: %s", thr$gamma_m, round(thr$threshold_hsc_loss, 2)), collapse = "; "))
  list(tables = list(map_ab = map_ab, map_reserve_clone = map_rf, map_link = map_link,
                     pseudo_second_wave = pseudo, pseudo_threshold = thr),
       map_ab = map_ab, map_rf = map_rf, map_link = map_link, pseudo = pseudo, thr = thr, fingerprint = fp,
       driver = drv)
}
