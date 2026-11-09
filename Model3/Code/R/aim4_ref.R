# =============================================================================
# aim4_ref.R
# Purpose : The Aim4 v1.3 lineage-tree model as a frozen reference inside Model1
#           (engine copied to tests/ref/aim4_v13.cpp, hashes in
#           SOURCE_aim4_v13_sha256.txt). Used for the K1 calibration targets
#           (typical ANC nadir, ALC course) and to generate the S4 cohorts
#           (Aim4 truth, Model1 fit). Nothing here is fitted.
# Inputs  : config (paths$aim4_ref), individual values, structural values.
# Outputs : Aim4 v1.3 parameter vectors and trajectories.
# =============================================================================

# Aim4 v1.3 structural defaults (Aim4/Code/config.R, 2026-11-04 09:43)
aim4_v13_struct <- c(
  alpha = 1, gmax = 3, gamma_m = 0.5,
  k_mpp = 0.3, sr_mpp = 0.21,
  k_cmp = 2, sr_cmp = 0.75, k_gmp = 2, sr_gmp = 0.75, k_mep = 2, sr_mep = 0.75,
  k_clp = 1, sr_clp = 0.5, k_ly = 0.05, sr_ly = 0.04,
  fb_mpp = 1, fb_cmp = 1, fb_gmp = 1, fb_mep = 1, fb_clp = 1, fb_ly = 1,
  r_a = 1, r_b = 1, delta = 1,
  mtt_n = 5, mtt_mo = 1.5, mtt_pl = 5,
  c_n = 2.3, c_mo = 0.5, c_pl = 0.12, c_l = 2,
  gam_n = 0.2, gam_mo = 0.2, gam_pl = 0.2,
  k_i = 0.2, emax = 0.9, ic50 = 0.7, hill = 2, e_mo = 0.7, e_pl = 0.7, e_l = 0.7, eta_syn = 0,
  e_sr = 0.5, th_al = 1.05,
  k_e = 0.8, beta_renal = 0.7, u_br_rgo = 0.05, u_br_cy = 0.01,
  gcsf_fold = 2.5, gcsf_rel = 2, k_g = 1, tpo_fold = 1.8,
  d_n = 4, d_mo = 2, d_pl = 0,
  car_t0 = 0.001, k_ag = 0.1, kappa = 10.65, k_t = 0.2, beta_a = 1, k_z = 0.5, z50 = 0.0140,
  th_zs = 2.64, e_z = 0.5)
aim4_v13_th_ap_centre <- 2.31

load_aim4_reference <- function(config) {
  if (exists(".aim4ref", envir = globalenv())) return(invisible(get(".aim4ref", envir = globalenv())))
  e <- new.env()
  dir.create(config$paths$rcpp_cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(config$paths$aim4_ref, cacheDir = config$paths$rcpp_cache, env = e)
  assign(".aim4ref", e, envir = globalenv())
  invisible(e)
}

# Aim4 v1.3 build_par (Aim4/Code/R/model.R) for a Model1 model spec (same ids and nesting)
aim4_build_par <- function(indiv, struct, spec) {
  s <- struct
  on <- function(k) if (k %in% spec$re) indiv[[k]] else 0
  c(r_s = indiv[["r_s"]], f0 = on("f0"), th_as = on("th_as"), th_ap = indiv[["th_ap"]], th_bs = on("th_bs"),
    alpha = s[["alpha"]],
    k_mpp = s[["k_mpp"]], sr_mpp = s[["sr_mpp"]], k_cmp = s[["k_cmp"]], sr_cmp = s[["sr_cmp"]],
    k_gmp = s[["k_gmp"]], sr_gmp = s[["sr_gmp"]], k_mep = s[["k_mep"]], sr_mep = s[["sr_mep"]],
    k_clp = s[["k_clp"]], sr_clp = s[["sr_clp"]], k_ly = s[["k_ly"]], sr_ly = s[["sr_ly"]],
    gmax = s[["gmax"]], gamma_m = s[["gamma_m"]],
    fb_mpp = s[["fb_mpp"]], fb_cmp = s[["fb_cmp"]], fb_gmp = s[["fb_gmp"]], fb_mep = s[["fb_mep"]],
    fb_clp = s[["fb_clp"]], fb_ly = s[["fb_ly"]],
    r_a = s[["r_a"]], r_b = s[["r_b"]], delta = s[["delta"]],
    ktr_n = 3 / s[["mtt_n"]], ktr_mo = 3 / s[["mtt_mo"]], ktr_pl = 3 / s[["mtt_pl"]],
    c_n = s[["c_n"]], c_mo = s[["c_mo"]], c_pl = s[["c_pl"]], c_l = s[["c_l"]],
    gam_n = s[["gam_n"]], gam_mo = s[["gam_mo"]], gam_pl = s[["gam_pl"]],
    k_i = s[["k_i"]], emax = s[["emax"]], ic50 = s[["ic50"]], hill = s[["hill"]],
    e_mo = s[["e_mo"]], e_pl = s[["e_pl"]], e_l = s[["e_l"]],
    eta_syn = if (spec$q2 == "Mabx") s[["eta_syn"]] else 0,
    e_sr = s[["e_sr"]], th_al = s[["th_al"]],
    k_e = s[["k_e"]], gcsf_fold = s[["gcsf_fold"]], tpo_fold = s[["tpo_fold"]],
    d_n = s[["d_n"]], d_mo = s[["d_mo"]], d_pl = s[["d_pl"]],
    rho = if ("rho" %in% names(indiv)) indiv[["rho"]] else 0,
    a0 = if ("a0" %in% names(indiv)) indiv[["a0"]] else 1,
    k_t = s[["k_t"]], k_ag = s[["k_ag"]], kappa = s[["kappa"]], k_z = s[["k_z"]], z50 = s[["z50"]],
    th_zs = if (isTRUE(spec$z_hsc)) s[["th_zs"]] else 0,
    e_z = if (isTRUE(spec$z_mpp)) s[["e_z"]] else 0,
    k_g = s[["k_g"]], gcsf_rel = s[["gcsf_rel"]])
}

# Aim4 v1.3 trajectory matrix (same segment and driver conventions as Model1)
aim4_sim <- function(indiv, struct, spec, seg, drv, times, config) {
  e <- load_aim4_reference(config)
  sv <- config$model$solver
  if (spec$clone != "none") struct[names(config$model$clone_types[[spec$clone]])] <- config$model$clone_types[[spec$clone]]
  e$a4_sim(aim4_build_par(indiv, struct, spec), seg, drv$v, drv$t0, drv$dt, times, sv$rtol, sv$atol, sv$hmax)
}

# Daily Aim4 v1.3 trajectory on the observed scale (as sim_daily; ALC = lymphocytes + CAR T)
aim4_daily <- function(indiv, struct, spec, design, drv, config, days, base) {
  seg <- make_segments(design, struct, max(days), config)
  out <- data.table::as.data.table(aim4_sim(indiv, struct, spec, seg, drv, days, config))
  data.table::setnames(out, "t", "day")
  out[, `:=`(anc = anc_rel * base$anc, plt = plt_rel * base$plt, amc = amc_rel * base$amc,
             alc = alc_rel * base$alc + cart)]
  out
}
