# =============================================================================
# model.R
# Purpose : R side of the marrow model: compile/load the C++ code, parameter
#           bookkeeping and the nested-model switches, chemotherapy dosing
#           segments (a), the CAR T-phase driver from ferritin (b), and
#           simulation wrappers.
# Inputs  : config$model, patient rows (io_read.R), inflammation series.
# Outputs : parameter vectors in C++ order; per-patient design objects ("pdes");
#           simulated trajectories.
#
# Nested models (HSC-level injury only; progenitor-level effects are always on
# and individual: th_aP for chemotherapy, q_b (suppression) for the CAR T phase).
# The accelerated-release term eps is off (0): with the ferritin driver it could
# not produce second waves; it is kept for steroid / G-CSF inputs once dosing
# dates are available. The model describes the drug-free trajectory.
#   M0   th_aS = th_bS = 0            progenitors only
#   Ma   th_bS = 0                    chemotherapy hits the HSC pool
#   Mb   th_aS = 0                    CAR T phase hits the HSC pool
#   Mab  both, eta_syn = 0            additive
#   Mabx both, eta_syn estimated      synergy (cycling HSC more sensitive to b)
# =============================================================================

load_hsc_model <- function(config) {
  if (exists("hsc_sim", mode = "function")) return(invisible(TRUE))
  dir.create(config$paths$rcpp_cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(file.path(config$paths$code, "src", "hsc_model.cpp"),
                  cacheDir = config$paths$rcpp_cache, env = globalenv())
  invisible(TRUE)
}

re_names_all <- c("th_ap", "th_as", "th_bs", "r_s", "q_b")

model_re_names <- function(model) {
  switch(model,
         M0 = c("th_ap", "q_b"),
         Ma = c("th_ap", "th_as", "r_s", "q_b"),
         Mb = c("th_ap", "th_bs", "r_s", "q_b"),
         Mab = re_names_all, Mabx = re_names_all,
         stop("unknown model: ", model))
}

# Structural parameters on the C++ scale (transit rates from mean transit times).
struct_to_cpp <- function(st) {
  c(alpha = st[["alpha"]], k_m = st[["k_m"]], phi = st[["phi"]], gmax = st[["gmax"]],
    gamma_m = st[["gamma_m"]],
    ktr_n = 4 / st[["mtt_n"]], ktr_pl = 4 / st[["mtt_pl"]], ktr_mo = 4 / st[["mtt_mo"]],
    c_n = st[["c_n"]], c_pl = st[["c_pl"]], c_mo = st[["c_mo"]],
    gam_n = st[["gam_n"]], gam_pl = st[["gam_pl"]], gam_mo = st[["gam_mo"]],
    k_e = st[["k_e"]], q_b = st[["q_b"]], eta_syn = st[["eta_syn"]], rho_m = st[["rho_m"]],
    eps = st[["eps"]], w_pl = st[["w_pl"]], w_mo = st[["w_mo"]])
}

# Full parameter vector in C++ order with the nesting switches applied.
# indiv: named natural-scale values for any of th_ap, th_as, th_bs, r_s
build_par <- function(indiv, struct, model) {
  nm <- hsc_param_names()
  par <- stats::setNames(numeric(length(nm)), nm)
  cp <- struct_to_cpp(struct)
  par[names(cp)] <- cp
  par[c("th_ap", "th_as", "th_bs", "r_s")] <- c(1, 0, 0, 0.05)
  par[names(indiv)] <- indiv
  if (model %in% c("M0", "Mb")) par["th_as"] <- 0
  if (model %in% c("M0", "Ma")) par["th_bs"] <- 0
  if (model != "Mabx") par["eta_syn"] <- 0
  par[nm]
}

# ---- chemotherapy dosing (a) -----------------------------------------------------------------------
# LD flu-cy on days -5, -4, -3 -> constant rate on [-5, -2), scaled by renal
# function; bridging chemotherapy as a low constant rate from apheresis + 2 to day -7.
patient_dosing <- function(prow, struct, config) {
  ld_rate <- (config$data$egfr_ref / prow$egfr)^struct[["beta_renal"]]
  br <- config$data$bridging
  br_rate <- switch(prow$bridge_regimen, rgo = struct[["u_br_rgo"]], oral_cy = struct[["u_br_cy"]], 0)
  br_start <- prow$apheresis_day + br$start_after_apheresis
  list(ld = c(start = -5, end = -2, rate = ld_rate),
       br = if (br_rate > 0) c(start = br_start, end = br$end_day, rate = br_rate) else NULL)
}

make_segments <- function(t0, t_end, dosing) {
  iv <- rbind(dosing$ld, dosing$br)
  bp <- sort(unique(c(t0, t_end, pmin(pmax(c(iv[, "start"], iv[, "end"]), t0), t_end))))
  seg <- cbind(start = bp[-length(bp)], end = bp[-1], dose = 0)
  mid <- (seg[, "start"] + seg[, "end"]) / 2
  for (i in seq_len(nrow(iv))) {
    on <- mid >= iv[i, "start"] & mid < iv[i, "end"]
    seg[on, "dose"] <- seg[on, "dose"] + iv[i, "rate"]
  }
  seg
}

# ---- CAR T-phase driver (b) --------------------------------------------------------------------------
# D(t) = max(0, s(t - delta_b)); s = smoothed log(ferritin / baseline) over the
# observed post-infusion days, exponential decay after the last observation,
# 0 before infusion.
ferritin_signal <- function(days, ferritin, base) {
  k <- !is.na(ferritin) & days >= 0
  x <- days[k]
  y <- log(ferritin[k] / base)
  if (length(x) < 4) return(list(fun = function(t) rep(0, length(t)), last = 0))
  df <- min(8, max(3, length(x) / 3))
  sp <- stats::smooth.spline(x, y, df = df)
  last <- max(x)
  list(fun = function(t) stats::predict(sp, t)$y, last = last)
}

make_driver <- function(sig, struct, config, t0 = config$model$t_start, t_end = config$model$t_end) {
  dt <- config$model$drv_dt
  grid <- seq(t0, t_end, by = dt)
  tt <- grid - struct[["delta_b"]]
  v <- numeric(length(grid))
  on <- tt >= 0 & tt <= sig$last
  if (any(on)) v[on] <- sig$fun(tt[on])
  after <- tt > sig$last
  if (any(after)) v[after] <- sig$fun(sig$last) * exp(-struct[["k_dec"]] * (tt[after] - sig$last))
  list(v = pmax(v, 0), t0 = t0, dt = dt)
}

# Typical driver: pointwise median of the patients' drivers (design steps).
median_driver <- function(drivers) {
  v <- apply(do.call(cbind, lapply(drivers, `[[`, "v")), 1, stats::median)
  list(v = v, t0 = drivers[[1]]$t0, dt = drivers[[1]]$dt)
}

# ---- per-patient design object -------------------------------------------------------------------
# Everything needed to simulate a patient except the individual parameters.
make_pdes <- function(prow, inflam_p, struct, config) {
  dosing <- patient_dosing(prow, struct, config)
  t0 <- min(config$model$t_start, if (is.null(dosing$br)) Inf else dosing$br[["start"]] - 1)
  sig <- ferritin_signal(inflam_p$day, inflam_p$ferritin, prow$ferr_base)
  list(id = prow$patient_id, t0 = t0,
       seg = make_segments(t0, config$model$t_end, dosing),
       drv = make_driver(sig, struct, config, t0 = t0),
       sig = sig,
       set = c(anc = prow$anc_set, plt = prow$plt_set, amc = prow$amc_set),
       x = list(th_ap = numeric(0), th_as = numeric(0), th_bs = numeric(0), q_b = numeric(0),
                r_s = c(prior_lines_c = prow$prior_lines_c, prior_asct = prow$prior_asct)))
}

# Rebuild the parts of a pdes that depend on structural parameters (dosing rates, driver lag).
refresh_pdes <- function(pd, prow, struct, config) {
  dosing <- patient_dosing(prow, struct, config)
  pd$seg <- make_segments(pd$t0, config$model$t_end, dosing)
  pd$drv <- make_driver(pd$sig, struct, config, t0 = pd$t0)
  pd
}

# ---- simulation ------------------------------------------------------------------------------------
sim_raw <- function(par, pd, times, config) {
  s <- config$model$solver
  hsc_sim(par, pd$seg, pd$drv$v, pd$drv$t0, pd$drv$dt, times, s$rtol, s$atol, s$hmax)
}

# Daily trajectory on the absolute scale.
sim_daily <- function(par, pd, config, days = seq(-10, 60, by = 1)) {
  m <- sim_raw(par, pd, days, config)
  data.table::data.table(day = m[, "t"], anc = m[, "anc_rel"] * pd$set[["anc"]],
                         plt = m[, "plt_rel"] * pd$set[["plt"]], amc = m[, "amc_rel"] * pd$set[["amc"]],
                         S = m[, "S"], M = m[, "M"], C = m[, "C"], D = m[, "D"])
}

# Phenotype of the latent (treatment-free) trajectory.
sim_phenotype <- function(par, pd, config, variant = NULL) {
  days <- 0:config$phenotype$primary$window_end
  m <- sim_raw(par, pd, days, config)
  if (anyNA(m[, 2])) return(NA_character_)
  classify_sim(days, m[, "anc_rel"] * pd$set[["anc"]], config, variant)
}

# Direct HSC kill fractions in each phase (continuous read-outs f_a, f_b).
kill_fractions <- function(par, pd, config) {
  days <- seq(pd$t0, config$model$t_end, by = 0.25)
  m <- sim_raw(par, pd, days, config)
  if (anyNA(m[, 2])) return(c(f_a = NA_real_, f_b = NA_real_))
  dt <- 0.25
  ka <- par[["th_as"]] * m[, "C"]
  kb <- par[["th_bs"]] * m[, "D"] * (1 + par[["eta_syn"]] * pmax(1 - m[, "S"], 0))
  c(f_a = 1 - exp(-sum(ka) * dt), f_b = 1 - exp(-sum(kb) * dt))
}

# Generic design patient: median set points, no bridging, eGFR = reference, typical driver.
design_pdes <- function(struct, config, drv, set = c(anc = 3.4, plt = 180, amc = 0.6)) {
  prow <- list(patient_id = "design", egfr = config$data$egfr_ref, bridge_regimen = "none",
               apheresis_day = -27)
  dosing <- patient_dosing(prow, struct, config)
  t0 <- config$model$t_start
  list(id = "design", t0 = t0, seg = make_segments(t0, config$model$t_end, dosing),
       drv = drv, set = set, x = NULL)
}
