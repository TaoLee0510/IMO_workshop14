# =============================================================================
# model.R
# Purpose : R side of the parent model and the Friberg reference model:
#           compile/load the C++ code, model specifications (nesting
#           switches), parameter vectors, piecewise inputs (segments) and the
#           simulation wrappers used by every step.
# Inputs  : config, individual parameters, patient design.
# Outputs : parameter vectors, segment matrices, simulated trajectories.
#
# Model ids
#   Q2 axis  M0 (no HSC injury) | Ma (chemo) | Mb (CAR T phase) | Mab | Mabx (synergy)
#   Q1 axis  suffix "C_ddr" / "C_dta" adds the clone layer with that clone type
#            e.g. "MabC_ddr"; "MF" is the Friberg (2002) reference model
# =============================================================================

load_aim4_model <- function(config) {
  if (exists("a4_sim", mode = "function")) return(invisible(TRUE))
  dir.create(config$paths$rcpp_cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(file.path(config$paths$code, "src", "aim4_model.cpp"),
                  cacheDir = config$paths$rcpp_cache, env = globalenv())
  invisible(TRUE)
}

# ---- model specifications ----------------------------------------------------------------
model_spec <- function(id) {
  if (id == "MF") return(list(id = id, family = "mf", q2 = NA_character_, clone = "none", route = "none",
                              re = c("slope", "mtt_scale")))
  # suffix _nc: clone model fitted without covariates on f0 (test T6)
  m <- regmatches(id, regexec("^(M0|Ma|Mb|Mab|Mabx)(C_(ddr|dta))?(_Z(h|m))?(_nc)?$", id))[[1]]
  if (!length(m)) stop("unknown model id: ", id)
  q2 <- m[2]
  clone <- if (nzchar(m[4])) m[4] else "none"
  # route of the CAR T - target signal Z: both (default), HSC only (_Zh), MPP only (_Zm);
  # the HSC route is part of the CAR T-phase injury (b models only)
  route <- switch(m[6], h = "hsc", m = "mpp", "both")
  re <- c("r_s", "th_ap")
  if (q2 %in% c("Ma", "Mab", "Mabx")) re <- c(re, "th_as")
  if (q2 %in% c("Mb", "Mab", "Mabx")) re <- c(re, "th_bs")
  if (clone != "none") re <- c(re, "f0")
  re <- c(re, "rho")
  list(id = id, family = "a4", q2 = q2, clone = clone, route = route, re = re,
       z_hsc = q2 %in% c("Mb", "Mab", "Mabx") && route != "mpp", z_mpp = route != "hsc",
       f0_nocov = nzchar(m[7]))
}

model_id <- function(q2, clone = "none", route = "both") {
  paste0(q2, if (clone != "none") paste0("C_", clone) else "",
         switch(route, hsc = "_Zh", mpp = "_Zm", ""))
}

model_re_names <- function(spec) spec$re

# Individual baselines (v1.3): one per lineage in the likelihood (MF has no ALC)
baseline_names <- function(spec, config) {
  v <- config$model$baseline$vars
  if (spec$family == "mf") v <- setdiff(v, "alc")
  paste0("b_", v)
}

# Structural parameters estimated for a model.
struct_est_names <- function(spec, config) {
  if (spec$family == "mf") return(config$fit$struct_est$MF)
  est <- config$fit$struct_est$base
  if (isTRUE(spec$z_hsc)) est <- c(est, config$fit$struct_est$route_hsc)
  if (isTRUE(spec$z_mpp)) est <- c(est, config$fit$struct_est$route_mpp)
  if (spec$q2 == "Mabx") est <- c(est, config$fit$struct_est$Mabx)
  est
}

# Default structural vector for a model (clone traits from the clone type).
default_struct <- function(spec, config) {
  if (spec$family == "mf") return(config$model$mf_struct)
  st <- config$model$struct
  if (spec$clone != "none") st[names(config$model$clone_types[[spec$clone]])] <- config$model$clone_types[[spec$clone]]
  if (spec$q2 == "Mabx") st[["eta_syn"]] <- config$fit$eta_syn_prior[["centre"]]
  st
}

# ---- working scale for structural parameters ----------------------------------------------
struct_to_work <- function(x, config) {
  u <- names(x) %in% config$model$unit_interval
  out <- x
  out[u] <- stats::qlogis(pmin(pmax(x[u], 1e-6), 1 - 1e-6))
  out[!u] <- log(pmax(x[!u], 1e-10))
  out
}
work_to_struct <- function(w, config) {
  u <- names(w) %in% config$model$unit_interval
  out <- w
  out[u] <- stats::plogis(w[u])
  out[!u] <- exp(w[!u])
  out
}

# ---- individual parameters -----------------------------------------------------------------
# eta (named, model REs) -> natural-scale individual values incl. switched-off ones
eta_to_indiv <- function(eta, spec) {
  if (spec$family == "mf") return(c(slope = exp(eta[["slope"]]), mtt_scale = exp(eta[["mtt_scale"]])))
  out <- c(r_s = exp(eta[["r_s"]]), f0 = 0, th_as = 0, th_ap = exp(eta[["th_ap"]]), th_bs = 0, rho = 0)
  if ("th_as" %in% names(eta)) out[["th_as"]] <- exp(eta[["th_as"]])
  if ("th_bs" %in% names(eta)) out[["th_bs"]] <- exp(eta[["th_bs"]])
  if ("f0" %in% names(eta)) out[["f0"]] <- stats::plogis(eta[["f0"]])
  if ("rho" %in% names(eta)) out[["rho"]] <- exp(eta[["rho"]])
  out
}

# Initial antigen load of a patient: a0 = exp(beta_a * ldh_c) (ldh_c: centred log LDH at LD)
patient_a0 <- function(struct, ldh_c) exp(struct[["beta_a"]] * (if (is.null(ldh_c) || is.na(ldh_c)) 0 else ldh_c))

# Full C++ parameter vector.
build_par <- function(indiv, struct, spec) {
  if (spec$family == "mf") {
    return(c(slope = indiv[["slope"]], mtt_scale = indiv[["mtt_scale"]],
             ktr_n = 4 / struct[["mtt_n"]], ktr_pl = 4 / struct[["mtt_pl"]], ktr_mo = 4 / struct[["mtt_mo"]],
             gam_n = struct[["gam_mf"]], gam_pl = struct[["gam_mf"]], gam_mo = struct[["gam_mf"]],
             slope_pl = struct[["slope_pl"]], slope_mo = struct[["slope_mo"]], k_e = struct[["k_e"]]))
  }
  s <- struct
  eta_syn <- if (spec$q2 == "Mabx") s[["eta_syn"]] else 0
  # nesting switches: terms a model does not contain are zero
  on <- function(k) if (k %in% spec$re) indiv[[k]] else 0
  prog <- c("mpp", "cmp", "gmp", "mep", "clp", "ly")
  bad <- prog[s[paste0("sr_", prog)] >= s[paste0("k_", prog)]]
  if (length(bad)) stop("self-renewal rate must be below the exit rate: ", paste(bad, collapse = ", "))
  # maturation chains have three compartments: ktr = 3 / MTT
  c(r_s = indiv[["r_s"]], f0 = on("f0"), th_as = on("th_as"), th_ap = indiv[["th_ap"]],
    th_bs = on("th_bs"),
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
    e_mo = s[["e_mo"]], e_pl = s[["e_pl"]], e_l = s[["e_l"]], eta_syn = eta_syn,
    e_sr = s[["e_sr"]],
    th_al = s[["th_al"]],
    k_e = s[["k_e"]], gcsf_fold = s[["gcsf_fold"]], tpo_fold = s[["tpo_fold"]],
    d_n = s[["d_n"]], d_mo = s[["d_mo"]], d_pl = s[["d_pl"]],
    # v1.3: CAR T - target interaction (Z routes switched by the model id) and G-CSF effect
    rho = if ("rho" %in% names(indiv)) indiv[["rho"]] else 0,
    a0 = if ("a0" %in% names(indiv)) indiv[["a0"]] else 1,
    k_t = s[["k_t"]], k_ag = s[["k_ag"]], kappa = s[["kappa"]], k_z = s[["k_z"]], z50 = s[["z50"]],
    th_zs = if (isTRUE(spec$z_hsc)) s[["th_zs"]] else 0,
    e_z = if (isTRUE(spec$z_mpp)) s[["e_z"]] else 0,
    k_g = s[["k_g"]], gcsf_rel = s[["gcsf_rel"]])
}

# ---- inputs ----------------------------------------------------------------------------------
# design: list(t0, egfr, bridging = "none" | "rgo" | "cy", apheresis_day, ld_scale,
#              gcsf = matrix(start, end) or NULL, tpo = matrix or NULL,
#              cd34 = data.frame(day, add) or NULL)
# Segments: start, end, dose, gcsf, tpo, s_add (s_add applied at the segment start).
make_segments <- function(design, struct, t_end, config) {
  bw <- config$data$bridging_window
  br_rate <- switch(design$bridging %||% "none", rgo = struct[["u_br_rgo"]], cy = struct[["u_br_cy"]], 0)
  br <- if (br_rate > 0) c(design$apheresis_day + bw[1], bw[2]) else NULL
  ld_rate <- (90 / design$egfr)^struct[["beta_renal"]] * (design$ld_scale %||% 1)
  gc <- design$gcsf
  tp <- design$tpo
  cd <- design$cd34
  bp <- c(design$t0, t_end, -5, -2, 0, br, if (!is.null(gc)) c(gc[, 1], gc[, 2]),
          if (!is.null(tp)) c(tp[, 1], tp[, 2]), if (!is.null(cd)) cd$day)
  bp <- sort(unique(bp[bp >= design$t0 & bp <= t_end]))
  st <- bp[-length(bp)]
  en <- bp[-1]
  mid <- (st + en) / 2
  in_win <- function(w, x) if (is.null(w) || !nrow(w)) rep(FALSE, length(x)) else
    vapply(x, function(z) any(z >= w[, 1] & z < w[, 2]), logical(1))
  dose <- ifelse(mid >= -5 & mid < -2, ld_rate, 0)
  if (!is.null(br)) dose <- dose + ifelse(mid >= br[1] & mid < br[2], br_rate, 0)
  s_add <- rep(0, length(st))
  if (!is.null(cd)) for (j in seq_len(nrow(cd))) s_add[st == cd$day[j]] <- s_add[st == cd$day[j]] + cd$add[j]
  # CAR T infusion on day 0 (design$no_car: counterfactual without the CAR T phase)
  car <- if (isTRUE(design$no_car) || is.na(struct["car_t0"])) rep(0, length(st)) else
    ifelse(st == 0, struct[["car_t0"]], 0)
  cbind(start = st, end = en, dose = dose, gcsf = as.numeric(in_win(gc, mid)),
        tpo = as.numeric(in_win(tp, mid)), s_add = s_add, car = car)
}

# ---- CAR T-phase driver -----------------------------------------------------------------------
drv_grid_t <- function(t0, config) seq(t0, config$model$drv_end, by = config$model$drv_dt)

# Gamma-shaped kernel on log(ferritin / baseline): peak A at day tp.
kernel_fun <- function(t, A, tp, h) ifelse(t > 0, A * (t / tp)^h * exp(h * (1 - t / tp)), 0)

driver_from_kernel <- function(t0, A, tp, config) {
  tt <- drv_grid_t(t0, config)
  list(t0 = t0, dt = config$model$drv_dt, v = kernel_fun(tt, A, tp, config$model$kernel$h))
}

# Smoothed log(ferritin / baseline) from observations (day, value); positive part,
# zero before day 0, exponential decay after the last observation.
driver_from_ferritin <- function(day, value, base, t0, config) {
  tt <- drv_grid_t(t0, config)
  w <- config$data$fer_smooth_window
  ok <- !is.na(value) & value > 0 & day >= w[1] - 3 & day <= w[2]
  v <- numeric(length(tt))
  if (sum(ok) >= 4) {
    x <- day[ok]
    y <- log(value[ok] / base)
    sp <- stats::smooth.spline(x, y, df = min(max(4, floor(sum(ok) / 3)), 8))
    last <- max(x)
    pos <- tt >= 0 & tt <= last
    v[pos] <- pmax(stats::predict(sp, tt[pos])$y, 0)
    after <- tt > last
    v_last <- max(stats::predict(sp, last)$y, 0)
    v[after] <- v_last * exp(-config$model$fer_decay * (tt[after] - last))
  }
  list(t0 = t0, dt = config$model$drv_dt, v = v)
}

zero_driver <- function(t0, config) {
  tt <- drv_grid_t(t0, config)
  list(t0 = t0, dt = config$model$drv_dt, v = numeric(length(tt)))
}

# Scale the driver after a given day (anti-inflammatory what-if).
scale_driver <- function(drv, after_day, factor) {
  tt <- drv$t0 + (seq_along(drv$v) - 1) * drv$dt
  drv$v[tt > after_day] <- drv$v[tt > after_day] * factor
  drv
}

# ---- simulation -------------------------------------------------------------------------------
# par: full vector (build_par); seg: make_segments; drv: driver list; times: ascending days
# track: also return the mean time (age_*) and number of divisions (gen_*) since
# MPP entry of the blood cells (parent model only)
sim_model <- function(par, seg, drv, times, spec, config, track = FALSE) {
  sv <- config$model$solver
  if (spec$family == "mf") return(mf_sim(par, seg, times, sv$rtol, sv$atol, sv$hmax))
  a4_sim(par, seg, drv$v, drv$t0, drv$dt, times, sv$rtol, sv$atol, sv$hmax, track = track)
}

# Daily trajectory as a data.table (counts on the observed scale when base given;
# ALC only when base$alc is given and the model has a lymphoid branch).
sim_daily <- function(indiv, struct, spec, design, drv, config, days = NULL, base = NULL, track = FALSE) {
  if (is.null(days)) days <- seq(ceiling(design$t0), config$model$sim_end)
  seg <- make_segments(design, struct, max(days), config)
  m <- sim_model(build_par(indiv, struct, spec), seg, drv, days, spec, config, track = track)
  out <- data.table::as.data.table(m)
  data.table::setnames(out, "t", "day")
  if (!is.null(base)) out[, `:=`(anc = anc_rel * base$anc, plt = plt_rel * base$plt, amc = amc_rel * base$amc)]
  # v1.3: observed ALC = lymphocytes (relative x baseline) + CAR T cells (absolute)
  if (!is.null(base$alc) && "alc_rel" %in% names(out)) out[, alc := alc_rel * base$alc + (if ("cart" %in% names(out)) cart else 0)]
  out
}

# Closed-form pre-treatment means of the time (days) and number of divisions since
# MPP entry, per compartment (independent of the C++ initialisation; K1 compares
# the two). Renewing compartment: age = a_in + 1/(k - sr), gen = g_in + inc + 2 sr/(k - sr).
steady_age_generation <- function(struct) {
  s <- struct
  ren <- function(k, sr, a_in, g_in, inc) c(age = a_in + 1 / (k - sr), gen = g_in + inc + 2 * sr / (k - sr))
  mpp <- ren(s[["k_mpp"]], s[["sr_mpp"]], 0, 0, 0)
  cmp <- ren(s[["k_cmp"]], s[["sr_cmp"]], mpp[["age"]], mpp[["gen"]], 1)
  gmp <- ren(s[["k_gmp"]], s[["sr_gmp"]], cmp[["age"]], cmp[["gen"]], 1)
  mep <- ren(s[["k_mep"]], s[["sr_mep"]], cmp[["age"]], cmp[["gen"]], 1)
  clp <- ren(s[["k_clp"]], s[["sr_clp"]], mpp[["age"]], mpp[["gen"]], 1)
  ly  <- ren(s[["k_ly"]],  s[["sr_ly"]],  clp[["age"]], clp[["gen"]], 1)
  blood <- function(up, mtt, c, d) c(age = up[["age"]] + mtt + 1 / c, gen = up[["gen"]] + 1 + d)
  out <- rbind(MPP = mpp, CMP = cmp, GMP = gmp, MEP = mep, CLP = clp, L = ly,
               ANC = blood(gmp, s[["mtt_n"]], s[["c_n"]], s[["d_n"]]),
               AMC = blood(gmp, s[["mtt_mo"]], s[["c_mo"]], s[["d_mo"]]),
               PLT = blood(mep, s[["mtt_pl"]], s[["c_pl"]], s[["d_pl"]]),
               ALC = c(age = ly[["age"]] + 1 / s[["c_l"]], gen = ly[["gen"]]))
  data.table::data.table(compartment = rownames(out), age_days = out[, "age"], generations = out[, "gen"])
}

# Typical design used by the simulation steps (no bridging, normal kidney function).
typical_design <- function(config, t0 = config$model$t_min) {
  list(t0 = t0, egfr = 90, bridging = "none", apheresis_day = -27, ld_scale = 1,
       gcsf = NULL, tpo = NULL, cd34 = NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
