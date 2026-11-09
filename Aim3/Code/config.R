# =============================================================================
# config.R
# Purpose : Single source of truth for paths, data-handling rules, phenotype
#           rule, model defaults, fitting controls, step settings and plotting
#           for Aim 3 (when is the stem cell pool injured: chemotherapy (a)
#           vs CAR T phase (b)). Every function receives the config object
#           explicitly; nothing reads globals.
# Inputs  : quick flag, cores.
# Outputs : a nested list.
# =============================================================================

make_config <- function(quick = FALSE, cores = NULL, root = NULL) {
  if (is.null(root)) root <- aim3_root()
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  cfg <- list(
    quick = quick,
    runtime = list(cores = cores, seed = 20261103L),

    # ---- paths ----------------------------------------------------------------
    paths = list(
      root       = root,
      code       = file.path(root, "Code"),
      raw_xlsx   = file.path(root, "..", "Aim2", "Data", "Copy of ICAHT_twin_dataset_deidentified.xlsx"),
      derived    = file.path(root, "Data", "derived"),
      results    = file.path(root, "Results", "real"),
      rcpp_cache = file.path(root, "Code", ".rcpp_cache")
    ),

    # ---- data-handling rules (see Data_review.md section 3) ----------------
    data = list(
      lineages = c("anc", "plt", "amc"),
      wbc_no_diff = 0.5,             # blank differential below this WBC -> value <= WBC
      manual_diff_cells = 100,       # zero on a manual differential -> value < WBC / 100
      gcsf_from_day = 14,            # need_gcsf_after_d14: ANC after this day is an upper bound
      # steroids (dates not exported): drug-free ANC <= observed (neutrophilia) and
      # drug-free AMC >= observed (monocytopenia) in this window. Decided 2026-11-03
      # after the second wave was found to track steroids (8/10).
      steroid_window = c(0, 21),
      plt_tx = list(plt_below = 30, rise = 30, hgb_rise_max = 0.5, drop_days = 3),
      setpoint_window = c(-30, -6),  # pre-treatment set point (bridging chemo: <= apheresis day)
      fit_window = c(-30, 120),
      egfr_ref = 90,
      bridging = list(start_after_apheresis = 2, end_day = -7,
                      regimen = c(rgo = "R-Gem-Ox", oral_cy = "cy"))
    ),

    # ---- phenotype rule (locked 2026-11-03; Data_review.md section 5.2) --------
    phenotype = list(
      primary = list(window_end = 30, nadir_window = c(0, 14), rec_thr = 1.5, rec_n = 1,
                     dip_thr = 1.0, dip_n = 1, dip_after = 0, drop_spikes = TRUE),
      # sensitivity variants: changes relative to the primary rule
      variants = list(
        win35     = list(window_end = 35),
        rejeski21 = list(dip_after = 21),
        sust2     = list(rec_n = 2, dip_n = 2),
        rec10     = list(rec_thr = 1.0),
        dip05     = list(dip_thr = 0.5)),
      classes = c("quick", "second_wave", "delayed")
    ),

    # ---- model ----------------------------------------------------------------------
    model = list(
      t_start = -30, t_end = 120,
      solver = list(rtol = 1e-6, atol = 1e-9, hmax = 0.5),
      drv_dt = 0.25,
      # structural parameters (natural scale; prior centres). Tuned 2026-11-03 so
      # that a typical patient matches the cohort's median shape (nadir days 5-9,
      # recovery days 11-14, monocytes first). phi * gmax < 1: multipotent
      # progenitors cannot maintain themselves without stem cell input; rho_m:
      # multipotent progenitors are half as sensitive to injury as committed ones.
      struct = c(
        alpha = 1, k_m = 0.2, phi = 0.5, gmax = 1.3, gamma_m = 0.5,
        mtt_n = 6, mtt_pl = 5, mtt_mo = 2,
        c_n = 2.3, c_pl = 0.15, c_mo = 1.5,
        gam_n = 0.25, gam_pl = 0.2, gam_mo = 0.2,
        k_e = 0.4, q_b = 0.1, eta_syn = 0, rho_m = 0.5,
        eps = 0, w_pl = 0.2, w_mo = 0.5,   # accelerated release: off by default (kept for drug inputs)
        beta_renal = 0.7, u_br_rgo = 0.05, u_br_cy = 0.01,
        delta_b = 1, k_dec = 0.07),
      # log-scale prior SD of estimated structural parameters (eta_syn centred at 1)
      struct_prior_sd = c(mtt_n = 0.3, mtt_pl = 0.3, mtt_mo = 0.4,
                          u_br_rgo = 1.0, eta_syn = 1.0),
      # random effects (log scale): prior centres and between-patient SD.
      # q_b (inflammatory suppression of progenitors in the CAR T phase) is
      # individual, so that deep reversible dips are not forced onto the HSC terms.
      re_prior = list(
        mu    = c(th_ap = log(4), th_as = log(0.1), th_bs = log(0.03), r_s = log(0.05), q_b = log(0.1)),
        omega = c(th_ap = 0.5, th_as = 1.0, th_bs = 1.0, r_s = 0.5, q_b = 1.0)),
      covariates = list(th_ap = character(0), th_as = character(0), th_bs = character(0),
                        r_s = c("prior_lines_c", "prior_asct"), q_b = character(0)),
      sigma_init = c(anc = 0.3, plt = 0.25, amc = 0.35),
      # design values for steps that use no counts (W2, W3)
      design = list(
        mu = c(th_ap = log(4), th_as = log(0.1), th_bs = log(0.03), r_s = log(0.05), q_b = log(0.1)),
        omega = c(th_ap = 0.4, th_as = 0.8, th_bs = 0.8, r_s = 0.4, q_b = 0.8),
        target = c(quick = 19, second_wave = 10, delayed = 6) / 35,
        sigma = c(anc = 0.25, plt = 0.2, amc = 0.3),
        gcsf_fold = c(2, 4), steroid_anc_fold = c(1.5, 3), steroid_amc_fold = c(0.3, 0.7))
    ),

    # ---- fitting -----------------------------------------------------------------
    fit = list(
      models = c("M0", "Ma", "Mb", "Mab", "Mabx"),
      struct_est = list(M0 = c("mtt_n", "mtt_pl", "mtt_mo", "u_br_rgo"),
                        Ma = c("mtt_n", "mtt_pl", "mtt_mo", "u_br_rgo"),
                        Mb = c("mtt_n", "mtt_pl", "mtt_mo", "u_br_rgo"),
                        Mab = c("mtt_n", "mtt_pl", "mtt_mo", "u_br_rgo"),
                        Mabx = c("mtt_n", "mtt_pl", "mtt_mo", "u_br_rgo", "eta_syn")),
      its_max_iter = 10L, its_tol = 0.02, struct_every = 2L, map_maxit = 200L,
      mcmc = list(n_iter = 2000L, burn = 500L, thin = 2L, scale = 2.38),
      attrib_draws = 300L
    ),

    # ---- step settings ---------------------------------------------------------------
    steps = list(
      W2 = list(f_grid = c(0, 0.3, 0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995),
                loss_demo = 0.9, eta_levels = c(0, 3),
                peak_shifts = c(-3, 0, 3, 6),
                # MPP regimes for the pseudo-second-wave scan (phi * gmax >= 1 lets
                # progenitors expand transiently beyond stem cell supply)
                regimes = list(default = c(k_m = 0.2, phi = 0.5, gmax = 1.3, rho_m = 0.5),
                               burst   = c(k_m = 0.3, phi = 0.7, gmax = 2.0, rho_m = 0.5),
                               strong  = c(k_m = 0.5, phi = 0.8, gmax = 3.0, rho_m = 0.3))),
      W3 = list(n_rep = 20L, scenarios = c("a_only", "b_only", "mix", "synergy"),
                mix = c(a = 0.4, b = 0.3, ab = 0.3), eta_syn = 3,
                calib_n = 1500L, n_large = 350L, n_rep_large = 2L, attrib_max = 70L,
                route = list(cohort_correct = 0.80, indiv_prob = 0.60)),
      W6 = list(landmarks = c(0, 7, 10, 14, 21))
    ),

    # ---- plotting -------------------------------------------------------------------
    plot = list(
      width = 7, height = 4.5, dpi = 300,
      col_phenotype = c(quick = "#3F7D5C", second_wave = "#C9A227", delayed = "#9E2A2B"),
      col_mech = c(none = "#3F7D5C", a = "#2E6B8A", b = "#C9A227", ab_joint = "#9E2A2B",
                   ab_redundant = "#6E3B6E", other = "#6B7A8F"),
      col_lineage = c(anc = "#2E6B8A", plt = "#9E2A2B", amc = "#C9A227"),
      col_accent = c(crimson = "#9E2A2B", steel = "#2E6B8A", slate = "#6B7A8F",
                     ink = "#1A2433", gold = "#C9A227", green = "#3F7D5C")
    )
  )

  if (quick) cfg <- apply_quick_mode(cfg)
  cfg
}

# Smaller grids and fewer draws for a smoke test.
apply_quick_mode <- function(cfg) {
  cfg$steps$W2$f_grid <- c(0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99)
  cfg$steps$W3$n_rep <- 1L
  cfg$steps$W3$n_rep_large <- 0L
  cfg$steps$W3$calib_n <- 300L
  cfg$fit$its_max_iter <- 3L
  cfg$fit$map_maxit <- 100L
  cfg$fit$mcmc <- list(n_iter = 600L, burn = 200L, thin = 2L, scale = 2.38)
  cfg$fit$attrib_draws <- 60L
  cfg
}

# Locate Aim3/ from this file's location (works with Rscript and source()).
aim3_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(normalizePath(file.path(dirname(f), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), "..")))
  }
  normalizePath(getwd())
}
