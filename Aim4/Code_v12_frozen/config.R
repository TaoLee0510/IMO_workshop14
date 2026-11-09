# =============================================================================
# config.R
# Purpose : Single source of truth for paths, data-handling rules, phenotype
#           rules, model defaults, fitting controls, step settings and plotting
#           for Aim 4 (one parent marrow model: who fails to recover - reserve
#           or clone - and when the stem cell pool is hit - chemotherapy (a) or
#           CAR T phase (b)). Every function receives the config object
#           explicitly; nothing reads globals.
# Inputs  : data source ("real" | "test"), quick flag, cores.
# Outputs : a nested list.
# =============================================================================

make_config <- function(data = "test", quick = FALSE, cores = NULL, root = NULL) {
  if (is.null(root)) root <- aim4_root()
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  stopifnot(data %in% c("real", "test"))

  cfg <- list(
    data_source = data,
    quick = quick,
    runtime = list(cores = cores, seed = 20261103L),

    # ---- paths ----------------------------------------------------------------
    paths = list(
      root       = root,
      code       = file.path(root, "Code"),
      # raw extract: read-only reference, never copied (SHA-256 checked)
      raw_xlsx   = normalizePath(file.path(root, "..", "Aim2", "Data",
                                           "Copy of ICAHT_twin_dataset_deidentified.xlsx"), mustWork = FALSE),
      derived    = file.path(root, "Data", "derived"),
      testdata   = file.path(root, "TestData"),
      results    = file.path(root, "Results", data),
      rcpp_cache = file.path(root, "Code", ".rcpp_cache"),
      greenlight = file.path(root, "GreenLight")
    ),

    # ---- data-handling rules (Aim 3 Data_review section 3, reused) ---------------
    data = list(
      obs_vars = c("anc", "plt", "amc", "alc"),   # alc: converted and plotted, not in the likelihood
      wbc_no_diff = 0.5,             # blank differential below this WBC -> value <= WBC
      manual_diff_cells = 100,       # zero on a manual differential -> value < WBC / 100
      gcsf_from_day = 14,            # need_gcsf_after_d14: ANC after this day is an upper bound
      plt_tx = list(plt_below = 30, rise = 30, hgb_rise_max = 0.5, drop_days = 3),
      setpoint_window = c(-30, -6),  # pre-treatment set point (bridging chemo: <= apheresis day)
      obs_max_day = 120,             # last day entering the likelihood (before progression)
      crp_mgdl_to_mgl = 10,
      bridging_window = c(2, -7),    # apheresis day + 2 .. day -7 (no recorded dates)
      fer_baseline = "ferritin_ld",
      fer_smooth_window = c(0, 35)
    ),

    # ---- phenotype rules ------------------------------------------------------------
    # primary rule locked 2026-11-03 (Aim 3 Data_review 5.2), applied unchanged here
    phenotype = list(
      primary = list(window_end = 30, nadir_window = c(0, 14), rec_thr = 1.5, rec_n = 1,
                     dip_thr = 1.0, dip_n = 1, dip_after = 0, drop_spikes = TRUE,
                     cens_value = 0.25),
      variants = list(
        win35     = list(window_end = 35),
        rejeski21 = list(dip_after = 21),
        sust2     = list(rec_n = 2, dip_n = 2),
        rec10     = list(rec_thr = 1.0),
        dip05     = list(dip_thr = 0.5)),
      # simulations / GreenLight only: Aim 2 late class (median ANC on days 45-75 < 0.5)
      persistent = list(window = c(45, 75), thr = 0.5),
      classes = c("quick", "second_wave", "delayed"),
      classes4 = c("quick", "second_wave", "delayed", "persistent")
    ),

    # ---- model ------------------------------------------------------------------------
    model = list(
      t_min = -31,                   # simulations start at min(t_min, apheresis + 1, first obs - 1)
      sim_end = 90,
      solver = list(rtol = 1e-6, atol = 1e-9, hmax = 0.5),
      drv_dt = 0.25,
      drv_end = 130,
      fer_decay = 0.07,              # D(t) after the last ferritin: exponential decay
      # parent model v1.2 (lineage tree Aim4/image (2).png), natural scale (prior centres).
      # k_* exit rate, sr_* self-renewal rate (1/day, sr < k); effective residence
      # 1/(k - sr). fb_* = 1: the compartment's self-renewal senses the blood counts
      # of its own progeny (0: no feedback). mtt_* are converted to ktr = 3 / mtt
      # (three maturation compartments); mtt_pl = megakaryocyte maturation.
      # d_*: divisions inside the maturation chain (bookkeeping only; endomitosis
      # of megakaryocytes is not a division).
      struct = c(
        alpha = 1, gmax = 3, gamma_m = 0.5,
        k_mpp = 0.3, sr_mpp = 0.21,      # = v1.1 k_M and phi * k_M
        # sr of CMP/GMP/MEP (K1 calibration): the value closest to the v1.1 ANC recovery day at
        # which one chemotherapy hit can still give quick, pseudo second wave and delayed
        # recovery (0.5: no pseudo second wave; 1.0: no delayed; K1 table progenitor_screen)
        k_cmp = 2, sr_cmp = 0.75, k_gmp = 2, sr_gmp = 0.75, k_mep = 2, sr_mep = 0.75,
        k_clp = 1, sr_clp = 0.5, k_ly = 0.05, sr_ly = 0.04,
        fb_mpp = 1, fb_cmp = 1, fb_gmp = 1, fb_mep = 1, fb_clp = 1, fb_ly = 1,
        r_a = 1, r_b = 1, delta = 1,
        mtt_n = 5, mtt_mo = 1.5, mtt_pl = 5,
        c_n = 2.3, c_mo = 0.5, c_pl = 0.12, c_l = 2,
        gam_n = 0.2, gam_mo = 0.2, gam_pl = 0.2,
        k_i = 0.2, emax = 0.9, ic50 = 0.7, hill = 2, e_mo = 0.7, e_pl = 0.7, e_l = 0.7, eta_syn = 0,
        e_sr = 0.5,                      # inflammation on MPP/CMP/GMP/MEP self-renewal; estimated in F1 (v1.2)
        th_al = 1.05,                    # lymphodepletion of lymphocytes: typical-design model ALC on days -2..0
                                         # = 4.2% of baseline (median per patient, real ALC, 35 patients; K1 calibrate_th_al)
        k_e = 0.8, beta_renal = 0.7, u_br_rgo = 0.05, u_br_cy = 0.01,
        gcsf_fold = 2.5, tpo_fold = 1.8,
        d_n = 4, d_mo = 2, d_pl = 0),
      # bookkeeping constants reported with the generation outputs (not in the ODEs)
      bookkeeping = list(endomitosis_rounds = 4, endomitosis_range = c(3, 5), ploidy_range = c(16, 64),
                         platelets_per_megakaryocyte = 1e3),
      # parameters on (0, 1) use a logit working scale, the rest log
      unit_interval = c("emax", "e_pl", "e_mo", "e_l", "e_sr"),
      struct_prior_sd = 0.5,
      # Friberg reference model MF (Friberg et al. 2002); kprol = kcirc = ktr
      mf_struct = c(mtt_n = 5, mtt_pl = 6, mtt_mo = 3.5, gam_mf = 0.2,
                    slope_pl = 1, slope_mo = 1, k_e = 0.8, beta_renal = 0.7,
                    u_br_rgo = 0.05, u_br_cy = 0.01),
      # clone types (not estimable on the current data: fixed scenarios)
      clone_types = list(
        ddr = c(r_a = 0.3, r_b = 1, delta = 0.3),   # PPM1D / TP53-like: chemo-tolerant, low output
        dta = c(r_a = 1, r_b = 0.2, delta = 1)),    # DNMT3A / TET2-like: inflammation-tolerant
      # random effects (eta scale: log, logit for f0): prior centres and spread
      re_prior = list(
        # th_ap centre 2.31: the typical-design ANC nadir equals v1.1 at th_ap = 1 (0.221; K1 calibrate_th_ap)
        mu    = c(r_s = log(0.08), f0 = qlogis(0.02), th_as = log(0.12), th_ap = log(2.31),
                  th_bs = log(0.02), slope = log(1.5), mtt_scale = 0),
        omega = c(r_s = 0.5, f0 = 1.5, th_as = 0.8, th_ap = 0.5, th_bs = 1.0,
                  slope = 0.5, mtt_scale = 0.2)),
      covariates = list(r_s = c("prior_lines_c", "prior_asct", "age_c"),
                        f0 = c("age_c", "prior_lines_c", "prior_asct", "log_rdw_c", "mcv_c"),
                        th_as = character(0), th_ap = character(0), th_bs = character(0),
                        slope = character(0), mtt_scale = character(0)),
      sigma_init = c(anc = 0.3, plt = 0.25, amc = 0.35),   # ALC is not in the likelihood
      # forecast driver: Gamma-shaped kernel log(fer/base) = A (t/tp)^h exp(h (1 - t/tp))
      kernel = list(h = 2),
      # design values for steps that use no counts (K2, K3, TestData)
      design = list(
        mu = c(r_s = log(0.08), f0 = qlogis(0.15), th_as = log(0.12), th_ap = log(2.31),
               th_bs = log(0.02)),
        omega = c(r_s = 0.4, f0 = 0.8, th_as = 0.6, th_ap = 0.4, th_bs = 0.8),
        target = c(quick = 19, second_wave = 10, delayed = 6) / 35,
        sigma = c(anc = 0.25, plt = 0.2, amc = 0.3),
        gcsf = list(p = 0.7, fold = c(2, 4), days = 5, anc_thr = 1.0, from = 15),
        plt_tx_below = 15,
        cens = list(anc = 0.2, anc_bound = 0.3, amc = 0.02, amc_bound = 0.05),
        eta_syn = 3,
        p_clone = 0.35)
    ),

    # ---- fitting -----------------------------------------------------------------------
    fit = list(
      q2_models = c("M0", "Ma", "Mb", "Mab", "Mabx"),
      struct_est = list(base = c("mtt_n", "mtt_pl", "mtt_mo", "ic50", "k_i", "gamma_m", "u_br_rgo", "e_sr"),
                        Mabx = "eta_syn",
                        MF = c("mtt_n", "mtt_pl", "mtt_mo", "gam_mf", "slope_pl", "slope_mo", "u_br_rgo")),
      eta_syn_prior = c(centre = 1, sd = 1.2),   # working scale log(eta_syn)
      e_sr_prior_sd = 1,                         # working scale logit(e_sr); centre = struct value
      omega_nu = 4, omega_max = 2.5,     # shrinkage of random-effect SDs towards the prior spread
      its_max_iter = 12L, its_tol = 0.02, struct_every = 2L, map_maxit = 150L,
      mcmc = list(n_chains = 2L, n_iter = 1500L, burn = 500L, thin = 2L, scale = 2.38),
      attrib_draws = 400L
    ),

    # ---- step settings ---------------------------------------------------------------------
    steps = list(
      K2 = list(n_grid = 25, th_as_range = c(0.005, 1.2), th_bs_range = c(0.0005, 0.15),
                r_s_range = c(0.01, 0.4), f0_range = c(0, 0.8), eta_levels = c(0, 3)),
      K3 = list(n_rep = 20L, mech = c("a_only", "b_only", "mix", "synergy"),
                clone = c(none = 0, clone = 0.35), mix = c(a = 0.4, b = 0.3, ab = 0.3),
                calib_n = 1200L, calib_iter = 8L,
                n350 = list(n = 350L, n_rep = 5L, n_vaf = 15L, vaf_sd = 0.5),
                route = list(cohort_correct = 0.80, indiv_prob = 0.60)),
      F3 = list(landmarks = c(0, 7, 10, 14, 21), n_draws = 300L, n_kernel = 200L)
    ),

    # ---- plotting -------------------------------------------------------------------
    plot = list(
      width = 7, height = 4.5, dpi = 300,
      col_phenotype = c(quick = "#3F7D5C", second_wave = "#C9A227", delayed = "#9E2A2B",
                        persistent = "#6E3B6E"),
      col_mech = c(none = "#3F7D5C", a = "#2E6B8A", b = "#C9A227", ab_joint = "#9E2A2B",
                   ab_redundant = "#6E3B6E", other = "#6B7A8F"),
      col_who = c(clone = "#9E2A2B", reserve = "#2E6B8A", both = "#6E3B6E", other = "#6B7A8F"),
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
  cfg$steps$K2$n_grid <- 10
  cfg$steps$K3$n_rep <- 1L
  cfg$steps$K3$calib_n <- 300L
  cfg$steps$K3$calib_iter <- 4L
  cfg$steps$K3$n350$n_rep <- 1L
  cfg$steps$K3$n350$n <- 120L
  cfg$fit$its_max_iter <- 4L
  cfg$fit$map_maxit <- 80L
  cfg$fit$mcmc <- list(n_chains = 1L, n_iter = 500L, burn = 200L, thin = 2L, scale = 2.38)
  cfg$fit$attrib_draws <- 60L
  cfg$steps$F3$n_draws <- 60L
  cfg$steps$F3$n_kernel <- 60L
  cfg
}

# Locate Aim4/ from this file's location (works with Rscript and source()).
aim4_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(normalizePath(file.path(dirname(f), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), "..")))
  }
  normalizePath(getwd())
}
