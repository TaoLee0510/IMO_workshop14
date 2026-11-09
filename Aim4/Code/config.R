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
    # v1.3 (2026-11-03): CAR T - target interaction module, G-CSF as input, individual
    # baselines, DE-ITS estimation. v1.2 is frozen in Aim4/Code_v12_frozen with its
    # results in Results/<data> and Data/derived.
    version = "v1.3",
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
      derived    = file.path(root, "Data", "derived_v13"),
      testdata   = file.path(root, "TestData_v13"),
      results    = file.path(root, "Results", paste0(data, "_v13")),
      rcpp_cache = file.path(root, "Code", ".rcpp_cache"),
      greenlight = file.path(root, "GreenLight")
    ),

    # ---- data-handling rules (Aim 3 Data_review section 3, reused) ---------------
    data = list(
      obs_vars = c("anc", "plt", "amc", "alc"),   # v1.3: all four in the likelihood (ALC = lymphocytes + CAR T)
      wbc_no_diff = 0.5,             # blank differential below this WBC -> value <= WBC
      manual_diff_cells = 100,       # zero on a manual differential -> value < WBC / 100
      gcsf_from_day = 14,            # v1.2 only: ANC after this day was an upper bound (G-CSF)
      # v1.3 G-CSF windows (only "needed after day 14" is recorded): one dose the day before
      # each isolated ANC spike; flagged patients: 5-day course from the first day >= 15 with
      # ANC < 1.0, repeated while ANC falls below 1.0 again (up to day 100)
      gcsf_rule = list(from = 15, thr = 1.0, days = 5, until = 100),
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
        # G-CSF (v1.3: input with imputed windows; effect GC rises at k_g while given):
        # chain entry x (1 + (gcsf_fold - 1) GC), release of mature neutrophils x (1 + gcsf_rel GC)
        gcsf_fold = 2.5, gcsf_rel = 2, k_g = 1, tpo_fold = 1.8,
        d_n = 4, d_mo = 2, d_pl = 0,
        # CAR T - target interaction (v1.3). car_t0: CAR T in blood right after the infusion
        # (~1 cell/uL, assumption); k_ag: antigen half-saturation (relative load, median
        # patient a0 = 1); kappa: calibrated so that the typical design peaks on day 8 (K1
        # calibrate_cart); k_t: contraction (prior centre); beta_a: antigen load per log LDH;
        # k_z: rate of each stage of the two-stage lag of the signal Z; z50: typical peak of Z
        # (K1); th_zs: HSC loss per unit Z (centre: typical Z-route HSC loss = typical
        # ferritin-route loss at the prior centres, K1); e_z: maximal suppression of MPP
        # self-renewal by Z
        car_t0 = 0.001, k_ag = 0.1, kappa = 10.65, k_t = 0.2, beta_a = 1, k_z = 0.5, z50 = 0.0140,
        th_zs = 2.64, e_z = 0.5),
      # bookkeeping constants reported with the generation outputs (not in the ODEs)
      bookkeeping = list(endomitosis_rounds = 4, endomitosis_range = c(3, 5), ploidy_range = c(16, 64),
                         platelets_per_megakaryocyte = 1e3),
      # parameters on (0, 1) use a logit working scale, the rest log
      unit_interval = c("emax", "e_pl", "e_mo", "e_l", "e_sr", "e_z"),
      struct_prior_sd = 0.5,
      # prior SDs (working scale) that differ from struct_prior_sd
      struct_prior_sd_k = c(e_sr = 1, th_zs = 1, e_z = 1),
      # lower / upper limits (natural scale) of the structural search box that are
      # tighter than prior centre +- 4 prior SD
      struct_limits = list(gcsf_fold = c(1, 10)),
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
        # rho centre: typical CAR T peak 0.05 x 10^9/L on day 8 = 50-fold over car_t0 (K1 calibrate_cart)
        mu    = c(r_s = log(0.08), f0 = qlogis(0.02), th_as = log(0.12), th_ap = log(2.31),
                  th_bs = log(0.02), rho = log(0.811), slope = log(1.5), mtt_scale = 0),
        omega = c(r_s = 0.5, f0 = 1.5, th_as = 0.8, th_ap = 0.5, th_bs = 1.0, rho = 0.5,
                  slope = 0.5, mtt_scale = 0.2)),
      covariates = list(r_s = c("prior_lines_c", "prior_asct", "age_c"),
                        f0 = c("age_c", "prior_lines_c", "prior_asct", "log_rdw_c", "mcv_c"),
                        th_as = character(0), th_ap = character(0), th_bs = character(0),
                        rho = character(0), slope = character(0), mtt_scale = character(0)),
      sigma_init = c(anc = 0.3, plt = 0.25, amc = 0.35, alc = 0.4),
      # individual baselines (v1.3): log b = log(set-point median) + d, d ~ N(0, s^2/n + tau^2);
      # s from the data (pooled within-patient SD of the set-point values), tau estimated
      baseline = list(vars = c("anc", "plt", "amc", "alc"), tau_init = 0.3, tau_min = 0.02,
                      bound_sd = 4, bound_min = 1),
      # physiological search box of the individual parameters (natural scale; DE-ITS)
      indiv_limits = list(r_s = c(0.003, 1), th_ap = c(0.2, 30), th_as = c(0.001, 3),
                          th_bs = c(1e-4, 2), f0 = c(0.001, 0.95), rho = c(0.1, 3),
                          slope = c(0.05, 30), mtt_scale = c(0.3, 3)),
      # forecast driver: Gamma-shaped kernel log(fer/base) = A (t/tp)^h exp(h (1 - t/tp))
      kernel = list(h = 2),
      # design values for steps that use no counts (K2, K3, TestData)
      design = list(
        mu = c(r_s = log(0.08), f0 = qlogis(0.15), th_as = log(0.12), th_ap = log(2.31),
               th_bs = log(0.02), rho = log(0.811)),
        omega = c(r_s = 0.4, f0 = 0.8, th_as = 0.6, th_ap = 0.4, th_bs = 0.8, rho = 0.4),
        target = c(quick = 19, second_wave = 10, delayed = 6) / 35,
        sigma = c(anc = 0.25, plt = 0.2, amc = 0.3, alc = 0.3),
        tau = c(anc = 0.2, plt = 0.15, amc = 0.2, alc = 0.25),   # set-point drift (synthetic)
        gcsf = list(p = 0.7, fold = c(2, 4), days = 5, anc_thr = 1.0, from = 15),
        plt_tx_below = 15,
        cens = list(anc = 0.2, anc_bound = 0.3, amc = 0.02, amc_bound = 0.05),   # amc rule also used for ALC
        eta_syn = 3,
        p_clone = 0.35)
    ),

    # ---- fitting -----------------------------------------------------------------------
    fit = list(
      q2_models = c("M0", "Ma", "Mb", "Mab", "Mabx"),
      # v1.3 waves: wave 1 the Q2 set and MF; wave 2 for the selected Q2 model: the two clone
      # types and the two single-route variants (with and without the DDR clone)
      wave1 = c("M0", "Ma", "Mb", "Mab", "Mabx", "MF"),
      struct_est = list(base = c("mtt_n", "mtt_pl", "mtt_mo", "ic50", "k_i", "gamma_m", "u_br_rgo", "e_sr",
                                 "k_t", "beta_a", "gcsf_fold", "gcsf_rel"),
                        route_hsc = "th_zs",   # b models without the _Zm suffix
                        route_mpp = "e_z",     # parent models without the _Zh suffix
                        Mabx = "eta_syn",
                        MF = c("mtt_n", "mtt_pl", "mtt_mo", "gam_mf", "slope_pl", "slope_mo", "u_br_rgo")),
      eta_syn_prior = c(centre = 1, sd = 1.2),   # working scale log(eta_syn)
      e_sr_prior_sd = 1,                         # working scale logit(e_sr); centre = struct value
      omega_nu = 4, omega_max = 2.5,     # shrinkage of random-effect SDs towards the prior spread
      its_max_iter = 12L, its_tol = 0.02, struct_every = 2L, map_maxit = 150L,
      use_gcsf_upper = FALSE,            # v1.2 handling (G-CSF ANC as upper bounds); v1.3: ordinary rows
      # DE-ITS (v1.3). E step: DEoptim per patient (strategy 2, NP = np, F, CR), run in chunks;
      # converged when the best value improved by < tol_best * (|f| + 1) over `window`
      # generations and the population objective spread (90% quantile - min) < tol_spread
      de = list(np = 220L, F = 0.8, CR = 0.9, strategy = 2L, chunk = 50L, window = 100L,
                tol_best = 1e-6, tol_spread = 1e-3, max_gen = 20000L,
                struct_np_per_dim = 10L, struct_np_min = 60L, polish_maxit = 500L),
      # outer loop: converged when the max change of the population vector < tol_par and the
      # relative change of sum(nll) < tol_obj on n_ok consecutive rounds; oscillating
      # components are damped; max_iter = safety cap (flag, not a fit).
      # Deviation 2026-11-04 (stochastic E steps leave a noise floor, weak components drift
      # slowly): components whose change (> relax_min) keeps its direction for relax_after
      # rounds are over-relaxed (step x relax_grow per round, at most lam_max); when no
      # component shows a material (> tol_drift over the window) and significant (> 2 SE) trend
      # over `window` rounds, the averaging phase starts from the window mean (stochastic
      # approximation, step 1/n) and is converged when the change of the average < tol_par on
      # n_ok consecutive rounds.
      # Deviation 2026-11-04 09:40 (over-relaxed steps overshot; covariate effects of f0 drifted
      # without bound): over-relaxation is switched off for every component when sum(nll) rose by
      # more than reset_rise (Salakhutdinov & Roweis 2003), the extrapolated part of a step is at
      # most max_extra per component, and covariate effects are kept in their legal range
      outer = list(tol_par = 0.01, tol_obj = 1e-4, n_ok = 2L, damp = 0.5, damp_min = 0.25, max_iter = 120L,
                   relax_min = 0.002, relax_after = 3L, relax_grow = 1.5, lam_max = 10, window = 10L,
                   tol_drift = 0.02, reset_rise = 3, max_extra = 0.5),
      # importance-sampling marginal likelihood (final E step)
      is = list(n = 2000L, df = 4, scale = 1.5, ess_min = 30),
      mcmc = list(n_chains = 2L, n_iter = 1500L, burn = 500L, thin = 2L, scale = 2.38,
                  rhat_tol = 1.05, extend = 1500L, max_iter = 9000L),
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
      # v1.3 F3: approximate landmark prediction (population of the full-data fit)
      F3 = list(landmarks = c(0, 7, 10, 14, 21), n_draws = 300L, n_kernel = 200L, refit_population = FALSE),
      # v1.3 recovery test (replaces K3): synthetic real35 cohorts from the final model with
      # structural values off the prior centres (multiplier on the natural scale)
      R = list(n_cohorts = 10L, model = "MabC_ddr",
               truth_mult = c(th_zs = 2, e_z = 0.6, k_t = 1.5, gcsf_rel = 0.7, beta_a = 0.7, e_sr = 0.8))
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
