# =============================================================================
# config.R
# Purpose : Single source of truth for paths, data-handling rules, phenotype
#           rules, model defaults, fitting controls, step settings and plotting
#           for Model1 (Aim4 v1.3 with a Friberg-like haematopoietic layer:
#           stem cells -> four MPP components -> maturation chains -> blood;
#           recovery phenotype in the likelihood). Every function receives the
#           config object explicitly; nothing reads globals.
# Inputs  : data source ("real" | "test"), quick flag, cores.
# Outputs : a nested list.
# =============================================================================

make_config <- function(data = "test", quick = FALSE, cores = NULL, root = NULL) {
  if (is.null(root)) root <- model1_root()
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  stopifnot(data %in% c("real", "test"))

  cfg <- list(
    # Model1 v1.0 (2026-11-04): copied from Aim4/Code v1.3 (hashes in SOURCE_aim4_v13_sha256.txt),
    # plan Model1/Model1_plan.md v0.3
    version = "Model1 v1.0",
    data_source = data,
    quick = quick,
    runtime = list(cores = cores, seed = 20261104L),

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
      aim4_ref   = file.path(root, "Code", "tests", "ref", "aim4_v13.cpp"),
      rcpp_cache = file.path(root, "Code", ".rcpp_cache"),
      greenlight = file.path(root, "GreenLight")
    ),

    # ---- data-handling rules (Aim 3 Data_review section 3, reused) ---------------
    data = list(
      obs_vars = c("anc", "plt", "amc", "alc"),   # all four in the likelihood (ALC = lymphocytes + CAR T)
      wbc_no_diff = 0.5,             # blank differential below this WBC -> value <= WBC
      manual_diff_cells = 100,       # zero on a manual differential -> value < WBC / 100
      gcsf_from_day = 14,            # Aim4 v1.2 only: ANC after this day was an upper bound (G-CSF)
      # G-CSF windows (only "needed after day 14" is recorded): one dose the day before
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
      # Model1 (plan section 2), natural scale (prior centres). Stem cells as Aim4 v1.3 without
      # the MPP-depletion activation alpha (decision D1). The four MPP components (N, Mo, Pl, L)
      # share the exit rate k_p and the self-renewal rate sr_p (= Aim4 MPP; effective residence
      # 1/(k_p - sr_p) = 11 days; D3) and one feedback exponent gamma (D4) with cap gmax; their
      # self-renewal is suppressed by marrow inflammation E(I) (maximal effect emax, which
      # replaces Aim4 e_sr) and by the CAR T - target signal (e_z). mtt_* are converted to
      # ktr = 3 / mtt (three transit compartments); mtt_pl = megakaryocyte maturation;
      # lymphocytes have no chain (D2): blood turnover c_l only.
      struct = c(
        r_a = 1, r_b = 1, delta = 1, eta_syn = 0,
        th_zs = 2.64,                    # HSC loss per unit Z (Aim4 v1.3 calibration, independent of the layers below)
        k_p = 0.3, sr_p = 0.21, gmax = 3, gamma = 0.5,
        mtt_n = 5, mtt_mo = 1.5, mtt_pl = 5,
        c_n = 2.3, c_mo = 0.5, c_pl = 0.12,
        c_l = 0.156, th_al = 0.921,      # K1 calibrate_lymph: typical ALC on days -2..0 = 0.042 (real data) and
                                         # on day 30 = Aim4 v1.3 (0.848)
        # emax centre 0.7: the smallest value of the K1 phenotype_screen grid (0.3-0.9) at which the CAR T
        # phase alone (driver A 2.5, th_bS 0.06) gives a second wave
        k_i = 0.2, emax = 0.7, ic50 = 0.7, hill = 2,
        k_e = 0.8, beta_renal = 0.7, u_br_rgo = 0.05, u_br_cy = 0.01,
        # G-CSF (input with imputed windows; effect GC rises at k_g while given):
        # chain entry x (1 + (gcsf_fold - 1) GC), release of mature neutrophils x (1 + gcsf_rel GC)
        gcsf_fold = 2.5, gcsf_rel = 2, k_g = 1, tpo_fold = 1.8,
        # CAR T - target interaction (Aim4 v1.3, calibrations unchanged: they do not involve the
        # haematopoietic layers). car_t0: CAR T in blood right after the infusion (~1 cell/uL);
        # k_ag: antigen half-saturation (relative load, median patient a0 = 1); kappa: typical
        # design peaks on day 8; k_t: contraction (prior centre); beta_a: antigen load per log
        # LDH; k_z: rate of each stage of the two-stage lag of Z; z50: typical peak of Z; e_z:
        # maximal suppression of MPP self-renewal by Z
        car_t0 = 0.001, k_ag = 0.1, kappa = 10.65, k_t = 0.2, beta_a = 1, k_z = 0.5, z50 = 0.0140,
        e_z = 0.5),
      # parameters on (0, 1) use a logit working scale, the rest log
      unit_interval = c("emax", "e_z"),
      struct_prior_sd = 0.5,
      # prior SDs (working scale) that differ from struct_prior_sd
      struct_prior_sd_k = c(emax = 1, th_zs = 1, e_z = 1),
      # lower / upper limits (natural scale) of the structural search box that are
      # tighter than prior centre +- 4 prior SD; gamma: the pre-treatment steady state loses
      # stability (Hopf) at gamma = 2.37 with sr_p = 0.21 (K1 2026-11-04), so the box stops at 2.2
      struct_limits = list(gcsf_fold = c(1, 10), gamma = c(0.05, 2.2)),
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
        # th_ap centre 2.71: the typical-design ANC nadir equals Aim4 v1.3 at its prior centre 2.31
        # (K1 calibrate_th_ap); rho centre: typical CAR T peak 0.05 x 10^9/L on day 8 = 50-fold over
        # car_t0 (Aim4 v1.3 K1)
        mu    = c(r_s = log(0.08), f0 = qlogis(0.02), th_as = log(0.12), th_ap = log(2.71),
                  th_bs = log(0.02), rho = log(0.811), slope = log(1.5), mtt_scale = 0),
        omega = c(r_s = 0.5, f0 = 1.5, th_as = 0.8, th_ap = 0.5, th_bs = 1.0, rho = 0.5,
                  slope = 0.5, mtt_scale = 0.2)),
      covariates = list(r_s = c("prior_lines_c", "prior_asct", "age_c"),
                        f0 = c("age_c", "prior_lines_c", "prior_asct", "log_rdw_c", "mcv_c"),
                        th_as = character(0), th_ap = character(0), th_bs = character(0),
                        rho = character(0), slope = character(0), mtt_scale = character(0)),
      sigma_init = c(anc = 0.3, plt = 0.25, amc = 0.35, alc = 0.4),
      # individual baselines: log b = log(set-point median) + d, d ~ N(0, s^2/n + tau^2);
      # s from the data (pooled within-patient SD of the set-point values), tau estimated
      baseline = list(vars = c("anc", "plt", "amc", "alc"), tau_init = 0.3, tau_min = 0.02,
                      bound_sd = 4, bound_min = 1),
      # physiological search box of the individual parameters (natural scale; DE-ITS)
      indiv_limits = list(r_s = c(0.003, 1), th_ap = c(0.2, 30), th_as = c(0.001, 3),
                          th_bs = c(1e-4, 2), f0 = c(0.001, 0.95), rho = c(0.1, 3),
                          slope = c(0.05, 30), mtt_scale = c(0.3, 3)),
      # forecast driver: Gamma-shaped kernel log(fer/base) = A (t/tp)^h exp(h (1 - t/tp))
      kernel = list(h = 2),
      # design values for steps that use no counts (K2, synthetic cohorts, TestData)
      design = list(
        mu = c(r_s = log(0.08), f0 = qlogis(0.15), th_as = log(0.12), th_ap = log(2.71),
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
      # waves: wave 1 the Q2 set and MF; wave 2 for the selected Q2 model: the two clone
      # types and the two single-route variants (with and without the DDR clone)
      wave1 = c("M0", "Ma", "Mb", "Mab", "Mabx", "MF"),
      struct_est = list(base = c("mtt_n", "mtt_pl", "mtt_mo", "ic50", "k_i", "gamma", "u_br_rgo", "emax",
                                 "k_t", "beta_a", "gcsf_fold", "gcsf_rel"),
                        route_hsc = "th_zs",   # b models without the _Zm suffix
                        route_mpp = "e_z",     # models without the _Zh suffix
                        Mabx = "eta_syn",
                        MF = c("mtt_n", "mtt_pl", "mtt_mo", "gam_mf", "slope_pl", "slope_mo", "u_br_rgo")),
      eta_syn_prior = c(centre = 1, sd = 1.2),   # working scale log(eta_syn)
      omega_nu = 4, omega_max = 2.5,     # shrinkage of random-effect SDs towards the prior spread
      its_max_iter = 12L, its_tol = 0.02, struct_every = 2L, map_maxit = 150L,
      use_gcsf_upper = FALSE,            # G-CSF is an input; ANC under G-CSF are ordinary rows
      # recovery phenotype in the likelihood (plan section 5): - w log P'(c | psi) with the
      # margin-based soft classification of the model ANC on the patient's label days
      # (smoothing scale s on the log scale) and label error eps: P' = (1 - 3 eps) P + eps.
      # Composite likelihood (the day 0-30 ANC enter twice); not used at landmarks (F3).
      pheno = list(w = 1, s = 0.3, eps = 0.02, sens_s = c(0.15, 0.6), sens_w = c(0, 5)),
      # round 0 (plan section 6): every patient fitted on its own with the prior SDs times
      # sd_mult (search box unchanged); the MAPs and Hessians initialise mu, beta, omega, tau
      # and sigma by the M-step formulas. Runs warm-started from another fit skip round 0.
      round0 = list(sd_mult = 4),
      # DE-ITS. E step: DEoptim per patient (strategy 2, NP = np, F, CR), run in chunks;
      # converged when the best value improved by < tol_best * (|f| + 1) over `window`
      # generations and the population objective spread (90% quantile - min) < tol_spread
      de = list(np = 220L, F = 0.8, CR = 0.9, strategy = 2L, chunk = 50L, window = 100L,
                tol_best = 1e-6, tol_spread = 1e-3, max_gen = 20000L,
                struct_np_per_dim = 10L, struct_np_min = 60L, polish_maxit = 500L),
      # outer loop (Aim4 v1.3 settings incl. its two deviations of 2026-11-04): converged when
      # the max change of the population vector < tol_par and the relative change of sum(nll)
      # < tol_obj on n_ok consecutive rounds; oscillating components are damped; components
      # whose change (> relax_min) keeps its direction for relax_after rounds are over-relaxed
      # (step x relax_grow per round, at most lam_max; extrapolated part at most max_extra;
      # switched off when sum(nll) rose by more than reset_rise); when no component shows a
      # material (> tol_drift over the window) and significant (> 2 SE) trend over `window`
      # rounds, the averaging phase starts from the window mean (stochastic approximation, step
      # 1/n) and is converged when the change of the average < tol_par on n_ok consecutive
      # rounds; covariate effects are kept in their legal range; max_iter = safety cap.
      outer = list(tol_par = 0.01, tol_obj = 1e-4, n_ok = 2L, damp = 0.5, damp_min = 0.25, max_iter = 120L,
                   relax_min = 0.002, relax_after = 3L, relax_grow = 1.5, lam_max = 10, window = 10L,
                   tol_drift = 0.02, reset_rise = 3, max_extra = 0.5),
      # importance-sampling marginal likelihood (final E step)
      is = list(n = 2000L, df = 4, scale = 1.5, ess_min = 30),
      mcmc = list(n_chains = 2L, n_iter = 1500L, burn = 500L, thin = 2L, scale = 2.38,
                  rhat_tol = 1.05, extend = 1500L, max_iter = 9000L),
      attrib_draws = 400L
    ),

    # ---- validation (plan section 8) --------------------------------------------------------
    validation = list(
      # V1: real data, 5 folds stratified by the primary phenotype (fixed at PR); every fold
      # refit is warm-started from the full-data fit of the final model
      V1 = list(k = 5L, seed = 20261104L),
      # S1: synthetic cohorts with known truth from Model1 on the real designs (recovery of the
      # parameters and of the attribution); structural truth off the prior centres
      S1 = list(n_cohorts = 10L, model = "MabC_ddr",
                truth_mult = c(th_zs = 2, e_z = 0.6, k_t = 1.5, gcsf_rel = 0.7, beta_a = 0.7, emax = 0.8)),
      # S2: the S1 cohorts refitted with w = 0 (bias of the phenotype term, coverage)
      S2 = list(w = 0),
      # S3: V1 on the first n_cohorts S1 cohorts (prediction ceiling)
      S3 = list(n_cohorts = 2L),
      # S4: cohorts generated with the Aim4 v1.3 lineage-tree model (truth: Aim4 fitted population
      # when available, else its prior centres), fitted with Model1
      S4 = list(n_cohorts = 5L, aim4_model = "MabC_ddr"),
      # S5: visual predictive check from the fitted population on the real designs
      S5 = list(n_rep = 500L)
    ),

    # ---- step settings ---------------------------------------------------------------------
    steps = list(
      K2 = list(n_grid = 25, th_as_range = c(0.005, 1.2), th_bs_range = c(0.0005, 0.15),
                r_s_range = c(0.01, 0.4), f0_range = c(0, 0.8), eta_levels = c(0, 3)),
      # F3: landmark prediction, out of sample through the V1 folds (no phenotype term)
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
      col_lineage = c(anc = "#2E6B8A", plt = "#9E2A2B", amc = "#C9A227", alc = "#6E3B6E"),
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
  cfg$fit$its_max_iter <- 4L
  cfg$fit$map_maxit <- 80L
  cfg$fit$mcmc <- list(n_chains = 1L, n_iter = 500L, burn = 200L, thin = 2L, scale = 2.38)
  cfg$fit$attrib_draws <- 60L
  cfg$steps$F3$n_draws <- 60L
  cfg$steps$F3$n_kernel <- 60L
  cfg$validation$S5$n_rep <- 50L
  cfg
}

# Locate Model1/ from this file's location (works with Rscript and source()).
model1_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(normalizePath(file.path(dirname(f), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), "..")))
  }
  normalizePath(getwd())
}
