# =============================================================================
# config.R
# Purpose : Single source of truth for paths, thresholds, model defaults,
#           fitting controls, step settings and plotting. Every function
#           receives the config object explicitly; nothing reads globals.
# Inputs  : data source ("test", "test_stress", "real"), quick flag, cores.
# Outputs : a nested list.
# =============================================================================

make_config <- function(data_source = "test", quick = FALSE, cores = NULL,
                        root = NULL) {
  if (is.null(root)) root <- aim2_root()
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  data_dir <- switch(data_source,
                     test        = file.path(root, "TestData"),
                     test_stress = file.path(root, "TestData", "stress"),
                     real        = file.path(root, "Data"),
                     stop("unknown data source: ", data_source))

  cfg <- list(
    data_source = data_source,
    quick = quick,
    runtime = list(cores = cores, seed = 20261102L),

    # ---- paths ----------------------------------------------------------------
    paths = list(
      root        = root,
      code        = file.path(root, "Code"),
      data        = data_dir,
      results     = file.path(root, "Results", data_source),
      rcpp_cache  = file.path(root, "Code", ".rcpp_cache")
    ),

    # ---- data contract ---------------------------------------------------------
    data = list(
      files = c(patients = "patients.csv", labs = "labs_long.csv",
                inflam = "inflam_long.csv", events = "events.csv",
                split = "split.csv", mut_tested = "mutation_tested.csv",
                mut_pre = "mutations_pre.csv",
                post_planned = "post_samples_planned.csv",
                mut_post = file.path("later", "mutations_post.csv")),
      lod = list(anc = 0.01, plt = 1, crp = 0.5, vaf = 0.01),
      censor = list(
        plt_after_transfusion_days = 3,     # drop PLT for N days after a PLT transfusion
        crp_after_tocilizumab = "censor",   # "censor" or "keep"
        drop_after = c("cd34_boost", "relapse", "subsequent_tx", "death")
      ),
      exclude_flags = c("marrow_involvement", "iec_hs"),
      # baseline counts: median over these days (pre-lymphodepletion); the
      # earliest value on or before day 0 when there is none
      baseline_window = c(-Inf, -5),
      # treatment days when only yes/no is recorded (R/events_approx.R)
      approx_events = list(toci_no_crs_day = 3, steroid_days = 4, steroid_no_crs_day = 5,
                           gcsf_from = 14, gcsf_anc = 1.0, gcsf_days = 5),
      ch_genes = list(DDR = c("PPM1D", "TP53", "CHEK2"),
                      DTA = c("DNMT3A", "TET2", "ASXL1"))
    ),

    # ---- phenotype rule (PLACEHOLDER thresholds - lock with clinicians) ------
    phenotype = list(
      min_followup_day = 45,
      early_window   = c(0, 28),
      recover_thr    = 1.0, recover_run = 3,
      dip_window     = c(21, 60), dip_thr = 1.0, dip_run = 2,
      persist_window = c(45, 75), persist_thr = 0.5, persist_min_obs = 3
    ),

    # ---- model -------------------------------------------------------------------
    model = list(
      time = list(t0 = -6, tmax = 90, sim_days = seq(-5, 90, by = 1)),
      solver = list(rtol = 1e-6, atol = 1e-9, hmax = 0.5),
      obs = list(crp_scale = 100, fer_scale = 2000),
      # structural parameters (prior centres; some are estimated, see fit$struct_est)
      struct = c(
        a = 0.15, kel_ld = 0.45, klds = 0.04, kldp = 12,
        r_ld = 0.5, r_inf = 0.5, delta = 0.3,
        q = 0.10, emax = 0.9, ic50 = 0.7, hill = 2,
        k_out = 0.2, k_m = 0.2, t_pk = 7, h_shape = 2.5,
        ktr_n = 0.67, kcirc_n = 2.0, gamma_n = 0.5,
        ktr_p = 0.5, kcirc_p = 0.12, gamma_p = 0.4, e_plt = 0.7,
        gcsf_fold = 2.5, tpo_fold = 1.8,
        toci_red = 0.5, ster_red = 0.6, ster_boost = 1.0
      ),
      # prior SD (working scale) for structural parameters that are estimated
      struct_prior_sd = 0.7,
      # population distribution of random effects (eta scale): prior centres
      re_prior = list(
        mu    = c(rho = 0, g = log(0.4), f0 = qlogis(0.02)),
        omega = c(rho = 0.6, g = 0.7, f0 = 2.0)
      ),
      # covariate models on random-effect means (columns built in io_read.R);
      # log_rdw_c and mcv_c are used only when the CBC has RDW / MCV
      covariates = list(
        rho = c("prior_lines_c", "prior_asct", "prior_bendamustine", "log_anc_base_c"),
        g   = c("costim_cd28"),
        f0  = c("age_c", "prior_lines_c", "prior_asct", "log_rdw_c", "mcv_c")
      ),
      vaf_prior_sd = list(pre_LD = 0.5, post_bridging = 0.7, apheresis = 1.2,
                          negative = 1.5),
      vaf_noise_sd = 0.15,
      # log-scale SD of a VAF measured on the day the model predicts it
      # (assay noise only; no timing mismatch as for the pre-treatment sample)
      vaf_post_sd = 0.25,
      # design values used before any data are seen (steps D1-D3): clone
      # traits and population spread considered plausible for planning
      design = list(
        struct = c(r_ld = 0.3, r_inf = 0.2, delta = 0.03, q = 0.12),
        pop = list(mu = c(rho = -0.3, g = log(0.4), f0 = -0.9),
                   omega = c(rho = 0.7, g = 0.6, f0 = 0.9),
                   p_ch = 0.45, f0_non_ch = 0.001),
        base = list(anc = 3.0, plt = 185, crp = 8, ferritin = 600),
        # D3/D3b recalibrate mu_rho and mu_g so that, under the clinical
        # process below, the phenotype mix matches the target
        target = c(quick = 0.50, second_wave = 0.40, persistent = 0.10),
        calib = list(n = 1500, max_iter = 8, tol = 0.02),
        # supportive care and sampling applied to virtual patients (D3/D3b);
        # same rules and rates as the main TestData scenario
        clinical = list(
          events = list(p_gcsf = 0.6, p_gcsf_second = 0.5, p_toci_g2 = 0.8, p_toci_g1 = 0.2,
                        p_steroid_g3 = 0.7, p_steroid_g2 = 0.15, p_anakinra_g3 = 0.2,
                        p_tpo_persistent = 0.3, p_cd34_persistent = 0.3,
                        p_relapse = 0.15, p_death = 0.05),
          crs_cut = c(0.3, 1.2, 2.4, 3.6),
          sampling = list(sparse_frac = 0.30, sparse_after = 30, missing = 0.08, inflam_missing = 0.05))
      )
    ),

    # ---- fitting -----------------------------------------------------------------
    fit = list(
      models = c("M0", "M1", "M2"),
      struct_est = list(M0 = character(0),
                        M1 = c("q", "emax"),
                        M2 = c("q", "emax", "r_ld")),
      # clone traits that trade off with f0 are estimated by profile likelihood
      # (whole population refitted at every grid point)
      profile = list(delta = c(0.03, 0.1, 0.3, 0.8), r_inf = c(0.2, 0.5, 1.0)),
      mech_draws = 200L,           # posterior draws per patient for P(mechanism)
      its_max_iter = 10L, its_tol = 0.02,
      struct_every = 2L,           # update structural parameters every k ITS iterations
      map_maxit = 200L,
      sigma_init = c(anc = 0.25, plt = 0.2, crp = 0.35, ferritin = 0.3)
    ),

    # ---- step settings -------------------------------------------------------------
    steps = list(
      D2 = list(n_rho = 30, n_f0 = 30, rho_range = c(0.15, 2.5), f0_range = c(0, 0.7),
                g_levels = c(low = 0.2, mid = 0.6, high = 1.3),
                delta_levels = c(low_output = 0.03, functional = 0.9)),
      D3 = list(n_virtual = 500, scenarios = c("counts", "counts_pre_vaf", "counts_pre_post_vaf"),
                post_vaf_day = 30, vaf_noise_sd = 0.15),
      # D3b: clone traits unknown - profile likelihood over (delta, r_inf) in
      # virtual cohorts under four information scenarios
      D3b = list(n_virtual = 160, n_rep = 2, post_vaf_day = 30, n_pre_subset = 15,
                 scenarios = c("counts", "pre_subset", "pre_all", "paired"),
                 post_rule = list(anc_thr = 1.0, random_frac = 0.25),
                 grid = list(delta = c(0.03, 0.1, 0.3, 0.8), r_inf = c(0.2, 0.5, 1.0)),
                 low_output_max = 0.1, region_ll = 3.0),
      # D3c: what the current data can say - cohorts at the TestData scale
      # re-simulated and refitted with every parameter estimated (D4 procedure),
      # with and without RDW / MCV in the f0 prior
      D3c = list(n_rep = 10, calib_n = 2000, variants = c("indices", "no_indices"),
                 low_output_max = 0.1, region_ll = 3.0),
      D6 = list(n_virtual = 5000, ch_prevalence = c(0.20, 0.35, 0.50),
                n_enrol = c(50, 100, 150, 200, 250, 300, 400), n_rep = 500,
                post_vaf_day = 30, quick_sample_frac = 0.25),
      H  = list(landmark_day = 10, n_draws = 200, n_boot = 500),
      PR = list(n_draws = 1000, intervals = c(0.80, 0.95)),
      P1 = list(window = c(45, 90), min_points = 4, slope_flat = 0.01, level_low = 0.5),
      P2 = list(nadir_window = c(21, 60), peak_window = c(0, 14)),
      P3 = list(window = c(45, 75), plt_low = 50),
      D4 = list(n_examples_per_phenotype = 2),
      D5 = list(n_boot = 1000)
    ),

    # ---- plotting -------------------------------------------------------------------
    plot = list(
      width = 7, height = 4.5, dpi = 300,
      col_phenotype = c(quick = "#3F7D5C", second_wave = "#C9A227", persistent = "#9E2A2B"),
      col_model = c(M0 = "#6B7A8F", M1 = "#2E6B8A", M2 = "#9E2A2B"),
      col_accent = c(crimson = "#9E2A2B", steel = "#2E6B8A", slate = "#6B7A8F",
                     ink = "#1A2433", gold = "#C9A227", green = "#3F7D5C")
    )
  )

  # steps that do not apply to a data source (data-free design steps are run
  # once, under test)
  cfg$skip_steps <- "D3d"
  if (data_source == "real") cfg <- apply_real_data_mode(cfg)
  if (quick) cfg <- apply_quick_mode(cfg)
  cfg
}

# Real-data mode (Aim2_plan.md v6). The extract has 35 patients, counts to
# about day 30, no mutation data, and treatment flags without dates, so:
#   - phenotype: day 0-30 scheme (quick / intermittent / no recovery by d30)
#   - models M0 and M1 only (no clone layer), fitted to data up to day 30
#   - covariates that exist and vary (all patients received a CD28 product)
#   - early prediction evaluated by stratified cross-validation (too few
#     patients for a locked hold-out)
#   - steps that need mutation data or a fitted clone layer are skipped;
#     D3d simulates cohorts at this scale
apply_real_data_mode <- function(cfg) {
  cfg$phenotype <- list(
    scheme = "early", max_day = 30, min_followup_day = 21,
    recover_thr = 1.5, recover_run = 1,
    dip_thr = 1.0, dip_run = 1, dip_window = c(21, 30),
    aplastic_thr = 0.5, aplastic_days = 14,
    labels = list(quick = "Quick", second_wave = "Intermittent", persistent = "No recovery by d30"))
  cfg$data$baseline_window <- c(-7, -5)
  cfg$fit$models <- c("M0", "M1")
  cfg$fit$obs_max_day <- 30
  # lymphodepletion kill and clearance, maturation rates and the timing of the
  # CAR T inflammatory signal are design values in the synthetic data; the
  # real nadir is deeper and longer and CRS peaks earlier, so they are
  # estimated as population parameters (dense day 0-21 counts and CRP), with
  # a wider prior. The marrow lag of inflammation (k_m) stays at its
  # literature value: day 0-30 data cannot bound it (estimated freely it
  # drifts to a ~45-day delay that predicts a collapse after day 30 not seen
  # in the late counts), and P2 then tests the mechanism's lag
  kin <- c("klds", "kldp", "kel_ld", "ktr_n", "ktr_p", "t_pk", "k_out")
  cfg$fit$struct_est <- list(M0 = kin, M1 = c(kin, "q", "emax"))
  cfg$model$struct_prior_sd <- 1.5
  cfg$fit$its_max_iter <- 30L
  cfg$model$covariates <- list(rho = c("prior_lines_c", "prior_asct", "log_anc_base_c"),
                               g = character(0), f0 = character(0))
  cfg$steps$P2$nadir_window <- c(21, 30)
  cfg$steps$P3$window <- c(21, 30)
  cfg$steps$H$mode <- "cv"
  cfg$steps$H$cv_folds <- 5L
  # D3d: cohorts simulated at the real-data scale, fitted as the real data
  cfg$steps$D3d <- list(
    n_rep = 50L, variants = c("inflam", "no_inflam"), n = 35L,
    lab_days = c(-7:21, seq(23, 29, by = 2)), lab_missing = 0.12,
    late_frac = 0.7, late_days = 85:90,
    inflam_days = c(-7:21, 24, 28), inflam_missing = 0.15,
    p_gcsf = 0.35)
  cfg$skip_steps <- c("D2", "D3", "D3b", "D3c", "D5", "PR", "D5b", "D6", "V")
  cfg
}

# Smaller grids and fewer draws for a smoke test.
apply_quick_mode <- function(cfg) {
  cfg$steps$D2$n_rho <- 14
  cfg$steps$D2$n_f0 <- 14
  cfg$steps$D3$n_virtual <- 60
  cfg$model$design$calib$n <- 500
  cfg$steps$D3b$n_virtual <- 80
  cfg$steps$D3b$n_rep <- 1L
  cfg$steps$D3b$grid <- list(delta = c(0.03, 0.8), r_inf = c(0.2, 0.5))
  cfg$steps$D3c$n_rep <- 1L
  cfg$steps$D3c$calib_n <- 500L
  if (!is.null(cfg$steps$D3d)) cfg$steps$D3d$n_rep <- 2L
  cfg$steps$D6$n_virtual <- 800
  cfg$steps$D6$n_rep <- 100
  cfg$steps$D6$n_enrol <- c(100, 200, 300)
  cfg$steps$H$n_draws <- 50
  cfg$steps$H$n_boot <- 100
  cfg$steps$PR$n_draws <- 300
  cfg$fit$its_max_iter <- 3L
  cfg$fit$map_maxit <- 100L
  cfg$fit$profile <- list(delta = c(0.03, 0.3), r_inf = c(0.2, 1.0))
  cfg$fit$mech_draws <- 50L
  cfg
}

# Locate Aim2/ from this file's location (works with Rscript and source()).
aim2_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f)) return(normalizePath(file.path(dirname(f), "..")))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), "..")))
  }
  normalizePath(getwd())
}
