# =============================================================================
# an_G_export.R
# Purpose : Step G - hand the fitted parent model to GreenLight (Aim4/GreenLight):
#           population constants for the in-browser module (data/aim4_pop.js),
#           the real cohort as a third roster tab (data/real_cohort.js, local
#           use only, de-identified) with the R engine's per-patient
#           attribution, and reference trajectories from the C++ integrator for
#           the JavaScript consistency test (tools/ref_trajectories.json).
#           Writes only inside Aim4/GreenLight/.
# Inputs  : ctx (config, data, results F1 and F2).
# Outputs : list(tables = list(export))
#
# data/aim4_pop.js:  window.GL_AIM4_POP = {
#   source: "prior" | "fitted", data_source, fitted_at, model,
#   struct: {natural-scale structural parameters, incl. mtt_*, u_br_*, beta_renal},
#   mu: {re: eta-scale mean}, omega: {re: sd}, beta: {re: {covariate: coef}},
#   sigma: {anc, plt, amc},
#   covariates: {age_ref, age_scale, prior_lines_median, log_rdw_mean, mcv_mean},
#   kernel_pop: {mu: [log A, log tp], cov: [[..], [..]], h}, cd34_add }
# =============================================================================

# ---- population object ------------------------------------------------------------------
g_covariate_centres <- function(patients) {
  list(age_ref = 60, age_scale = 10, prior_lines_median = stats::median(patients$prior_lines),
       log_rdw_mean = mean(log(patients$rdw_base)), mcv_mean = mean(patients$mcv_base))
}

g_kernel_pop <- function(drivers, cfg) {
  kp <- fit_kernel_population(drivers, cfg)
  list(mu = unname(kp$mu), cov = unname(lapply(seq_len(2), function(i) unname(kp$cov[i, ]))), h = kp$h)
}

g_pop_object <- function(model, struct, mu, omega, beta, sigma, d, cfg, source, data_source = NA) {
  list(source = source, data_source = data_source, fitted_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
       model = model, struct = as.list(struct), mu = as.list(mu), omega = as.list(omega),
       beta = stats::setNames(lapply(beta[vapply(beta, length, integer(1)) > 0], as.list),
                              names(beta)[vapply(beta, length, integer(1)) > 0]) %||% structure(list(), names = character(0)),
       sigma = as.list(sigma),
       covariates = g_covariate_centres(d$patients), kernel_pop = g_kernel_pop(d$drivers, cfg),
       cd34_add = 0.3)
}

g_prior_pop <- function(d, cfg, model = "MabC_ddr") {
  sp <- model_spec(model)
  g_pop_object(model, default_struct(sp, cfg), cfg$model$re_prior$mu[sp$re], cfg$model$re_prior$omega[sp$re],
               list(), cfg$model$sigma_init, d, cfg, source = "prior")
}

g_fitted_pop <- function(fit, d, cfg) {
  g_pop_object(fit$model, fit$pop$struct, fit$pop$mu, fit$pop$omega, fit$pop$beta, fit$pop$sigma, d, cfg,
               source = "fitted", data_source = cfg$data_source)
}

# ---- real cohort --------------------------------------------------------------------------
# d: read_aim4_data() on the real extract; attrib: optional F2 summary (per patient)
g_real_cohort <- function(d, cfg, attrib = NULL) {
  raw <- suppressWarnings(data.table::as.data.table(
    readxl::read_excel(cfg$paths$raw_xlsx, sheet = "daily_labs", guess_max = 2000)))
  wh <- raw[, .(patient_id = study_id, day = as.numeric(day), wbc, hgb)]
  ph <- d$phenotype
  obs_rows <- function(pid) {
    o <- d$obs[patient_id == pid]
    a <- lapply(seq_len(nrow(o)), function(i) {
      r <- list(day = o$day[i], variable = o$var[i], value = signif(o$value[i], 4), source = "real")
      if (isTRUE(o$cens[i])) r$cens <- TRUE
      if (isTRUE(o$upper[i])) r$upper <- TRUE
      if (!isTRUE(o$keep[i])) r$excluded <- TRUE
      r
    })
    f <- d$inflam[patient_id == pid]
    b <- c(lapply(which(!is.na(f$ferritin)), function(i) list(day = f$day[i], variable = "ferritin", value = f$ferritin[i], source = "real")),
           lapply(which(!is.na(f$crp)), function(i) list(day = f$day[i], variable = "crp", value = signif(f$crp[i], 4), source = "real")))
    w <- wh[patient_id == pid]
    cc <- c(lapply(which(!is.na(w$wbc)), function(i) list(day = w$day[i], variable = "wbc", value = w$wbc[i], source = "real")),
            lapply(which(!is.na(w$hgb)), function(i) list(day = w$day[i], variable = "hgb", value = w$hgb[i], source = "real")))
    out <- c(a, b, cc)
    out[order(vapply(out, `[[`, numeric(1), "day"))]
  }
  lapply(d$patients$patient_id, function(pid) {
    p <- d$patients[patient_id == pid]
    at <- if (!is.null(attrib) && pid %in% attrib$patient_id) as.list(attrib[patient_id == pid]) else NULL
    list(
      patient_id = paste0("RC-", pid), label = pid,
      demographics = list(age = p$age, sex = if (p$female == 1) "F" else "M"),
      product = "Yescarta", costim = "CD28",
      history = list(prior_lines = p$prior_lines, prior_asct = p$prior_asct, prior_bendamustine = 0,
                     months_since_last_chemo = NULL, marrow_involvement = 0),
      bridging = if (p$bridging == "none") "none" else "chemo",
      ld = list(flu_total = 90, cy_total = 1500, dose_reduced = 0), car_dose_frac = 1,
      crs = NULL, crs_max_grade = p$max_crs,
      mutations = list(tested = FALSE, variants = list()),
      observations = obs_rows(pid),
      events = list(),   # no dated drug events in the extract; LD (D-5..-3) and infusion are drawn by the UI
      recorded = list(phenotype_clinical = ph[patient_id == pid]$primary, phenotype_rule = "aim4_primary_d30",
                      last_followup_day = p$last_lab_day, gcsf_after_d14 = p$gcsf_after_d14,
                      tocilizumab = p$tocilizumab, steroids = p$steroids, max_crs = p$max_crs,
                      relapse_day = if (is.na(p$progression_day)) NULL else p$progression_day),
      aim4 = list(egfr = round(p$egfr, 1), bridging_class = p$bridging, bridging_regimen = p$bridging_regimen,
                  apheresis_day = p$apheresis_day, t0 = p$t0,
                  base = list(anc = p$anc_base, plt = p$plt_base, amc = p$amc_base,
                              alc = if ("alc_base" %in% names(p)) p$alc_base else NULL), ferritin_base = p$ferritin_ld,
                  cov = list(age_c = p$age_c, prior_lines_c = p$prior_lines_c, prior_asct = p$prior_asct,
                             log_rdw_c = p$log_rdw_c, mcv_c = p$mcv_c),
                  r_engine = at),
      demo = list(status = "real", checkpoints = c(-5, 0, 7, 10, 14, 21, 28), initialCheckpoint = 6,
                  deidentified = TRUE))
  })
}

# ---- reference trajectories (C++) for the JavaScript port -----------------------------------
g_ref_trajectories <- function(d, cfg, pop) {
  st <- unlist(pop$struct)
  pid <- names(d$drivers)[which.max(vapply(d$drivers, function(x) max(x$v), numeric(1)))]
  cases <- list(
    list(name = "typical_Mab_kernel", model = "Mab", indiv = c(r_s = 0.08, f0 = 0, th_as = 0.2, th_ap = 1.5, th_bs = 0.03, rho = 0.81, a0 = 1),
         design = list(t0 = -31, egfr = 90, bridging = "none", apheresis_day = -27, ld_scale = 1),
         drv = driver_from_kernel(-31, 1.2, 8, cfg)),
    list(name = "clone_ddr_rgo_low_egfr_real_driver", model = "MabC_ddr",
         indiv = c(r_s = 0.05, f0 = 0.3, th_as = 0.5, th_ap = 2.5, th_bs = 0.05, rho = 1.1, a0 = 2),
         design = list(t0 = -31, egfr = 55, bridging = "rgo", apheresis_day = -28, ld_scale = 1),
         drv = d$drivers[[pid]]),
    list(name = "synergy_Mabx_deep", model = "Mabx", indiv = c(r_s = 0.03, f0 = 0, th_as = 1.2, th_ap = 3, th_bs = 0.08, rho = 0.6, a0 = 0.5),
         design = list(t0 = -31, egfr = 120, bridging = "cy", apheresis_day = -25, ld_scale = 0.75),
         drv = driver_from_kernel(-31, 2.0, 10, cfg)),
    list(name = "supportive_gcsf_tpo_cd34", model = "MabC_dta", indiv = c(r_s = 0.06, f0 = 0.2, th_as = 0.3, th_ap = 2, th_bs = 0.04, rho = 0.9, a0 = 1),
         design = list(t0 = -31, egfr = 90, bridging = "none", apheresis_day = -27, ld_scale = 1,
                       gcsf = matrix(c(15, 20), 1), tpo = matrix(c(21, 60), 1), cd34 = data.frame(day = 45, add = 0.3)),
         drv = driver_from_kernel(-31, 1.5, 9, cfg)))
  days <- -31:90
  lapply(cases, function(cs) {
    sp <- model_spec(cs$model)
    s2 <- st
    s2[names(cfg$model$clone_types$ddr)] <- if (sp$clone == "none") c(1, 1, 1) else cfg$model$clone_types[[sp$clone]]
    if (sp$q2 == "Mabx") s2[["eta_syn"]] <- 3
    par <- build_par(cs$indiv, s2, sp)
    seg <- make_segments(cs$design, s2, max(days), cfg)
    m <- a4_sim(par, seg, cs$drv$v, cs$drv$t0, cs$drv$dt, days, track = TRUE)
    des <- cs$design
    des$gcsf <- if (is.null(des$gcsf)) list() else lapply(seq_len(nrow(des$gcsf)), function(i) des$gcsf[i, ])
    des$tpo <- if (is.null(des$tpo)) list() else lapply(seq_len(nrow(des$tpo)), function(i) des$tpo[i, ])
    des$cd34 <- if (is.null(des$cd34)) list() else lapply(seq_len(nrow(des$cd34)), function(i) as.list(des$cd34[i, ]))
    list(name = cs$name, model = cs$model, par = as.list(par), design = des,
         struct = as.list(s2[c("u_br_rgo", "u_br_cy", "beta_renal")]),
         drv = list(t0 = cs$drv$t0, dt = cs$drv$dt, v = cs$drv$v), days = days,
         out = lapply(as.data.frame(m[, c("anc_rel", "plt_rel", "amc_rel", "alc_rel", "S", "f_clone", "M", "I",
                                           "cart", "ag", "z", "gcsf",
                                           paste0("age_", c("anc", "plt", "amc", "alc")),
                                           paste0("gen_", c("anc", "plt", "amc", "alc")))]), identity))
  })
}

# ---- writers -------------------------------------------------------------------------------------
g_write_js <- function(obj, var, path, comment) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  json <- jsonlite::toJSON(obj, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  writeLines(c(paste0("/* ", comment, " */"), paste0("window.", var, " = ", json, ";")), path)
  file.size(path)
}

g_export_all <- function(gl_dir, pop, cohort, refs) {
  out <- data.table::data.table(file = character(0), bytes = numeric(0), note = character(0))
  if (!is.null(pop)) {
    f <- file.path(gl_dir, "data", if (identical(pop$data_source, "test")) "aim4_pop_test.js" else "aim4_pop.js")
    out <- rbind(out, data.table::data.table(file = f, note = paste("population:", pop$source, pop$model),
      bytes = g_write_js(pop, "GL_AIM4_POP", f, paste0("GreenLight · Aim 4 population constants (", pop$source,
                                                       ", ", pop$fitted_at, "). Generated by Aim4/Code step G."))))
  }
  if (!is.null(cohort)) {
    f <- file.path(gl_dir, "data", "real_cohort.js")
    obj <- list(generated = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), source = "ICAHT extract (de-identified), Aim4/Data/derived",
                note = "Real cohort - local use only. Do not publish.", patients = cohort)
    out <- rbind(out, data.table::data.table(file = f, note = paste(length(cohort), "real patients"),
      bytes = g_write_js(obj, "GL_REAL_DATA", f, "GreenLight · real cohort (de-identified, local use only). Generated by Aim4/Code step G.")))
  }
  if (!is.null(refs)) {
    f <- file.path(gl_dir, "tools", "ref_trajectories.json")
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(refs, f, auto_unbox = TRUE, digits = NA)
    out <- rbind(out, data.table::data.table(file = f, bytes = file.size(f), note = paste(length(refs), "C++ reference cases")))
  }
  out
}

an_G_export <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  f2 <- ctx$results$F2
  fit <- f1$fits[[f1$selected$final]]
  pop <- g_fitted_pop(fit, d, cfg)
  # the real cohort and the kernel population always come from the real extract
  dr <- if (cfg$data_source == "real") d else ctx$design()
  cohort <- g_real_cohort(dr, make_config("real", root = cfg$paths$root),
                          attrib = if (cfg$data_source == "real" && !is.null(f2)) f2$summary else NULL)
  # v1.3: written next to the results; GreenLight/data keeps the v1.2 files until its JavaScript
  # port is updated to v1.3 (CAR T - target module, G-CSF release, individual baselines)
  pop$tau <- as.list(fit$pop$tau)
  pop$version <- cfg$version
  tab <- g_export_all(file.path(cfg$paths$results, "greenlight_export"), pop, cohort, g_ref_trajectories(dr, cfg, pop))
  ctx$log("  G: exported ", paste(basename(tab$file), collapse = ", "), " (model ", pop$model, ")")
  list(tables = list(export = tab))
}
