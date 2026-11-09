# =============================================================================
# io_convert.R
# Purpose : Convert the de-identified real-data extract (xlsx, read-only, in
#           Aim2/Data/) into Aim 4's standard tables. Applies the data rules of
#           Aim 3's Data_review (section 3): censoring, manual-differential
#           zeros, G-CSF upper bounds and spikes, platelet-transfusion
#           exclusions, CKD-EPI 2021 eGFR from creatinine_ld, pre-treatment
#           set points, bridging regimen classes.
# Inputs  : config$paths$raw_xlsx.
# Outputs : list(patients, obs, inflam, meta) written to Data/derived/ with the
#           SHA-256 of the source file.
#
# obs columns: patient_id, day, var (anc | plt | amc | alc), value, cens (value is an
#   upper bound: left-censored), upper (G-CSF: value bounds the untreated
#   count from above), keep (enters the likelihood), reason
# =============================================================================

# CKD-EPI 2021 (race-free), creatinine in mg/dL.
egfr_ckd_epi_2021 <- function(scr, age, female) {
  k <- ifelse(female, 0.7, 0.9)
  a <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / k, 1)^a * pmax(scr / k, 1)^-1.2 * 0.9938^age * ifelse(female, 1.012, 1)
}

bridging_class <- function(chemo, regimen) {
  r <- tolower(ifelse(is.na(regimen), "", regimen))
  ifelse(chemo != 1, "none",
         ifelse(grepl("gem", r), "rgo",
                ifelse(grepl("cy|cytox|cyclophos", r), "cy", "cy")))
}

yes_no <- function(x) as.integer(tolower(as.character(x)) %in% c("yes", "1", "true"))

convert_real_data <- function(config) {
  f <- config$paths$raw_xlsx
  if (!file.exists(f)) stop("raw extract not found: ", f)
  p <- data.table::as.data.table(readxl::read_excel(f, sheet = "patient_level"))
  d <- suppressWarnings(data.table::as.data.table(
    readxl::read_excel(f, sheet = "daily_labs", guess_max = 2000)))
  dr <- config$data

  # ---- patients --------------------------------------------------------------------
  pt <- p[, .(patient_id = study_id, age, female = as.integer(sex == "Female"),
              prior_lines = as.numeric(prior_lines), prior_asct = yes_no(prior_auto_sct),
              bridging = bridging_class(as.integer(bridge_chemo), bridging_regimen),
              bridging_regimen, apheresis_day = apheresis_date_day,
              creatinine_ld = as.numeric(creatinine_ld),
              ferritin_ld = as.numeric(ferritin_ld), crp_ld = as.numeric(crp_ld),
              anc_ld = as.numeric(anc_ld), hgb_ld = as.numeric(hgb_ld), plt_ld = as.numeric(plt_ld),
              gcsf_after_d14 = yes_no(need_gcsf_after_d14),
              tocilizumab = yes_no(tocilizumab), steroids = yes_no(steroids),
              max_crs = as.numeric(max_crs), peak_expansion = as.numeric(peak_expansion),
              progression_day = as.numeric(progression_date_day), death_day = as.numeric(death_date_day),
              last_lab_day = as.numeric(last_lab_day))]
  pt[, egfr := egfr_ckd_epi_2021(creatinine_ld, age, female == 1)]
  data.table::setkey(pt, patient_id)

  # ---- daily labs -> long table ---------------------------------------------------------
  lab <- d[, .(patient_id = study_id, day = as.numeric(day), wbc, anc, amc, plt, hgb, alc = as.numeric(alc),
               ferritin, crp = crp * dr$crp_mgdl_to_mgl, transfused = as.integer(transfused),
               mcv, rdw)]
  data.table::setkey(lab, patient_id, day)
  lab <- lab[!is.na(day)]

  long <- function(v) {
    x <- lab[, .(patient_id, day, var = v, value = get(v), wbc)]
    x[, `:=`(cens = FALSE, upper = FALSE, keep = TRUE, reason = "")]
    x
  }
  obs <- data.table::rbindlist(lapply(dr$obs_vars, long))

  # rule 1: blank differential with WBC < 0.5 -> value <= WBC (left-censored)
  obs[var %in% c("anc", "amc", "alc") & is.na(value) & !is.na(wbc) & wbc < dr$wbc_no_diff,
      `:=`(value = pmax(wbc, 0.01), cens = TRUE, reason = "no_diff_low_wbc")]
  # rule 3: zero on a (manual) differential -> value < WBC / 100
  obs[var %in% c("anc", "amc", "alc") & !is.na(value) & value == 0 & !is.na(wbc),
      `:=`(value = pmax(wbc / dr$manual_diff_cells, 0.001), cens = TRUE, reason = "manual_diff_zero")]
  obs <- obs[!is.na(value)]

  # rule 11: isolated G-CSF spikes and G-CSF use after day 14 -> upper bounds for ANC
  data.table::setorder(obs, patient_id, var, day)
  obs[var == "anc", spike := {
    x <- value
    s <- rep(FALSE, .N)
    if (.N >= 3) { i <- gcsf_spikes(x); s[i] <- TRUE }
    s
  }, by = patient_id]
  obs[var == "anc" & spike == TRUE, `:=`(upper = TRUE, reason = "gcsf_spike")]
  g_ids <- pt[gcsf_after_d14 == 1, patient_id]
  obs[var == "anc" & patient_id %in% g_ids & day > dr$gcsf_from_day & upper == FALSE,
      `:=`(upper = TRUE, reason = "gcsf_after_d14")]
  obs[, spike := NULL]

  # rule 8: platelet transfusions (type not recorded): PLT < 30 on a transfusion day,
  # or next-day PLT rise > 30 without Hb rise > 0.5 -> drop PLT on the next 3 days
  tx <- lab[transfused == 1, .(patient_id, day, plt, hgb)]
  nxt <- lab[, .(day, day_next = data.table::shift(day, -1), plt_next = data.table::shift(plt, -1),
                 hgb_next = data.table::shift(hgb, -1)), by = patient_id]
  tx <- nxt[tx, on = .(patient_id, day)]
  tx[, plt_tx := (!is.na(plt) & plt < dr$plt_tx$plt_below) |
                 (!is.na(plt_next) & !is.na(plt) & plt_next - plt > dr$plt_tx$rise &
                    (is.na(hgb_next) | is.na(hgb) | hgb_next - hgb <= dr$plt_tx$hgb_rise_max))]
  ptx <- tx[plt_tx == TRUE, .(patient_id, tx_day = day)]
  for (j in seq_len(nrow(ptx))) {
    obs[var == "plt" & patient_id == ptx$patient_id[j] & day > ptx$tx_day[j] &
          day <= ptx$tx_day[j] + dr$plt_tx$drop_days, `:=`(keep = FALSE, reason = "after_plt_tx")]
  }

  # follow-up limits: after progression / death, beyond obs_max_day
  obs[pt, prog := i.progression_day, on = "patient_id"]
  obs[!is.na(prog) & day >= prog, `:=`(keep = FALSE, reason = "after_progression")]
  obs[day > dr$obs_max_day, `:=`(keep = FALSE, reason = "after_obs_max_day")]
  obs[, prog := NULL]

  # ---- set points, red-cell indices, covariates ---------------------------------------------
  sp_window <- function(pid) {
    r <- pt[pid]
    if (r$bridging != "none") c(-Inf, r$apheresis_day) else dr$setpoint_window
  }
  base <- obs[cens == FALSE, {
    w <- sp_window(patient_id[1])
    sel <- day >= w[1] & day <= w[2]
    if (!any(sel)) sel <- day <= -6
    if (!any(sel)) sel <- day == min(day)
    .(base = stats::median(value[sel]), n_base = sum(sel))
  }, by = .(patient_id, var)]
  bw <- data.table::dcast(base, patient_id ~ var, value.var = "base")
  data.table::setnames(bw, dr$obs_vars, paste0(dr$obs_vars, "_base"))
  pt <- bw[pt, on = "patient_id"]
  ri <- lab[day <= -6, .(mcv_base = stats::median(mcv, na.rm = TRUE),
                         rdw_base = stats::median(rdw, na.rm = TRUE)), by = patient_id]
  pt <- ri[pt, on = "patient_id"]
  for (v in c("mcv_base", "rdw_base")) pt[is.na(get(v)), (v) := mean(pt[[v]], na.rm = TRUE)]
  pt[, `:=`(age_c = (age - 60) / 10, prior_lines_c = prior_lines - stats::median(prior_lines),
            log_rdw_c = log(rdw_base) - mean(log(rdw_base)), mcv_c = (mcv_base - mean(mcv_base)) / 10)]
  first_day <- obs[, .(first_day = min(day)), by = patient_id]
  pt[first_day, first_day := i.first_day, on = "patient_id"]
  pt[, t0 := pmin(config$model$t_min, apheresis_day + 1, first_day - 1)]
  pt[, first_day := NULL]

  # ---- inflammation markers ------------------------------------------------------------------
  inflam <- lab[, .(patient_id, day, ferritin, crp)]
  inflam <- inflam[!(is.na(ferritin) & is.na(crp))]

  data.table::setcolorder(pt, "patient_id")
  data.table::setorder(obs, patient_id, var, day)
  list(patients = pt[], obs = obs[, .(patient_id, day, var, value, cens, upper, keep, reason)],
       inflam = inflam[],
       meta = list(source = basename(f), sha256 = digest::digest(file = f, algo = "sha256"),
                   converted = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   n_patients = nrow(pt)))
}

write_derived <- function(d, config) {
  dir.create(config$paths$derived, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(d$patients, file.path(config$paths$derived, "patients.csv"))
  data.table::fwrite(d$obs, file.path(config$paths$derived, "obs.csv"))
  data.table::fwrite(d$inflam, file.path(config$paths$derived, "inflam.csv"))
  jsonlite::write_json(d$meta, file.path(config$paths$derived, "meta.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(TRUE)
}
