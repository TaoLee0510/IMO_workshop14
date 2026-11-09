# =============================================================================
# io_read.R
# Purpose : Read the ICAHT workbook and convert it into the standardised Aim 3
#           tables, applying every data-handling rule of Data_review.md:
#           censored blanks, zeros on the differential, G-CSF upper bounds,
#           platelet-transfusion exclusions, the shifted chemistry columns,
#           pre-treatment set points and the CAR-HEMATOTOX score.
# Inputs  : config$paths$raw_xlsx (sheets daily_labs, patient_level).
# Outputs : list(obs, labs, inflam, patients) as data.tables; write_derived()
#           saves them as CSV in Data/derived/.
# =============================================================================

read_raw_xlsx <- function(config) {
  f <- config$paths$raw_xlsx
  if (!file.exists(f)) stop("raw workbook not found: ", f)
  d <- suppressWarnings(readxl::read_excel(f, sheet = "daily_labs", guess_max = 5000))
  p <- suppressWarnings(readxl::read_excel(f, sheet = "patient_level", guess_max = 5000))
  list(daily = data.table::as.data.table(d), patient = data.table::as.data.table(p))
}

# CKD-EPI 2021 creatinine equation (mL/min/1.73 m2), no race term.
egfr_ckdepi2021 <- function(creat, age, female) {
  k <- ifelse(female, 0.7, 0.9)
  a <- ifelse(female, -0.241, -0.302)
  142 * pmin(creat / k, 1)^a * pmax(creat / k, 1)^(-1.2) * 0.9938^age * ifelse(female, 1.012, 1)
}

# CAR-HEMATOTOX (Rejeski 2021) from values at lymphodepletion. CRP in mg/dL.
car_hematotox <- function(anc, hgb, plt, crp_mgdl, ferritin) {
  (anc < 1.2) + (hgb < 9.0) + ifelse(plt < 76, 2, ifelse(plt <= 175, 1, 0)) +
    (crp_mgdl > 3.0) + ifelse(ferritin > 2000, 2, ifelse(ferritin >= 650, 1, 0))
}

yes_no <- function(x) as.integer(tolower(trimws(as.character(x))) == "yes")

# ---- patient table ------------------------------------------------------------------------------
build_patients <- function(raw, config) {
  p <- raw$patient
  regimen <- tolower(as.character(p$bridging_regimen))
  bridge_chemo <- as.integer(p$bridge_chemo == 1)
  bridge_regimen <- ifelse(bridge_chemo == 0, "none", ifelse(grepl("gem", regimen), "rgo", "oral_cy"))
  female <- p$sex == "Female"
  out <- data.table::data.table(
    patient_id = p$study_id,
    age = as.numeric(p$age), female = as.integer(female),
    prior_lines = as.numeric(p$prior_lines), prior_asct = yes_no(p$prior_auto_sct),
    bridge_chemo = bridge_chemo, bridge_regimen = bridge_regimen,
    apheresis_day = as.numeric(p$apheresis_date_day),
    creat_ld = as.numeric(p$creatinine_ld),
    egfr = egfr_ckdepi2021(as.numeric(p$creatinine_ld), as.numeric(p$age), female),
    ferr_base = as.numeric(p$ferritin_ld), crp_base_mgdl = as.numeric(p$crp_ld),
    anc_ld = as.numeric(p$anc_ld), hgb_ld = as.numeric(p$hgb_ld), plt_ld = as.numeric(p$plt_ld),
    gcsf14 = yes_no(p$need_gcsf_after_d14), toci = yes_no(p$tocilizumab),
    steroids = yes_no(p$steroids), max_crs = as.numeric(p$max_crs),
    day_max_crs = as.numeric(p$days_to_max_crs), peak_expansion = as.numeric(p$peak_expansion),
    label_clinician = as.character(p$cohort_phenotype_label))
  out[, hematotox := car_hematotox(anc_ld, hgb_ld, plt_ld, crp_base_mgdl, ferr_base)]
  out[, prior_lines_c := prior_lines - stats::median(prior_lines)]
  data.table::setorder(out, patient_id)
  out
}

# ---- daily labs: censoring and exclusion rules ---------------------------------------------------
build_labs <- function(raw, patients, config) {
  d <- raw$daily
  rule <- config$data
  labs <- data.table::data.table(patient_id = d$study_id, day = as.numeric(d$day),
                                 wbc = d$wbc, anc = d$anc, alc = d$alc, amc = d$amc,
                                 plt = d$plt, hgb = d$hgb, transfused = as.integer(d$transfused == 1),
                                 crp_mgdl = d$crp, ferritin = d$ferritin)
  data.table::setorder(labs, patient_id, day)
  # blank differential at low WBC: value <= WBC
  labs[, no_diff := is.na(anc) & !is.na(wbc) & wbc < rule$wbc_no_diff]
  # zero on the differential: below the resolution of a 100-cell count
  for (v in c("anc", "alc", "amc")) {
    labs[, (paste0(v, "_zero")) := !is.na(get(v)) & get(v) == 0]
  }
  # likely platelet transfusion
  labs[, plt_next := data.table::shift(plt, type = "lead"), by = patient_id]
  labs[, hgb_next := data.table::shift(hgb, type = "lead"), by = patient_id]
  labs[, day_next := data.table::shift(day, type = "lead"), by = patient_id]
  tx <- rule$plt_tx
  labs[, plt_tx := transfused == 1 & (
    (!is.na(plt) & plt < tx$plt_below) |
      (!is.na(plt_next) & day_next - day <= 2 & plt_next - plt > tx$rise &
         (is.na(hgb_next) | hgb_next - hgb <= tx$hgb_rise_max)))]
  labs[, c("plt_next", "hgb_next", "day_next") := NULL]
  labs[, plt_drop := FALSE]
  for (r in which(labs$plt_tx)) {
    labs[patient_id == labs$patient_id[r] & day > labs$day[r] & day <= labs$day[r] + tx$drop_days,
         plt_drop := TRUE]
  }
  labs[]
}

# Long observation table for the likelihood: one row per (patient, day, lineage)
# value: observed value or bound; ub: latent (drug-free) value <= value;
# lb: latent value >= value
build_obs <- function(labs, patients, config) {
  rule <- config$data
  win <- rule$fit_window
  gc <- patients[gcsf14 == 1, patient_id]
  ster <- patients[steroids == 1, patient_id]
  sw <- rule$steroid_window
  long <- list()
  for (v in rule$lineages) {
    o <- labs[day >= win[1] & day <= win[2],
              .(patient_id, day, value = get(v), ub = FALSE, lb = FALSE, reason = "observed", wbc, no_diff,
                zero = if (v == "plt") FALSE else get(paste0(v, "_zero")), plt_drop)]
    o[, var := v]
    if (v != "plt") {
      o[no_diff == TRUE, `:=`(value = wbc, ub = TRUE, reason = "no_differential")]
      o[zero == TRUE, `:=`(value = pmax(wbc / rule$manual_diff_cells, 0.005), ub = TRUE, reason = "zero_on_differential")]
    }
    if (v == "anc") {
      o[patient_id %in% ster & day >= sw[1] & day <= sw[2] & ub == FALSE,
        `:=`(ub = TRUE, reason = "steroid_upper_bound")]
      o[patient_id %in% gc & day >= rule$gcsf_from_day & ub == FALSE,
        `:=`(ub = TRUE, reason = "gcsf_upper_bound")]
    }
    o[, keep := !is.na(value) & value > 0]
    if (v == "amc") {
      in_sw <- o$patient_id %in% ster & o$day >= sw[1] & o$day <= sw[2]
      o[in_sw & ub == TRUE, `:=`(keep = FALSE, reason = "steroid_uninformative")]
      o[in_sw & ub == FALSE, `:=`(lb = TRUE, reason = "steroid_lower_bound")]
    }
    if (v == "plt") o[plt_drop == TRUE, `:=`(keep = FALSE, reason = "after_plt_transfusion")]
    long[[v]] <- o[, .(patient_id, day, var, value, ub, lb, reason, keep)]
  }
  out <- data.table::rbindlist(long)
  data.table::setorder(out, patient_id, var, day)
  out
}

# Pre-treatment set points (medians of uncensored values). Bridging chemotherapy:
# only values up to the apheresis day + 1; fallback to the whole window.
build_setpoints <- function(obs, patients, config) {
  w <- config$data$setpoint_window
  sp <- lapply(patients$patient_id, function(pid) {
    pr <- patients[patient_id == pid]
    o <- obs[patient_id == pid & keep & !ub & !lb & day >= w[1] & day <= -5]
    if (pr$bridge_chemo == 1) {
      o2 <- o[day <= pr$apheresis_day + 1]
      if (all(config$data$lineages %in% o2$var)) o <- o2
    }
    v <- o[, .(set = stats::median(value)), by = var]
    data.table::data.table(patient_id = pid,
                           anc_set = v[var == "anc", set], plt_set = v[var == "plt", set],
                           amc_set = v[var == "amc", set])
  })
  data.table::rbindlist(sp, fill = TRUE)
}

# ---- one call: everything W0 needs ----------------------------------------------------------------
convert_icaht <- function(config) {
  raw <- read_raw_xlsx(config)
  patients <- build_patients(raw, config)
  labs <- build_labs(raw, patients, config)
  obs <- build_obs(labs, patients, config)
  patients <- merge(patients, build_setpoints(obs, patients, config), by = "patient_id")
  last <- labs[!is.na(plt), .(last_lab_day = max(day)), by = patient_id]
  patients <- merge(patients, last, by = "patient_id")
  inflam <- labs[!is.na(crp_mgdl) | !is.na(ferritin),
                 .(patient_id, day, crp = crp_mgdl * 10, ferritin)]
  if (anyNA(patients[, .(anc_set, plt_set, amc_set, egfr, ferr_base)]))
    stop("missing set point, eGFR or baseline ferritin for some patients")
  list(patients = patients, labs = labs, obs = obs, inflam = inflam)
}

write_derived <- function(d, config) {
  dir.create(config$paths$derived, recursive = TRUE, showWarnings = FALSE)
  for (nm in c("patients", "labs", "obs", "inflam", "phenotype")) {
    if (!is.null(d[[nm]])) data.table::fwrite(d[[nm]], file.path(config$paths$derived, paste0(nm, ".csv")))
  }
  invisible(TRUE)
}

data_hash <- function(config) digest::digest(file = config$paths$raw_xlsx, algo = "sha256")
