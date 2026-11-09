# =============================================================================
# io_read.R
# Purpose : Read the Aim 2 data files, validate them against the data
#           contract, and return analysis-ready, standardised tables.
#           This layer reads and checks; it does not analyse.
# Inputs  : config$paths$data, config$data (file names, detection limits,
#           censoring rules, exclusion flags, CH gene classes).
# Outputs : list(patients, obs, base, pheno, vaf, events, post_planned, checks)
#
# Access rules (enforced here)
#   - later/mutations_post.csv is read only by read_post_mutations(), which
#     step D5b calls after verifying the pre-registration.
#   - truth/ is read only by read_truth(), which only step V calls.
# =============================================================================

# labs may also carry anc_upper: upper bound of an ANC that was not reported
# because it was too low to measure (left-censored), with anc left empty.
required_columns <- list(
  patients = c("patient_id", "age", "costim", "prior_lines", "prior_asct",
               "prior_bendamustine", "bridging", "ld_dose_reduced", "crs_grade",
               "marrow_involvement", "iec_hs", "relapse_day", "subsequent_tx_day",
               "death_day", "last_followup_day"),
  labs = c("patient_id", "day", "anc", "plt"),
  inflam = c("patient_id", "day", "crp", "ferritin"),
  events = c("patient_id", "event", "start_day", "end_day"),
  split = c("patient_id", "split"),
  mut_tested = c("patient_id", "sample_day", "sample_timepoint", "lod_vaf"),
  mut_pre = c("patient_id", "gene", "vaf"),
  post_planned = c("patient_id", "sample_day")
)

read_csv_dt <- function(path) {
  data.table::fread(path, na.strings = c("", "NA"), showProgress = FALSE)
}

read_raw <- function(config) {
  f <- config$data$files
  rd <- function(key, required = TRUE) {
    path <- file.path(config$paths$data, f[[key]])
    if (!file.exists(path)) {
      if (required) stop("missing data file: ", path)
      return(NULL)
    }
    read_csv_dt(path)
  }
  list(patients = rd("patients"), labs = rd("labs"), inflam = rd("inflam"),
       events = rd("events"), split = rd("split"),
       mut_tested = rd("mut_tested", FALSE), mut_pre = rd("mut_pre", FALSE),
       post_planned = rd("post_planned", FALSE))
}

# ---- validation ----------------------------------------------------------------------------
validate_data <- function(raw) {
  out <- list()
  add <- function(check, ok, detail = "", severity = "error") {
    out[[length(out) + 1]] <<- data.table::data.table(
      check = check, status = if (ok) "ok" else severity, detail = detail)
  }
  for (k in names(required_columns)) {
    x <- raw[[k]]
    if (is.null(x)) {
      add(paste0(k, ": present"), FALSE, "file absent", "warn")
      next
    }
    miss <- setdiff(required_columns[[k]], names(x))
    add(paste0(k, ": required columns"), !length(miss), paste(miss, collapse = ", "))
  }
  ids <- raw$patients$patient_id
  add("patients: unique patient_id", !anyDuplicated(ids))
  for (k in c("labs", "inflam", "events", "split")) {
    unk <- setdiff(unique(raw[[k]]$patient_id), ids)
    add(paste0(k, ": ids known"), !length(unk), paste(head(unk, 5), collapse = ", "))
  }
  add("split: covers all patients", setequal(raw$split$patient_id, ids))
  add("split: values train/test", all(raw$split$split %in% c("train", "test")))
  add("labs: day range -30..400", all(between(raw$labs$day, -30, 400)), severity = "warn")
  add("labs: non-negative counts", all(raw$labs$anc >= 0 & raw$labs$plt >= 0, na.rm = TRUE))
  add("labs: every patient has a pre-LD value",
      all(ids %in% raw$labs[day <= -3, unique(patient_id)]), severity = "warn")
  if (!is.null(raw$mut_pre) && !is.null(raw$mut_tested))
    add("mutations: variants only in tested patients",
        all(raw$mut_pre$patient_id %in% raw$mut_tested$patient_id))
  if (!is.null(raw$post_planned))
    add("post samples: ids known", all(raw$post_planned$patient_id %in% ids))
  data.table::rbindlist(out)
}

# ---- standardisation -------------------------------------------------------------------------
baseline_values <- function(labs, inflam, window = c(-Inf, -5)) {
  pre <- function(dt, var) {
    x <- dt[!is.na(get(var))]
    b <- x[day >= window[1] & day <= window[2], .(v = stats::median(get(var))), by = patient_id]
    alt <- x[day <= 0 & !patient_id %in% b$patient_id, .(v = get(var)[which.min(day)]), by = patient_id]
    data.table::setnames(rbind(b, alt), "v", var)
  }
  # red-cell indices are optional CBC columns
  red <- intersect(c("mcv", "rdw"), names(labs))
  b <- Reduce(function(a, b) merge(a, b, by = "patient_id", all = TRUE),
              c(list(pre(labs, "anc"), pre(labs, "plt"), pre(labs, "hgb")),
                lapply(red, function(v) pre(labs, v)),
                list(pre(inflam, "crp"), pre(inflam, "ferritin"))))
  b[, anc := pmax(anc, 0.05)]
  b[, plt := pmax(plt, 5)]
  b
}

censor_days <- function(patients, events, config) {
  cd <- patients[, .(patient_id, relapse_day, subsequent_tx_day, death_day)]
  boost <- events[event == "cd34_boost", .(cd34_day = min(start_day)), by = patient_id]
  cd <- merge(cd, boost, by = "patient_id", all.x = TRUE)
  keys <- c(cd34_boost = "cd34_day", relapse = "relapse_day",
            subsequent_tx = "subsequent_tx_day", death = "death_day")
  use <- keys[config$data$censor$drop_after]
  cd[, censor_day := do.call(pmin, c(lapply(use, function(k) data.table::fifelse(is.na(get(k)), Inf, as.numeric(get(k))))))]
  cd[, censor_reason := {
    m <- as.matrix(.SD)
    m[is.na(m)] <- Inf
    r <- names(use)[max.col(-m, ties.method = "first")]
    ifelse(is.finite(censor_day), r, NA_character_)
  }, .SDcols = unname(use)]
  cd[, .(patient_id, censor_day, censor_reason)]
}

build_obs <- function(labs, inflam, events, cens, config) {
  lod <- config$data$lod
  long <- rbind(
    labs[, .(patient_id, day, var = "anc", value = anc)],
    labs[, .(patient_id, day, var = "plt", value = plt)],
    inflam[, .(patient_id, day, var = "crp", value = crp)],
    inflam[, .(patient_id, day, var = "ferritin", value = ferritin)])
  long <- long[!is.na(value)][, upper := FALSE]
  # ANC not reported because the count was too low for a differential (e.g.
  # WBC < 0.5): left-censored at the bound given in the optional anc_upper
  if ("anc_upper" %in% names(labs))
    long <- rbind(long, labs[is.na(anc) & !is.na(anc_upper),
                             .(patient_id, day, var = "anc", value = anc_upper, upper = TRUE)])
  long <- merge(long, cens, by = "patient_id", all.x = TRUE)
  long[, keep := day <= censor_day]
  long[, reason := data.table::fifelse(keep, NA_character_, paste0("after_", censor_reason))]

  # platelet transfusion: drop PLT for N days after each transfusion
  ntx <- config$data$censor$plt_after_transfusion_days
  tx <- events[event == "plt_tx", .(patient_id, tx_day = start_day)]
  if (nrow(tx)) {
    hit <- tx[long[var == "plt"], on = "patient_id", allow.cartesian = TRUE, nomatch = NULL][
      day > tx_day & day <= tx_day + ntx, unique(paste(patient_id, day))]
    long[var == "plt" & paste(patient_id, day) %in% hit & keep == TRUE,
         `:=`(keep = FALSE, reason = "plt_transfusion")]
  }
  # tocilizumab: CRP no longer reflects inflammation
  if (identical(config$data$censor$crp_after_tocilizumab, "censor")) {
    toci <- events[event == "tocilizumab", .(toci_day = min(start_day)), by = patient_id]
    long <- merge(long, toci, by = "patient_id", all.x = TRUE)
    long[var == "crp" & !is.na(toci_day) & day > toci_day & keep == TRUE,
         `:=`(keep = FALSE, reason = "after_tocilizumab")]
    long[, toci_day := NULL]
  }
  lodv <- c(anc = lod$anc, plt = lod$plt, crp = lod$crp, ferritin = 1)
  long[, cens := upper | value <= lodv[var]]
  long[cens == TRUE, value := pmax(value, lodv[var])]
  data.table::setorder(long, patient_id, var, day)
  long[, .(patient_id, day, var, value, cens, keep, reason)]
}

build_vaf <- function(raw, patients, config) {
  if (is.null(raw$mut_tested) || !nrow(raw$mut_tested)) return(NULL)
  cls <- config$data$ch_genes
  gene_class <- function(g) ifelse(g %in% cls$DDR, "DDR", ifelse(g %in% cls$DTA, "DTA", "other"))
  tested <- merge(raw$mut_tested, patients[, .(patient_id, bridging)], by = "patient_id", all.x = TRUE)
  mp <- if (is.null(raw$mut_pre)) data.table::data.table(patient_id = character(0), gene = character(0),
                                                         vaf = numeric(0)) else raw$mut_pre
  dom <- mp[order(-vaf), .SD[1], by = patient_id][, .(patient_id, gene, vaf, n_variants = 0L)]
  nv <- mp[, .N, by = patient_id]
  dom[nv, n_variants := i.N, on = "patient_id"]
  v <- merge(tested, dom, by = "patient_id", all.x = TRUE)
  v[, status := ifelse(is.na(vaf), "negative", "positive")]
  v[, gene_class := ifelse(status == "positive", gene_class(gene), "none")]
  v[, f_measured := ifelse(status == "positive", pmin(2 * vaf, 1), NA_real_)]
  sds <- config$model$vaf_prior_sd
  v[, timepoint_sd := vapply(sample_timepoint, function(tp)
    if (!is.null(sds[[tp]])) sds[[tp]] else sds$apheresis, numeric(1))]
  # apheresis samples followed by bridging chemotherapy are less informative about f0
  v[sample_timepoint == "apheresis" & bridging == "chemo", timepoint_sd := timepoint_sd * 1.5]
  v[status == "negative", timepoint_sd := config$model$vaf_prior_sd$negative]
  v[]
}

standardise <- function(raw, config) {
  p <- data.table::copy(raw$patients)
  p <- merge(p, raw$split[, .(patient_id, split)], by = "patient_id", all.x = TRUE)
  events <- if (is.null(raw$events)) data.table::data.table(patient_id = character(0), event = character(0),
                                                             start_day = numeric(0), end_day = numeric(0))
            else raw$events
  base <- baseline_values(raw$labs, raw$inflam, config$data$baseline_window)
  p <- merge(p, base[, .(patient_id, anc_base = anc, plt_base = plt, hgb_base = hgb,
                         crp_base = crp, fer_base = ferritin)], by = "patient_id", all.x = TRUE)
  # a patient without any pre-treatment value of a marker gets the cohort
  # median (recorded in base_imputed)
  bcols <- c("anc_base", "plt_base", "hgb_base", "crp_base", "fer_base")
  p[, base_imputed := apply(is.na(as.matrix(.SD)), 1, function(m) paste(sub("_base$", "", bcols[m]), collapse = ",")),
    .SDcols = bcols]
  for (v in bcols) {
    miss <- which(is.na(p[[v]]))
    if (length(miss)) data.table::set(p, miss, v, stats::median(p[[v]], na.rm = TRUE))
  }
  # centred covariates (centring uses all patients and no outcomes)
  p[, `:=`(age_c = age - mean(age, na.rm = TRUE),
           prior_lines_c = prior_lines - 2,
           costim_cd28 = as.integer(costim == "CD28"),
           log_anc_base_c = log(anc_base) - mean(log(anc_base), na.rm = TRUE),
           ld_scale = ifelse(ld_dose_reduced == 1, 0.75, 1))]
  # red-cell indices: proxies of clone burden for the f0 prior (log RDW; MCV
  # per 10 fL); patients without a baseline value are imputed at the mean
  centre_impute <- function(x) { x <- x - mean(x, na.rm = TRUE); x[is.na(x)] <- 0; x }
  for (v in intersect(c("mcv", "rdw"), names(base))) {
    p[, (paste0(v, "_base")) := base[[v]][match(p$patient_id, base$patient_id)]]
  }
  if (any(!is.na(p$rdw_base))) p[, log_rdw_c := centre_impute(log(rdw_base))]
  if (any(!is.na(p$mcv_base))) p[, mcv_c := centre_impute(mcv_base / 10)]
  flags <- intersect(config$data$exclude_flags, names(p))
  p[, excluded := Reduce(`|`, lapply(flags, function(f) get(f) %in% 1))]
  cens <- censor_days(p, events, config)
  p <- merge(p, cens, by = "patient_id")
  obs <- build_obs(raw$labs, raw$inflam, events, cens, config)
  pheno <- obs[var == "anc" & keep == TRUE,
               .(phenotype = classify_phenotype(day, value, config$phenotype)), by = patient_id]
  pheno <- merge(p[, .(patient_id, split, excluded)], pheno, by = "patient_id", all.x = TRUE)
  if ("phenotype_clinical" %in% names(p))
    pheno <- merge(pheno, p[, .(patient_id, phenotype_clinical)], by = "patient_id")
  list(patients = p, obs = obs, pheno = pheno, events = events,
       vaf = build_vaf(raw, p, config), post_planned = raw$post_planned)
}

# ---- public API ------------------------------------------------------------------------------------
read_aim2_data <- function(config) {
  raw <- read_raw(config)
  checks <- validate_data(raw)
  if (any(checks$status == "error"))
    stop("data validation failed:\n", paste(checks[status == "error", paste(check, detail)], collapse = "\n"))
  d <- standardise(raw, config)
  d$checks <- checks
  d$data_hash <- data_files_hash(config)
  d
}

# Restrict every table to one split (hold-out protection).
subset_split <- function(d, which = "train") {
  split_ids <- d$patients[split == which, patient_id]
  out <- d
  out$patients <- d$patients[patient_id %in% split_ids]
  out$obs <- d$obs[patient_id %in% split_ids]
  out$pheno <- d$pheno[patient_id %in% split_ids]
  out$events <- d$events[patient_id %in% split_ids]
  if (!is.null(d$vaf)) out$vaf <- d$vaf[patient_id %in% split_ids]
  out$split_used <- which
  out
}

# Restrict every table to the patients `ids` (cross-validation folds).
subset_ids <- function(d, ids) {
  out <- d
  out$patients <- d$patients[patient_id %in% ids]
  out$obs <- d$obs[patient_id %in% ids]
  out$pheno <- d$pheno[patient_id %in% ids]
  out$events <- d$events[patient_id %in% ids]
  if (!is.null(d$vaf)) out$vaf <- d$vaf[patient_id %in% ids]
  out
}

# Analysis set: one split, excluded patients removed.
analysis_set <- function(d, which = "train") {
  s <- subset_split(d, which)
  keep_ids <- s$patients[excluded == FALSE, patient_id]
  s$patients <- s$patients[patient_id %in% keep_ids]
  s$obs <- s$obs[patient_id %in% keep_ids]
  s$pheno <- s$pheno[patient_id %in% keep_ids]
  s$events <- s$events[patient_id %in% keep_ids]
  if (!is.null(s$vaf)) s$vaf <- s$vaf[patient_id %in% keep_ids]
  s
}

# Post-treatment mutations: only for step D5b, after pre-registration is verified.
read_post_mutations <- function(config) {
  path <- file.path(config$paths$data, config$data$files[["mut_post"]])
  if (!file.exists(path)) return(NULL)
  read_csv_dt(path)
}

# Ground truth: only for validation step V on synthetic data.
read_truth <- function(config) {
  if (config$data_source == "real") stop("truth is not available for real data")
  dir <- file.path(config$paths$data, "truth")
  list(patient = read_csv_dt(file.path(dir, "true_params_patient.csv")),
       traj = read_csv_dt(file.path(dir, "true_trajectories.csv")),
       model = jsonlite::read_json(file.path(dir, "true_model.json"), simplifyVector = TRUE))
}

# Hash of the data files used by the analysis (truth/ and later/ excluded).
data_files_hash <- function(config) {
  f <- config$data$files
  f <- f[names(f) != "mut_post"]
  paths <- file.path(config$paths$data, f)
  paths <- paths[file.exists(paths)]
  digest::digest(vapply(paths, function(x) digest::digest(file = x, algo = "md5"), character(1)))
}
