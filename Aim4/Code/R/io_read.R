# =============================================================================
# io_read.R
# Purpose : Read Aim 4's standard tables (real: Data/derived, converted from
#           the read-only xlsx; test: TestData, synthetic with truth), check
#           that the real extract has not changed since conversion, and attach
#           per-patient CAR T-phase drivers and observed phenotypes.
# Inputs  : config.
# Outputs : list(patients, obs, inflam, drivers, phenotype, truth (test only))
# =============================================================================

data_dir <- function(config) if (config$data_source == "real") config$paths$derived else config$paths$testdata

# Convert the real extract when the derived tables are missing or stale.
ensure_real_derived <- function(config, log_fun = message) {
  meta_p <- file.path(config$paths$derived, "meta.json")
  sha <- digest::digest(file = config$paths$raw_xlsx, algo = "sha256")
  if (file.exists(meta_p)) {
    meta <- jsonlite::read_json(meta_p)
    if (identical(meta$sha256, sha)) return(invisible(FALSE))
    log_fun("  source extract changed (SHA-256) - re-converting")
  }
  d <- suppressWarnings(convert_real_data(config))
  write_derived(d, config)
  log_fun("  real extract converted: ", nrow(d$patients), " patients")
  invisible(TRUE)
}

read_aim4_data <- function(config) {
  dd <- data_dir(config)
  rd <- function(f) data.table::fread(file.path(dd, f))
  patients <- rd("patients.csv")
  obs <- rd("obs.csv")
  inflam <- rd("inflam.csv")
  truth <- if (file.exists(file.path(dd, "truth.csv"))) rd("truth.csv") else NULL
  for (v in c("cens", "upper", "keep")) obs[, (v) := as.logical(get(v))]
  data.table::setkey(patients, patient_id)
  drivers <- lapply(stats::setNames(patients$patient_id, patients$patient_id), function(pid) {
    f <- inflam[patient_id == pid & !is.na(ferritin)]
    driver_from_ferritin(f$day, f$ferritin, patients[pid]$ferritin_ld, patients[pid]$t0, config)
  })
  gcsf <- if (file.exists(file.path(dd, "gcsf.csv"))) rd("gcsf.csv") else NULL
  if (!"setpoint" %in% names(obs)) obs[, setpoint := FALSE]
  obs[, setpoint := as.logical(setpoint)]
  list(patients = patients, obs = obs, inflam = inflam, drivers = drivers,
       phenotype = phenotype_all_rules(obs, config),
       truth = truth, gcsf = gcsf, base_sd = baseline_within_sd(obs, config))
}

# Pooled within-patient SD of the log set-point values per lineage (v1.3 baseline prior)
baseline_within_sd <- function(obs, config) {
  v <- config$model$baseline$vars
  out <- vapply(v, function(k) {
    x <- obs[var == k & setpoint == TRUE & cens == FALSE, .(lv = log(pmax(value, 1e-3))), by = patient_id]
    g <- x[, .(ss = sum((lv - mean(lv))^2), df = .N - 1L), by = patient_id][df > 0]
    if (!nrow(g) || sum(g$df) < 3) return(0.2)
    sqrt(sum(g$ss) / sum(g$df))
  }, numeric(1))
  stats::setNames(out, v)
}

data_files_hash <- function(config) {
  dd <- data_dir(config)
  digest::digest(vapply(c("patients.csv", "obs.csv", "inflam.csv"),
                        function(f) digest::digest(file = file.path(dd, f)), character(1)))
}

# Design list for one patient (inputs of make_segments).
# gcsf: G-CSF window table (patient_id, start, end) or NULL; until: keep windows that start
# before this day (landmark forecasts use no future G-CSF)
patient_design <- function(prow, gcsf = NULL, until = Inf) {
  g <- NULL
  if (!is.null(gcsf) && nrow(gcsf)) {
    w <- gcsf[patient_id == prow$patient_id & start < until]
    if (nrow(w)) g <- cbind(w$start, pmin(w$end, until))
  }
  list(t0 = prow$t0, egfr = prow$egfr, bridging = prow$bridging, apheresis_day = prow$apheresis_day,
       ld_scale = 1, gcsf = g, tpo = NULL, cd34 = NULL)
}
