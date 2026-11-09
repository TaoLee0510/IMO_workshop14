# =============================================================================
# phenotype.R
# Purpose : The Aim 3 recovery phenotype (quick / second wave / delayed) as one
#           rule that is applied identically to observed ANC series and to model
#           trajectories, plus its sensitivity variants.
# Inputs  : day and ANC vectors (censored values already mapped below 0.5),
#           a rule list (config$phenotype$primary, optionally modified).
# Outputs : class labels; per-patient phenotype table for the observed data.
#
# Primary rule (locked 2026-11-03, Data_review.md 5.2)
#   window     day 0..30
#   spikes     isolated G-CSF jumps (previous < 1.0, value >= 1.5 and >= 2x,
#              next < 1.0 or none) are removed before classification
#   nadir      lowest value in days 0..14
#   recovery   first value >= 1.5 after the nadir (rec_n consecutive values)
#   second wave  any value < 1.0 after recovery (dip_n consecutive, day >= dip_after)
#   delayed    no recovery by the end of the window
# =============================================================================

pheno_rule <- function(config, variant = NULL) {
  r <- config$phenotype$primary
  if (!is.null(variant)) r[names(config$phenotype$variants[[variant]])] <- config$phenotype$variants[[variant]]
  r
}

gcsf_spikes <- function(anc) {
  n <- length(anc)
  sp <- rep(FALSE, n)
  if (n < 2) return(sp)
  for (j in 2:n) {
    if (anc[j] >= 1.5 && anc[j - 1] < 1.0 && anc[j] >= 2 * anc[j - 1] && (j == n || anc[j + 1] < 1.0))
      sp[j] <- TRUE
  }
  sp
}

first_run <- function(ok, n) {
  if (length(ok) < n) return(NA_integer_)
  for (j in seq_len(length(ok) - n + 1)) if (all(ok[j:(j + n - 1)])) return(j)
  NA_integer_
}

classify_anc <- function(day, anc, rule) {
  k <- !is.na(anc) & day >= 0 & day <= rule$window_end
  day <- day[k]
  anc <- anc[k]
  if (!length(anc)) return(NA_character_)
  if (isTRUE(rule$drop_spikes)) {
    sp <- gcsf_spikes(anc)
    day <- day[!sp]
    anc <- anc[!sp]
  }
  early <- day >= rule$nadir_window[1] & day <= rule$nadir_window[2]
  if (!any(early)) return(NA_character_)
  nd <- day[early][which.min(anc[early])]
  post <- day > nd
  r <- first_run(anc[post] >= rule$rec_thr, rule$rec_n)
  if (is.na(r)) return("delayed")
  rec_day <- day[post][r + rule$rec_n - 1]
  aft <- day > rec_day & day >= rule$dip_after
  if (!is.na(first_run(anc[aft] < rule$dip_thr, rule$dip_n))) "second_wave" else "quick"
}

# ANC series for classification: censored blanks (no differential) count as < 0.5.
anc_for_phenotype <- function(labs, pid) {
  x <- labs[patient_id == pid & day >= 0]
  a <- x$anc
  a[x$no_diff] <- pmin(x$wbc[x$no_diff], 0.499)
  list(day = x$day, anc = a)
}

# Observed phenotype of every patient: primary rule and all variants.
phenotype_table <- function(labs, patients, config) {
  vars <- c("primary", names(config$phenotype$variants))
  out <- data.table::data.table(patient_id = patients$patient_id)
  for (v in vars) {
    rule <- pheno_rule(config, if (v == "primary") NULL else v)
    out[, (v) := vapply(patient_id, function(pid) {
      s <- anc_for_phenotype(labs, pid)
      classify_anc(s$day, s$anc, rule)
    }, character(1))]
  }
  out[, stable := apply(.SD, 1, function(r) length(unique(r)) == 1), .SDcols = vars]
  out[, longest_lt05 := vapply(patient_id, function(pid) {
    s <- anc_for_phenotype(labs, pid)
    k <- !is.na(s$anc) & s$day <= 30
    longest_below(s$day[k], s$anc[k], 0.5)
  }, numeric(1))]
  out[, spike_days := vapply(patient_id, function(pid) {
    s <- anc_for_phenotype(labs, pid)
    k <- !is.na(s$anc) & s$day <= 30
    paste(s$day[k][gcsf_spikes(s$anc[k])], collapse = ",")
  }, character(1))]
  out
}

# Longest span (days, first to last draw inclusive) of consecutive draws below thr.
longest_below <- function(day, anc, thr) {
  best <- 0
  st <- NA
  for (j in seq_along(anc)) {
    if (anc[j] < thr) {
      if (is.na(st)) st <- day[j]
      best <- max(best, day[j] - st + 1)
    } else st <- NA
  }
  best
}

# Phenotype of a simulated daily ANC trajectory (absolute scale, x10^9/L).
classify_sim <- function(day, anc, config, variant = NULL) {
  classify_anc(day, anc, pheno_rule(config, variant))
}

pheno_levels <- function(config) config$phenotype$classes
