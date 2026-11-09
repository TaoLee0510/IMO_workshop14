# =============================================================================
# phenotype.R
# Purpose : Classify an ANC trajectory into quick / second_wave / persistent.
#           Thresholds come from config$phenotype (placeholders until locked
#           with the clinicians); the same rule is used for observed data,
#           simulated data and model predictions.
# Inputs  : day and ANC vectors for one patient (x10^9/L).
# Outputs : character label, or NA when follow-up is too short.
#
# Two schemes (config$phenotype$scheme)
#   "late"  : the planned three classes, which need counts to day 75
#   "early" : the day 0-30 version for data without late counts (Rejeski et
#             al. 2021 neutrophil recovery patterns). The codes stay the same,
#             with these meanings (display names in config$phenotype$labels):
#               quick       = sustained recovery, no later dip
#               second_wave = "intermittent": recovery, then a dip after day 21
#               persistent  = "no recovery by day 30", or aplastic (ANC < 0.5
#                             over >= aplastic_days)
#
# Late rule (applied to observations with day >= 0)
#   unclassifiable : last observation before min_followup_day, or fewer than
#                    persist_min_obs values inside persist_window
#   persistent     : median ANC inside persist_window (late, after the usual
#                    second-wave recovery) below persist_thr
#   early recovery : run of >= recover_run values >= recover_thr starting in early_window
#   second_wave    : early recovery, then a run of >= dip_run values < dip_thr
#                    starting after that recovery and inside dip_window
#   quick          : everything else
#
# Early rule (observations with 0 <= day <= max_day)
#   unclassifiable : last observation before min_followup_day
#   quick          : ANC never below recover_thr
#   persistent     : no run of >= recover_run values >= recover_thr after the
#                    nadir (lowest value), or a run of values < aplastic_thr
#                    spanning >= aplastic_days days
#   second_wave    : after that recovery, a run of >= dip_run values < dip_thr
#                    inside dip_window
#   quick          : everything else
#   Recovery is counted from the nadir, because the day 0-2 counts are often
#   still high before lymphodepletion takes effect.
# =============================================================================

phenotype_levels <- c("quick", "second_wave", "persistent")

# Start indices of runs of TRUE with length >= k.
run_starts <- function(x, k) {
  if (length(x) < k) return(integer(0))
  r <- rle(x)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  starts[r$values & r$lengths >= k]
}

classify_phenotype <- function(day, anc, rule) {
  if (identical(rule$scheme, "early")) return(classify_phenotype_early(day, anc, rule))
  keep <- !is.na(anc) & day >= 0
  day <- day[keep]
  anc <- anc[keep]
  if (!length(day) || max(day) < rule$min_followup_day) return(NA_character_)
  o <- order(day)
  day <- day[o]
  anc <- anc[o]

  # persistently low late (checked first)
  inw <- day >= rule$persist_window[1] & day <= rule$persist_window[2]
  if (sum(inw) < rule$persist_min_obs) return(NA_character_)
  if (stats::median(anc[inw]) < rule$persist_thr) return("persistent")

  # early recovery followed by a later dip
  rec <- anc >= rule$recover_thr
  rec_starts <- run_starts(rec, rule$recover_run)
  rec_starts <- rec_starts[day[rec_starts] >= rule$early_window[1] &
                           day[rec_starts] <= rule$early_window[2]]
  if (length(rec_starts) > 0) {
    r <- rle(rec[rec_starts[1]:length(rec)])
    rec_end <- rec_starts[1] + r$lengths[1] - 1L
    after <- seq_along(day) > rec_end
    dip <- anc < rule$dip_thr & after &
           day >= rule$dip_window[1] & day <= rule$dip_window[2]
    if (length(run_starts(dip, rule$dip_run)) > 0) return("second_wave")
  }
  "quick"
}

classify_phenotype_early <- function(day, anc, rule) {
  keep <- !is.na(anc) & day >= 0 & day <= rule$max_day
  day <- day[keep]
  anc <- anc[keep]
  if (!length(day) || max(day) < rule$min_followup_day) return(NA_character_)
  o <- order(day)
  day <- day[o]
  anc <- anc[o]

  # aplastic: values below aplastic_thr spanning aplastic_days or more
  low <- rle(anc < rule$aplastic_thr)
  ends <- cumsum(low$lengths)
  starts <- ends - low$lengths + 1L
  span <- day[ends] - day[starts] + 1
  if (any(low$values & span >= rule$aplastic_days)) return("persistent")

  nadir <- which.min(anc)
  if (anc[nadir] >= rule$recover_thr) return("quick")
  rec_starts <- run_starts(anc >= rule$recover_thr & seq_along(day) > nadir, rule$recover_run)
  if (!length(rec_starts)) return("persistent")
  after <- seq_along(day) >= rec_starts[1] + rule$recover_run
  dip <- after & anc < rule$dip_thr & day >= rule$dip_window[1] & day <= rule$dip_window[2]
  if (length(run_starts(dip, rule$dip_run)) > 0) "second_wave" else "quick"
}

# Display name of each phenotype code under the configured scheme.
phenotype_display <- function(rule) {
  if (is.null(rule$labels)) c(quick = "Quick", second_wave = "Second wave", persistent = "Persistent")
  else unlist(rule$labels)
}

# Vectorised helper over a long data.table (patient_id, day, anc).
classify_phenotype_dt <- function(dt, rule, id_col = "patient_id") {
  dt[, .(phenotype = classify_phenotype(day, anc, rule)), by = id_col]
}

# Ordinal severity used for "buffer vs aggravate" comparisons.
phenotype_severity <- function(x) {
  match(x, phenotype_levels) - 1L
}
