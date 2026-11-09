# =============================================================================
# events_approx.R
# Purpose : Treatment days when only "given yes/no" is known. The real-data
#           extract records tocilizumab, steroids and G-CSF after day 14
#           without dates, but the model needs days. The same rules are
#           applied to the simulated cohorts at the real-data scale (D3d), so
#           the design study carries the same approximation.
# Inputs  : flags per patient (tocilizumab, steroid, gcsf_after_d14 as 0/1,
#           crs_max_day), observed ANC (patient_id, day, anc), rule =
#           config$data$approx_events.
# Outputs : events table (patient_id, event, start_day, end_day, dose).
#
# Rules
#   tocilizumab : one dose on the day of maximal CRS (toci_no_crs_day without CRS)
#   steroid     : steroid_days days from the day of maximal CRS
#                 (steroid_no_crs_day without CRS)
#   gcsf        : gcsf_days days from the first observed ANC < gcsf_anc on or
#                 after day gcsf_from (day gcsf_from when none is below)
# =============================================================================

approx_treatment_events <- function(flags, labs, rule) {
  ev <- list()
  for (i in seq_len(nrow(flags))) {
    f <- flags[i]
    crs_day <- if (is.na(f$crs_max_day)) NA_real_ else f$crs_max_day
    if (f$tocilizumab %in% 1) {
      d <- if (is.na(crs_day)) rule$toci_no_crs_day else crs_day
      ev[[length(ev) + 1]] <- data.table::data.table(patient_id = f$patient_id, event = "tocilizumab",
                                                     start_day = d, end_day = d, dose = NA_real_)
    }
    if (f$steroid %in% 1) {
      d <- if (is.na(crs_day)) rule$steroid_no_crs_day else crs_day
      ev[[length(ev) + 1]] <- data.table::data.table(patient_id = f$patient_id, event = "steroid",
                                                     start_day = d, end_day = d + rule$steroid_days - 1,
                                                     dose = NA_real_)
    }
    if (f$gcsf_after_d14 %in% 1) {
      a <- labs[patient_id == f$patient_id & day >= rule$gcsf_from & !is.na(anc) & anc < rule$gcsf_anc]
      d <- if (nrow(a)) min(a$day) else rule$gcsf_from
      ev[[length(ev) + 1]] <- data.table::data.table(patient_id = f$patient_id, event = "gcsf",
                                                     start_day = d, end_day = d + rule$gcsf_days - 1,
                                                     dose = NA_real_)
    }
  }
  if (!length(ev)) return(data.table::data.table(patient_id = character(0), event = character(0),
                                                 start_day = numeric(0), end_day = numeric(0),
                                                 dose = numeric(0)))
  data.table::rbindlist(ev)
}
