# =============================================================================
# model_inputs.R
# Purpose : Turn a patient's treatment record into the piecewise-constant input
#           segments consumed by marrow_sim().
# Inputs  : lymphodepletion dose scale, the patient's rows of events.csv.
# Outputs : numeric matrix with columns start, end, ld, gcsf, toci, ster, tpo.
#
# Conventions
#   - lymphodepletion runs on days -5, -4, -3  -> interval [-5, -2)
#   - an event on days start..end is active on [start, end + 1)
#   - tocilizumab stays active from the first dose to the end of the horizon
# =============================================================================

input_columns <- c("ld", "gcsf", "toci", "ster", "tpo")

event_to_input <- c(gcsf = "gcsf", tocilizumab = "toci", steroid = "ster", tpo_ra = "tpo")

make_segments <- function(t0, tmax, ld_scale = 1, events = NULL, ld_days = c(-5, -2)) {
  iv <- data.table::data.table(start = ld_days[1], end = ld_days[2], col = "ld", val = ld_scale)

  if (!is.null(events) && nrow(events) > 0) {
    ev <- events[events$event %in% names(event_to_input), ]
    if (nrow(ev) > 0) {
      col <- unname(event_to_input[ev$event])
      end <- ifelse(col == "toci", tmax, ev$end_day + 1)
      iv <- rbind(iv, data.table::data.table(start = ev$start_day, end = end, col = col, val = 1))
    }
  }

  iv <- iv[iv$end > t0 & iv$start < tmax, ]
  bp <- sort(unique(c(t0, tmax, pmin(pmax(c(iv$start, iv$end), t0), tmax))))
  nseg <- length(bp) - 1L
  seg <- matrix(0, nrow = nseg, ncol = 7,
                dimnames = list(NULL, c("start", "end", input_columns)))
  seg[, "start"] <- bp[-length(bp)]
  seg[, "end"]   <- bp[-1]
  mid <- (seg[, "start"] + seg[, "end"]) / 2
  for (i in seq_len(nrow(iv))) {
    on <- mid >= iv$start[i] & mid < iv$end[i]
    seg[on, iv$col[i]] <- pmax(seg[on, iv$col[i]], iv$val[i])
  }
  seg
}

# Inputs for natural-history simulation (lymphodepletion only).
make_segments_ld_only <- function(t0, tmax, ld_scale = 1) {
  make_segments(t0, tmax, ld_scale = ld_scale, events = NULL)
}
