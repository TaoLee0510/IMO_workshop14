# =============================================================================
# an_S0_data.R
# Purpose : Step S0 - data overview and phenotype labels. Checks the converted
#           real extract against the Aim 3 data review (19/10/6 and the six rule
#           variants). Does not relate inflammation to outcomes (that is a
#           pre-registered test, T1).
# Inputs  : ctx (config, data).
# Outputs : list(tables, phenotype, counts, acceptance, n)
# =============================================================================

an_S0_data <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  ph <- d$phenotype
  rules <- setdiff(names(ph), "patient_id")
  counts <- data.table::rbindlist(lapply(rules, function(r) {
    tb <- table(factor(ph[[r]], levels = cfg$phenotype$classes))
    data.table::data.table(rule = r, quick = tb[["quick"]], second_wave = tb[["second_wave"]],
                           delayed = tb[["delayed"]])
  }))
  stable <- sum(apply(ph[, rules, with = FALSE], 1, function(x) length(unique(x)) == 1))
  acceptance <- NULL
  if (cfg$data_source == "real") {
    ref <- data.table::data.table(rule = c("primary", "win35", "rejeski21", "sust2", "rec10", "dip05"),
                                  ref_quick = c(19, 17, 19, 18, 19, 23), ref_second_wave = c(10, 12, 10, 8, 14, 6),
                                  ref_delayed = c(6, 6, 6, 9, 2, 6))
    acceptance <- ref[counts, on = "rule"]
    acceptance[, match := quick == ref_quick & second_wave == ref_second_wave & delayed == ref_delayed]
  }
  pts <- d$patients
  dens <- d$obs[var == "anc", .(n_d0_14 = sum(day >= 0 & day <= 14), n_d15_30 = sum(day >= 15 & day <= 30),
                                n_after30 = sum(day > 30), last_day = max(day)), by = patient_id]
  drv <- data.table::rbindlist(lapply(names(d$drivers), function(pid) {
    dr <- d$drivers[[pid]]
    tt <- dr$t0 + (seq_along(dr$v) - 1) * dr$dt
    data.table::data.table(patient_id = pid, driver_peak = max(dr$v), driver_peak_day = tt[which.max(dr$v)],
                           driver_auc_0_30 = trapz(tt[tt >= 0 & tt <= 30], dr$v[tt >= 0 & tt <= 30]))
  }))
  overview <- pts[, .(patient_id, age, female, prior_lines, prior_asct, bridging, egfr = round(egfr, 1),
                      anc_base, plt_base, amc_base, gcsf_after_d14, max_crs)]
  overview <- dens[drv[overview, on = "patient_id"], on = "patient_id"]
  overview[ph, phenotype := i.primary, on = "patient_id"]
  ctx$log("  S0: ", nrow(pts), " patients; primary phenotype ",
          paste(counts[rule == "primary", c(quick, second_wave, delayed)], collapse = "/"),
          "; stable across rules: ", stable)
  list(tables = list(overview = overview, phenotype = ph, phenotype_counts = counts,
                     acceptance = acceptance, censoring = d$obs[, .N, by = .(var, cens, upper, keep, reason)]),
       phenotype = ph, counts = counts, stable = stable, acceptance = acceptance,
       anc = d$obs[var == "anc"], drivers = drv, n = nrow(pts),
       all_match = if (is.null(acceptance)) NA else all(acceptance$match))
}
