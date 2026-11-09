# =============================================================================
# an_F2_attrib.R
# Purpose : Step F2 - individual posteriors (adaptive Metropolis) under the
#           final model (selected Q2 model + DDR-type clone layer) and the two
#           counterfactual readouts per patient: when (a / b / a+b, HSC level;
#           plus the whole-CAR T-phase readout) and who (reserve / clone).
#           Cohort-level expected counts are sums of probabilities.
# Inputs  : ctx (config, data, results F1).
# Outputs : list(summary, draws_tab, cohort, tables)
# =============================================================================

an_F2_attrib <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  fit <- f1$fits[[f1$selected$final]]
  pdats <- build_pdats(d, cfg)
  mc <- mcmc_all(pdats, fit, cfg, seed = cfg$runtime$seed)
  pds <- prep_segments(pdats, fit$pop$struct, cfg)
  tabs <- par_lapply(names(pds), function(pid)
    attribute_patient(mc[[pid]]$draws, pds[[pid]], fit, cfg, seed = cfg$runtime$seed), cfg$runtime$cores)
  draws_tab <- data.table::rbindlist(tabs)
  summ <- data.table::rbindlist(lapply(tabs, summarise_attribution, config = cfg), fill = TRUE)
  summ[d$phenotype, observed := i.primary, on = "patient_id"]
  diag <- data.table::rbindlist(lapply(names(mc), function(pid) data.table::data.table(
    patient_id = pid, accept = mean(mc[[pid]]$accept), max_rhat = suppressWarnings(max(mc[[pid]]$rhat, na.rm = TRUE)))))
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  hc <- c("clone", "reserve", "both", "other")
  cohort <- rbind(
    data.table::data.table(readout = "when (HSC level)", class = wc,
                           expected_n = vapply(wc, function(k) sum(summ[[paste0("p_when_", k)]], na.rm = TRUE), numeric(1))),
    data.table::data.table(readout = "when (whole CAR T phase)", class = wc,
                           expected_n = vapply(wc, function(k) sum(summ[[paste0("p_whentot_", k)]], na.rm = TRUE), numeric(1))),
    data.table::data.table(readout = "who (non-quick)", class = hc,
                           expected_n = vapply(hc, function(k) sum(summ[[paste0("p_who_", k)]] *
                                                                     (1 - summ$p_quick), na.rm = TRUE), numeric(1))))
  # calibration of the posterior phenotype probabilities against the observed labels
  cls <- cfg$phenotype$classes
  P <- as.matrix(summ[, paste0("p_", cls), with = FALSE])
  Y <- vapply(cls, function(k) as.numeric(summ$observed == k), numeric(nrow(summ)))
  brier <- mean(rowSums((P - Y)^2))
  ex <- summ[, .SD[order(-get(paste0("p_", observed[1])))][seq_len(min(2, .N))], by = observed]
  n_na <- sum(is.na(summ$p_when_a))
  ctx$log(sprintf("  F2: %s, %d patients; expected a/b/ab_joint/ab_redundant (HSC level) %.1f/%.1f/%.1f/%.1f%s; Brier %.2f",
                  fit$model, nrow(summ), sum(summ$p_when_a, na.rm = TRUE), sum(summ$p_when_b, na.rm = TRUE),
                  sum(summ$p_when_ab_joint, na.rm = TRUE), sum(summ$p_when_ab_redundant, na.rm = TRUE),
                  if (n_na) sprintf(" (%d patient(s) with NA draws)", n_na) else "", brier))
  list(model = fit$model, summary = summ, draws_tab = draws_tab, cohort = cohort, brier = brier,
       diagnostics = diag, examples = ex$patient_id,
       mcmc = lapply(mc, function(x) x$draws),
       tables = list(patients = summ, cohort = cohort, mcmc_diagnostics = diag,
                     calibration = data.table::data.table(brier_3class = brier)))
}
