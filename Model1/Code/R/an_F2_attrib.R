# =============================================================================
# an_F2_attrib.R
# Purpose : Step F2 - individual posteriors (adaptive Metropolis inside the
#           search box, chains extended until R-hat < tolerance) under the final
#           model, and the counterfactual readouts per patient: when (a / b /
#           a+b, HSC level; whole-CAR T-phase readout), the route of b
#           (ferritin-driven inflammation vs CAR T - target signal) and who
#           (reserve / clone). Cohort-level expected counts are sums of
#           probabilities. (v1.3)
# Inputs  : ctx (config, data, results F1).
# Outputs : list(summary, draws_tab, cohort, tables)
# =============================================================================

an_F2_attrib <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  fit <- f1$fits[[f1$selected$final]]
  pdats <- prep_segments(build_pdats(d, cfg), fit$pop$struct, cfg)
  mc <- mcmc_all13(pdats, fit, cfg, seed = cfg$runtime$seed)
  tabs <- par_lapply(names(pdats), function(pid)
    attribute_patient13(mc[[pid]]$draws, pdats[[pid]], fit, cfg, seed = cfg$runtime$seed), cfg$runtime$cores)
  draws_tab <- data.table::rbindlist(tabs)
  summ <- data.table::rbindlist(lapply(tabs, summarise_attribution13, config = cfg), fill = TRUE)
  summ[d$phenotype, observed := i.primary, on = "patient_id"]
  diag <- data.table::rbindlist(lapply(names(mc), function(pid) data.table::data.table(
    patient_id = pid, accept = mean(mc[[pid]]$accept), max_rhat = suppressWarnings(max(mc[[pid]]$rhat, na.rm = TRUE)),
    n_iter = mc[[pid]]$n_iter, converged = mc[[pid]]$converged)))
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  rc <- c("ferritin", "cart_target", "both", "joint")
  hc <- c("clone", "reserve", "both", "other")
  s0 <- function(x) sum(x, na.rm = TRUE)
  cohort <- rbind(
    data.table::data.table(readout = "when (HSC level)", class = wc,
                           expected_n = vapply(wc, function(k) s0(summ[[paste0("p_when_", k)]]), numeric(1))),
    data.table::data.table(readout = "when (whole CAR T phase)", class = wc,
                           expected_n = vapply(wc, function(k) s0(summ[[paste0("p_whentot_", k)]]), numeric(1))),
    data.table::data.table(readout = "route of b (b-involved)", class = rc,
                           expected_n = vapply(rc, function(k) s0(summ[[paste0("p_route_", k)]] * summ$p_b_involved),
                                               numeric(1))),
    data.table::data.table(readout = "who (non-quick)", class = hc,
                           expected_n = vapply(hc, function(k) s0(summ[[paste0("p_who_", k)]] * (1 - summ$p_quick)),
                                               numeric(1))))
  # calibration of the posterior phenotype probabilities against the observed labels
  cls <- cfg$phenotype$classes
  P <- as.matrix(summ[, paste0("p_", cls), with = FALSE])
  Y <- vapply(cls, function(k) as.numeric(summ$observed == k), numeric(nrow(summ)))
  brier <- mean(rowSums((P - Y)^2))
  ex <- summ[, .SD[order(-get(paste0("p_", observed[1])))][seq_len(min(2, .N))], by = observed]
  ctx$log(sprintf(paste("  F2: %s, %d patients; when a/b/ab_joint/ab_redundant %.1f/%.1f/%.1f/%.1f; route of b",
                        "ferritin/cart_target/both/joint %.1f/%.1f/%.1f/%.1f; Brier %.2f; MCMC R-hat <= %.2f in %d/%d"),
                  fit$model, nrow(summ), s0(summ$p_when_a), s0(summ$p_when_b), s0(summ$p_when_ab_joint),
                  s0(summ$p_when_ab_redundant), cohort[readout == "route of b (b-involved)"]$expected_n[1],
                  cohort[readout == "route of b (b-involved)"]$expected_n[2], cohort[readout == "route of b (b-involved)"]$expected_n[3],
                  cohort[readout == "route of b (b-involved)"]$expected_n[4], brier, cfg$fit$mcmc$rhat_tol,
                  sum(diag$converged), nrow(diag)))
  list(model = fit$model, summary = summ, draws_tab = draws_tab, cohort = cohort, brier = brier,
       diagnostics = diag, examples = ex$patient_id,
       mcmc = lapply(mc, function(x) x$draws),
       tables = list(patients = summ, cohort = cohort, mcmc_diagnostics = diag,
                     calibration = data.table::data.table(brier_3class = brier)))
}
