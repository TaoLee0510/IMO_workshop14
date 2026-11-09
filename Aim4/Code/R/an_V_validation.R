# =============================================================================
# an_V_validation.R
# Purpose : Step V (TestData only) - compare the pipeline's answers with the
#           known truth: model selection, individual parameters, and the two
#           attribution readouts (truth classes from the same counterfactuals
#           at the true parameters).
# Inputs  : ctx (config, data with truth, results F1, F2).
# Outputs : list(tables)
# =============================================================================

an_V_validation <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  tr <- d$truth
  if (is.null(tr)) return(list(tables = list()))
  f1 <- ctx$results$F1
  f2 <- ctx$results$F2
  fit <- f1$fits[[f1$selected$final]]
  e <- fit_eta(fit)
  rec <- data.table::rbindlist(lapply(intersect(colnames(e), c("r_s", "th_ap", "th_as", "th_bs")), function(k) {
    tv <- log(pmax(tr[match(rownames(e), patient_id)][[k]], 1e-6))
    data.table::data.table(parameter = k, cor_with_truth = suppressWarnings(stats::cor(e[, k], tv)))
  }))
  mech <- unique(tr$mech)[1]
  tf <- truth_fit(mech, tr, cfg)
  pdats <- prep_segments(build_pdats(d, cfg), tf$pop$struct, cfg)
  truth_cls <- data.table::rbindlist(lapply(names(pdats), function(pid) {
    cf <- counterfactual_draw(truth_eta(tr[patient_id == pid], tf$spec), pdats[[pid]], tf, cfg, 0:30)
    data.table::data.table(patient_id = pid, true_pheno = cf[["actual"]],
                           true_when = when_class(cf[["actual"]], cf[["no_a"]], cf[["no_bS"]], cf[["no_abS"]]),
                           true_who = who_class(cf[["actual"]], cf[["no_clone"]], cf[["reserve"]]))
  }))
  s <- f2$summary[truth_cls, on = "patient_id"]
  s[, p_true_when := vapply(seq_len(.N), function(i) s[[paste0("p_when_", true_when[i])]][i], numeric(1))]
  s[, p_true_who := vapply(seq_len(.N), function(i) if (is.na(true_who[i])) NA_real_ else
    s[[paste0("p_who_", true_who[i])]][i], numeric(1))]
  summ <- data.table::data.table(
    selected_q2 = f1$selected$q2, true_mechanism = mech,
    mean_p_true_when = mean(s$p_true_when), mean_p_true_who = mean(s$p_true_who, na.rm = TRUE),
    carriers = sum(tr$carrier))
  ctx$log(sprintf("  V: selected %s (truth: %s); mean P(true when) %.2f; mean P(true who) %.2f",
                  f1$selected$q2, mech, summ$mean_p_true_when, summ$mean_p_true_who))
  list(summary = summ, recovery = rec, patients = s,
       tables = list(summary = summ, parameter_recovery = rec, patients = s))
}
