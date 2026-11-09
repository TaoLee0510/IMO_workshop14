# =============================================================================
# an_V_validation.R
# Purpose : V - pipeline validation on synthetic data only. Compares what the
#           pipeline estimated (training split) with the generating truth.
#           This is the only step allowed to read truth/.
# Inputs  : read_truth(config), ctx$data(), ctx$results$D4.
# Outputs : list(recovery, struct, phenotype_confusion, mechanism_confusion,
#           mechanism_probability, tables)
# =============================================================================

an_V_validation <- function(ctx) {
  cfg <- ctx$config
  tru <- read_truth(cfg)
  d4 <- ctx$results$D4
  tr <- analysis_set(ctx$data(), "train")

  x <- merge(d4$fits$M2$indiv, tru$patient, by = "patient_id")
  rec <- rbind(
    x[, .(parameter = "log rho", n = .N, cor = stats::cor(eta_rho, log(rho)))],
    x[, .(parameter = "log g", n = .N, cor = stats::cor(eta_g, log(g)))],
    x[ch == 1, .(parameter = "logit f0 (CH carriers)", n = .N, cor = stats::cor(eta_f0, stats::qlogis(f0)))],
    x[ch == 1 & g > stats::median(g), .(parameter = "logit f0 (CH, high g)", n = .N,
                                        cor = stats::cor(eta_f0, stats::qlogis(f0)))])
  rec[, cor := round(cor, 3)]

  est <- intersect(c(unique(unlist(cfg$fit$struct_est)), names(cfg$fit$profile)), names(tru$model$struct))
  struct <- data.table::data.table(parameter = est, truth = round(unlist(tru$model$struct[est]), 4),
                                   prior_centre = round(cfg$model$struct[est], 4),
                                   estimate_M2 = round(d4$fits$M2$pop$struct[est], 4))

  ph <- merge(tr$pheno[, .(patient_id, phenotype)], tru$patient[, .(patient_id, phenotype_true, mechanism_true)],
              by = "patient_id")
  conf_ph <- ph[, .N, by = .(phenotype_true, phenotype)]
  acc_ph <- ph[!is.na(phenotype) & !is.na(phenotype_true), round(mean(phenotype == phenotype_true), 3)]

  mech <- merge(d4$inferred[, .(patient_id, mechanism)], ph, by = "patient_id")
  conf_mech <- mech[phenotype_true == "persistent", .N, by = .(mechanism_true, mechanism)]
  # posterior P(H3) (observed-persistent patients) against the true mechanism
  mp <- merge(d4$mech_prob, ph[, .(patient_id, mechanism_true)], by = "patient_id")
  prob_mech <- mp[mechanism_true %in% c("H1", "H3"),
                  .(n = .N, map_called_h3 = sum(mechanism_map == "H3", na.rm = TRUE),
                    mean_p_h3 = round(mean(p_h3), 3), sum_p_h3 = round(sum(p_h3), 1)), by = mechanism_true]

  model_sel <- data.table::data.table(truth_model = "M2", chosen_by_BIC = d4$chosen,
                                      correct = d4$chosen == "M2")
  p4 <- d4$P4[term != "(Intercept)", .(term, estimate, ci_low, ci_high,
                                       truth_sign = "+", sign_recovered = estimate > 0)]

  list(recovery = rec, struct = struct, phenotype_confusion = conf_ph, phenotype_accuracy = acc_ph,
       mechanism_confusion = conf_mech, mechanism_probability = prob_mech, scatter = x,
       tables = list(parameter_recovery = rec, structural = struct, phenotype_confusion = conf_ph,
                     persistent_mechanism_confusion = conf_mech, mechanism_probability = prob_mech,
                     model_selection = model_sel,
                     P4_signs = p4),
       key_numbers = list(recovery = rec, phenotype_accuracy = acc_ph, model_sel = model_sel,
                          struct = struct, p4 = p4))
}
