# =============================================================================
# plot_D5_D6_H_V.R
# Purpose : Figures for D5 (CH case series), PR / D5b (pre-registered
#           predictions), D6 (virtual cohort), H (hold-out) and V (validation).
# Inputs  : step result lists, config.
# Outputs : named lists of ggplots.
# =============================================================================

plot_D5 <- function(res, cfg, ...) {
  if (is.null(res$table) || !nrow(res$table)) return(list())
  x <- data.table::copy(res$table)
  lod <- cfg$data$lod$vaf
  x[, meas := ifelse(status == "positive", f_measured, lod)]
  p <- ggplot2::ggplot(x, ggplot2::aes(meas, f0_without_vaf, colour = gene_class)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = f0_without_lo, ymax = f0_without_hi), width = 0) +
    ggplot2::geom_point(ggplot2::aes(shape = phenotype), size = 2.5) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = c(DDR = cfg$plot$col_accent[["crimson"]], DTA = cfg$plot$col_accent[["steel"]],
                                            other = cfg$plot$col_accent[["gold"]], none = "grey55"), name = "Gene class") +
    ggplot2::labs(title = "CH case series: model-inferred vs measured clone fraction",
                  subtitle = paste0("Leave-VAF-out estimates with 80% intervals; negatives at the detection limit. Spearman ",
                                    res$concordance$spearman, " (n = ", res$concordance$n, ", descriptive)"),
                  x = "Measured clone fraction (2 x VAF, log)", y = "Inferred f0 without VAF (log)",
                  caption = synthetic_caption(cfg)) + theme_aim2()
  list(case_series = p)
}

plot_PR <- function(res, cfg, ...) {
  if (is.null(res$prediction)) return(list())
  x <- res$prediction
  p <- ggplot2::ggplot(x, ggplot2::aes(y = patient_id)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = vaf_q2.5, xmax = vaf_q97.5), width = 0, orientation = "y", colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = vaf_q10, xmax = vaf_q90), width = 0, orientation = "y", linewidth = 2,
                            colour = cfg$plot$col_accent[["steel"]]) +
    ggplot2::geom_point(ggplot2::aes(x = vaf_q50), size = 2.5, colour = cfg$plot$col_accent[["ink"]]) +
    ggplot2::geom_point(ggplot2::aes(x = pre_vaf), shape = 4, size = 3, colour = cfg$plot$col_accent[["crimson"]]) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = "Pre-registered prediction of post-treatment VAF",
                  subtitle = "Bars: 80% and 95% prediction intervals; x: pre-treatment VAF",
                  x = "VAF (log)", y = NULL, caption = synthetic_caption(cfg)) + theme_aim2()
  list(prediction = list(plot = p, width = 7, height = 3))
}

plot_D5b <- function(res, cfg, ...) {
  x <- res$table
  p <- ggplot2::ggplot(x, ggplot2::aes(y = patient_id)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = pred_lo95, xmax = pred_hi95), width = 0, orientation = "y", colour = "grey60") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = pred_lo80, xmax = pred_hi80), width = 0, orientation = "y", linewidth = 2,
                            colour = cfg$plot$col_accent[["steel"]]) +
    ggplot2::geom_point(ggplot2::aes(x = measured_vaf, colour = verdict), size = 3.5, shape = 18) +
    ggplot2::geom_point(ggplot2::aes(x = pre_vaf), shape = 4, size = 3) + ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = c(supports = cfg$plot$col_accent[["green"]],
                                            inconclusive = cfg$plot$col_accent[["gold"]],
                                            contradicts = cfg$plot$col_accent[["crimson"]]), name = NULL) +
    ggplot2::labs(title = "Pre-registered prediction vs later measurement",
                  subtitle = "Selected patients (marrow work-up for low counts): tests direction and size, not frequency",
                  x = "VAF (log)", y = NULL, caption = synthetic_caption(cfg)) + theme_aim2()
  list(comparison = list(plot = p, width = 7, height = 3))
}

plot_D6 <- function(res, cfg, ...) {
  af <- res$tables$attributable_fraction
  a <- ggplot2::ggplot(af, ggplot2::aes(scenario, attributable_fraction)) +
    ggplot2::geom_col(fill = cfg$plot$col_accent[["crimson"]], width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(100 * attributable_fraction), "%")), vjust = -0.3) +
    ggplot2::labs(title = "Share of persistent cytopenia attributable to clones",
                  subtitle = "Virtual cohort; counterfactual removes the clone (f0 = 0)",
                  x = "CH prevalence scenario", y = "Attributable fraction", caption = synthetic_caption(cfg)) +
    theme_aim2()
  pw <- ggplot2::ggplot(res$tables$power, ggplot2::aes(n_enrolled, power, colour = scenario)) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dotted") + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(title = "Paired design: power to detect clonal expansion by day +30",
                  subtitle = "Stratified sampling: day +30 VAF if ANC < 1.0, plus a random 25% of others",
                  x = "Patients enrolled", y = "Power", colour = NULL, caption = synthetic_caption(cfg)) + theme_aim2()
  list(attributable_fraction = a, power = pw)
}

plot_H <- function(res, cfg, ...) {
  m <- data.table::copy(res$metrics)
  m[, method := factor(method, c("twin", "baseline_score", "ml_elastic_net"),
                       c("Digital twin", "Baseline score", "Elastic net"))]
  lab <- pheno_lab(cfg)
  m[, metric := factor(metric, c("auc_persistent", "auc_second_wave", "brier"),
                       c(paste0("AUC: ", lab[["persistent"]]), paste0("AUC: ", lab[["second_wave"]]),
                         "Brier score (multiclass)"))]
  f <- ggplot2::ggplot(m, ggplot2::aes(estimate, method)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y") +
    ggplot2::geom_point(size = 2.5, colour = cfg$plot$col_accent[["crimson"]]) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::labs(title = if (identical(res$mode, "cv"))
                    paste0(cfg$steps$H$cv_folds, "-fold cross-validation, prediction at day +", cfg$steps$H$landmark_day,
                           " (n = ", res$n_test, ")")
                  else paste0("Hold-out at day +", cfg$steps$H$landmark_day, " (n = ", res$n_test, ")"),
                  subtitle = "Brier: lower is better; AUC: higher is better; 95% bootstrap CI",
                  x = NULL, y = NULL, caption = synthetic_caption(cfg)) + theme_aim2()
  cal <- ggplot2::ggplot(res$calibration, ggplot2::aes(mean_pred, obs_rate, colour = method)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~outcome) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "Calibration (risk quintiles)", x = "Mean predicted probability",
                  y = "Observed rate", colour = NULL) + theme_aim2()
  list(metrics = list(plot = f, width = 10, height = 3.5), calibration = list(plot = cal, width = 8, height = 4.5))
}

plot_V <- function(res, cfg, ...) {
  x <- res$scatter
  long <- rbind(x[, .(parameter = "log rho", truth = log(rho), estimate = eta_rho, ch)],
                x[, .(parameter = "log g", truth = log(g), estimate = eta_g, ch)],
                x[ch == 1, .(parameter = "logit f0 (CH)", truth = stats::qlogis(f0), estimate = eta_f0, ch)])
  r <- ggplot2::ggplot(long, ggplot2::aes(truth, estimate)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(alpha = 0.5, colour = cfg$plot$col_accent[["steel"]]) +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::labs(title = "Validation: individual parameter recovery (training split, M2)",
                  x = "True value", y = "Estimate (MAP)", caption = synthetic_caption(cfg)) + theme_aim2()
  cm <- ggplot2::ggplot(res$phenotype_confusion[!is.na(phenotype) & !is.na(phenotype_true)],
                        ggplot2::aes(phenotype_true, phenotype, fill = N)) +
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = N), colour = "white") +
    ggplot2::scale_fill_gradient(low = cfg$plot$col_accent[["slate"]], high = cfg$plot$col_accent[["ink"]]) +
    ggplot2::labs(title = paste0("Observed vs true phenotype (accuracy ", res$phenotype_accuracy, ")"),
                  x = "True (noise-free)", y = "Observed rule") + theme_aim2()
  list(parameter_recovery = list(plot = r, width = 10, height = 3.8),
       phenotype_confusion = list(plot = cm, width = 5, height = 4.5))
}
