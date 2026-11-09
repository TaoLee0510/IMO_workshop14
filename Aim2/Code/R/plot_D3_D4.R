# =============================================================================
# plot_D3_D4.R
# Purpose : Figures for D3 (identifiability), D3b (clone traits unknown) and D4 (population fit and
#           predictions P1-P5). Read result objects only.
# Inputs  : step result lists, config; D4 also uses the D2 grid if present.
# Outputs : named lists of ggplots.
# =============================================================================

scenario_labels <- c(counts = "Counts + CRP", counts_pre_vaf = "+ pre-treatment VAF",
                     counts_pre_post_vaf = "+ day +30 VAF")

plot_D3 <- function(res, cfg, ...) {
  a <- data.table::copy(res$tables$mechanism_accuracy)[, what := "All mechanisms"]
  b <- res$tables$persistent_h1_vs_h3[, .(scenario, n, mechanism_accuracy = h1_vs_h3_accuracy)][, what := "Persistent: H1 vs H3"]
  x <- rbind(a, b)
  x[, scenario := factor(scenario, names(scenario_labels), scenario_labels)]
  acc <- ggplot2::ggplot(x, ggplot2::aes(scenario, mechanism_accuracy, fill = what)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(0.7), width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mechanism_accuracy)),
                       position = ggplot2::position_dodge(0.7), vjust = -0.3, size = 3) +
    ggplot2::scale_fill_manual(values = c(cfg$plot$col_accent[["steel"]], cfg$plot$col_accent[["crimson"]]), name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
    ggplot2::labs(title = "How well can the mechanism be identified?",
                  subtitle = "Virtual patients refitted under three information scenarios",
                  x = NULL, y = "Accuracy") + theme_aim2()
  p <- res$patients[model == "M2" & ch == TRUE]
  p[, scenario := factor(scenario, names(scenario_labels), scenario_labels)]
  sc <- ggplot2::ggplot(p, ggplot2::aes(f0, f0_est, colour = mechanism_true)) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey60") + ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::facet_wrap(~scenario) + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = mech_colours(cfg), name = "True mechanism") +
    ggplot2::labs(title = "Clone fraction: truth vs estimate (clone carriers)", x = "True f0", y = "Estimated f0") +
    theme_aim2()
  list(accuracy = acc, f0_recovery = list(plot = sc, width = 10, height = 4))
}

plot_D3b <- function(res, cfg, ...) {
  sc_lab <- d3b_labels[res$settings$scenarios]
  sm <- data.table::copy(res$surface_mean)
  sm[, scenario := factor(sc_lab[scenario], sc_lab)]
  sm[, shown := pmax(delta_logLik, -15)]
  surf <- ggplot2::ggplot(sm, ggplot2::aes(factor(delta), factor(r_inf), fill = shown)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(delta_logLik, 1)), size = 2.8,
                       colour = ifelse(sm$shown > -6, "white", "#1A2433")) +
    ggplot2::geom_point(data = sm[is_truth == TRUE], shape = 0, size = 11, stroke = 1.1,
                        colour = cfg$plot$col_accent[["gold"]]) +
    ggplot2::facet_wrap(~scenario, nrow = 1) +
    ggplot2::scale_fill_gradient(low = "#F2F4F7", high = cfg$plot$col_accent[["crimson"]],
                                 limits = c(-15, 0), name = "logLik - max") +
    ggplot2::labs(title = "Clone traits unknown: profile likelihood over clone output and inflammation sensitivity",
                  subtitle = "Mean over virtual cohorts; gold square = true traits (low-output clone)",
                  x = "Clone output delta", y = "Clone sensitivity to inflammation r_inf") +
    theme_aim2() + ggplot2::theme(legend.position = "right")

  pc <- data.table::copy(res$per_cohort)
  pc[, scenario := factor(sc_lab[scenario], sc_lab)]
  ev <- ggplot2::ggplot(pc, ggplot2::aes(scenario, llr_low_vs_functional)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_hline(yintercept = res$settings$region_ll, linetype = "dotted", colour = "grey50") +
    ggplot2::stat_summary(fun = mean, geom = "col", fill = cfg$plot$col_accent[["steel"]], width = 0.6) +
    ggplot2::geom_point(size = 2, colour = cfg$plot$col_accent[["ink"]]) +
    ggplot2::labs(title = "Can the data tell a low-output clone from a functional one?",
                  subtitle = "Log-likelihood ratio, low-output (TP53-like) vs functional (DNMT3A-like) clone; dots = cohorts; dotted = 3",
                  x = NULL, y = "Log-likelihood ratio") + theme_aim2()

  acc <- data.table::melt(pc, id.vars = "scenario", measure.vars = c("h1_vs_h3_estimated", "h1_vs_h3_oracle"),
                          variable.name = "traits", value.name = "accuracy")
  acc[, traits := factor(traits, c("h1_vs_h3_estimated", "h1_vs_h3_oracle"),
                         c("traits estimated from the data", "traits known (oracle)"))]
  ac <- ggplot2::ggplot(acc, ggplot2::aes(scenario, accuracy, fill = traits)) +
    ggplot2::stat_summary(fun = mean, geom = "col", position = ggplot2::position_dodge(0.7), width = 0.65) +
    ggplot2::geom_point(position = ggplot2::position_dodge(0.7), size = 1.6, colour = cfg$plot$col_accent[["ink"]]) +
    ggplot2::scale_fill_manual(values = c(cfg$plot$col_accent[["crimson"]], "grey70"), name = NULL) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Persistent patients: low reserve (H1) vs clone (H3)",
                  subtitle = sprintf("Mechanism accuracy when the clone traits must be estimated (%.0f persistent patients per cohort, %.0f clonal)",
                                     mean(res$per_cohort$n_persistent), mean(res$per_cohort$n_h3)),
                  x = NULL, y = "Accuracy") + theme_aim2()
  list(profile_surface = list(plot = surf, width = 12, height = 3.8),
       clone_type_evidence = list(plot = ev, width = 7.5, height = 4.2),
       mechanism_accuracy = list(plot = ac, width = 7.5, height = 4.2))
}

plot_D3c <- function(res, cfg, ...) {
  lab <- d3c_labels[res$settings$variants]
  pr <- data.table::copy(res$per_rep)
  pr[, variant := factor(lab[variant], lab)]
  ev <- ggplot2::ggplot(pr, ggplot2::aes(variant, llr_low_vs_functional)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_hline(yintercept = 2, linetype = "dotted", colour = "grey50") +
    ggplot2::geom_line(ggplot2::aes(group = rep), colour = "grey80") +
    ggplot2::stat_summary(fun = stats::median, fun.min = stats::median, fun.max = stats::median,
                          geom = "crossbar", width = 0.4, colour = cfg$plot$col_accent[["crimson"]]) +
    ggplot2::geom_point(size = 2, colour = cfg$plot$col_accent[["ink"]]) +
    ggplot2::labs(title = "At the planned data scale: can the fit tell the clone type?",
                  subtitle = sprintf("%d simulated cohorts at the current scale (%d training patients, ~%.0f sequenced), every parameter estimated; dotted = 2",
                                     max(pr$rep), round(mean(pr$n_train)), mean(pr$n_sequenced)),
                  x = NULL, y = "Log-likelihood ratio\nlow-output vs functional clone") + theme_aim2()

  pr[, p_h3_clonal := if (n_h3 > 0) sum_p_h3_h3 / n_h3 else NA_real_, by = .(rep, variant)]
  mv <- c("cor_logit_f0", "h3_rate", "p_h3_clonal", "h1_fp_rate")
  m <- data.table::melt(pr, id.vars = c("rep", "variant"), measure.vars = mv, variable.name = "metric")
  m <- m[!is.na(value)]
  m[, metric := factor(metric, mv,
                       c("Clone fraction recovery\n(correlation, CH carriers)",
                         "Clonal persistent cases\ncalled clonal (MAP)",
                         "Clonal persistent cases\nmean P(H3)",
                         "Low-reserve cases\nwrongly called clonal"))]
  mt <- ggplot2::ggplot(m, ggplot2::aes(variant, value)) +
    ggplot2::stat_summary(fun = mean, geom = "col", fill = cfg$plot$col_accent[["steel"]], width = 0.6) +
    ggplot2::geom_jitter(width = 0.08, height = 0, size = 1.5, colour = cfg$plot$col_accent[["ink"]]) +
    ggplot2::facet_wrap(~metric, nrow = 1) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(title = "Per-patient conclusions at the planned data scale",
                  subtitle = sprintf("Bars = mean over cohorts; ~%.0f persistent patients per cohort, ~%.0f of them clonal",
                                     mean(pr$n_h1 + pr$n_h3), mean(pr$n_h3)),
                  x = NULL, y = NULL) + theme_aim2() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
  list(clone_type_evidence = list(plot = ev, width = 7, height = 4.5),
       patient_level = list(plot = mt, width = 12, height = 4.2))
}

plot_D3d <- function(res, cfg, ...) {
  pr <- data.table::copy(res$per_rep)
  pr[, variant := factor(variant, names(d3d_labels), d3d_labels)]
  dl <- data.table::melt(pr, id.vars = c("rep", "variant"), measure.vars = c("dbic_m1", "dbic_m1_exact"),
                         variable.name = "days", value.name = "dbic")
  dl[, days := factor(days, c("dbic_m1", "dbic_m1_exact"), c("Treatment days approximated (as real data)",
                                                              "True treatment days"))]
  choice <- ggplot2::ggplot(dl, ggplot2::aes(variant, dbic, colour = days)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
    ggplot2::geom_hline(yintercept = 10, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(position = ggplot2::position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
                        size = 1.3, alpha = 0.7) +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1), breaks = c(-100, -10, 0, 10, 100, 1000)) +
    ggplot2::scale_colour_manual(values = unname(cfg$plot$col_accent[c("steel", "gold")]), name = NULL) +
    ggplot2::labs(title = paste0("Can ", res$settings$n, " patients tell M1 (inflammation) from M0?"),
                  subtitle = "Each point = one simulated cohort; dashed line: strong evidence (dBIC = 10)",
                  x = NULL, y = "BIC(M0) - BIC(M1)") + theme_aim2()
  rl <- data.table::melt(pr, id.vars = c("rep", "variant"),
                         measure.vars = c("cor_log_rho", "cor_log_g", "cor_log_rho_exact", "cor_log_g_exact"),
                         variable.name = "what", value.name = "r")
  rl[, days := ifelse(grepl("exact", what), "True treatment days", "Treatment days approximated (as real data)")]
  rl[, param := ifelse(grepl("rho", what), "Reserve rho", "Inflammatory gain g")]
  rec <- ggplot2::ggplot(rl[, .(r = mean(r, na.rm = TRUE)), by = .(variant, days, param)],
                         ggplot2::aes(param, r, fill = days)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7), width = 0.6) +
    ggplot2::facet_wrap(~variant) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::scale_fill_manual(values = unname(cfg$plot$col_accent[c("steel", "gold")]), name = NULL) +
    ggplot2::labs(title = "Recovery of individual parameters (correlation with truth, log scale)",
                  subtitle = "Mean over simulated cohorts", x = NULL, y = "Correlation") + theme_aim2()
  list(model_choice = list(plot = choice, width = 8, height = 4.5), recovery = list(plot = rec, width = 9, height = 4))
}

plot_D4 <- function(res, cfg, res_D2 = NULL, ...) {
  out <- list()
  labs_ph <- pheno_lab(cfg)
  comps <- c(sum_gain_M1_over_M0 = "M1 vs M0 (inflammation)", sum_gain_M2_over_M1 = "M2 vs M1 (clone)")
  comps <- comps[names(comps) %in% names(res$P5$table)]
  t5 <- data.table::melt(res$P5$table[!is.na(phenotype)], id.vars = c("phenotype", "n"),
                         measure.vars = names(comps), variable.name = "comparison", value.name = "gain")
  t5[, comparison := factor(comparison, names(comps), comps)]
  out$P5_gain_by_phenotype <- ggplot2::ggplot(t5, ggplot2::aes(factor(phenotype, phenotype_levels), gain, fill = phenotype)) +
    ggplot2::geom_col(width = 0.6) + ggplot2::facet_wrap(~comparison, scales = "free_y") +
    scale_fill_pheno(cfg, guide = "none") + ggplot2::scale_x_discrete(labels = labs_ph) +
    ggplot2::labs(title = "Where does each model layer earn its keep? (P5)",
                  subtitle = "Gain in Laplace log-likelihood summed over training patients",
                  x = NULL, y = "Log-likelihood gain", caption = synthetic_caption(cfg)) + theme_aim2()

  ex <- res$examples
  cur <- merge(ex$curves, ex$patients, by = "patient_id")
  ob <- merge(ex$obs[var == "anc" & keep == TRUE], ex$patients, by = "patient_id")
  cur[, lab := paste0(patient_id, " (", labs_ph[phenotype], ")")]
  ob[, lab := paste0(patient_id, " (", labs_ph[phenotype], ")")]
  last <- cfg$fit$obs_max_day %||% Inf
  if (is.finite(last)) {
    cur <- cur[day <= max(ob$day)]
    ob[, used := day <= last]
  } else ob[, used := TRUE]
  out$example_fits <- list(plot = ggplot2::ggplot() +
    (if (is.finite(last)) ggplot2::geom_vline(xintercept = last, linetype = "dotted", colour = "grey50")) +
    ggplot2::geom_point(data = ob, ggplot2::aes(day, value, shape = used), colour = "grey30", size = 0.9) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    ggplot2::geom_line(data = cur, ggplot2::aes(day, pmax(anc, 0.01), colour = model), linewidth = 0.7) +
    ggplot2::facet_wrap(~lab, ncol = 3) + ggplot2::scale_y_log10() + scale_colour_model(cfg) +
    ggplot2::labs(title = "Individual fits (training split)",
                  subtitle = if (is.finite(last)) paste0("Fitted to day ", last, "; open circles after that are not used in the fit"),
                  x = "Day after infusion",
                  y = "ANC (x10^9/L, log scale)", caption = synthetic_caption(cfg)) + theme_aim2(),
    width = 10, height = 6)

  if (!is.null(res$P1$table) && nrow(res$P1$table)) {
    p1 <- res$P1$table
    out$P1_late_shape <- ggplot2::ggplot(p1, ggplot2::aes(level_rel, slope_log_per_day, colour = mechanism)) +
      ggplot2::geom_hline(yintercept = c(-1, 1) * cfg$steps$P1$slope_flat, linetype = "dotted") +
      ggplot2::geom_vline(xintercept = cfg$steps$P1$level_low, linetype = "dotted") +
      ggplot2::geom_point(size = 2.2) + ggplot2::scale_x_log10() +
      ggplot2::scale_colour_manual(values = mech_colours(cfg), name = "Model-attributed mechanism") +
      ggplot2::labs(title = "P1: late plateau or slow climb? (persistent patients, day 45-90)",
                    x = "Late ANC level / baseline (log)", y = "Slope of log ANC per day",
                    caption = synthetic_caption(cfg)) + theme_aim2()
  }
  if (nrow(res$P2$table)) {
    out$P2_timing <- ggplot2::ggplot(res$P2$table, ggplot2::aes(peak_day, nadir_day)) +
      ggplot2::geom_jitter(width = 0.15, height = 0.15, colour = cfg$plot$col_phenotype[["second_wave"]]) +
      ggplot2::geom_abline(intercept = res$P2$summary$model_implied_lag, slope = 1, linetype = "dashed") +
      ggplot2::labs(title = paste0("P2: ", tolower(labs_ph[["second_wave"]]), " nadir follows the inflammatory peak"),
                    subtitle = "Dashed: model-implied lag", x = "Ferritin peak day (0-14)",
                    y = paste0(labs_ph[["second_wave"]], " ANC nadir day"), caption = synthetic_caption(cfg)) + theme_aim2()
  }
  pts <- res$inferred[!is.na(phenotype)]
  out$reserve_vs_gain <- ggplot2::ggplot(pts, ggplot2::aes(rho, g, colour = phenotype)) +
    ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + scale_colour_pheno(cfg) +
    ggplot2::labs(title = "Estimated reserve and inflammatory gain by phenotype",
                  subtitle = paste0("MAP estimates, model ", utils::tail(cfg$fit$models %||% "M2", 1)),
                  x = "Marrow reserve rho (log)", y = "Inflammatory gain g (log)",
                  caption = synthetic_caption(cfg)) + theme_aim2()
  if (!is.null(res$late$table)) {
    lt <- res$late$table
    out$late_check <- ggplot2::ggplot(lt, ggplot2::aes(anc_pred, anc_obs, colour = window)) +
      ggplot2::geom_abline(linetype = "dashed", colour = "grey60") +
      ggplot2::geom_point(size = 2) + ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(values = unname(cfg$plot$col_accent[c("steel", "gold", "crimson")]), name = "Day") +
      ggplot2::labs(title = paste0("Counts after day ", cfg$fit$obs_max_day, " vs forward prediction"),
                    subtitle = "Prediction from the fit to the first 30 days (MAP parameters)",
                    x = "Predicted ANC (log)", y = "Observed ANC (log)", caption = synthetic_caption(cfg)) +
      theme_aim2()
  }
  if (!is.null(res_D2) && any(!is.na(res$inferred$f0))) {
    # each patient is shown in the facet of the D2 gain level closest to its estimated gain
    lev <- cfg$steps$D2$g_levels
    lab <- stats::setNames(paste0("inflammatory gain: ", names(lev)), names(lev))
    bg <- data.table::copy(res_D2$grid[delta_level == "low_output"])
    pts <- data.table::copy(res$inferred)[!is.na(phenotype)]
    pts[, g_level := names(lev)[vapply(g, function(x) which.min(abs(log(x) - log(lev))), integer(1))]]
    pts[, f0_plot := pmin(f0, max(bg$f0))]
    bg[, facet := factor(lab[g_level], lab)]
    pts[, facet := factor(lab[g_level], lab)]
    out$patients_on_phasemap <- list(plot = ggplot2::ggplot() +
      ggplot2::geom_raster(data = bg, ggplot2::aes(rho, f0, fill = phenotype), alpha = 0.3) +
      ggplot2::geom_point(data = pts, ggplot2::aes(rho, f0_plot, colour = phenotype), size = 1.5) +
      ggplot2::facet_wrap(~facet, nrow = 1, drop = FALSE) +
      ggplot2::scale_x_log10(limits = range(bg$rho), oob = scales::squish) +
      scale_fill_pheno(cfg, guide = "none") + scale_colour_pheno(cfg) +
      ggplot2::labs(title = "Training patients placed on the phase map",
                    subtitle = "Points: observed phenotype at the patient's estimated parameters. Background: model phenotype for a low-output clone (design values)",
                    x = "Estimated reserve (log)", y = "Estimated clone fraction f0", caption = synthetic_caption(cfg)) +
      theme_aim2(), width = 11, height = 4.5)
  }
  out
}
