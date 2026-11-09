# =============================================================================
# plot_steps.R
# Purpose : Figures of every step. Plotting functions only read step results.
# Inputs  : step result lists, config.
# Outputs : named lists of ggplots (or list(plot, width, height)).
# =============================================================================

plot_S0 <- function(res, cfg, ...) {
  a <- res$anc[day >= -6 & day <= 35]
  a[res$phenotype, phenotype := i.primary, on = "patient_id"]
  a[, value_plot := pmax(value, 0.02)]
  p1 <- ggplot2::ggplot(a, ggplot2::aes(day, value_plot, group = patient_id, colour = phenotype)) +
    ggplot2::geom_hline(yintercept = c(0.5, 1, 1.5), linetype = 3, colour = "grey60") +
    ggplot2::geom_line(alpha = 0.6) + ggplot2::geom_point(ggplot2::aes(shape = cens), size = 0.9) +
    ggplot2::scale_y_log10() + ggplot2::facet_wrap(~phenotype, labeller = ggplot2::as_labeller(phenotype_display())) +
    scale_colour_pheno(cfg) + ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
    ggplot2::labs(title = "Daily ANC by recovery phenotype (primary rule)", x = "Day from infusion",
                  y = "ANC (x10^9/L, log)", caption = data_caption(cfg)) +
    theme_aim4() + ggplot2::theme(legend.position = "none")
  list(anc_by_phenotype = list(plot = p1, width = 9, height = 3.8))
}

plot_K1 <- function(res, cfg, ...) {
  tr <- res$traj[day >= -6 & day <= 60]
  long <- data.table::melt(tr, id.vars = c("day", "model"), measure.vars = intersect(c("anc", "plt", "amc", "alc"), names(tr)))
  long <- long[!is.na(value)]   # MF has no lymphoid branch
  p <- ggplot2::ggplot(long, ggplot2::aes(day, value, colour = model)) + ggplot2::geom_line() +
    ggplot2::facet_wrap(~variable, scales = "free_y") +
    ggplot2::labs(title = "Typical trajectories of the nested models and the Friberg reference (MF)",
                  x = "Day from infusion", y = "Count", colour = NULL) + theme_aim4()
  out <- list(typical = list(plot = p, width = 10, height = 4))
  # mean time since MPP entry of the circulating cells (parent model)
  ag <- tr[!is.na(age_anc) & model %in% c("M0", "Mab")]
  if (nrow(ag)) {
    la <- data.table::melt(ag, id.vars = c("day", "model"), measure.vars = c("age_anc", "age_plt", "age_amc", "age_alc"))
    la[, variable := factor(variable, labels = c("ANC", "PLT", "AMC", "ALC"))]
    p2 <- ggplot2::ggplot(la, ggplot2::aes(day, value, colour = model)) + ggplot2::geom_line() +
      ggplot2::facet_wrap(~variable, scales = "free_y", nrow = 1) +
      ggplot2::labs(title = "Mean time since MPP entry of the circulating cells", x = "Day from infusion",
                    y = "Days", colour = NULL) + theme_aim4()
    out$age_since_mpp <- list(plot = p2, width = 10, height = 3.5)
  }
  out
}

plot_K2 <- function(res, cfg, ...) {
  m <- res$map_ab
  p1 <- ggplot2::ggplot(m, ggplot2::aes(f_a, f_b, colour = phenotype)) + ggplot2::geom_point(shape = 15, size = 2.2) +
    ggplot2::facet_wrap(~eta, labeller = ggplot2::labeller(eta = function(x) paste0("a x b synergy eta = ", x))) +
    scale_colour_pheno(cfg) + ggplot2::coord_equal() +
    ggplot2::labs(title = "Recovery phenotype by HSC loss from chemotherapy (f_a) and the CAR T phase (f_b)",
                  x = "HSC loss attributable to chemotherapy, f_a", y = "HSC loss attributable to CAR T phase, f_b") + theme_aim4()
  p2 <- ggplot2::ggplot(res$map_rf, ggplot2::aes(r_s, f0, fill = phenotype)) + ggplot2::geom_tile() +
    ggplot2::scale_x_log10() + ggplot2::facet_wrap(~clone, labeller = ggplot2::labeller(clone = c(ddr = "DDR-type clone (PPM1D/TP53-like)", dta = "DTA-type clone (DNMT3A/TET2-like)"))) +
    scale_fill_pheno(cfg) + ggplot2::labs(title = "Who: marrow reserve x clone fraction", x = "HSC regrowth rate r_S (reserve)",
                                          y = "Clone fraction f0") + theme_aim4()
  p3 <- ggplot2::ggplot(res$map_link, ggplot2::aes(th_as, th_bs, fill = phenotype)) + ggplot2::geom_tile() +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~clone, labeller = ggplot2::labeller(clone = c(none = "No clone", ddr = "30% DDR-type clone", dta = "30% DTA-type clone"))) +
    scale_fill_pheno(cfg) + ggplot2::labs(title = "The link: the same a / b injuries with and without a clone",
                                          x = "Chemotherapy injury to HSC, theta_aS", y = "CAR T-phase injury to HSC, theta_bS") + theme_aim4()
  p4 <- ggplot2::ggplot(res$pseudo, ggplot2::aes(hsc_loss, factor(gamma_m), colour = phenotype)) +
    ggplot2::geom_point(size = 2.5) + scale_colour_pheno(cfg) +
    ggplot2::labs(title = "A single chemotherapy hit: phenotype vs HSC loss", x = "HSC loss from chemotherapy",
                  y = "Self-renewal feedback gain gamma_M") + theme_aim4()
  fp <- data.table::melt(res$fingerprint, id.vars = c("day", "route"), measure.vars = c("anc_rel", "plt_rel", "amc_rel", "S"))
  p5 <- ggplot2::ggplot(fp, ggplot2::aes(day, value, colour = route)) + ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~variable, nrow = 1, labeller = ggplot2::labeller(variable = c(anc_rel = "ANC", plt_rel = "PLT", amc_rel = "AMC", S = "HSC pool"))) +
    ggplot2::scale_colour_manual(values = c(cfg$plot$col_accent[["steel"]], cfg$plot$col_accent[["gold"]])) +
    ggplot2::labs(title = "Fingerprints: the same 70% HSC loss placed in a vs b", x = "Day from infusion",
                  y = "Relative to pre-treatment", colour = NULL) + theme_aim4()
  list(map_ab = list(plot = p1, width = 9, height = 4.8), map_reserve_clone = list(plot = p2, width = 9, height = 4.5),
       map_link = list(plot = p3, width = 11, height = 4.5), pseudo_second_wave = p4,
       fingerprint = list(plot = p5, width = 11, height = 3.5))
}

plot_K3 <- function(res, cfg, ...) {
  s <- res$summary
  p1 <- ggplot2::ggplot(s$q2, ggplot2::aes(mech, q2_correct, fill = clone)) +
    ggplot2::geom_col(position = "dodge") + ggplot2::geom_hline(yintercept = cfg$steps$K3$route$cohort_correct, linetype = 2) +
    ggplot2::scale_fill_manual(values = c(none = cfg$plot$col_accent[["steel"]], clone = cfg$plot$col_accent[["crimson"]])) +
    ggplot2::labs(title = "Q2 (when): share of synthetic cohorts where BIC picks the true model", x = "True mechanism",
                  y = "Correct model selected", fill = NULL, caption = "Real sampling design (35 patients); counts simulated.") + theme_aim4()
  att <- s$attribution[design == "real35"]
  p2 <- ggplot2::ggplot(att[true_pheno != "quick"], ggplot2::aes(p_true_when, fill = mech)) +
    ggplot2::geom_histogram(bins = 20, position = "stack") +
    ggplot2::labs(title = "Individual: posterior probability of the true timing class (non-quick patients)",
                  x = "P(true class)", y = "Patients", fill = NULL) + theme_aim4()
  list(q2_selection = p1, individual_when = p2)
}

plot_F1 <- function(res, cfg, ...) {
  cmp <- res$comparison
  p1 <- ggplot2::ggplot(cmp, ggplot2::aes(stats::reorder(model, -delta_BIC), delta_BIC, fill = family)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(a4 = cfg$plot$col_accent[["steel"]], mf = cfg$plot$col_accent[["slate"]])) +
    ggplot2::labs(title = "Model comparison (BIC; lower is better)", x = NULL, y = "Delta BIC vs best",
                  caption = data_caption(cfg)) + theme_aim4() + ggplot2::theme(legend.position = "none")
  fv <- res$fitted
  ph <- data.table::as.data.table(res$tests$T1$data)
  ex <- unique(fv$patient_id)[seq_len(min(6, data.table::uniqueN(fv$patient_id)))]
  fv <- fv[patient_id %in% ex & day >= -10 & day <= 40]
  p2 <- ggplot2::ggplot(fv, ggplot2::aes(day, anc, colour = model)) + ggplot2::geom_line() +
    ggplot2::facet_wrap(~patient_id, scales = "free_y") + ggplot2::scale_y_log10() +
    ggplot2::labs(title = "Fitted ANC (MAP): final model vs Friberg reference", x = "Day", y = "ANC", colour = NULL) + theme_aim4()
  out <- list(model_comparison = p1, fitted_examples = list(plot = p2, width = 9, height = 5.5))
  # ALC: model (lymphoid branch, not fitted) next to the observations (include CAR T cells)
  fa <- fv[model == res$selected$final & !is.na(alc)]
  if (nrow(fa) && !is.null(res$alc_obs)) {
    oa <- res$alc_obs[patient_id %in% ex & day >= -10 & day <= 40]
    p3 <- ggplot2::ggplot(fa, ggplot2::aes(day, alc)) + ggplot2::geom_line(colour = cfg$plot$col_accent[["steel"]]) +
      ggplot2::geom_point(data = oa, ggplot2::aes(day, pmax(value, 0.005), shape = cens), size = 1, colour = "grey30") +
      ggplot2::facet_wrap(~patient_id, scales = "free_y") + ggplot2::scale_y_log10() +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), guide = "none") +
      ggplot2::labs(title = if (identical(cfg$version, "v1.3")) "ALC: fitted lymphocytes + CAR T cells (line) vs observed (points)" else
                      "ALC: model lymphoid branch (line, not fitted) vs observed (points, include CAR T cells)",
                    x = "Day", y = "ALC") + theme_aim4()
    out$alc_examples <- list(plot = p3, width = 9, height = 5.5)
  }
  out
}

plot_F2 <- function(res, cfg, results = list(), ...) {
  s <- data.table::copy(res$summary)
  s[, clone_likely := ifelse(!is.na(p_who_clone) & p_who_clone >= 0.5, "P(clone) >= 0.5", "P(clone) < 0.5")]
  k2 <- results$K2
  p1 <- ggplot2::ggplot() +
    { if (!is.null(k2)) ggplot2::geom_point(data = k2$map_ab[eta == 0], ggplot2::aes(f_a, f_b, colour = phenotype),
                                            shape = 15, size = 2.2, alpha = 0.25) } +
    ggplot2::geom_errorbar(data = s, ggplot2::aes(x = f_a, ymin = f_b_lo, ymax = f_b_hi), colour = "grey40", width = 0, linewidth = 0.3) +
    ggplot2::geom_errorbarh(data = s, ggplot2::aes(y = f_b, xmin = f_a_lo, xmax = f_a_hi), colour = "grey40", height = 0, linewidth = 0.3) +
    ggplot2::geom_point(data = s, ggplot2::aes(f_a, f_b, colour = observed, shape = clone_likely), size = 2.6) +
    scale_colour_pheno(cfg) + ggplot2::scale_shape_manual(values = c(16, 17), name = NULL) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(title = "When was the stem cell pool hit? Each patient on the phenotype map",
                  subtitle = "Points: posterior median HSC loss from chemotherapy (f_a) and the CAR T phase (f_b), 90% intervals; colour = observed phenotype",
                  x = "HSC loss attributable to chemotherapy, f_a", y = "HSC loss attributable to the CAR T phase, f_b",
                  caption = data_caption(cfg)) + theme_aim4()
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  lw <- data.table::melt(s, id.vars = c("patient_id", "observed"), measure.vars = paste0("p_when_", wc))
  lw[, class := factor(sub("p_when_", "", variable), levels = wc)]
  p2 <- ggplot2::ggplot(lw, ggplot2::aes(patient_id, value, fill = class)) + ggplot2::geom_col(width = 0.9) +
    ggplot2::facet_grid(~observed, scales = "free_x", space = "free_x", labeller = ggplot2::as_labeller(phenotype_display())) +
    ggplot2::scale_fill_manual(values = cfg$plot$col_mech, labels = when_labels, name = NULL) +
    ggplot2::labs(title = "When (HSC level): posterior probability of each timing class", x = NULL, y = "Probability") +
    theme_aim4() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6))
  hc <- c("clone", "reserve", "both", "other")
  lh <- data.table::melt(s[observed != "quick" | p_quick < 0.5], id.vars = c("patient_id", "observed"), measure.vars = paste0("p_who_", hc))
  lh[, class := factor(sub("p_who_", "", variable), levels = hc)]
  p3 <- ggplot2::ggplot(lh, ggplot2::aes(patient_id, value, fill = class)) + ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_fill_manual(values = cfg$plot$col_who, labels = who_labels, name = NULL) +
    ggplot2::labs(title = "Who (non-quick patients): reserve or clone", x = NULL, y = "Probability",
                  caption = "Clone fraction is latent on these data (no sequencing): probabilities are prior-dominated.") +
    theme_aim4() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 6))
  list(main_fa_fb = list(plot = p1, width = 7.5, height = 7.5), when_by_patient = list(plot = p2, width = 11, height = 4),
       who_by_patient = list(plot = p3, width = 9, height = 4))
}

plot_F3 <- function(res, cfg, ...) {
  m <- res$metrics
  lng <- data.table::melt(m, id.vars = c("L", "group"), measure.vars = c("auroc_model", "auroc_hematotox", "auroc_logistic"))
  lng[, variable := factor(variable, labels = c("Aim 4 model", "CAR-HEMATOTOX (at LD)", "Counts at L (logistic)"))]
  p <- ggplot2::ggplot(lng, ggplot2::aes(L, value, colour = variable)) + ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_ribbon(data = m, ggplot2::aes(L, ymin = auroc_lo, ymax = auroc_hi), inherit.aes = FALSE, alpha = 0.12) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 3) + ggplot2::facet_wrap(~group) +
    ggplot2::scale_colour_manual(values = unname(cfg$plot$col_accent[c("crimson", "slate", "gold")])) +
    ggplot2::labs(title = if (identical(cfg$version, "v1.3")) "Landmark prediction of non-quick recovery (approximate: population of the full-data fit)" else
                    "Leave-one-out prediction of non-quick recovery by landmark day", x = "Landmark day L",
                  y = "AUROC", colour = NULL, caption = data_caption(cfg)) + theme_aim4()
  list(auroc_by_landmark = list(plot = p, width = 9, height = 4.2))
}

plot_F4 <- function(res, cfg, ...) {
  b <- res$by_class
  p <- ggplot2::ggplot(b, ggplot2::aes(intervention, when_map, fill = mean_delta_nonquick)) + ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%+.2f\nn=%d", mean_delta_nonquick, n)), size = 2.6) +
    ggplot2::scale_fill_gradient2(low = cfg$plot$col_accent[["green"]], mid = "white", high = cfg$plot$col_accent[["crimson"]]) +
    ggplot2::scale_y_discrete(labels = when_labels) +
    ggplot2::labs(title = "What-if: change in P(non-quick recovery) by the patient's most likely timing class",
                  x = NULL, y = NULL, fill = "Delta P") + theme_aim4() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
  list(by_class = list(plot = p, width = 9, height = 4.2))
}

plot_V <- function(res, cfg, ...) {
  if (is.null(res$patients)) return(list())
  p <- ggplot2::ggplot(res$patients, ggplot2::aes(true_when, p_true_when)) + ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 1) +
    ggplot2::labs(title = "TestData: posterior probability of the true timing class", x = "True class", y = "P(true class)") + theme_aim4()
  list(when_vs_truth = p)
}
