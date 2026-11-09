# =============================================================================
# plot_steps.R
# Purpose : Figures for W0-W3. Every function only reads a saved step result
#           (res$tables) and returns a named list of ggplots.
# Inputs  : step result, config.
# Outputs : list(name = ggplot | list(plot, width, height)).
# =============================================================================

log_breaks_anc <- c(0.01, 0.05, 0.1, 0.5, 1, 1.5, 5)

plot_W0 <- function(res, cfg, ...) {
  d <- res$derived
  ph <- d$phenotype[, .(patient_id, primary)]
  o <- merge(d$obs[keep == TRUE & day >= -5 & day <= 35 & var %in% c("anc", "plt")], ph, by = "patient_id")
  o[, primary := factor(primary, levels = cfg$phenotype$classes)]
  o[, var := factor(var, levels = c("anc", "plt"), labels = c("ANC (x10^9/L)", "Platelets (x10^9/L)"))]
  traj <- ggplot2::ggplot(o, ggplot2::aes(day, value, group = patient_id, colour = primary)) +
    ggplot2::geom_line(alpha = 0.55, linewidth = 0.35) +
    ggplot2::geom_point(data = o[ub == TRUE], shape = 6, size = 0.9, alpha = 0.7) +
    ggplot2::geom_hline(data = data.table::data.table(var = factor("ANC (x10^9/L)", levels = levels(o$var)),
                                                      y = c(1.0, 1.5)),
                        ggplot2::aes(yintercept = y), linetype = 3, colour = "#6B7A8F") +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(var ~ primary, scales = "free_y", labeller = ggplot2::labeller(primary = pheno_labels)) +
    scale_colour_pheno(cfg) + theme_aim3() + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = "Recovery phenotypes in the 35 patients (locked rule)",
                  subtitle = "Triangles = upper bounds (no differential, zero on the differential, ANC under G-CSF)",
                  x = "Day after infusion", y = NULL)
  cnt <- data.table::melt(res$tables$phenotype_counts, id.vars = "rule", variable.name = "class", value.name = "n")
  cnt[, rule := factor(rule, levels = rev(res$tables$phenotype_counts$rule))]
  bars <- ggplot2::ggplot(cnt, ggplot2::aes(n, rule, fill = class)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::geom_text(ggplot2::aes(label = n), position = ggplot2::position_stack(vjust = 0.5),
                                                        colour = "white", size = 3) +
    scale_fill_pheno(cfg) + theme_aim3() +
    ggplot2::labs(title = "Phenotype counts under the primary rule and its variants", x = "Patients", y = NULL)
  list(trajectories = list(plot = traj, width = 9, height = 5.5), rule_variants = bars)
}

plot_W1 <- function(res, cfg, ...) {
  tr <- data.table::melt(res$tables$typical[, .(day, anc, plt = plt / 100, amc, case, pheno)],
                         id.vars = c("day", "case", "pheno"), variable.name = "lineage")
  tr[, case := factor(case, levels = c("none", "a", "b", "ab"),
                      labels = c("Progenitors only", "HSC hit in a", "HSC hit in b", "HSC hit split a + b"))]
  tr[, lineage := factor(lineage, levels = c("anc", "plt", "amc"),
                         labels = c("ANC", "Platelets (x100)", "Monocytes"))]
  p <- ggplot2::ggplot(tr, ggplot2::aes(day, value, colour = case)) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~lineage, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = unname(cfg$plot$col_accent[c("slate", "steel", "gold", "crimson")]), name = NULL) +
    theme_aim3() + ggplot2::labs(title = "Typical design patient: same total HSC loss, different timing",
                                 x = "Day after infusion", y = "Count (x10^9/L)")
  tests <- res$tables$unit_tests
  tt <- ggplot2::ggplot(tests, ggplot2::aes(1, factor(test, levels = rev(test)), fill = passed)) +
    ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = detail), size = 2.6, hjust = 0.5) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#CFE3D6", `FALSE` = "#E9C3C3")) +
    theme_aim3() + ggplot2::theme(axis.text.x = ggplot2::element_blank(), legend.position = "none") +
    ggplot2::labs(title = "Unit tests", x = NULL, y = NULL)
  list(typical = list(plot = p, width = 9, height = 3.8), unit_tests = list(plot = tt, width = 9, height = 3.5))
}

plot_W2 <- function(res, cfg, ...) {
  tb <- res$tables
  tr <- data.table::melt(tb$trajectories[day >= -5 & day <= 50, .(day, anc, plt = plt / 100, amc, case)],
                         id.vars = c("day", "case"), variable.name = "lineage")
  tr[, case := factor(case, levels = c("none", "a", "b", "ab"),
                      labels = c("Progenitors only", "90% HSC loss in a", "90% HSC loss in b", "90% split a + b"))]
  tr[, lineage := factor(lineage, levels = c("anc", "plt", "amc"), labels = c("ANC", "Platelets (x100)", "Monocytes"))]
  fp <- ggplot2::ggplot(tr, ggplot2::aes(day, value, colour = case)) + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_hline(data = data.table::data.table(lineage = factor("ANC", levels = levels(tr$lineage)), y = c(1, 1.5)),
                        ggplot2::aes(yintercept = y), linetype = 3, colour = "#6B7A8F") +
    ggplot2::facet_wrap(~lineage, scales = "free_y", nrow = 1) +
    ggplot2::scale_colour_manual(values = unname(cfg$plot$col_accent[c("slate", "steel", "gold", "crimson")]), name = NULL) +
    theme_aim3() + ggplot2::labs(title = "Fingerprints: the same stem cell loss, delivered before or after infusion",
                                 x = "Day after infusion", y = "Count (x10^9/L)")
  pm <- data.table::copy(tb$phase_map)
  pm[, `:=`(fa = factor(f_a), fb = factor(f_b), eta = paste0("synergy eta = ", eta),
            set = paste0("ANC set point: ", set))]
  phase <- ggplot2::ggplot(pm, ggplot2::aes(fa, fb, fill = pheno)) + ggplot2::geom_tile(colour = "white", linewidth = 0.2) +
    ggplot2::facet_grid(eta ~ set) + scale_fill_pheno(cfg) + theme_aim3(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
    ggplot2::labs(title = "Phenotype map over stem cell loss in each phase",
                  x = "HSC fraction killed by chemotherapy (f_a)", y = "HSC fraction killed in the CAR T phase (f_b)")
  pw <- tb$pseudo_wave
  pseudo <- ggplot2::ggplot(pw, ggplot2::aes(f_a, min_after, colour = pheno)) +
    ggplot2::geom_line(ggplot2::aes(group = regime), colour = "#B8C0CC") + ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_hline(yintercept = 1, linetype = 3) + ggplot2::facet_wrap(~regime, nrow = 1) +
    scale_colour_pheno(cfg) + theme_aim3() +
    ggplot2::labs(title = "Can chemotherapy alone produce a second wave?",
                  subtitle = "Lowest ANC after the first recovery to 1.5, by multipotent-progenitor regime",
                  x = "HSC fraction killed by chemotherapy (f_a)", y = "Lowest ANC after recovery")
  tm <- data.table::copy(tb$timing)
  tm[, lineage := factor(lineage, levels = c("amc", "anc", "plt"), labels = c("Monocytes", "ANC", "Platelets"))]
  tm[, phase := factor(phase, levels = c("a", "b"), labels = c("HSC hit in a", "HSC hit in b"))]
  timing <- ggplot2::ggplot(tm, ggplot2::aes(shift, footprint_day, colour = lineage)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_wrap(~phase) +
    ggplot2::scale_colour_manual(values = unname(cfg$plot$col_lineage[c("amc", "anc", "plt")]), name = NULL) +
    theme_aim3() + ggplot2::labs(title = "Timing fingerprint: the footprint follows the inflammatory peak only in b",
                                 x = "Shift of the patient's ferritin peak (days)", y = "Day of largest HSC footprint")
  list(fingerprint = list(plot = fp, width = 9, height = 3.8), phase_map = list(plot = phase, width = 9, height = 7),
       pseudo_wave = list(plot = pseudo, width = 9, height = 3.6), timing = list(plot = timing, width = 7, height = 3.6))
}

plot_W3 <- function(res, cfg, ...) {
  tb <- res$tables
  best <- data.table::copy(tb$best)
  best[, model := factor(model, levels = cfg$fit$models)]
  sel <- ggplot2::ggplot(best, ggplot2::aes(scenario, fill = model)) + ggplot2::geom_bar(position = "fill") +
    ggplot2::facet_wrap(~paste0("n = ", n)) +
    ggplot2::scale_fill_manual(values = c(M0 = "#B8C0CC", Ma = "#2E6B8A", Mb = "#C9A227", Mab = "#9E2A2B", Mabx = "#6E3B6E"),
                               labels = model_labels, name = NULL) +
    theme_aim3() + ggplot2::labs(title = "Which model wins (lowest BIC) in virtual cohorts at the real design",
                                 x = "True scenario", y = "Share of cohorts")
  cf <- tb$confusion[n == min(n)]
  cf[, share := N / sum(N), by = truth]
  conf <- ggplot2::ggplot(cf, ggplot2::aes(map, truth, fill = share)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = N), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F4F6F8", high = "#2E6B8A", limits = c(0, 1)) + theme_aim3() +
    ggplot2::labs(title = "Per-patient attribution (Mab fit), n = 35 cohorts", x = "Most probable class", y = "True class")
  cal <- tb$calibration[N > 0]
  calp <- ggplot2::ggplot(cal, ggplot2::aes(mean_p, observed, size = N, colour = class)) +
    ggplot2::geom_abline(linetype = 3) + ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_colour_manual(values = cfg$plot$col_mech, labels = mech_labels, name = NULL) +
    theme_aim3() + ggplot2::labs(title = "Calibration of P(class)", x = "Predicted probability", y = "Observed frequency")
  list(model_selection = list(plot = sel, width = 8, height = 4), confusion = list(plot = conf, width = 6, height = 4.5),
       calibration = list(plot = calp, width = 6, height = 4.5))
}

plot_PR <- function(res, cfg, ...) list()
