# =============================================================================
# plot_S0_D1_D2.R
# Purpose : Figures for S0 (data overview), D1 (model archetypes) and
#           D2 (phenotype phase map). Read result objects only.
# Inputs  : step result lists, config.
# Outputs : named lists of ggplots (saved by io_write.R).
# =============================================================================

plot_S0 <- function(res, cfg, ...) {
  smp <- ggplot2::ggplot(res$sampling, ggplot2::aes(day, N)) +
    ggplot2::geom_col(fill = cfg$plot$col_accent[["steel"]], width = 0.8) +
    ggplot2::labs(title = "ANC observations per day (all patients)", x = "Day after infusion",
                  y = "Observations", caption = synthetic_caption(cfg)) + theme_aim2()

  tr <- res$traj_train[var == "anc" & !is.na(phenotype)]
  band <- tr[, .(med = stats::median(value), lo = stats::quantile(value, 0.25),
                 hi = stats::quantile(value, 0.75), n = .N), by = .(phenotype, day)][n >= 3]
  traj <- ggplot2::ggplot(band, ggplot2::aes(day, med, colour = phenotype, fill = phenotype)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = c(0.5, 1), linetype = "dotted", colour = "grey50") +
    ggplot2::scale_y_log10() + scale_colour_pheno(cfg) + scale_fill_pheno(cfg) +
    ggplot2::labs(title = "ANC by recovery phenotype (training split)", subtitle = "Median and IQR",
                  x = "Day after infusion", y = "ANC (x10^9/L, log scale)", caption = synthetic_caption(cfg)) +
    theme_aim2()

  ph <- res$tables$phenotype_train[!is.na(phenotype)]
  bar <- ggplot2::ggplot(ph, ggplot2::aes(factor(phenotype, phenotype_levels), N, fill = phenotype)) +
    ggplot2::geom_col(width = 0.6) + ggplot2::geom_text(ggplot2::aes(label = paste0(pct, "%")), vjust = -0.3) +
    scale_fill_pheno(cfg, guide = "none") + ggplot2::scale_x_discrete(labels = pheno_lab(cfg)) +
    ggplot2::labs(title = "Phenotypes (training split)", x = NULL, y = "Patients") + theme_aim2()
  list(sampling = smp, anc_by_phenotype = traj, phenotype_bar = list(plot = bar, width = 4.5, height = 4))
}

plot_D1 <- function(res, cfg, ...) {
  tr <- data.table::copy(res$traj)[day >= -5]
  tr[, archetype := factor(archetype, c("quick", "second_wave", "persistent_H1", "persistent_H3"),
                           c("Quick", "Second wave (H2)", "Persistent: low reserve (H1)", "Persistent: clone (H3)"))]
  p <- ggplot2::ggplot(tr, ggplot2::aes(day, pmax(anc, 0.01), colour = model)) +
    ggplot2::geom_hline(yintercept = c(0.5, 1), linetype = "dotted", colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::facet_wrap(~archetype, nrow = 1) +
    ggplot2::scale_y_log10() + scale_colour_model(cfg) +
    ggplot2::labs(title = "One model, three switches: archetypal ANC trajectories",
                  subtitle = "M0 reserve only; M1 + marrow inflammation; M2 + clonal layer (design values). M1 and M2 coincide when f0 = 0",
                  x = "Day after infusion", y = "ANC (x10^9/L, log scale)") + theme_aim2()
  fc <- ggplot2::ggplot(tr[model == "M2"], ggplot2::aes(day, f_clone, colour = archetype)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = unname(c(cfg$plot$col_phenotype, cfg$plot$col_accent[["ink"]])), name = NULL) +
    ggplot2::labs(title = "Clone fraction of the stem pool (M2)", x = "Day after infusion", y = "Clone fraction") +
    theme_aim2()
  list(archetypes = list(plot = p, width = 11, height = 3.8), clone_fraction = fc)
}

plot_D2 <- function(res, cfg, ...) {
  g <- data.table::copy(res$grid)
  g[, g_level := factor(g_level, names(cfg$steps$D2$g_levels),
                        paste0("inflammatory gain: ", names(cfg$steps$D2$g_levels)))]
  g[, delta_level := factor(delta_level, names(cfg$steps$D2$delta_levels),
                            c("low-output clone (TP53-like)", "functional clone (DNMT3A/TET2-like)"))]
  p <- ggplot2::ggplot(g, ggplot2::aes(rho, f0, fill = phenotype)) +
    ggplot2::geom_raster() + ggplot2::facet_grid(delta_level ~ g_level) +
    ggplot2::scale_x_log10() + scale_fill_pheno(cfg) +
    ggplot2::labs(title = "Phenotype phase map: reserve x clone fraction",
                  subtitle = "Natural history after lymphodepletion; design parameter values",
                  x = "Marrow reserve (regenerative capacity, log scale)", y = "Baseline clone fraction f0") +
    theme_aim2()
  e <- ggplot2::ggplot(g[f0 > min(f0)], ggplot2::aes(rho, f0, fill = effect)) +
    ggplot2::geom_raster() + ggplot2::facet_grid(delta_level ~ g_level) + ggplot2::scale_x_log10() +
    ggplot2::scale_fill_manual(values = c(buffer = cfg$plot$col_accent[["green"]], `no change` = "grey88",
                                          aggravate = cfg$plot$col_accent[["crimson"]]), name = NULL) +
    ggplot2::labs(title = "Does the clone buffer or aggravate?",
                  subtitle = "Phenotype severity versus the same patient without a clone",
                  x = "Marrow reserve (log scale)", y = "Baseline clone fraction f0") + theme_aim2()
  list(phasemap = list(plot = p, width = 10, height = 6.5), buffer_aggravate = list(plot = e, width = 10, height = 6.5))
}
