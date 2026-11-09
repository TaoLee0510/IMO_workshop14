# =============================================================================
# plot_theme.R
# Purpose : Shared look for every Aim 3 figure (palette matches the team
#           slides: crimson 9E2A2B, steel blue 2E6B8A, slate, gold, green).
# Inputs  : config$plot.
# Outputs : ggplot theme, labels and scale helpers.
# =============================================================================

theme_aim3 <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#1A2433", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(colour = "#4A5568"),
      plot.caption = ggplot2::element_text(colour = "#6B7A8F", hjust = 0, size = base_size - 2),
      strip.text = ggplot2::element_text(face = "bold", colour = "#1A2433"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")
}

pheno_labels <- c(quick = "Quick", second_wave = "Second wave", delayed = "Delayed")
mech_labels <- c(none = "Quick (no HSC role)", a = "Chemotherapy (a)", b = "CAR T phase (b)",
                 ab_joint = "a + b, both needed", ab_redundant = "a + b, either suffices",
                 other = "Other")
lineage_labels <- c(anc = "Neutrophils (ANC)", plt = "Platelets", amc = "Monocytes (AMC)")
model_labels <- c(M0 = "M0 progenitors only", Ma = "Ma chemotherapy", Mb = "Mb CAR T phase",
                  Mab = "Mab additive", Mabx = "Mabx synergy")

scale_colour_pheno <- function(cfg, ...) ggplot2::scale_colour_manual(values = cfg$plot$col_phenotype,
                                                                      labels = pheno_labels, name = NULL, ...)
scale_fill_pheno <- function(cfg, ...) ggplot2::scale_fill_manual(values = cfg$plot$col_phenotype,
                                                                  labels = pheno_labels, name = NULL, ...)
scale_fill_mech <- function(cfg, ...) ggplot2::scale_fill_manual(values = cfg$plot$col_mech,
                                                                 labels = mech_labels, name = NULL, ...)
scale_colour_lineage <- function(cfg, ...) ggplot2::scale_colour_manual(values = cfg$plot$col_lineage,
                                                                        labels = lineage_labels, name = NULL, ...)
