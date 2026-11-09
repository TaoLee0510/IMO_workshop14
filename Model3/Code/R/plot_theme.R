# =============================================================================
# plot_theme.R
# Purpose : Shared look for every Aim 4 figure (palette of the team slides:
#           crimson 9E2A2B, steel blue 2E6B8A, slate, gold, green).
# Inputs  : config$plot.
# Outputs : ggplot theme and scale helpers.
# =============================================================================

theme_model3 <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#1A2433", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(colour = "#4A5568"),
      plot.caption = ggplot2::element_text(colour = "#6B7A8F", hjust = 0, size = base_size - 2),
      strip.text = ggplot2::element_text(face = "bold", colour = "#1A2433"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom")
}

scale_colour_pheno <- function(cfg, ...) {
  ggplot2::scale_colour_manual(values = cfg$plot$col_phenotype, labels = phenotype_display(), name = NULL, ...)
}
scale_fill_pheno <- function(cfg, ...) {
  ggplot2::scale_fill_manual(values = cfg$plot$col_phenotype, labels = phenotype_display(), name = NULL, ...)
}

when_labels <- c(none = "Quick (HSC not limiting)", a = "Chemotherapy (a)", b = "CAR T phase (b)",
                 ab_joint = "a + b, both needed", ab_redundant = "a + b, either suffices", other = "Other")
who_labels <- c(clone = "Clone", reserve = "Low reserve", both = "Both", other = "Other")

data_caption <- function(cfg) {
  if (cfg$data_source == "real") "Real cohort (n = 35, de-identified)." else
    "Synthetic TestData (real sampling design, simulated counts) - illustrates the pipeline, not patient results."
}

plot_none <- function(res, cfg, ...) list()
