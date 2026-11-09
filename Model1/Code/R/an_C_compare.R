# =============================================================================
# an_C_compare.R
# Purpose : Step C - one comparison table with the same metrics for Aim 4, its
#           sub-models and the Friberg reference, so the later choice between
#           Aim 4 and Aim 2 + Aim 3 rests on common numbers. Results of Aim 2
#           and Aim 3 on the real data are read (never written) when present.
# Inputs  : ctx (config, results F1, F3).
# Outputs : list(table, external, tables)
# =============================================================================

an_C_compare <- function(ctx) {
  cfg <- ctx$config
  f1 <- ctx$results$F1
  f3 <- ctx$results$F3
  cmp <- data.table::copy(f1$comparison)
  cmp[, role := data.table::fcase(model == f1$selected$final, "final (Aim 4)",
                                  model == f1$selected$q2, "selected Q2",
                                  model == "MF", "reference (Friberg 2002)",
                                  default = "")]
  m3 <- f3$metrics[group == "all", .(L, auroc_model, auroc_lo, auroc_hi, auroc_hematotox, auroc_logistic, brier_3class)]
  # other aims (read-only; absent until they are fitted to the real data)
  ext_paths <- c(Aim2 = file.path(cfg$paths$root, "..", "Aim2", "Results", "real", "tables"),
                 Aim3 = file.path(cfg$paths$root, "..", "Aim3", "Results", "real", "tables"))
  external <- data.table::rbindlist(lapply(names(ext_paths), function(a) {
    p <- ext_paths[[a]]
    files <- if (dir.exists(p)) list.files(p, pattern = "comparison|metrics|BIC", full.names = TRUE) else character(0)
    data.table::data.table(aim = a, available = length(files) > 0,
                           files = if (length(files)) paste(basename(files), collapse = "; ") else "not available yet")
  }))
  list(table = cmp, landmark = m3, external = external,
       tables = list(models = cmp, landmark = m3, other_aims = external))
}
