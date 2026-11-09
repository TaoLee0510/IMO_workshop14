# =============================================================================
# report_summary.R
# Purpose : Compose Results/<source>/summary.md from the saved step results:
#           key numbers per deliverable plus the list of figures and tables.
# Inputs  : named list of step results (loaded from objects/), config.
# Outputs : character vector (markdown lines); run_all.R writes the file.
# =============================================================================

md_table <- function(x) {
  if (is.null(x) || !nrow(x)) return("_(none)_")
  knitr::kable(x, format = "pipe")
}

report_summary <- function(results, config) {
  r <- results
  L <- c(paste0("# Aim 2 results summary (", config$data_source, ")"), "",
         paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), ". Quick mode: ", config$quick, "."))
  if (config$data_source != "real")
    L <- c(L, "", "> Synthetic TestData: these numbers check the pipeline, they are not patient findings.")
  if (identical(config$phenotype$scheme, "early"))
    L <- c(L, "", paste0("Phenotype scheme: day 0-", config$phenotype$max_day, " (codes quick / second_wave / persistent = ",
                         paste(phenotype_display(config$phenotype), collapse = " / "), "). Models: ",
                         paste(config$fit$models, collapse = ", "), ", fitted to day ", config$fit$obs_max_day, "."))

  if (!is.null(r$S0)) L <- c(L, "", "## S0 Data", "",
    "Late follow-up density (decides how P1 can be tested):", "", md_table(r$S0$tables$late_density), "",
    "Phenotypes, training split:", "", md_table(r$S0$tables$phenotype_train),
    if (!is.null(r$S0$tables$phenotype_vs_clinical))
      c("", "Rule phenotype vs the phenotype recorded in the data:", "", md_table(r$S0$tables$phenotype_vs_clinical)))

  if (!is.null(r$D1)) L <- c(L, "", "## D1 Model tests", "",
    sprintf("Unit tests passed: %d / %d.", r$D1$key_numbers$tests_passed, r$D1$key_numbers$tests_total), "",
    md_table(r$D1$tables$archetypes))

  if (!is.null(r$D2)) L <- c(L, "", "## D2 Phase map", "", "Clone effect by facet (share of grid cells with f0 > 0):", "",
    md_table(r$D2$tables$clone_effect))

  if (!is.null(r$D3)) L <- c(L, "", "## D3 Identifiability", "",
    md_table(r$D3$tables$mechanism_accuracy), "", "Persistent patients, H1 vs H3:", "",
    md_table(r$D3$tables$persistent_h1_vs_h3), "", "Clone-fraction recovery among carriers:", "",
    md_table(r$D3$tables$f0_recovery))

  if (!is.null(r$D3b)) L <- c(L, "", "## D3b Clone traits unknown", "",
    "Profile likelihood over (delta, r_inf) in virtual cohorts; truth = low-output clone. ",
    "`llr_low_vs_functional` > 0 means the data favour the true clone type.", "",
    md_table(r$D3b$summary))

  if (!is.null(r$D3c)) L <- c(L, "", "## D3c What the current data can say", "",
    "TestData-scale cohorts re-simulated and refitted with every parameter estimated (D4 procedure). ",
    "`llr_*` = log-likelihood ratio, low-output vs functional clone (> 0 = true clone type favoured).",
    if (!is.null(r$D4$fits$M2$profile)) {
      pt <- r$D4$fits$M2$profile
      low <- pt$delta <= r$D3c$settings$low_output_max
      sprintf("For comparison, D4 on this dataset: LLR %.1f.", max(pt$logLik[low]) - max(pt$logLik[!low]))
    }, "",
    md_table(r$D3c$summary))

  if (!is.null(r$D3d)) L <- c(L, "", "## D3d Cohorts at the real-data scale", "",
    paste0("Cohorts of ", r$D3d$settings$n, " patients simulated with the real sampling (counts to day 30, one late draw), ",
           "treatment days hidden and approximated as in the real data, fitted with M0 and M1. ",
           "`dbic_m1` = BIC(M0) - BIC(M1) (> 0 favours the inflammation layer); variant no_inflam has no marrow ",
           "suppression in the truth."), "", md_table(r$D3d$summary))

  if (!is.null(r$D4)) L <- c(L, "", "## D4 Population fit (training split)", "",
    md_table(r$D4$compare), "", paste0("Model chosen for the hold-out: **", r$D4$chosen, "**."),
    if (!is.null(r$D4$P1)) c("", "P1 late shape x attributed mechanism:", "", md_table(r$D4$P1$cross)), "",
    "P2 timing:", "", md_table(r$D4$P2$summary), "",
    "P3 platelets:", "", md_table(r$D4$P3$table), "", md_table(r$D4$P3$test),
    if (!is.null(r$D4$P4)) c("", "P4 risk factors for clone fraction (covariate-free fit):", "", md_table(r$D4$P4)), "",
    "P5 log-likelihood gain by phenotype:", "", md_table(r$D4$P5$table),
    if (!is.null(r$D4$mech_summary)) c("",
      "Mechanism of persistent patients: MAP calls and expected counts from posterior P(H3) (a lower bound):", "",
      md_table(r$D4$mech_summary)),
    if (!is.null(r$D4$event_sensitivity)) c("",
      "Treatment-day sensitivity (M1 refitted; correlations with the main estimates):", "",
      md_table(r$D4$event_sensitivity)),
    if (!is.null(r$D4$late$summary)) c("",
      "Counts after the fitting window vs forward prediction:", "", md_table(r$D4$late$summary)))

  if (!is.null(r$D5) && !is.null(r$D5$table)) L <- c(L, "", "## D5 CH case series (descriptive)", "",
    md_table(r$D5$concordance), "", md_table(r$D5$table[, .(patient_id, gene, gene_class, vaf,
                                                            f0_without_vaf = round(f0_without_vaf, 3),
                                                            f0_with_vaf = round(f0_with_vaf, 3),
                                                            phenotype, mechanism_with_vaf)]))

  if (!is.null(r$PR) && !is.null(r$PR$prediction)) L <- c(L, "", "## PR Pre-registered predictions", "",
    md_table(r$PR$prediction))
  if (!is.null(r$D5b)) L <- c(L, "", "## D5b Prediction vs measurement", "",
    paste0("Pre-registration check: ", r$D5b$verification$msg), "", md_table(r$D5b$table))

  if (!is.null(r$D6)) L <- c(L, "", "## D6 Virtual cohort", "", md_table(r$D6$tables$attributable_fraction), "",
    "Effect size of clonal expansion (persistent, clone detectable at baseline):", "",
    md_table(r$D6$tables$effect_size), "", "Paired-design power:", "", md_table(r$D6$tables$power))

  if (!is.null(r$H) && identical(r$H$mode, "cv")) L <- c(L, "",
    paste0("## H Early prediction (", config$steps$H$cv_folds, "-fold cross-validation, day +",
           config$steps$H$landmark_day, ")"), "",
    paste0("Patients: ", r$H$n_test, "; model: ", r$H$plan$chosen_model,
           " (population refitted in every fold, started from the D4 estimates)."), "",
    md_table(r$H$metrics), "", md_table(r$H$diffs), "", "Accuracy (most probable class):", "", md_table(r$H$accuracy))
  else if (!is.null(r$H)) L <- c(L, "", "## H Hold-out (single pre-registered evaluation)", "",
    paste0("Test patients evaluated: ", r$H$n_test, "; model: ", r$H$plan$chosen_model, "."), "",
    md_table(r$H$metrics), "", md_table(r$H$diffs))

  if (!is.null(r$V)) L <- c(L, "", "## V Validation against the generating truth", "",
    md_table(r$V$tables$parameter_recovery), "", md_table(r$V$tables$structural), "",
    paste0("Observed-rule phenotype accuracy vs truth: ", r$V$phenotype_accuracy), "",
    md_table(r$V$tables$persistent_mechanism_confusion), "",
    "Posterior P(H3) by true mechanism (observed-persistent patients):", "",
    md_table(r$V$tables$mechanism_probability), "", md_table(r$V$tables$model_selection))

  figs <- list.files(file.path(config$paths$results, "figures"), pattern = "\\.pdf$")
  tabs <- list.files(file.path(config$paths$results, "tables"), pattern = "\\.csv$")
  c(L, "", "## Files", "", paste0("Figures (", length(figs), "): ", paste0("`", figs, "`", collapse = ", ")), "",
    paste0("Tables (", length(tabs), "): ", paste0("`", tabs, "`", collapse = ", ")))
}
