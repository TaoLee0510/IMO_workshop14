# =============================================================================
# report_summary.R
# Purpose : summary.md with the key numbers of every step that has results.
# Inputs  : named list of step results (NULL when not run), config.
# Outputs : character vector (markdown lines).
# =============================================================================

report_summary <- function(r, cfg) {
  f <- function(x, d = 2) formatC(x, format = "f", digits = d)
  out <- c(paste0("# Aim 4 results summary (", cfg$data_source, ")"), "",
           paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"), ". ", data_caption(cfg)), "")
  if (!is.null(r$S0)) {
    pc <- r$S0$counts[rule == "primary"]
    out <- c(out, "## S0 data", paste0("- ", r$S0$n, " patients; phenotype quick / second wave / delayed = ",
                                       pc$quick, " / ", pc$second_wave, " / ", pc$delayed,
                                       "; stable across six rules: ", r$S0$stable),
             if (!is.na(r$S0$all_match)) paste0("- Matches the Aim 3 data review for all six rules: ", r$S0$all_match), "")
  }
  if (!is.null(r$K1)) out <- c(out, "## K1 model", paste0("- Unit tests passed: ", sum(r$K1$tests$passed), "/", nrow(r$K1$tests)), "")
  if (!is.null(r$K2)) out <- c(out, "## K2 maps", paste0("- Pseudo second wave (single chemotherapy hit) needs HSC loss >= ",
                                                         paste(sprintf("%s (gamma_M %s)", f(r$K2$thr$threshold_hsc_loss), r$K2$thr$gamma_m), collapse = "; ")), "")
  if (!is.null(r$K3)) {
    ro <- r$K3$routes
    out <- c(out, "## K3 identifiability", paste0("- ", r$K3$n_parts, " synthetic cohorts"),
             paste0("- ", ro$question, ": cohort-level correct ", f(ro$cohort_correct), ", individual mean P(true) ",
                    f(ro$indiv_mean_p_true), " -> route ", ro$route), "")
  }
  if (!is.null(r$F1)) {
    s <- r$F1$selected
    out <- c(out, "## F1 model comparison",
             paste0("- Q2 selected: **", s$q2, "** (runner-up ", s$q2_runner_up, ", +", f(s$q2_margin, 1), " BIC; supported: ", s$q2_supported, ")"),
             paste0("- Q1: ", s$q1_best_clone, " vs ", s$q2, ": dBIC ", f(s$q1_delta_bic, 1), " (supported: ", s$q1_supported, ")"),
             paste0("- Friberg reference MF: +", f(s$mf_delta_bic, 1), " BIC vs ", s$q2),
             paste0("- T1 slope (second-wave nadir vs ferritin peak): ", f(r$F1$tests$T1$fit$slope), " [",
                    f(r$F1$tests$T1$fit$lo), ", ", f(r$F1$tests$T1$fit$hi), "], n = ", r$F1$tests$T1$fit$n), "")
  }
  if (!is.null(r$F2)) {
    co <- r$F2$cohort
    out <- c(out, "## F2 attribution (expected numbers of patients)",
             paste0("- ", co$readout, " / ", co$class, ": ", f(co$expected_n, 1)), "",
             paste0("- Brier (3-class) of the posterior phenotype: ", f(r$F2$brier)), "")
  }
  if (!is.null(r$F3)) {
    m <- r$F3$metrics[group == "all"]
    out <- c(out, "## F3 landmark prediction (LOO, all patients)",
             paste0("- L = ", m$L, ": AUROC ", f(m$auroc_model), " [", f(m$auroc_lo), ", ", f(m$auroc_hi), "]; CAR-HEMATOTOX ",
                    f(m$auroc_hematotox), "; counts-logistic ", f(m$auroc_logistic), "; Brier ", f(m$brier_3class)), "")
  }
  if (!is.null(r$F4)) out <- c(out, "## F4 what-if (mean change in P(non-quick))",
                               paste0("- ", r$F4$overall$intervention, ": ", sprintf("%+.3f", r$F4$overall$mean_delta_nonquick)), "")
  if (!is.null(r$V) && !is.null(r$V$summary)) out <- c(out, "## V validation (TestData)",
                                                         paste0("- Selected ", r$V$summary$selected_q2, " (truth: ", r$V$summary$true_mechanism,
                                                                "); mean P(true when) ", f(r$V$summary$mean_p_true_when),
                                                                "; mean P(true who) ", f(r$V$summary$mean_p_true_who)), "")
  out
}
