# =============================================================================
# an_PR_prespecify.R
# Purpose : PR - pre-register predictions of the post-treatment clone size for
#           the patients with a planned post-treatment sample, BEFORE those
#           results are read. run_all.R writes the files and their SHA-256.
#           D5b later verifies the hashes and compares with the measurements.
# Inputs  : ctx$data() (all splits: the planned patients are known clinically),
#           ctx$results$D4.
# Outputs : list(prediction, prereg = list(files, hash))
# =============================================================================

an_PR_prespecify <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  pp <- d$post_planned
  if (is.null(pp) || !nrow(pp)) return(list(tables = list(), note = "no planned post-treatment samples"))
  fit <- ctx$results$D4$fits$M2
  s <- cfg$steps$PR
  sub <- d
  sub$patients <- d$patients[patient_id %in% pp$patient_id]
  pds <- build_pdats(sub, cfg, use_vaf = TRUE)
  fits <- fit_given_population(pds, fit, cfg, use_vaf = TRUE)
  set.seed(cfg$runtime$seed + 7L)
  probs <- sort(unique(c((1 - s$intervals) / 2, 0.5, 1 - (1 - s$intervals) / 2)))
  pred <- data.table::rbindlist(lapply(seq_len(nrow(pp)), function(k) {
    pid <- pp$patient_id[k]
    sday <- pp$sample_day[k]
    f <- fits[[pid]]
    prow <- d$patients[patient_id == pid]
    seg <- make_segments(cfg$model$time$t0, max(cfg$model$time$tmax, sday), prow$ld_scale,
                         d$events[patient_id == pid])
    draws <- MASS::mvrnorm(s$n_draws, f$eta, f$cov)
    if (is.null(dim(draws))) draws <- matrix(draws, nrow = 1, dimnames = list(NULL, names(f$eta)))
    f_t <- apply(draws, 1, function(e) {
      e <- stats::setNames(e, names(f$eta))
      m <- marrow_sim(build_par(eta_named_to_indiv(e), fit$pop$struct, "M2"), seg, c(-5, sday))
      m[2, "f_clone"]
    })
    vaf_pred <- f_t / 2 * exp(stats::rnorm(length(f_t), 0, cfg$model$vaf_noise_sd))
    q <- stats::quantile(vaf_pred, probs, na.rm = TRUE)
    pre <- if (!is.null(d$vaf) && pid %in% d$vaf$patient_id) d$vaf[patient_id == pid, vaf] else NA_real_
    out <- data.table::data.table(patient_id = pid, split = prow$split, sample_day = sday,
                                  pre_vaf = pre, f0_map = stats::plogis(f$eta[["f0"]]),
                                  p_vaf_above_pre = if (is.na(pre)) NA_real_ else round(mean(vaf_pred > pre), 3))
    out[, paste0("vaf_q", formatC(100 * probs, format = "g")) := as.list(round(q, 4))]
    out
  }))

  int_txt <- paste0(100 * s$intervals, "%", collapse = " and ")
  md <- c(
    "# Pre-registered prediction: post-treatment clone size",
    "",
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " (data source: ", cfg$data_source, ")"),
    "",
    "Model: M2 population fit from D4 (training split); each patient fitted individually with",
    "all available counts, inflammatory markers and the pre-treatment VAF.",
    paste0("Prediction intervals (", int_txt, ") include VAF measurement noise (log-SD ",
           cfg$model$vaf_noise_sd, ")."),
    "",
    "## Decision rule (fixed now)",
    "- Supports the model: measured VAF inside the 80% interval, or inside the 95% interval",
    "  with the predicted direction relative to the pre-treatment VAF.",
    "- Contradicts the model: measured VAF outside the 95% interval in the opposite direction.",
    "",
    "## Selection bias (stated in advance)",
    "These patients have a post-treatment sample because their counts stayed low and they had a",
    "marrow work-up. They are therefore enriched for the clonal hypothesis (H3) and cannot be",
    "used to estimate how often clones expand; they test only whether the predicted direction",
    "and size are right in this kind of patient.",
    "",
    "## Predictions",
    knitr::kable(pred, format = "pipe"))
  list(prediction = pred,
       tables = list(prediction = pred),
       prereg = list(files = list(PR_prediction.csv = pred, PR_prediction.md = md), hash = "PR_prediction.sha256"),
       key_numbers = list(prediction = pred))
}

# D5b - compare the pre-registered predictions with the later measurements.
an_D5b_compare <- function(ctx) {
  cfg <- ctx$config
  pred <- ctx$prereg_prediction
  post <- ctx$post_mutations
  d <- ctx$data()
  lo80 <- "vaf_q10"; hi80 <- "vaf_q90"; lo95 <- "vaf_q2.5"; hi95 <- "vaf_q97.5"
  tab <- data.table::rbindlist(lapply(seq_len(nrow(pred)), function(k) {
    pid <- pred$patient_id[k]
    pre_gene <- if (!is.null(d$vaf) && pid %in% d$vaf$patient_id) d$vaf[patient_id == pid, gene] else NA_character_
    m <- post[patient_id == pid]
    meas <- if (!nrow(m)) 0 else if (!is.na(pre_gene) && pre_gene %in% m$gene) m[gene == pre_gene, max(vaf)] else max(m$vaf)
    p <- pred[k]
    in80 <- meas >= p[[lo80]] & meas <= p[[hi80]]
    in95 <- meas >= p[[lo95]] & meas <= p[[hi95]]
    dir_pred <- if (is.na(p$pre_vaf)) NA else p$vaf_q50 > p$pre_vaf
    dir_obs <- if (is.na(p$pre_vaf)) NA else meas > p$pre_vaf
    verdict <- if (in80 || (in95 && isTRUE(dir_pred == dir_obs))) "supports"
               else if (!in95 && isTRUE(dir_pred != dir_obs)) "contradicts" else "inconclusive"
    data.table::data.table(patient_id = pid, sample_day = p$sample_day, pre_vaf = p$pre_vaf,
                           measured_vaf = meas, pred_median = p$vaf_q50, pred_lo80 = p[[lo80]],
                           pred_hi80 = p[[hi80]], pred_lo95 = p[[lo95]], pred_hi95 = p[[hi95]],
                           in_80 = in80, in_95 = in95, direction_correct = dir_pred == dir_obs,
                           verdict = verdict)
  }))
  list(table = tab, verification = ctx$prereg_verification,
       tables = list(comparison = tab),
       key_numbers = list(comparison = tab))
}
