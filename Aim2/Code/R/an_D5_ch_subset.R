# =============================================================================
# an_D5_ch_subset.R
# Purpose : D5 - CH case series (training split). Each sequenced patient is
#           fitted twice with the D4 population: with the pre-treatment VAF
#           (best fit) and without it (leave-VAF-out). Concordance between the
#           leave-VAF-out estimate and the measured clone fraction avoids
#           circularity. Descriptive only: n is 10-20.
# Inputs  : ctx$data(), ctx$results$D4.
# Outputs : list(table, concordance, tables)
# =============================================================================

an_D5_ch_subset <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  tr <- analysis_set(d, "train")
  if (is.null(tr$vaf) || !nrow(tr$vaf)) return(list(tables = list(), table = NULL, note = "no sequenced patients"))
  fit <- ctx$results$D4$fits$M2
  ids <- tr$vaf$patient_id
  sub <- tr
  sub$patients <- tr$patients[patient_id %in% ids]
  pd_with <- build_pdats(sub, cfg, use_vaf = TRUE)
  pd_without <- build_pdats(sub, cfg, use_vaf = FALSE)
  f_with <- fit_given_population(pd_with, fit, cfg, use_vaf = TRUE)
  f_without <- fit_given_population(pd_without, fit, cfg, use_vaf = FALSE)
  z <- stats::qnorm(0.9)
  est <- function(f) {
    e <- f$eta[["f0"]]
    s <- sqrt(f$cov["f0", "f0"])
    c(est = stats::plogis(e), lo = stats::plogis(e - z * s), hi = stats::plogis(e + z * s))
  }
  tab <- data.table::rbindlist(lapply(ids, function(pid) {
    a <- est(f_without[[pid]])
    b <- est(f_with[[pid]])
    indiv <- eta_named_to_indiv(f_with[[pid]]$eta)
    prow <- tr$patients[patient_id == pid]
    mech <- infer_mechanism(indiv, fit$pop$struct, "M2", prow$anc_base, prow$ld_scale, cfg)
    v <- tr$vaf[patient_id == pid]
    data.table::data.table(patient_id = pid, status = v$status, gene = v$gene, gene_class = v$gene_class,
                           sample_timepoint = v$sample_timepoint, vaf = v$vaf, f_measured = v$f_measured,
                           f0_without_vaf = a[["est"]], f0_without_lo = a[["lo"]], f0_without_hi = a[["hi"]],
                           f0_with_vaf = b[["est"]], f0_with_lo = b[["lo"]], f0_with_hi = b[["hi"]],
                           mechanism_with_vaf = mech$mechanism)
  }))
  tab <- merge(tab, tr$pheno[, .(patient_id, phenotype)], by = "patient_id", all.x = TRUE)

  # concordance (leave-VAF-out vs measured); negatives placed at the detection limit
  x <- tab[, .(est = f0_without_vaf, meas = ifelse(status == "positive", f_measured, 2 * cfg$data$lod$vaf / 2))]
  conc <- data.table::data.table(n = nrow(x), spearman = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  if (nrow(x) >= 5) {
    set.seed(cfg$runtime$seed + 5L)
    r <- suppressWarnings(stats::cor(x$est, x$meas, method = "spearman"))
    bs <- replicate(cfg$steps$D5$n_boot, {
      i <- sample(nrow(x), replace = TRUE)
      suppressWarnings(stats::cor(x$est[i], x$meas[i], method = "spearman"))
    })
    conc <- data.table::data.table(n = nrow(x), spearman = round(r, 3),
                                   ci_low = round(stats::quantile(bs, 0.025, na.rm = TRUE), 3),
                                   ci_high = round(stats::quantile(bs, 0.975, na.rm = TRUE), 3))
  }
  list(table = tab, concordance = conc,
       tables = list(case_series = tab, concordance = conc),
       key_numbers = list(concordance = conc, n = nrow(tab)))
}
