# =============================================================================
# an_D3b_traits.R
# Purpose : D3b - identifiability when the clone traits are NOT known.
#           D3 fixes the clone traits at their true values; in real data they
#           are population parameters to be estimated, and D4 showed a ridge
#           between clone output (delta) and clone sensitivity to inflammation
#           (r_inf). Here virtual cohorts are drawn from the design population
#           and the profile likelihood over a (delta, r_inf) grid is computed
#           under four information scenarios:
#             counts      blood counts + CRP/ferritin only
#             pre_subset  + pre-treatment VAF in n_pre_subset patients (today)
#             pre_all     + pre-treatment VAF in every patient
#             paired      + pre-treatment VAF in every patient and a day +30
#                         VAF by the stratified rule (ANC < thr on that day,
#                         plus a random fraction of the others)
#           All scenarios use the same simulated cohorts; only the information
#           passed to the fit differs.
# Inputs  : ctx$config (design values; no patient data).
# Outputs : list(surface, summary, tables, key_numbers)
#
# Metrics per cohort x scenario
#   llr_low_vs_functional : max logLik over low-output points (delta <=
#                           low_output_max) minus max over the others; the
#                           design truth is low-output, so > 0 is correct
#   region_size           : grid points within region_ll of the maximum
#   argmax_is_truth       : the profile maximum sits on the true traits
#   h1_vs_h3_*            : mechanism accuracy for truly persistent patients
#                           with traits at the profile maximum (estimated)
#                           versus at the true traits (oracle)
# =============================================================================

d3b_labels <- c(counts = "Counts + CRP", pre_subset = "+ pre VAF (subset)",
                pre_all = "+ pre VAF (all)", paired = "+ pre VAF + day +30 VAF")

# Mechanism call for every patient given one eta matrix and one set of traits.
d3b_mechanisms <- function(eta, vp, struct, cfg) {
  unlist(par_lapply(seq_along(vp), function(i) {
    e <- eta[i, ]
    infer_mechanism(eta_named_to_indiv(e), struct, "M2", vp[[i]]$base$anc, 1, cfg, vp[[i]]$events)$mechanism
  }, cfg$runtime$cores))
}

d3b_cohort <- function(r, cfg, st, s, dp) {
  n <- s$n_virtual
  truth <- d3_truth_cohort(n, cfg, st, dp, prefix = sprintf("R%dV", r))
  vp <- lapply(seq_len(n), function(i)
    d3_virtual_patient(truth$id[i], c(rho = truth$rho[i], g = truth$g[i], f0 = truth$f0[i]), st, cfg,
                       d3_schedule(cfg$model$design$clinical)))
  truth <- d3_add_truth_mechanism(truth, vp, st, cfg)
  # day +30 sampling is decided from the most recent ANC in the week before
  anc_d <- vapply(vp, function(v) {
    o <- v$obs[var == "anc" & day <= s$post_vaf_day & day > s$post_vaf_day - 7]
    if (nrow(o)) o$value[which.max(o$day)] else NA_real_
  }, numeric(1))
  truth[, `:=`(take_post = (!is.na(anc_d) & anc_d < s$post_rule$anc_thr) | stats::runif(n) < s$post_rule$random_frac,
               pre_subset = seq_len(n) %in% sample(n, s$n_pre_subset))]
  truth[, post_detected := take_post & !vapply(vp, function(v) is.null(v$vaf_post), logical(1))]
  list(truth = truth, vp = vp)
}

d3b_scenario_pdats <- function(coh, scenario, cfg) {
  tr <- coh$truth
  pds <- lapply(seq_along(coh$vp), function(i) {
    v <- coh$vp[[i]]
    has_pre <- scenario %in% c("pre_all", "paired") || (scenario == "pre_subset" && tr$pre_subset[i])
    d3_pdat(v, cfg, vaf = if (has_pre) v$vaf_pre else NULL,
            vaf_post = if (scenario == "paired" && tr$take_post[i]) v$vaf_post else NULL)
  })
  stats::setNames(pds, tr$id)
}

an_D3b_traits <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$D3b
  cfg$s3_post_day <- s$post_vaf_day
  st <- design_struct(cfg)
  cal <- calibrate_design_pop(cfg, st)
  pop <- d3_design_pop(cfg, st, cal$pop)
  grid <- do.call(data.table::CJ, s$grid)
  truth_traits <- st[names(s$grid)]
  is_truth <- Reduce(`&`, lapply(names(s$grid), function(k) abs(grid[[k]] - truth_traits[[k]]) < 1e-9))
  set.seed(cfg$runtime$seed + 31L)

  surfaces <- list()
  rows <- list()
  for (r in seq_len(s$n_rep)) {
    coh <- d3b_cohort(r, cfg, st, s, cal$pop)
    tr <- coh$truth
    pers <- which(tr$mechanism_true %in% c("H1", "H3"))
    for (sc in s$scenarios) {
      t0 <- Sys.time()
      pds <- d3b_scenario_pdats(coh, sc, cfg)
      pr <- profile_grid(pds, "M2", pop, grid, cfg, use_vaf = sc != "counts", keep_eta = TRUE)
      tab <- pr$table
      j_max <- which.max(tab$logLik)
      low <- tab$delta <= s$low_output_max
      st_hat <- st
      st_hat[names(pr$best)] <- pr$best
      mech_est <- d3b_mechanisms(pr$eta[[j_max]][pers, , drop = FALSE], coh$vp[pers], st_hat, cfg)
      mech_orc <- if (any(is_truth)) d3b_mechanisms(pr$eta[[which(is_truth)]][pers, , drop = FALSE], coh$vp[pers], st, cfg)
                  else rep(NA_character_, length(pers))
      ch <- which(tr$ch)
      f0_est <- stats::plogis(pr$eta[[j_max]][ch, "f0"])
      rows[[length(rows) + 1]] <- data.table::data.table(
        rep = r, scenario = sc, n = nrow(tr),
        n_pre_vaf = if (sc == "counts") 0L else if (sc == "pre_subset") sum(tr$pre_subset) else nrow(tr),
        n_post_assays = if (sc == "paired") sum(tr$take_post) else 0L,
        n_post_detected = if (sc == "paired") sum(tr$post_detected) else 0L,
        argmax_delta = pr$best[["delta"]], argmax_r_inf = pr$best[["r_inf"]],
        argmax_is_truth = is_truth[j_max],
        llr_low_vs_functional = max(tab$logLik[low]) - max(tab$logLik[!low]),
        region_size = sum(tab$logLik >= max(tab$logLik) - s$region_ll),
        truth_in_region = any(is_truth & tab$logLik >= max(tab$logLik) - s$region_ll),
        n_persistent = length(pers), n_h3 = sum(tr$mechanism_true[pers] == "H3"),
        h1_vs_h3_estimated = mean(mech_est == tr$mechanism_true[pers]),
        h1_vs_h3_oracle = mean(mech_orc == tr$mechanism_true[pers]),
        h3_sensitivity_estimated = mean(mech_est[tr$mechanism_true[pers] == "H3"] == "H3"),
        cor_logit_f0 = stats::cor(stats::qlogis(f0_est), stats::qlogis(tr$f0[ch])))
      surfaces[[length(surfaces) + 1]] <- tab[, `:=`(rep = r, scenario = sc, is_truth = is_truth)]
      ctx$log(sprintf("  D3b: cohort %d, %-10s argmax delta=%g r_inf=%g, LLR(low vs functional) %.1f (%.0fs)",
                      r, sc, pr$best[["delta"]], pr$best[["r_inf"]], rows[[length(rows)]]$llr_low_vs_functional,
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }
  per_cohort <- data.table::rbindlist(rows)
  surface <- data.table::rbindlist(surfaces)
  num <- c("llr_low_vs_functional", "region_size", "h1_vs_h3_estimated", "h1_vs_h3_oracle",
           "h3_sensitivity_estimated", "cor_logit_f0", "n_post_assays", "n_post_detected")
  summ <- per_cohort[, c(list(cohorts = .N, n_per_cohort = n[1], n_pre_vaf = n_pre_vaf[1],
                              clone_type_correct = mean(llr_low_vs_functional > 0),
                              argmax_is_truth = mean(argmax_is_truth),
                              truth_in_region = mean(truth_in_region)),
                         lapply(.SD, function(x) round(mean(x, na.rm = TRUE), 3))),
                     by = scenario, .SDcols = num]
  summ[, scenario := factor(scenario, s$scenarios)]
  data.table::setorder(summ, scenario)
  summ[, scenario := as.character(scenario)]
  surf_mean <- surface[, .(delta_logLik = mean(delta_logLik)), by = .(scenario, delta, r_inf, is_truth)]

  list(surface = surface, surface_mean = surf_mean, per_cohort = per_cohort, summary = summ,
       truth_traits = truth_traits, settings = s,
       tables = list(design_calibration = cal$trace, summary = summ, per_cohort = per_cohort, profile_surface = surface),
       key_numbers = list(summary = summ))
}
