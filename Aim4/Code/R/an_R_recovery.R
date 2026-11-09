# =============================================================================
# an_R_recovery.R
# Purpose : Step R (v1.3, replaces K3) - recovery of the final structure from
#           synthetic real35 cohorts (sim_cohort.R simulate_cohort13) fitted
#           with DE-ITS on the cluster: structural parameters vs truth, random
#           effects and baselines vs truth (correlations), and the HSC loss
#           through the two CAR T-phase routes (ferritin-driven I vs CAR T -
#           target signal Z) at the MAP vs truth.
# Inputs  : ctx (config); Results/<data>_v13/recovery/cohorts/synNN.rds and the
#           copied-back fits Results/<data>_v13/hpc_runs/synNN/fit.rds.
# Outputs : list(struct, indiv, routes, summary, tables)
# =============================================================================

an_R_recovery <- function(ctx) {
  cfg <- ctx$config
  cdir <- file.path(cfg$paths$results, "recovery", "cohorts")
  rdir <- file.path(cfg$paths$results, "hpc_runs")
  cohorts <- sub("\\.rds$", "", list.files(cdir, pattern = "^syn[0-9]+\\.rds$"))
  have <- cohorts[file.exists(file.path(rdir, cohorts, "fit.rds"))]
  if (!length(have)) stop("no recovery fits found under ", rdir)
  per <- lapply(have, function(cn) {
    co <- readRDS(file.path(cdir, paste0(cn, ".rds")))
    fit <- readRDS(file.path(rdir, cn, "fit.rds"))
    est <- struct_est_names(fit$spec, cfg)
    st <- data.table::data.table(cohort = cn, parameter = est, truth = unname(co$truth_struct[est]),
                                 estimate = unname(fit$pop$struct[est]),
                                 prior_centre = unname(default_struct(fit$spec, cfg)[est]))
    psi <- fit_psi13(fit)
    tr <- co$truth
    ind <- data.table::data.table(cohort = cn, patient_id = rownames(psi))
    for (k in fit$re) {
      tv <- tr[match(ind$patient_id, patient_id)][[paste0("true_", k)]]
      ind[, (paste0("true_", k)) := if (k == "f0") stats::qlogis(tv) else log(tv)]
      ind[, (paste0("est_", k)) := psi[, k]]
    }
    pts <- co$patients
    for (v in cfg$model$baseline$vars) {
      ind[, (paste0("true_b_", v)) := log(tr[match(ind$patient_id, patient_id)][[v]])]
      ind[, (paste0("est_b_", v)) := log(pts[match(ind$patient_id, patient_id)][[paste0(v, "_base")]]) +
            psi[, paste0("b_", v)]]
    }
    # HSC loss through the two b routes at the MAP (estimated) vs truth
    pdats <- prep_segments(build_pdats(co, cfg), fit$pop$struct, cfg)
    loss <- data.table::rbindlist(lapply(names(pdats), function(pid) {
      pd <- pdats[[pid]]
      x <- tryCatch(counterfactual_draw13(psi[pid, ], pd, fit, cfg), error = function(e) NULL)
      if (is.null(x)) return(NULL)
      data.table::data.table(patient_id = pid, est_f_bi = as.numeric(x[["f_bi"]]), est_f_bz = as.numeric(x[["f_bz"]]),
                             est_f_a = as.numeric(x[["f_a"]]), est_peak_cart = as.numeric(x[["peak_cart"]]))
    }))
    loss[tr, `:=`(f_bi = i.f_bi, f_bz = i.f_bz, f_a = i.f_a, peak_cart = i.peak_cart), on = "patient_id"]
    loss[, cohort := cn]
    conv <- data.table::data.table(cohort = cn, rounds = fit$rounds, converged = fit$converged, BIC = fit$BIC,
                                   laplace_fallback = sum(fit$indiv$logL_source == "Laplace"))
    list(struct = st, indiv = ind, loss = loss, conv = conv)
  })
  struct <- data.table::rbindlist(lapply(per, `[[`, "struct"))
  indiv <- data.table::rbindlist(lapply(per, `[[`, "indiv"))
  loss <- data.table::rbindlist(lapply(per, `[[`, "loss"))
  conv <- data.table::rbindlist(lapply(per, `[[`, "conv"))
  struct[, `:=`(ratio = estimate / truth, prior_ratio = prior_centre / truth)]
  ssum <- struct[, .(n = .N, median_ratio = stats::median(ratio), q10 = stats::quantile(ratio, 0.1),
                     q90 = stats::quantile(ratio, 0.9), prior_over_truth = prior_ratio[1],
                     moved_towards_truth = mean(abs(log(ratio)) < abs(log(prior_ratio)))), by = parameter]
  keys <- sub("^true_", "", grep("^true_", names(indiv), value = TRUE))
  isum <- data.table::rbindlist(lapply(keys, function(k) {
    x <- indiv[[paste0("true_", k)]]; y <- indiv[[paste0("est_", k)]]
    ok <- is.finite(x) & is.finite(y)
    data.table::data.table(parameter = k, n = sum(ok), pearson_r = stats::cor(x[ok], y[ok]),
                           bias = mean(y[ok] - x[ok]), rmse = sqrt(mean((y[ok] - x[ok])^2)))
  }))
  rsum <- loss[, .(n = .N, r_f_bi = stats::cor(f_bi, est_f_bi), r_f_bz = stats::cor(f_bz, est_f_bz),
                   r_f_a = stats::cor(f_a, est_f_a),
                   dominant_route_agreement = mean((f_bz > f_bi) == (est_f_bz > est_f_bi)),
                   r_log_peak_cart = stats::cor(log(peak_cart), log(est_peak_cart)))]
  ctx$log(sprintf("  R: %d/%d cohorts fitted (%d converged); struct median est/true: %s; route agreement %.2f",
                  length(have), length(cohorts), sum(conv$converged),
                  paste(sprintf("%s %.2f", ssum$parameter, ssum$median_ratio), collapse = ", "),
                  rsum$dominant_route_agreement))
  list(struct = struct, indiv = indiv, loss = loss, conv = conv, summary = list(struct = ssum, indiv = isum, routes = rsum),
       tables = list(struct_recovery = ssum, indiv_recovery = isum, route_recovery = rsum, convergence = conv))
}
