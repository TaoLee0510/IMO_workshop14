# =============================================================================
# an_F4_whatif.R
# Purpose : Step F4 - what-if interventions on each patient's posterior (final
#           model), summarised by the patient's most likely timing class. These
#           are hypothesis-generating; the efficacy trade-off of dampening
#           inflammation lies outside this model (GreenLight tumour module).
#   pre-LD (decided at day -6): no bridging chemotherapy; LD at 75% dose
#   day 7:  anti-inflammatory (driver x 0.5 after day 7)
#   day 14: CD34+ boost (S_W + 0.3); TPO-RA from day 14; G-CSF days 14-21
# Inputs  : ctx (config, data, results F1, F2).
# Outputs : list(table, by_class, tables)
# =============================================================================

f4_interventions <- function() list(
  none = list(),
  no_bridging = list(bridging = "none"),
  ld_75 = list(ld_scale = 0.75),
  anti_inflammatory_d7 = list(drv_after = 7, drv_factor = 0.5),
  cd34_d14 = list(cd34 = data.frame(day = 14, add = 0.3)),
  tpo_d14 = list(tpo = cbind(14, 60)),
  gcsf_d14 = list(gcsf = cbind(14, 21)))

f4_apply <- function(pd, iv) {
  des <- pd$design
  drv <- pd$drv
  if (!is.null(iv$bridging)) des$bridging <- iv$bridging
  if (!is.null(iv$ld_scale)) des$ld_scale <- iv$ld_scale
  if (!is.null(iv$cd34)) des$cd34 <- iv$cd34
  if (!is.null(iv$tpo)) des$tpo <- iv$tpo
  if (!is.null(iv$gcsf)) des$gcsf <- iv$gcsf
  if (!is.null(iv$drv_after)) drv <- scale_driver(drv, iv$drv_after, iv$drv_factor)
  list(design = des, drv = drv)
}

an_F4_whatif <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  f1 <- ctx$results$F1
  f2 <- ctx$results$F2
  fit <- f1$fits[[f1$selected$final]]
  pdats <- build_pdats(d, cfg)
  ivs <- f4_interventions()
  days <- 0:90
  n_draw <- min(100L, cfg$fit$attrib_draws)
  res <- data.table::rbindlist(par_lapply(names(pdats), function(pid) {
    pd <- pdats[[pid]]
    dr <- f2$mcmc[[pid]]
    set.seed(cfg$runtime$seed)
    ix <- sample(nrow(dr), min(n_draw, nrow(dr)))
    data.table::rbindlist(lapply(names(ivs), function(nm) {
      a <- f4_apply(pd, ivs[[nm]])
      seg <- make_segments(a$design, fit$pop$struct, max(days), cfg)
      ph <- vapply(ix, function(j) {
        m <- sim_model(build_par(eta_to_indiv(dr[j, ], fit$spec), fit$pop$struct, fit$spec), seg, a$drv,
                       days, fit$spec, cfg)
        if (anyNA(m[, 2])) return(NA_character_)
        phenotype_sim(m[, "t"], m[, "anc_rel"] * pd$base$anc, cfg, with_persistent = TRUE)
      }, character(1))
      data.table::data.table(patient_id = pid, intervention = nm, p_nonquick = mean(ph != "quick", na.rm = TRUE),
                             p_delayed_or_persistent = mean(ph %in% c("delayed", "persistent"), na.rm = TRUE))
    }))
  }, cfg$runtime$cores))
  base <- res[intervention == "none", .(patient_id, base_nonquick = p_nonquick)]
  res <- base[res, on = "patient_id"]
  res[, delta_nonquick := p_nonquick - base_nonquick]
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  cls <- f2$summary[, .(patient_id, when_map = wc[max.col(as.matrix(.SD), ties.method = "first")]),
                    .SDcols = paste0("p_when_", wc)]
  res[cls, when_map := i.when_map, on = "patient_id"]
  by_class <- res[intervention != "none", .(n = .N, mean_delta_nonquick = mean(delta_nonquick, na.rm = TRUE)),
                  by = .(intervention, when_map)]
  overall <- res[intervention != "none", .(mean_delta_nonquick = mean(delta_nonquick, na.rm = TRUE),
                                          expected_nonquick = sum(p_nonquick)), by = intervention]
  ctx$log("  F4: mean change in P(non-quick): ", paste(sprintf("%s %+.2f", overall$intervention,
                                                              overall$mean_delta_nonquick), collapse = ", "))
  list(table = res, by_class = by_class, overall = overall,
       tables = list(patients = res, by_class = by_class, overall = overall))
}
