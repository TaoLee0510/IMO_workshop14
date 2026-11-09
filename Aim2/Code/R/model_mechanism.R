# =============================================================================
# model_mechanism.R
# Purpose : Shared helpers that turn a parameter set into a predicted
#           phenotype and a dominant mechanism, via counterfactual simulation.
# Inputs  : natural-scale individual parameters, structural parameters.
# Outputs : phenotype label and mechanism label (none / H1 / H2 / H3).
#
# Mechanism rule
#   quick        -> "none"
#   second_wave  -> "H2" (inflammation-driven dip)
#   persistent   -> "H3" if removing the clone (f0 = 0) lifts the patient out
#                   of the persistent class, otherwise "H1" (low reserve)
# =============================================================================

mechanism_levels <- c("none", "H1", "H2", "H3")

simulate_anc_nat <- function(indiv, struct, model, anc_base, ld_scale, config, events = NULL) {
  days <- config$model$time$sim_days
  par <- build_par(indiv, struct, model)
  seg <- make_segments(config$model$time$t0, max(days), ld_scale = ld_scale, events = events)
  m <- marrow_sim(par, seg, days, config$model$solver$rtol, config$model$solver$atol,
                  config$model$solver$hmax)
  list(day = days, anc = m[, "anc_rel"] * anc_base, plt_rel = m[, "plt_rel"],
       f_clone = m[, "f_clone"], is = m[, "is"])
}

infer_mechanism <- function(indiv, struct, model, anc_base, ld_scale, config, events = NULL) {
  s <- simulate_anc_nat(indiv, struct, model, anc_base, ld_scale, config, events)
  ph <- classify_phenotype(s$day, s$anc, config$phenotype)
  mech <- if (is.na(ph)) NA_character_
          else if (ph == "quick") "none"
          else if (ph == "second_wave") "H2"
          else {
            if (model == "M2" && "f0" %in% names(indiv) && indiv[["f0"]] > 1e-4) {
              cf <- indiv
              cf["f0"] <- 0
              s0 <- simulate_anc_nat(cf, struct, model, anc_base, ld_scale, config, events)
              if (identical(classify_phenotype(s0$day, s0$anc, config$phenotype), "persistent")) "H1" else "H3"
            } else "H1"
          }
  list(phenotype = ph, mechanism = mech)
}

# Posterior probability of each mechanism for the patients `ids`: eta is drawn
# from each patient's Laplace posterior N(eta_hat, cov) and every draw is
# attributed as in infer_mechanism(). Summed over patients, p_h3 gives the
# expected number of clonal cases.
# Note: when the clone fraction is large the counts carry little information on
# reserve, and the Laplace posterior is too narrow there, so these
# probabilities are over-confident and clonal cases tend to be missed (TestData
# diagnostic, BUILD_STATUS.md); the expected count is a lower bound.
mechanism_probability <- function(fit, d, ids, config, n_draw = config$fit$mech_draws,
                                  cores = config$runtime$cores) {
  st <- fit$pop$struct
  out <- par_lapply(ids, function(id) {
    i <- match(id, fit$indiv$patient_id)
    set.seed(config$runtime$seed + i)
    mu <- unlist(fit$indiv[i, paste0("eta_", fit$re), with = FALSE])
    names(mu) <- fit$re
    draws <- MASS::mvrnorm(n_draw, mu, fit$cov[[id]])
    prow <- d$patients[patient_id == id]
    mech <- vapply(seq_len(n_draw), function(k) {
      m <- infer_mechanism(eta_named_to_indiv(stats::setNames(draws[k, ], fit$re)), st, fit$model,
                           prow$anc_base, prow$ld_scale, config)$mechanism
      if (is.na(m)) "unclassified" else m
    }, character(1))
    data.table::data.table(patient_id = id, p_h1 = mean(mech == "H1"), p_h2 = mean(mech == "H2"),
                           p_h3 = mean(mech == "H3"), p_none = mean(mech == "none"))
  }, cores)
  data.table::rbindlist(out)
}

# Pre-treatment ANC of a patient carrying a clone: a low-output clone lowers the
# baseline (same relation as in the test-data generator).
baseline_with_clone <- function(anc, f0, struct) {
  anc * ((1 - f0) + struct[["delta"]] * f0)^(1 / (1 + struct[["gamma_n"]]))
}

# Design structural parameters (pre-data steps D1-D3), see config$model$design.
design_struct <- function(config) {
  st <- config$model$struct
  st[names(config$model$design$struct)] <- config$model$design$struct
  st
}
