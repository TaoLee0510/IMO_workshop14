# =============================================================================
# an_D2_phasemap.R
# Purpose : D2 - phenotype phase map. Simulates the natural history over a
#           grid of reserve x clone fraction, faceted by inflammatory gain and
#           clone output, and marks where a clone buffers or aggravates.
# Inputs  : ctx$config (design values; no patient data).
# Outputs : list(grid, tables, settings)
# =============================================================================

an_D2_phasemap <- function(ctx) {
  cfg <- ctx$config
  s <- cfg$steps$D2
  st0 <- design_struct(cfg)
  anc0 <- cfg$model$design$base$anc
  rho <- exp(seq(log(s$rho_range[1]), log(s$rho_range[2]), length.out = s$n_rho))
  f0 <- seq(s$f0_range[1], s$f0_range[2], length.out = s$n_f0)
  grid <- data.table::CJ(g_level = names(s$g_levels), delta_level = names(s$delta_levels),
                         rho = rho, f0 = f0)
  grid[, `:=`(g = s$g_levels[g_level], delta = s$delta_levels[delta_level])]

  sim_one <- function(rho, g, f0, delta) {
    st <- st0
    st["delta"] <- delta
    a0 <- baseline_with_clone(anc0, f0, st)
    x <- simulate_anc_nat(c(rho = rho, g = g, f0 = f0), st, "M2", a0, 1, cfg)
    late <- x$day >= 75
    c(phenotype = classify_phenotype(x$day, x$anc, cfg$phenotype),
      anc_late_rel = mean(x$anc[late]) / anc0, f_clone_90 = utils::tail(x$f_clone, 1))
  }
  out <- par_lapply(seq_len(nrow(grid)), function(i)
    sim_one(grid$rho[i], grid$g[i], grid$f0[i], grid$delta[i]), ctx$config$runtime$cores)
  out <- do.call(rbind, out)
  grid[, `:=`(phenotype = out[, "phenotype"], anc_late_rel = as.numeric(out[, "anc_late_rel"]),
              f_clone_90 = as.numeric(out[, "f_clone_90"]))]

  # buffer vs aggravate: severity relative to the same patient without a clone
  ref <- grid[f0 == min(f0), .(g_level, delta_level, rho, sev0 = phenotype_severity(phenotype))]
  grid <- merge(grid, ref, by = c("g_level", "delta_level", "rho"))
  grid[, effect := data.table::fcase(phenotype_severity(phenotype) < sev0, "buffer",
                                     phenotype_severity(phenotype) > sev0, "aggravate",
                                     default = "no change")]

  summ <- grid[, .N, by = .(g_level, delta_level, phenotype)][, pct := round(100 * N / sum(N), 1),
                                                             by = .(g_level, delta_level)]
  eff <- grid[f0 > min(f0), .N, by = .(g_level, delta_level, effect)][, pct := round(100 * N / sum(N), 1),
                                                                      by = .(g_level, delta_level)]
  list(grid = grid, tables = list(phenotype_share = summ, clone_effect = eff),
       settings = list(struct = st0, anc0 = anc0),
       key_numbers = list(clone_effect = eff))
}
