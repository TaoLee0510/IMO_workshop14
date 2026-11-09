# =============================================================================
# sim_clinical.R
# Purpose : Clinical process applied to a simulated patient: CRS grade from the
#           systemic inflammation, tocilizumab / steroids / anakinra by grade,
#           then up to two G-CSF courses when ANC < 0.5 after day 15; for
#           persistently low patients TPO-RA and a CD34+ boost. Shared by
#           the TestData generator and the virtual cohorts of D3/D3b, so that
#           design studies see the same supportive care as the test data.
# Inputs  : individual parameters (+ anc_base), structural parameters, config,
#           `spec` = list(events = probabilities, crs_cut = Is thresholds for
#           CRS grades 1-4) - truth_spec() in the generator or
#           config$model$design$clinical in D3/D3b.
# Outputs : list(m = daily simulation matrix, events = data.table, crs = list)
# =============================================================================

crs_from_is <- function(is, days, cut) {
  pk <- max(is)
  grade <- sum(pk >= cut)
  on <- which(is > 0.2 & days >= 0)
  if (grade == 0 || !length(on)) return(list(grade = 0L, onset = NA_integer_, end = NA_integer_))
  list(grade = as.integer(grade), onset = as.integer(days[min(on)]),
       end = as.integer(days[max(on[days[on] <= 30])]))
}

simulate_truth_patient <- function(indiv, struct, cfg, spec, ld_scale) {
  # indiv: rho, g, f0 and anc_base (true pre-treatment ANC, used for G-CSF decisions)
  days <- cfg$model$time$sim_days
  t0 <- cfg$model$time$t0
  par <- build_par(indiv, struct, "M2")
  sim <- function(ev) {
    seg <- make_segments(t0, max(days), ld_scale = ld_scale, events = ev)
    marrow_sim(par, seg, days)
  }
  ev <- data.table::data.table(event = character(0), start_day = numeric(0),
                               end_day = numeric(0), dose = numeric(0))
  # pass 1: natural history -> CRS grade -> tocilizumab / steroids / anakinra
  m <- sim(ev)
  crs <- crs_from_is(m[, "is"], days, spec$crs_cut)
  e <- spec$events
  if (crs$grade >= 1) {
    p_toci <- if (crs$grade >= 2) e$p_toci_g2 else e$p_toci_g1
    if (stats::runif(1) < p_toci)
      ev <- rbind(ev, list("tocilizumab", crs$onset + 1, crs$onset + 1, 8))
    p_ster <- if (crs$grade >= 3) e$p_steroid_g3 else if (crs$grade == 2) e$p_steroid_g2 else 0
    if (stats::runif(1) < p_ster)
      ev <- rbind(ev, list("steroid", crs$onset + 2, crs$onset + 5, 10))
    if (crs$grade >= 3 && stats::runif(1) < e$p_anakinra_g3)
      ev <- rbind(ev, list("anakinra", crs$onset + 2, crs$onset + 6, 100))
  }
  # pass 2: with CRS treatment -> G-CSF courses (up to 2)
  m <- sim(ev)
  anc <- m[, "anc_rel"] * indiv[["anc_base"]]
  last_end <- -Inf
  for (course in 1:2) {
    cand <- which(days >= 15 & days > last_end + 7 & anc < 0.5)
    if (!length(cand)) break
    p <- if (course == 1) e$p_gcsf else e$p_gcsf_second
    if (stats::runif(1) >= p) break
    d <- days[cand[1]]
    ev <- rbind(ev, list("gcsf", d, d + 4, 5))
    last_end <- d + 4
    m <- sim(ev)
    anc <- m[, "anc_rel"] * indiv[["anc_base"]]
  }
  list(m = m, events = ev, crs = crs)
}

# Support given to persistently low patients: TPO receptor agonist from day 50
# and a CD34+ boost between day 55 and 75 (data after the boost are censored).
# `e` holds p_tpo_persistent and p_cd34_persistent.
add_persistent_support <- function(ev, e) {
  if (stats::runif(1) < e$p_tpo_persistent) ev <- rbind(ev, list("tpo_ra", 50, 90, 50))
  if (stats::runif(1) < e$p_cd34_persistent) {
    d34 <- sample(55:75, 1)
    ev <- rbind(ev, list("cd34_boost", d34, d34, 2))
  }
  ev
}
