# =============================================================================
# phenotype.R
# Purpose : Recovery phenotype from an ANC series. The same function serves
#           observed series (irregular days, censored values, G-CSF spikes)
#           and model trajectories (daily grid), so predicted probabilities and
#           observed labels are defined identically.
# Inputs  : day, anc, cens (logical), rule (config$phenotype$primary, possibly
#           modified by a variant).
# Outputs : class code ("quick" | "second_wave" | "delayed" [| "persistent"]);
#           phenotype_prob(): model probabilities of the three classes for the
#           phenotype term of the likelihood (Model1 plan section 5).
#
# Primary rule (day 0-30; thresholds of Rejeski 2021; locked 2026-11-03)
#   censored values count as cens_value (< 0.5); an isolated G-CSF spike
#   (previous < 1.0, value >= 1.5 and at least twice the previous, next < 1.0)
#   is not a recovery; early nadir = first minimum on days 0-14; recovery = the
#   first rec_n consecutive values >= rec_thr after the nadir; second wave =
#   dip_n consecutive values < dip_thr after the recovery (day >= dip_after);
#   delayed = no recovery by window_end; quick = the rest.
# Persistent (simulations only): median ANC on days 45-75 < 0.5.
# =============================================================================

phenotype_rule <- function(config, variant = NULL) {
  r <- config$phenotype$primary
  if (!is.null(variant)) r[names(config$phenotype$variants[[variant]])] <- config$phenotype$variants[[variant]]
  r
}

# indices of isolated G-CSF spikes in an ordered series; terminal = TRUE also
# flags a jump at the last value (no later value to confirm a recovery), as in
# the Aim 3 data review (P11 and P19 on day 28)
gcsf_spikes <- function(x, terminal = FALSE) {
  n <- length(x)
  if (n < 2) return(integer(0))
  i <- 2:n
  nxt <- c(x[-1], if (terminal) -Inf else Inf)[i]
  i[x[i - 1] < 1.0 & x[i] >= 1.5 & x[i] >= 2 * x[i - 1] & nxt < 1.0]
}

# first index i such that cond[i..i+k-1] are all TRUE (NA if none)
first_run <- function(cond, k, from = 1L) {
  n <- length(cond)
  if (from > n) return(NA_integer_)
  for (i in from:n) {
    if (i + k - 1 > n) return(NA_integer_)
    if (all(cond[i:(i + k - 1)])) return(i)
  }
  NA_integer_
}

# return_days: also return the days the rule used (the label days of the phenotype term)
phenotype_series <- function(day, anc, cens = NULL, rule, persistent = NULL, return_days = FALSE) {
  if (is.null(cens)) cens <- rep(FALSE, length(anc))
  o <- order(day)
  day <- day[o]; anc <- anc[o]; cens <- cens[o]
  ok <- !is.na(anc)
  day <- day[ok]; anc <- anc[ok]; cens <- cens[ok]
  late <- NULL
  if (!is.null(persistent)) {
    w <- day >= persistent$window[1] & day <= persistent$window[2]
    if (sum(w) >= 3 && stats::median(anc[w]) < persistent$thr) return("persistent")
  }
  keep <- day >= 0 & day <= rule$window_end
  day <- day[keep]; x <- anc[keep]; cens <- cens[keep]
  if (!length(x)) return(NA_character_)
  x[cens] <- rule$cens_value
  if (isTRUE(rule$drop_spikes)) {
    sp <- gcsf_spikes(x, terminal = TRUE)
    if (length(sp)) { day <- day[-sp]; x <- x[-sp] }
  }
  nw <- which(day >= rule$nadir_window[1] & day <= rule$nadir_window[2])
  cls <- if (!length(nw)) NA_character_ else {
    nad <- nw[which.min(x[nw])]
    rec <- first_run(x >= rule$rec_thr, rule$rec_n, nad + 1L)
    if (is.na(rec)) "delayed" else {
      dip <- first_run(x < rule$dip_thr & day >= rule$dip_after, rule$dip_n, rec + rule$rec_n)
      if (!is.na(dip)) "second_wave" else "quick"
    }
  }
  if (return_days) list(class = cls, days = day) else cls
}

# ---- phenotype term (Model1) ---------------------------------------------------------------
# Probabilities of quick / second wave / delayed for a model ANC series mu (absolute counts) on
# the label days, by a margin-based soft version of the rule:
#   j*  = day of the minimum on the nadir window
#   m_R = max over runs after j* of log(min of rec_n consecutive values) - log(rec_thr)
#   r*  = first run with a positive recovery margin, else the run with the largest margin
#   m_D = max over runs from r* + rec_n (day >= dip_after) of log(dip_thr) - log(max of dip_n values)
#   P(delayed) = 1 - Phi(m_R / s), P(second wave) = Phi(m_R / s) Phi(m_D / s),
#   P(quick) = Phi(m_R / s) (1 - Phi(m_D / s));  P' = (1 - 3 eps) P + eps
# s -> 0 gives the rule applied to the model series. Returns NULL without a nadir-window day.
phenotype_prob <- function(day, mu, rule, s, eps = 0) {
  o <- order(day)
  day <- day[o]
  lx <- log(pmax(mu[o], 1e-8))
  n <- length(lx)
  nw <- which(day >= rule$nadir_window[1] & day <= rule$nadir_window[2])
  if (!length(nw)) return(NULL)
  nad <- nw[which.min(lx[nw])]
  run_stat <- function(from, k, f) {
    if (from > n - k + 1) return(numeric(0))
    vapply(from:(n - k + 1), function(i) f(lx[i:(i + k - 1)]), numeric(1))
  }
  rr <- run_stat(nad + 1L, rule$rec_n, min) - log(rule$rec_thr)
  if (!length(rr)) {
    p <- c(quick = 0, second_wave = 0, delayed = 1)
  } else {
    m_r <- max(rr)
    r_star <- nad + if (any(rr >= 0)) which(rr >= 0)[1] else which.max(rr)
    from <- r_star + rule$rec_n
    dd <- run_stat(from, rule$dip_n, max)
    if (length(dd)) dd <- (log(rule$dip_thr) - dd)[day[from:(n - rule$dip_n + 1)] >= rule$dip_after]
    p_r <- stats::pnorm(m_r / s)
    p_d <- if (length(dd)) stats::pnorm(max(dd) / s) else 0
    p <- c(quick = p_r * (1 - p_d), second_wave = p_r * p_d, delayed = 1 - p_r)
  }
  (1 - 3 * eps) * p + eps
}

# Phenotype of every patient in an observation table (var == "anc").
# obs: data.table(patient_id, day, var, value, cens, ...)
phenotype_observed <- function(obs, config, variant = NULL) {
  rule <- phenotype_rule(config, variant)
  a <- obs[var == "anc" & !is.na(value)]
  a[, .(phenotype = phenotype_series(day, value, cens, rule)), by = patient_id]
}

# Phenotype table under the primary rule and every variant.
phenotype_all_rules <- function(obs, config) {
  out <- phenotype_observed(obs, config)
  data.table::setnames(out, "phenotype", "primary")
  for (v in names(config$phenotype$variants)) {
    pv <- phenotype_observed(obs, config, v)
    out[pv, (v) := i.phenotype, on = "patient_id"]
  }
  out
}

# Phenotype of a simulated daily trajectory (model latent ANC, no G-CSF spikes).
# with_persistent: also assign the late class (needs follow-up to day 75).
phenotype_sim <- function(day, anc, config, variant = NULL, with_persistent = FALSE) {
  phenotype_series(day, anc, NULL, phenotype_rule(config, variant),
                   if (with_persistent) config$phenotype$persistent else NULL)
}

pheno_rank <- function(x) match(x, c("quick", "second_wave", "delayed", "persistent"))

phenotype_display <- function(x) {
  lab <- c(quick = "Quick recovery", second_wave = "Second wave", delayed = "Delayed recovery",
           persistent = "Persistent cytopenia")
  if (missing(x)) lab else unname(lab[x])
}
