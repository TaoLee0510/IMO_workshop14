# =============================================================================
# model_simulate.R
# Purpose : R wrappers around the C++ integrator: compile/load the model and
#           map normalised model output to observable quantities.
# Inputs  : parameter vectors (model_params.R), segments (model_inputs.R).
# Outputs : matrices / data.tables of simulated trajectories.
# =============================================================================

# Compile (or load the cached build of) the C++ model once per session.
load_marrow_model <- function(config) {
  if (exists("marrow_sim", mode = "function")) return(invisible(TRUE))
  dir.create(config$paths$rcpp_cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(file.path(config$paths$code, "src", "marrow_model.cpp"),
                  cacheDir = config$paths$rcpp_cache, env = globalenv())
  invisible(TRUE)
}

# Simulate the normalised model at the requested times.
sim_raw <- function(par, seg, obs_t, solver) {
  marrow_sim(par, seg, obs_t, rtol = solver$rtol, atol = solver$atol, hmax = solver$hmax)
}

# Convert normalised output to observable scales.
# base : list(anc, plt, crp, ferritin) - patient's pre-treatment values
# obsm : config$model$obs (crp_scale, fer_scale)
sim_observables <- function(raw, base, obsm) {
  data.table::data.table(
    day      = raw[, "t"],
    anc      = raw[, "anc_rel"] * base$anc,
    plt      = raw[, "plt_rel"] * base$plt,
    crp      = base$crp + obsm$crp_scale * raw[, "is"],
    ferritin = base$ferritin + obsm$fer_scale * raw[, "is"],
    is       = raw[, "is"],
    im       = raw[, "im"],
    f_clone  = raw[, "f_clone"],
    x_tot    = raw[, "x_tot"],
    f_rel    = raw[, "f_rel"]
  )
}

# Convenience: natural-history daily trajectory for one parameter set.
sim_daily <- function(par, config, base = NULL, ld_scale = 1, events = NULL,
                      days = config$model$time$sim_days) {
  seg <- make_segments(config$model$time$t0, max(days), ld_scale = ld_scale, events = events)
  raw <- sim_raw(par, seg, days, config$model$solver)
  if (is.null(base)) base <- list(anc = 1, plt = 1, crp = 0, ferritin = 0)
  sim_observables(raw, base, config$model$obs)
}
