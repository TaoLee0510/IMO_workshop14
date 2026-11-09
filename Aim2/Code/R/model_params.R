# =============================================================================
# model_params.R
# Purpose : Parameter bookkeeping for the marrow twin: C++ parameter layout,
#           random-effect transforms, and the nested-model switches M0/M1/M2.
# Inputs  : named parameter lists from config (config$model).
# Outputs : numeric vectors in the exact order expected by marrow_sim().
# =============================================================================

# Individual (random-effect) parameters and their transforms.
# eta = transform(natural); natural = inverse(eta)
re_transform <- list(
  rho = list(to_eta = log,    to_nat = exp),
  g   = list(to_eta = log,    to_nat = exp),
  f0  = list(to_eta = qlogis, to_nat = plogis)
)

# Which random effects each nested model estimates.
#   M0: inflammation is observed (CRP) but has no effect on the marrow
#   M1: inflammation suppresses the marrow (second wave possible)
#   M2: M1 + clonal layer (new low set point possible)
model_re_names <- function(model) {
  switch(model,
         M0 = c("rho", "g"),
         M1 = c("rho", "g"),
         M2 = c("rho", "g", "f0"),
         stop("unknown model: ", model))
}

# Structural parameters that are estimated at the population level per model.
model_struct_est <- function(model, config) {
  config$fit$struct_est[[model]]
}

# Apply the nesting switches to a full named parameter vector.
apply_model_switch <- function(par, model) {
  if (model %in% c("M0", "M1")) par["f0"] <- 0
  if (model == "M0") {
    par["q"] <- 0
    par["emax"] <- 0
  }
  par
}

# Assemble the full parameter vector in C++ order.
# indiv  : named numeric (natural scale) with any of rho, g, f0
# struct : named numeric with all structural parameters
build_par <- function(indiv, struct, model = "M2") {
  nm <- marrow_param_names()
  par <- setNames(numeric(length(nm)), nm)
  par[names(struct)] <- struct
  par[c("rho", "g", "f0")] <- c(1, 0, 0)
  par[names(indiv)] <- indiv
  par <- apply_model_switch(par, model)
  missing <- setdiff(nm, c(names(struct), "rho", "g", "f0"))
  if (length(missing)) stop("structural parameters missing: ", paste(missing, collapse = ", "))
  par[nm]
}

# eta (named, transformed) -> natural-scale individual parameters
eta_to_indiv <- function(eta) {
  vapply(names(eta), function(k) re_transform[[k]]$to_nat(eta[[k]]), numeric(1))
}

# natural-scale individual parameters -> eta
indiv_to_eta <- function(indiv) {
  vapply(names(indiv), function(k) re_transform[[k]]$to_eta(indiv[[k]]), numeric(1))
}

# Structural parameters live on log scale during estimation, except those
# bounded in (0, 1), which use logit.
struct_unit_interval <- c("r_ld", "r_inf", "delta", "emax", "e_plt",
                          "toci_red", "ster_red")

struct_to_work <- function(x) {
  out <- x
  for (k in names(x)) out[k] <- if (k %in% struct_unit_interval) qlogis(x[[k]]) else log(x[[k]])
  out
}

work_to_struct <- function(w) {
  out <- w
  for (k in names(w)) out[k] <- if (k %in% struct_unit_interval) plogis(w[[k]]) else exp(w[[k]])
  out
}
