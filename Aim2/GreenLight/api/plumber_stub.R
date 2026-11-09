# GreenLight · R plumber stub implementing the module contract (docs/MODULE_API.md).
# Replace the bodies of `fit` and `forecast` with calls into the Aim 2 model (Code/R, src/marrow_model.cpp).
# Run:  Rscript -e 'plumber::pr_run(plumber::pr("api/plumber_stub.R"), port = 8788)'
# Then register http://localhost:8788/module on the Integration page.

#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Authorization")
  if (req$REQUEST_METHOD == "OPTIONS") { res$status <- 200; return(list()) }
  plumber::forward()
}

#* @get /module/meta
function() {
  list(
    meta = list(id = "marrow_r", nameKey = "Marrow module (Aim 2, R)", owner = "Aim 2 · Team Green", version = "0.1", kind = "remote"),
    consumes = list(variables = c("anc", "plt", "crp"), planKeys = c("ld.flu_frac", "ld.cy_frac", "cytopenia.*"), signals = c("car_expansion")),
    provides = list(
      signals = c("inflammation"), variables = c("anc", "plt", "crp"),
      params = list(list(id = "R", labelKey = "param.R", lo = 0, hi = 1.2), list(id = "g", labelKey = "param.g", log = TRUE, lo = 2, hi = 600), list(id = "f0", labelKey = "param.f0", lo = 0, hi = 0.9)),
      outcomes = list(list(id = "p_quick", labelKey = "out.p_quick", kind = "prob", direction = "higher"),
                      list(id = "p_second_wave", labelKey = "out.p_second_wave", kind = "prob", direction = "lower"),
                      list(id = "p_persistent", labelKey = "out.p_persistent", kind = "prob", direction = "lower"))
    )
  )
}

#* @post /module/fit
#* @serializer unboxedJSON
function(req) {
  body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
  snapshot <- body$snapshot; ctx <- body$ctx
  obs <- snapshot$observations
  # TODO: call the Aim 2 population fit / individual posterior here.
  list(moduleId = "marrow_r", asOfDay = ctx$asOfDay, nObs = length(obs), ess = 100,
       params = list(R = list(median = 0.6, lo = 0.4, hi = 0.8), g = list(median = 40, lo = 20, hi = 80), f0 = list(median = 0.02, lo = 0, hi = 0.1)),
       opaque = list(seed = ctx$seed))
}

#* @post /module/forecast
#* @serializer unboxedJSON
function(req) {
  body <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
  state <- body$state; plan <- body$plan; ctx <- body$ctx
  grid <- seq(-6, if (is.null(ctx$horizon)) 120 else ctx$horizon)
  # TODO: simulate the Aim 2 ODE under `plan` for posterior draws in `state`; return predictive quantiles.
  base <- 3 * (1 - 0.8 * exp(-pmax(grid, 0) / 10))
  q <- function(f) as.list(base * f)
  list(grid = as.list(grid),
       variables = list(anc = list(q05 = q(0.5), q25 = q(0.8), q50 = q(1), q75 = q(1.2), q95 = q(1.6)),
                        plt = list(q05 = q(40), q25 = q(55), q50 = q(70), q75 = q(85), q95 = q(110)),
                        crp = list(q05 = q(2), q25 = q(4), q50 = q(6), q75 = q(9), q95 = q(14))),
       outcomes = list(p_quick = list(value = 0.6), p_second_wave = list(value = 0.3), p_persistent = list(value = 0.1)),
       signals = list(inflammation = list(grid = as.list(grid), values = q(6))))
}
