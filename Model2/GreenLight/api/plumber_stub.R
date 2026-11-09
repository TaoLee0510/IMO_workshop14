# GreenLight (MarrowTwin) · plumber stub of a remote module (docs/MODULE_API.md, docs/ML_API.md).
# GET /card returns the model card, POST /predict takes a patient snapshot as of a day.
# Replace the body of /predict with your model (e.g. a MarrowTwin variant or an R ML model).
# Run:  plumber::pr("api/plumber_stub.R") |> plumber::pr_run(port = 8789)
# Then add http://localhost:8789 on Researcher -> Interfaces.

#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type")
  if (req$REQUEST_METHOD == "OPTIONS") { res$status <- 200; return(list()) }
  plumber::forward()
}

#* Model card
#* @get /card
function() {
  list(id = "anc_trend_r", version = "0.1", type = "ml",
       outputs = list(list(id = "p_nonquick", labelKey = "out.p_nonquick", unit = "%", direction = "lower")))
}

#* Predict from a snapshot (only data up to snapshot$asOf)
#* @post /predict
function(req) {
  snap <- jsonlite::fromJSON(req$postBody, simplifyVector = FALSE)
  anc <- Filter(function(o) o[[1]] >= 0, snap$obs$anc)
  if (!length(anc)) return(list(p_nonquick = NULL))
  last <- max(anc[[length(anc)]][[2]], 0.01)
  eta <- -0.5 - 1.2 * log(last)
  list(p_nonquick = 1 / (1 + exp(-eta)), landmark = snap$asOf)
}
