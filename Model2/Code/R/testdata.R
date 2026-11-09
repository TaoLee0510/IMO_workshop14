# =============================================================================
# testdata.R
# Purpose : Write the synthetic TestData cohort: the 35 real patients' design
#           (sampling days, ferritin, eGFR, bridging, set points) with counts
#           simulated from the parent model (mixed a/b mechanisms, 35% DDR-type
#           clone carriers). Lets the whole pipeline run, and be validated
#           against the truth, before any real-data fit.
# Inputs  : config (test), the converted real extract (design only).
# Outputs : TestData/{patients, obs, inflam, truth}.csv
# =============================================================================

make_testdata <- function(config, seed = 20261103L) {
  creal <- make_config("real", root = config$paths$root)
  d <- read_model2_data(creal)
  dz <- design_real35(d)
  ids <- sprintf("S%02d", seq_along(dz))
  co <- simulate_cohort(dz, "mix", config$model$design$p_clone, design_population(config), config,
                        seed, ids = ids)
  dir.create(config$paths$testdata, recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(co$patients, file.path(config$paths$testdata, "patients.csv"))
  data.table::fwrite(co$obs, file.path(config$paths$testdata, "obs.csv"))
  data.table::fwrite(co$inflam, file.path(config$paths$testdata, "inflam.csv"))
  data.table::fwrite(co$truth, file.path(config$paths$testdata, "truth.csv"))
  writeLines(c("# TestData (synthetic)",
               "",
               "Design of the 35 real patients (sampling days, ferritin, eGFR, bridging, set points);",
               "counts simulated from the Aim 4 parent model, scenario `mix` (a 40% / b 30% / a+b 30%),",
               "35% DDR-type clone carriers. `truth.csv` holds the true parameters and categories.",
               "Regenerate with `Rscript Code/run_all.R --make-testdata`."),
             file.path(config$paths$testdata, "README.md"))
  invisible(co)
}
