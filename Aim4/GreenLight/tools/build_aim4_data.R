# =============================================================================
# tools/build_aim4_data.R
# Purpose : Build GreenLight's Aim 4 data files without a model fit:
#           data/aim4_pop.js (prior centres from Aim4/Code/config.R, source
#           "prior"; kept if a fitted file already exists unless --force-pop),
#           data/real_cohort.js (35 real patients, de-identified, local use only)
#           and tools/ref_trajectories.json (C++ reference for the JS port).
#           Step G of Aim4/Code/run_all.R writes the same files after fitting.
# Usage   : (from Aim4/)  Rscript GreenLight/tools/build_aim4_data.R [--force-pop]
# =============================================================================

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".")
if (!file.exists(file.path(root, "Code", "config.R"))) stop("run from Aim4/")
source(file.path(root, "Code", "config.R"))
for (f in list.files(file.path(root, "Code", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- make_config("real", root = root)
load_aim4_model(cfg)
ensure_real_derived(cfg)
d <- read_aim4_data(cfg)
gl <- cfg$paths$greenlight
pop_file <- file.path(gl, "data", "aim4_pop.js")
keep_fitted <- file.exists(pop_file) && any(grepl('"source":"fitted"', readLines(pop_file, warn = FALSE), fixed = TRUE)) &&
  !"--force-pop" %in% args
pop <- g_prior_pop(d, cfg)
tab <- g_export_all(gl, if (keep_fitted) NULL else pop, g_real_cohort(d, cfg), g_ref_trajectories(d, cfg, pop))
print(tab)
if (keep_fitted) message("kept the fitted data/aim4_pop.js (use --force-pop to overwrite with priors)")
