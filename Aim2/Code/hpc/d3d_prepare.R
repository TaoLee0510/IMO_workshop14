# =============================================================================
# hpc/d3d_prepare.R
# Purpose : Write the D3d truth (population parameters of the real-data M1
#           fit, no patient-level data) to Results/real/objects/D3d_truth.rds,
#           so the array tasks on the cluster can run without the data.
# Usage   : Rscript Code/hpc/d3d_prepare.R   (after D4 on the real data)
# =============================================================================

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))
suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- make_config("real", root = normalizePath(file.path(code_dir, "..")))
d4 <- load_step_result("D4", cfg)
if (is.null(d4)) stop("D4 result for the real data not found - run D4 first")
truth <- d3d_truth(cfg, d4)
saveRDS(truth, d3d_truth_path(cfg))
cat("written", d3d_truth_path(cfg), "\n")
str(truth)
