# =============================================================================
# hpc/d3c_part.R
# Purpose : Compute one D3c part (cohort x variant) as a Slurm array task.
#           Sources the analysis code without the plotting packages, so it
#           runs in a container that has only the computing packages
#           (data.table, digest, MASS, Rcpp).
# Usage   : Rscript Code/hpc/d3c_part.R --task k [--cores n] [--quick] [--results dir]
#           task k -> cohort ceiling(k / n_variants), variant in the order of
#           config$steps$D3c$variants (k = 1 .. n_rep * n_variants)
#           --results overrides Results/test (e.g. for a test run)
# Output  : <results>/objects/D3c_parts/repNN_<variant>.rds, which the D3c step
#           of run_all.R combines once all parts are copied back.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))

suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- make_config("test", quick = "--quick" %in% args, cores = as.integer(get_arg("--cores", "1")),
                   root = normalizePath(file.path(code_dir, "..")))
if (!is.null(get_arg("--results"))) cfg$paths$results <- get_arg("--results")
cfg$paths$rcpp_cache <- file.path(tempdir(), "rcpp")   # node-local build, no clashes between tasks
load_marrow_model(cfg)

st <- d3c_setup(cfg)
task <- as.integer(get_arg("--task"))
nv <- length(st$s$variants)
if (is.na(task) || task < 1 || task > st$s$n_rep * nv) stop("--task must be 1..", st$s$n_rep * nv)
r <- (task - 1L) %/% nv + 1L
variant <- st$s$variants[(task - 1L) %% nv + 1L]
out <- d3c_part_file(cfg, r, variant)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

log_fun <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
log_fun("D3c part: cohort ", r, ", ", variant, ", cores ", cfg$runtime$cores, ", quick ", cfg$quick)
cal <- d3c_calibrate(cfg, st)
log_fun("calibrated")
part <- d3c_part(cfg, r, variant, cal, st, log_fun = log_fun)
saveRDS(part, out)
log_fun("written ", out, " (key ", part$key, ")")
