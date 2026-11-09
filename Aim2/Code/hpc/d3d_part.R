# =============================================================================
# hpc/d3d_part.R
# Purpose : Compute one D3d part (cohort x variant) as a Slurm array task.
#           Needs only the code and objects/D3d_truth.rds (population
#           parameters; written by hpc/d3d_prepare.R), not the patient data.
# Usage   : Rscript Code/hpc/d3d_part.R --task k [--cores n] [--quick]
#           task k -> cohort ceiling(k / n_variants), variant in the order of
#           config$steps$D3d$variants (k = 1 .. n_rep * n_variants)
# Output  : Results/real/objects/D3d_parts/repNN_<variant>.rds, combined by
#           the D3d step of run_all.R once copied back.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))

suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- make_config("real", quick = "--quick" %in% args, cores = as.integer(get_arg("--cores", "1")),
                   root = normalizePath(file.path(code_dir, "..")))
cfg$paths$rcpp_cache <- file.path(tempdir(), "rcpp")   # node-local build, no clashes between tasks
load_marrow_model(cfg)

s <- cfg$steps$D3d
task <- as.integer(get_arg("--task"))
nv <- length(s$variants)
if (is.na(task) || task < 1 || task > s$n_rep * nv) stop("--task must be 1..", s$n_rep * nv)
r <- (task - 1L) %/% nv + 1L
variant <- s$variants[(task - 1L) %% nv + 1L]
out <- d3d_part_file(cfg, r, variant)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

log_fun <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
truth <- readRDS(d3d_truth_path(cfg))
log_fun("D3d part: cohort ", r, ", ", variant, ", cores ", cfg$runtime$cores, ", quick ", cfg$quick)
cal <- d3d_calibrate(cfg, truth)
log_fun("calibrated: mu_rho ", round(cal$mu_rho, 3), ", mu_g ", round(cal$mu_g, 3))
part <- d3d_part(cfg, truth, r, variant, cal, log_fun = log_fun)
saveRDS(part, out)
log_fun("written ", out, " (key ", part$key, ")")
