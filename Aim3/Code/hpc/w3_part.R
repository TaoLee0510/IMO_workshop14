# =============================================================================
# hpc/w3_part.R
# Purpose : Compute one W3 task (scenario x replicate x cohort size) as a Slurm
#           array task. Sources the analysis code without the plotting
#           packages, so it runs in a container that only has the computing
#           packages (data.table, Rcpp).
# Usage   : Rscript Code/hpc/w3_part.R --task k [--cores n] [--quick]
#           k indexes w3_tasks(config) (1 .. 88 with the full settings)
# Needs   : Results/real/objects/W0.rds (converted data, copied from local)
# Output  : Results/real/objects/W3_parts/W3_task_kkk.rds with the code and
#           settings fingerprint; the local W3 step merges all parts.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))

suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cfg <- make_config(quick = "--quick" %in% args, cores = as.integer(get_arg("--cores", "1")),
                   root = normalizePath(file.path(code_dir, "..")))
cfg$paths$rcpp_cache <- file.path(tempdir(), "rcpp")   # node-local build, no clashes between tasks
load_hsc_model(cfg)

k <- as.integer(get_arg("--task"))
tasks <- w3_tasks(cfg)
if (is.na(k) || !k %in% tasks$task_id) stop("--task must be 1..", nrow(tasks))
d <- readRDS(file.path(cfg$paths$results, "objects", "W0.rds"))$derived
out <- file.path(cfg$paths$results, "objects", "W3_parts", sprintf("W3_task_%03d.rds", k))
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)

log_fun <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
task <- tasks[task_id == k]
log_fun("W3 task ", k, ": ", task$scenario, ", rep ", task$rep, ", n ", task$n, ", cores ", cfg$runtime$cores)
t0 <- Sys.time()
part <- w3_run_task(task, d, cfg, cores = cfg$runtime$cores, log_fun = log_fun)
part$fingerprint <- w3_fingerprint(cfg)
part$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
saveRDS(part, out)
log_fun("written ", out, " (", round(part$runtime_min, 1), " min)")
