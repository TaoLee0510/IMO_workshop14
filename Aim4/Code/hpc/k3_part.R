# =============================================================================
# hpc/k3_part.R
# Purpose : Compute one K3 part (synthetic cohort x refit x attribution) as a
#           Slurm array task. Sources the analysis code without plotting
#           packages, so it runs in a container that has only the computing
#           packages (data.table, digest, jsonlite, MASS, Rcpp).
# Usage   : Rscript Code/hpc/k3_part.R --task k [--cores n] [--quick]
# Output  : Results/real/objects/K3_parts/taskNNN.rds; the K3 step of
#           run_all.R combines the parts once they are copied back.
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
cfg$paths$rcpp_cache <- file.path(tempdir(), "rcpp")   # node-local build
load_aim4_model(cfg)
d <- read_aim4_data(cfg)                                # design only: Data/derived (synced)

tasks <- k3_tasks(cfg)
k <- as.integer(get_arg("--task"))
if (is.na(k) || !k %in% tasks$task) stop("--task must be 1..", nrow(tasks))
out <- k3_part_file(cfg, k)
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
log_fun <- function(...) cat(format(Sys.time(), "%H:%M:%S"), " ", ..., "\n", sep = "")
part <- k3_part(cfg, tasks[task == k], d, log_fun = log_fun)
saveRDS(part, out)
log_fun("written ", out)
