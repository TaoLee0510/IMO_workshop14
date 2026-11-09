# =============================================================================
# hpc/de_task.R
# Purpose : One step of a Model2 DE-ITS fit on the cluster: init | round0 | estep |
#           efinal | mstep | assemble (R/fit_de.R). The run directory holds the data
#           snapshot, the population of every round, the per-patient E-step
#           results and the status file read by hpc/de_next.sh. Sources the
#           analysis code without plotting packages (computing image).
# Usage   : Rscript Code/hpc/de_task.R --what init --run DIR --model Mab [--data real|FILE.rds]
#           Rscript Code/hpc/de_task.R --what estep --run DIR --iter k --patient i --cores n   (k = 0: round 0)
#           Rscript Code/hpc/de_task.R --what round0 --run DIR
#           Rscript Code/hpc/de_task.R --what mstep --run DIR --iter k --cores n
#           Rscript Code/hpc/de_task.R --what efinal --run DIR --iter k --patient i --cores n
#           Rscript Code/hpc/de_task.R --what assemble --run DIR
#           Rscript Code/hpc/de_task.R --what continue --run DIR   (run stopped at the cap)
#           Rscript Code/hpc/de_task.R --what reconfig --run DIR   (running run -> current outer settings)
#           init options --warm-pop FILE --warm-maps DIR: start from another run (de_init);
#           --pheno-w W / --pheno-s S: weight / smoothing scale of the phenotype term for this run
#           (S2: w 0; sensitivity analyses); --fold K: training set
#           of V1 fold K (the patients of fold K are left out; folds by v1_folds(), for real data
#           checked against the pre-registered Results/real/prereg/V1_folds.csv)
# Real data: init refuses to start unless the Model2 pre-registration verifies.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))

suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)

cores <- as.integer(get_arg("--cores", "1"))
cfg <- make_config("real", cores = cores, root = normalizePath(file.path(code_dir, "..")))
cfg$paths$rcpp_cache <- file.path(tempdir(), "rcpp")   # node-local build
what <- get_arg("--what")
if (!what %in% c("continue", "reconfig")) load_model2(cfg)
run <- get_arg("--run")
log_fun <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ..., "\n", sep = "")
iter <- as.integer(get_arg("--iter", "1"))
patient <- as.integer(get_arg("--patient", "1"))

if (what == "init") {
  src <- get_arg("--data", "real")
  if (src == "real") {
    v <- verify_prereg(cfg, "PR.sha256")
    if (!isTRUE(v$ok)) stop("Model2 pre-registration does not verify: ", v$msg)
    d <- read_model2_data(cfg)
  } else d <- readRDS(src)
  pw <- get_arg("--pheno-w")
  if (!is.null(pw) && nzchar(pw)) cfg$fit$pheno$w <- as.numeric(pw)
  ps <- get_arg("--pheno-s")
  if (!is.null(ps) && nzchar(ps)) cfg$fit$pheno$s <- as.numeric(ps)
  fk <- get_arg("--fold")
  if (!is.null(fk) && nzchar(fk)) {
    folds <- v1_folds(d$phenotype, cfg)
    if (src == "real") {
      pre <- data.table::fread(file.path(cfg$paths$results, "prereg", "V1_folds.csv"))
      if (!identical(pre[order(patient_id)]$fold, folds[order(patient_id)]$fold)) stop("V1 folds differ from the pre-registration")
    }
    keep <- folds[fold != as.integer(fk)]$patient_id
    d <- subset_data(d, keep)
  }
  wp <- get_arg("--warm-pop")
  warm <- if (is.null(wp) || !nzchar(wp)) NULL else list(pop = wp, maps = get_arg("--warm-maps"))
  rp <- de_init(run, get_arg("--model"), d, cfg, warm)
  writeLines(as.character(readRDS(rp$meta)$n), file.path(run, "n_patients"))
  log_fun("initialised ", run, " (", get_arg("--model"), ", ", readRDS(rp$meta)$n, " patients)")
} else if (what %in% c("estep", "efinal")) {
  f <- de_estep_task(run, iter, patient, cfg, cores, final = what == "efinal")
  log_fun(sprintf("%s %d patient %d (%s): nll %.3f, generations %d, converged %s, at bound [%s] (%.0f s)%s",
                  what, iter, patient, f$id, f$nll, f$gen, f$conv == 0L, paste(f$at_bound, collapse = ","),
                  f$runtime_s, if (!is.null(f$is)) sprintf("; IS log L %.2f (ESS %.0f), Laplace %.2f",
                                                          f$is[["logml"]], f$is[["ess"]], f$logL) else ""))
} else if (what == "round0") {
  st <- de_round0_task(run, cfg)
  log_fun("round 0 initialisation: ", st)
} else if (what == "mstep") {
  st <- de_mstep_task(run, iter, cfg, cores)
  if (file.exists(file.path(run, "trace.csv"))) {
    tr <- utils::read.csv(file.path(run, "trace.csv"))
    log_fun(sprintf("mstep %d: sum nll %.3f, max change %.4f (%s), rel obj %.2e -> %s", iter,
                    utils::tail(tr$sum_nll, 1), utils::tail(tr$max_change, 1), utils::tail(tr$which_max, 1),
                    utils::tail(tr$rel_obj, 1), st))
  } else log_fun("mstep ", iter, ": ", st)
} else if (what == "assemble") {
  fit <- de_assemble_task(run, cfg)
  if (is.list(fit)) log_fun(sprintf("assembled %s: logLik %.2f, BIC %.2f, n_par %d, rounds %d, converged %s, Laplace fallback %d",
                                    fit$model, fit$logLik, fit$BIC, fit$n_par, fit$rounds, fit$converged,
                                    sum(fit$indiv$logL_source == "Laplace")))
} else if (what == "continue") {
  k <- de_continue_task(run, cfg)
  log_fun("continued ", run, " from round ", k, " (outer settings: ",
          paste(names(cfg$fit$outer), unlist(cfg$fit$outer), sep = "=", collapse = ", "), ")")
} else if (what == "reconfig") {
  rc <- de_reconfig_task(run, cfg)
  log_fun("reconfigured ", run, " at round ", utils::tail(rc, 1)[[1]]$at_round, " (outer settings: ",
          paste(names(cfg$fit$outer), unlist(cfg$fit$outer), sep = "=", collapse = ", "), ")")
} else stop("unknown --what: ", what)
