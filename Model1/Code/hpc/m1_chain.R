# =============================================================================
# hpc/m1_chain.R
# Purpose : Decide which Model1 runs come next once runs have finished (called
#           after every assemble by hpc/m1_chain.sh; idempotent through marker
#           files in RUNS/_chain). Writes RUNS/_chain/todo.txt with one launch per
#           line: RUN MODEL DATA WARM_POP WARM_MAPS PHENO_W PHENO_S FOLD ("-" = none).
#   1. wave 1 assembled -> Q2 selection (pre-registered rule) -> wave 2 for q: qC_ddr qC_dta q_Zh
#      q_Zm qC_ddr_Zh qC_ddr_Zm (warm from q). Deviation 2026-11-04 17:55 (approved by the user):
#      wave 2 starts once Ma, Mb, Mab and Mabx are assembled, without waiting for the slow M0;
#   1b. when M0 is assembled the Q2 selection is repeated with all five: unchanged -> logged;
#      changed -> wave 2 for the new q and the final selection is repeated after it
#   2. wave 2 all assembled -> final model (Q1 and route rules) -> V1: five fold refits
#      (warm from the final run) + phenotype-term sensitivity (s 0.15 / 0.6, w 0 / 5)
#   3. S1 cohort synNN assembled -> S2 refit with w = 0 (synNNw0, warm) and, for the
#      first validation$S3$n_cohorts cohorts, S3 fold refits (synNN_foldK, warm)
# Usage   : Rscript Code/hpc/m1_chain.R --runs Results/real/hpc_runs
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))
suppressPackageStartupMessages(library(data.table))
source(file.path(code_dir, "config.R"))
for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- make_config("real", root = normalizePath(file.path(code_dir, "..")))

runs <- get_arg("--runs")
chain <- file.path(runs, "_chain")
dir.create(chain, showWarnings = FALSE)
done <- function(r) file.exists(file.path(runs, r, "fit.rds"))
mark <- function(m) file.exists(file.path(chain, m))
todo <- character(0)
launch <- function(run, model, data = "real", warm = NULL, w = "-", s = "-", fold = "-") {
  if (dir.exists(file.path(runs, run))) return(invisible(NULL))
  wp <- wm <- "-"
  if (!is.null(warm)) {
    pf <- sort(list.files(file.path(runs, warm), pattern = "^pop_[0-9]+\\.rds$", full.names = TRUE))
    wp <- utils::tail(pf, 1)
    wm <- file.path(runs, warm, "final")
  }
  todo <<- c(todo, paste(file.path(runs, run), model, data, wp, wm, w, s, fold))
}
logm <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", ..., "\n", sep = "",
                          file = file.path(chain, "chain.log"), append = TRUE)

fits_of <- function(ids) {
  ids <- ids[vapply(ids, done, logical(1))]
  lapply(stats::setNames(ids, ids), function(m) readRDS(file.path(runs, m, "fit.rds")))
}

# 1. wave 2
q2_ids <- cfg$fit$q2_models
wave2_of <- function(q) {
  w2 <- c(model_id(q, "ddr"), model_id(q, "dta"), model_id(q, "none", "hsc"), model_id(q, "none", "mpp"),
          model_id(q, "ddr", "hsc"), model_id(q, "ddr", "mpp"))
  # q_Zh is not defined for models without b (the HSC route belongs to b): keep the defined ones
  w2[!(model_spec(q)$q2 %in% c("M0", "Ma") & grepl("_Zh$", w2))]
}
early_ids <- setdiff(q2_ids, "M0")
if (all(vapply(early_ids, done, logical(1))) && !mark("wave2")) {
  sel <- f1_select13(fits_of(q2_ids), cfg)
  q <- sel$q2
  w2 <- wave2_of(q)
  for (m in w2) launch(m, m, warm = q)
  writeLines(c(q, w2), file.path(chain, "wave2"))
  if (!done("M0")) invisible(file.create(file.path(chain, "wave2_without_M0")))
  logm("wave 2 for ", q, " (runner-up ", sel$q2_runner_up, ", margin ", round(sel$q2_margin, 2), "", 
       if (!done("M0")) "; M0 not yet assembled - provisional" else "", "): ", paste(w2, collapse = " "))
}

# 1b. M0 check of a wave 2 started without it
if (mark("wave2_without_M0") && done("M0") && !mark("m0_check")) {
  sel <- f1_select13(fits_of(q2_ids), cfg)
  q_old <- readLines(file.path(chain, "wave2"))[1]
  if (identical(sel$q2, q_old)) {
    logm("M0 check: Q2 selection unchanged (", q_old, "; M0 BIC ", round(fits_of("M0")$M0$BIC, 1), ")")
  } else {
    w2 <- wave2_of(sel$q2)
    for (m in w2) launch(m, m, warm = sel$q2)
    writeLines(c(sel$q2, w2), file.path(chain, "wave2"))
    unlink(file.path(chain, "final"))
    logm("M0 check: Q2 selection changed from ", q_old, " to ", sel$q2, " - wave 2 for it launched; final selection repeated after it")
  }
  invisible(file.create(file.path(chain, "m0_check")))
}

# 2. V1 folds and phenotype-term sensitivity of the final model
if (mark("wave2") && !mark("final")) {
  w2 <- readLines(file.path(chain, "wave2"))
  if (all(vapply(w2[-1], done, logical(1)))) {
    sel <- f1_select13(fits_of(unique(c(q2_ids, w2))), cfg)
    fin <- sel$final
    for (k in seq_len(cfg$validation$V1$k)) launch(paste0(fin, "_fold", k), fin, warm = fin, fold = k)
    ph <- cfg$fit$pheno
    for (s in ph$sens_s) launch(paste0(fin, "_s", s), fin, warm = fin, s = s)
    for (w in ph$sens_w) launch(paste0(fin, "_w", w), fin, warm = fin, w = w)
    writeLines(fin, file.path(chain, "final"))
    logm("final model ", fin, ": V1 folds and phenotype-term sensitivity launched")
  }
}

# 3. S2 and S3 after each S1 cohort
cdir <- file.path(dirname(runs), "recovery", "cohorts")
for (cn in sub("\\.rds$", "", list.files(cdir, pattern = "^syn[0-9]+\\.rds$"))) {
  if (!done(cn) || mark(paste0("S_", cn))) next
  co <- file.path(cdir, paste0(cn, ".rds"))
  model <- readRDS(file.path(runs, cn, "fit.rds"))$model
  launch(paste0(cn, "w0"), model, data = co, warm = cn, w = cfg$validation$S2$w)
  if (as.integer(sub("syn", "", cn)) <= cfg$validation$S3$n_cohorts)
    for (k in seq_len(cfg$validation$V1$k)) launch(paste0(cn, "_fold", k), model, data = co, warm = cn, fold = k)
  file.create(file.path(chain, paste0("S_", cn)))
  logm("S2", if (as.integer(sub("syn", "", cn)) <= cfg$validation$S3$n_cohorts) " + S3", " launched for ", cn)
}

if (length(todo)) cat(todo, sep = "\n", file = file.path(chain, "todo.txt"), append = TRUE)
