# =============================================================================
# run_all.R  -  single entry point for Aim 3
# Purpose : Orchestrate the four layers (read -> analyse -> write -> plot) for
#           every step, with dependency resolution, caching, logging and the
#           pre-registration guard. Contains no analysis logic.
#
# Usage (from Aim3/):
#   Rscript Code/run_all.R                       # all steps that are implemented
#   Rscript Code/run_all.R --quick               # smoke test
#   Rscript Code/run_all.R --steps W2            # selected steps (+ dependencies)
#   Rscript Code/run_all.R --plots-only          # redraw figures from saved results
#   Rscript Code/run_all.R --force               # ignore the cache
#   Rscript Code/run_all.R --cores 8
# W3 on the HPC: see Code/README.md (hpc/w3_array.sbatch); the local W3 step
# then only merges the parts in Results/real/objects/W3_parts/.
# =============================================================================

.aim3_code_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(dirname(f)) else {
    of <- NULL
    for (i in rev(seq_len(sys.nframe()))) if (!is.null(sys.frame(i)$ofile)) { of <- sys.frame(i)$ofile; break }
    if (is.null(of)) normalizePath("Code") else normalizePath(dirname(of))
  }
})

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source(file.path(.aim3_code_dir, "config.R"))
for (.f in list.files(file.path(.aim3_code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(.f)
source(file.path(.aim3_code_dir, "tests", "test_model.R"))
rm(.f)

# ---- step registry ---------------------------------------------------------------------------------
aim3_steps <- list(
  W0 = list(fun = "an_W0_data",            deps = character(0), uses_data = TRUE, plot = "plot_W0"),
  W1 = list(fun = "an_W1_model",           deps = "W0",         uses_data = TRUE, plot = "plot_W1"),
  W2 = list(fun = "an_W2_fingerprint",     deps = "W1",         uses_data = TRUE, plot = "plot_W2"),
  W3 = list(fun = "an_W3_identifiability", deps = "W0",         uses_data = TRUE, plot = "plot_W3"),
  PR = list(fun = "an_PR_prespecify",      deps = "W0",         uses_data = TRUE, plot = "plot_PR")
)

parse_args <- function(args) {
  get <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
  list(steps = if (is.null(get("--steps"))) NULL else strsplit(get("--steps"), ",")[[1]],
       quick = "--quick" %in% args, force = "--force" %in% args, plots_only = "--plots-only" %in% args,
       cores = if (is.null(get("--cores"))) NULL else as.integer(get("--cores")))
}

step_key <- function(step, cfg, dhash, upstream_keys) {
  files <- c(list.files(file.path(cfg$paths$code, "R"), full.names = TRUE),
             file.path(cfg$paths$code, "src", "hsc_model.cpp"),
             file.path(cfg$paths$code, "tests", "test_model.R"))
  digest::digest(list(step = step, quick = cfg$quick,
                      cfg = cfg[c("data", "phenotype", "model", "fit", "steps")],
                      data = if (aim3_steps[[step]]$uses_data) dhash else NULL,
                      code = vapply(sort(files), function(f) digest::digest(file = f), character(1)),
                      upstream = upstream_keys))
}

resolve_steps <- function(requested) {
  out <- character(0)
  add <- function(s) {
    for (d in aim3_steps[[s]]$deps) add(d)
    if (!s %in% out) out <<- c(out, s)
  }
  for (s in requested) add(s)
  names(aim3_steps)[names(aim3_steps) %in% out]
}

run_aim3 <- function(steps = NULL, quick = FALSE, force = FALSE, plots_only = FALSE, cores = NULL) {
  root <- normalizePath(file.path(.aim3_code_dir, ".."))
  cfg <- make_config(quick = quick, cores = cores, root = root)
  if (quick) cfg$paths$results <- file.path(root, "Results", "quick")
  ensure_result_dirs(cfg)
  log <- make_logger(cfg)
  write_session_info(cfg)
  load_hsc_model(cfg)
  log("Aim 3 run | quick: ", quick, " | cores: ", cfg$runtime$cores)

  todo <- resolve_steps(if (is.null(steps)) names(aim3_steps) else steps)
  data_cache <- NULL
  get_data <- function() {
    if (is.null(data_cache)) data_cache <<- convert_icaht(cfg)
    data_cache
  }
  dhash <- data_hash(cfg)
  results <- list()
  keys <- list()

  for (s in todo) {
    st <- aim3_steps[[s]]
    key <- step_key(s, cfg, dhash, unlist(keys[st$deps]))
    keys[[s]] <- key
    if (plots_only) {
      results[[s]] <- load_step_result(s, cfg)
      next
    }
    if (s == "PR" && prereg_exists(cfg, "PR_plan.sha256")) {
      v <- verify_prereg(cfg, "PR_plan.sha256")
      log("  PR kept: the pre-registration exists and is never regenerated (", v$msg, ")")
      results[[s]] <- load_step_result(s, cfg)
      next
    }
    if (!force && identical(cached_key(s, cfg), key)) {
      log("  ", s, ": cached")
      results[[s]] <- load_step_result(s, cfg)
      next
    }
    for (dd in st$deps) if (is.null(results[[dd]])) results[[dd]] <- load_step_result(dd, cfg)
    ctx <- list(config = cfg, data = get_data, results = results, log = log)
    log("> ", s, " running")
    t0 <- Sys.time()
    res <- do.call(st$fun, list(ctx))
    if (s == "W1" && !isTRUE(res$all_passed)) {
      save_step_result(res, s, key, cfg)
      write_step_tables(res$tables, s, cfg)
      stop("W1 unit tests failed - see tables/W1_unit_tests.csv")
    }
    save_step_result(res, s, key, cfg)
    write_step_tables(res$tables, s, cfg)
    if (s == "W0") write_derived(res$derived, cfg)
    if (!is.null(res$prereg)) {
      write_prereg(res$prereg$files, cfg, res$prereg$hash)
      log("  ", s, ": pre-registration written and hashed")
    }
    results[[s]] <- res
    log("  ", s, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  }

  for (s in names(results)) {
    res <- results[[s]]
    if (is.null(res)) next
    plots <- tryCatch(do.call(aim3_steps[[s]]$plot, list(res, cfg)),
                      error = function(e) { log("  plot ", s, " failed: ", conditionMessage(e)); list() })
    save_step_plots(plots, s, cfg)
  }
  invisible(results)
}

if (sys.nframe() == 0L) {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  do.call(run_aim3, a)
}
