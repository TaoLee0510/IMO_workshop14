# =============================================================================
# run_all.R  -  single entry point for Aim 2
# Purpose : Orchestrate the four layers (read -> analyse -> write -> plot) for
#           every step, with dependency resolution, caching, logging, and the
#           hold-out / pre-registration guards. Contains no analysis logic.
#
# Usage (from Aim2/):
#   Rscript Code/run_all.R                              # TestData, all steps
#   Rscript Code/run_all.R --quick                      # smoke test
#   Rscript Code/run_all.R --steps S0,D1,D2             # selected steps (+ dependencies)
#   Rscript Code/run_all.R --plots-only                 # redraw figures from saved results
#   Rscript Code/run_all.R --force                      # ignore the cache
#   Rscript Code/run_all.R --cores 8
#   Rscript Code/run_all.R --make-testdata [--scenario stress]
#   Rscript Code/run_all.R --data test_stress           # run on the stress scenario
#   Rscript Code/run_all.R --data real                  # real data in Data/
#   Rscript Code/run_all.R --data real --steps H --unlock-holdout   # the single hold-out run
# From R:  source("Code/run_all.R"); run_aim2(data = "test", steps = c("D1", "D2"))
# =============================================================================

.aim2_code_dir <- local({
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
source(file.path(.aim2_code_dir, "config.R"))
for (.f in list.files(file.path(.aim2_code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(.f)
rm(.f)

# ---- step registry ---------------------------------------------------------------------------------
# fun: analysis function; deps: upstream steps; file: analysis source (for the cache key);
# uses_data: whether the data hash enters the cache key; plot: plotting function
aim2_steps <- list(
  S0  = list(fun = "an_S0_data",            deps = character(0), file = "an_S0_data.R",            uses_data = TRUE,  plot = "plot_S0"),
  D1  = list(fun = "an_D1_model",           deps = character(0), file = "an_D1_model.R",           uses_data = FALSE, plot = "plot_D1"),
  D2  = list(fun = "an_D2_phasemap",        deps = "D1",         file = "an_D2_phasemap.R",        uses_data = FALSE, plot = "plot_D2"),
  D3  = list(fun = "an_D3_identifiability", deps = "D1",         file = "an_D3_identifiability.R", uses_data = FALSE, plot = "plot_D3"),
  D3b = list(fun = "an_D3b_traits",         deps = "D1",         file = c("an_D3_identifiability.R", "an_D3b_traits.R"), uses_data = FALSE, plot = "plot_D3b"),
  D3c = list(fun = "an_D3c_current",        deps = "D1",         file = c("sim_cohort.R", "an_D4_popfit.R", "an_D3c_current.R"), uses_data = FALSE, plot = "plot_D3c"),
  D4  = list(fun = "an_D4_popfit",          deps = c("S0", "D1"), file = "an_D4_popfit.R",         uses_data = TRUE,  plot = "plot_D4"),
  D3d = list(fun = "an_D3d_realscale",      deps = "D4",         file = c("sim_cohort.R", "events_approx.R", "an_D3d_realscale.R"), uses_data = FALSE, plot = "plot_D3d"),
  D5  = list(fun = "an_D5_ch_subset",       deps = "D4",         file = "an_D5_ch_subset.R",       uses_data = TRUE,  plot = "plot_D5"),
  PR  = list(fun = "an_PR_prespecify",      deps = "D4",         file = "an_PR_prespecify.R",      uses_data = TRUE,  plot = "plot_PR"),
  D5b = list(fun = "an_D5b_compare",        deps = "PR",         file = "an_PR_prespecify.R",      uses_data = TRUE,  plot = "plot_D5b"),
  D6  = list(fun = "an_D6_virtual_cohort",  deps = "D4",         file = "an_D6_virtual_cohort.R",  uses_data = TRUE,  plot = "plot_D6"),
  H   = list(fun = "an_H_holdout",          deps = "D4",         file = "an_H_holdout.R",          uses_data = TRUE,  plot = "plot_H"),
  V   = list(fun = "an_V_validation",       deps = "D4",         file = "an_V_validation.R",       uses_data = TRUE,  plot = "plot_V")
)
core_files <- c("R/model_params.R", "R/model_inputs.R", "R/model_simulate.R", "R/model_mechanism.R",
                "R/sim_clinical.R",
                "R/phenotype.R", "R/likelihood.R", "R/fit_its.R", "R/io_read.R", "src/marrow_model.cpp")

parse_args <- function(args) {
  get <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
  list(data = get("--data", "test"),
       steps = if (is.null(get("--steps"))) NULL else strsplit(get("--steps"), ",")[[1]],
       quick = "--quick" %in% args, force = "--force" %in% args,
       plots_only = "--plots-only" %in% args, cores = if (is.null(get("--cores"))) NULL else as.integer(get("--cores")),
       make_testdata = "--make-testdata" %in% args, scenario = get("--scenario", "main"),
       unlock_holdout = "--unlock-holdout" %in% args)
}

step_key <- function(step, cfg, data_hash, upstream_keys) {
  st <- aim2_steps[[step]]
  files <- file.path(cfg$paths$code, c(core_files, file.path("R", st$file)))
  if (step == "D1") files <- c(files, file.path(cfg$paths$code, "tests", "test_model.R"))
  digest::digest(list(step = step, quick = cfg$quick,
                      cfg = cfg[c("data", "phenotype", "model", "fit", "steps")],
                      data = if (st$uses_data) data_hash else NULL,
                      code = vapply(files, function(f) digest::digest(file = f), character(1)),
                      upstream = upstream_keys))
}

resolve_steps <- function(requested) {
  out <- character(0)
  add <- function(s) {
    for (d in aim2_steps[[s]]$deps) add(d)
    if (!s %in% out) out <<- c(out, s)
  }
  for (s in requested) add(s)
  names(aim2_steps)[names(aim2_steps) %in% out]
}

run_aim2 <- function(data = "test", steps = NULL, quick = FALSE, force = FALSE, plots_only = FALSE,
                     cores = NULL, make_testdata = FALSE, scenario = "main", unlock_holdout = FALSE) {
  root <- normalizePath(file.path(.aim2_code_dir, ".."))

  # ---- test data ----
  if (make_testdata || (data %in% c("test", "test_stress") &&
                        !file.exists(file.path(make_config(data, root = root)$paths$data, "patients.csv")))) {
    source(file.path(root, "Code", "testdata", "make_testdata.R"))
    sc <- if (make_testdata) scenario else if (data == "test_stress") "stress" else "main"
    cfg0 <- make_config("test", root = root)
    load_marrow_model(cfg0)
    make_testdata(cfg0, scenario = sc)
    if (make_testdata && is.null(steps)) return(invisible(NULL))
  }

  cfg <- make_config(data, quick = quick, cores = cores, root = root)
  ensure_result_dirs(cfg)
  log <- make_logger(cfg)
  write_session_info(cfg)
  load_marrow_model(cfg)
  log("Aim 2 run | data: ", data, " | quick: ", quick, " | cores: ", cfg$runtime$cores)

  requested <- if (is.null(steps)) names(aim2_steps) else steps
  requested <- setdiff(requested, cfg$skip_steps)
  todo <- resolve_steps(requested)

  data_cache <- NULL
  get_data <- function() {
    if (is.null(data_cache)) data_cache <<- read_aim2_data(cfg)
    data_cache
  }
  data_hash <- if (file.exists(file.path(cfg$paths$data, "patients.csv"))) data_files_hash(cfg) else NA
  results <- list()
  keys <- list()

  for (s in todo) {
    st <- aim2_steps[[s]]
    up_keys <- unlist(keys[st$deps])
    key <- step_key(s, cfg, data_hash, up_keys)
    keys[[s]] <- key

    if (plots_only) {
      results[[s]] <- load_step_result(s, cfg)
      next
    }
    # ---- guards ----
    if (s == "H" && !identical(cfg$steps$H$mode, "cv")) {
      if (cfg$data_source == "real" && !unlock_holdout) { log("  H skipped: needs --unlock-holdout"); next }
      if (cfg$data_source == "real" && file.exists(holdout_lock_path(cfg))) {
        log("  H refused: the hold-out has already been evaluated (", holdout_lock_path(cfg), ")"); next
      }
      v <- verify_prereg(cfg, "holdout_plan.sha256")
      if (!v$ok) { log("  H refused: ", v$msg); next }
    }
    if (s == "D5b") {
      if (is.null(read_post_mutations(cfg))) { log("  D5b skipped: post-treatment results not available yet"); next }
      v <- verify_prereg(cfg, "PR_prediction.sha256")
      if (!v$ok) { log("  D5b refused: ", v$msg); next }
    }
    if (s == "PR" && cfg$data_source == "real" && prereg_exists(cfg, "PR_prediction.sha256")) {
      log("  PR kept: pre-registration already exists and is never regenerated for real data")
      results[[s]] <- load_step_result(s, cfg)
      next
    }

    if (!force && identical(cached_key(s, cfg), key)) {
      log("  ", s, ": cached")
      results[[s]] <- load_step_result(s, cfg)
      next
    }

    # ---- run ----
    for (dd in st$deps) if (is.null(results[[dd]])) results[[dd]] <- load_step_result(dd, cfg)
    ctx <- list(config = cfg, data = get_data, results = results, log = log)
    if (s == "H") ctx$holdout_plan <- jsonlite::read_json(res_path(cfg, "prereg", "holdout_plan.json"),
                                                         simplifyVector = TRUE)
    if (s == "D5b") {
      ctx$prereg_prediction <- data.table::fread(res_path(cfg, "prereg", "PR_prediction.csv"))
      ctx$post_mutations <- read_post_mutations(cfg)
      ctx$prereg_verification <- verify_prereg(cfg, "PR_prediction.sha256")
    }
    log("> ", s, " running")
    t0 <- Sys.time()
    res <- do.call(st$fun, list(ctx))
    if (s == "D1" && !isTRUE(res$all_passed)) {
      save_step_result(res, s, key, cfg)
      write_step_tables(res$tables, s, cfg)
      stop("D1 unit tests failed - see tables/D1_unit_tests.csv")
    }
    save_step_result(res, s, key, cfg)
    write_step_tables(res$tables, s, cfg)
    if (!is.null(res$prereg)) {
      write_prereg(res$prereg$files, cfg, res$prereg$hash)
      log("  ", s, ": pre-registration written (", res$prereg$hash, ")")
    }
    if (s == "H" && cfg$data_source == "real" && !identical(cfg$steps$H$mode, "cv")) write_holdout_lock(cfg)
    results[[s]] <- res
    log("  ", s, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  }

  # ---- plots ----
  for (s in names(results)) {
    res <- results[[s]]
    if (is.null(res)) next
    pf <- aim2_steps[[s]]$plot
    plots <- tryCatch(do.call(pf, list(res, cfg, res_D2 = results$D2 %||% load_step_result("D2", cfg))),
                      error = function(e) { log("  plot ", s, " failed: ", conditionMessage(e)); list() })
    save_step_plots(plots, s, cfg)
  }

  # ---- summary ----
  all_res <- lapply(stats::setNames(names(aim2_steps), names(aim2_steps)), function(s) results[[s]] %||% load_step_result(s, cfg))
  writeLines(report_summary(all_res, cfg), file.path(cfg$paths$results, "summary.md"))
  log("Summary written: ", file.path(cfg$paths$results, "summary.md"))
  invisible(results)
}

if (sys.nframe() == 0L) {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  do.call(run_aim2, a)
}
