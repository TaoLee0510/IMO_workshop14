# =============================================================================
# run_all.R  -  single entry point for Model3
# Purpose : Orchestrate the layers (read -> analyse -> write -> plot) for every
#           step, with dependency resolution, caching, logging and the
#           pre-registration guard (no real-data fit before PR is archived).
#           Contains no analysis logic.
#
# Usage (from Model3/):
#   Rscript Code/run_all.R                         # TestData (synthetic), all steps
#   Rscript Code/run_all.R --quick                 # smoke test
#   Rscript Code/run_all.R --steps S0,K1,K2        # selected steps (+ dependencies)
#   Rscript Code/run_all.R --data real --steps S0,K1,K2,K3,PR
#   Rscript Code/run_all.R --data real             # real data (needs the archived PR)
#   Rscript Code/run_all.R --plots-only            # redraw figures from saved results
#   Rscript Code/run_all.R --force                 # ignore the cache
#   Rscript Code/run_all.R --make-testdata         # regenerate TestData
# From R:  source("Code/run_all.R"); run_model3(data = "test", steps = c("K1"))
# =============================================================================

.model3_code_dir <- local({
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
source(file.path(.model3_code_dir, "config.R"))
for (.f in list.files(file.path(.model3_code_dir, "R"), pattern = "\\.R$", full.names = TRUE)) source(.f)
rm(.f)

# ---- step registry -------------------------------------------------------------------------------
# fun: analysis function; deps: upstream steps; file: analysis sources (cache key);
# uses_data: whether the data hash enters the cache key; plot: plotting function
model3_steps <- list(
  S0 = list(fun = "an_S0_data",      deps = character(0), file = "an_S0_data.R",      uses_data = TRUE,  plot = "plot_S0"),
  K1 = list(fun = "an_K1_model",     deps = character(0), file = c("model_checks.R", "an_K1_model.R"), uses_data = FALSE, plot = "plot_K1"),
  K2 = list(fun = "an_K2_phasemap",  deps = "K1",         file = "an_K2_phasemap.R",  uses_data = FALSE, plot = "plot_K2"),
  PR = list(fun = "an_PR_prereg",    deps = "S0",         file = "an_PR_prereg.R",   uses_data = TRUE,  plot = "plot_none"),
  F1 = list(fun = "an_F1_fit",       deps = c("S0", "PR"), file = "an_F1_fit.R",      uses_data = TRUE,  plot = "plot_F1"),
  F2 = list(fun = "an_F2_attrib",    deps = "F1",         file = "an_F2_attrib.R",    uses_data = TRUE,  plot = "plot_F2"),
  F3 = list(fun = "an_F3_landmark",  deps = "F1",         file = "an_F3_landmark.R",  uses_data = TRUE,  plot = "plot_F3"),
  F4 = list(fun = "an_F4_whatif",    deps = "F2",         file = "an_F4_whatif.R",    uses_data = TRUE,  plot = "plot_F4"),
  C  = list(fun = "an_C_compare",    deps = c("F1", "F3"), file = "an_C_compare.R",   uses_data = TRUE,  plot = "plot_none"),
  G  = list(fun = "an_G_export",     deps = c("F1", "F2"), file = "an_G_export.R",    uses_data = TRUE,  plot = "plot_none"),
  V  = list(fun = "an_V_validation", deps = "F2",         file = "an_V_validation.R", uses_data = TRUE,  plot = "plot_V"),
  # S: synthetic validation S1-S4 (cohorts + cluster fits) and the VPC S5 of the final fit
  S  = list(fun = "an_S_validation", deps = "F1",         file = c("sim_cohort.R", "an_S_validation.R", "an_F3_landmark.R"),
            uses_data = TRUE, plot = "plot_none")
)
core_files <- c("R/model.R", "R/phenotype.R", "R/likelihood.R", "R/fit_its.R", "R/fit_de.R", "R/fit_mcmc.R",
                "R/fit_mcmc13.R", "R/attribution13.R",
                "R/attribution.R", "R/io_read.R", "R/io_convert.R", "src/model3.cpp")

parse_args <- function(args) {
  get <- function(flag, default = NULL) if (flag %in% args) args[which(args == flag) + 1] else default
  list(data = get("--data", "test"),
       steps = if (is.null(get("--steps"))) NULL else strsplit(get("--steps"), ",")[[1]],
       quick = "--quick" %in% args, force = "--force" %in% args,
       plots_only = "--plots-only" %in% args, cores = if (is.null(get("--cores"))) NULL else as.integer(get("--cores")),
       make_testdata = "--make-testdata" %in% args)
}

step_key <- function(step, cfg, data_hash, upstream_keys) {
  st <- model3_steps[[step]]
  files <- file.path(cfg$paths$code, c(core_files, file.path("R", st$file)))
  if (step == "K1") files <- c(files, file.path(cfg$paths$code, c("tests/test_model.R", "tests/ref/aim4_v13.cpp", "R/aim4_ref.R")))
  digest::digest(list(step = step, quick = cfg$quick,
                      cfg = cfg[c("data", "phenotype", "model", "fit", "steps")],
                      data = if (st$uses_data) data_hash else NULL,
                      code = vapply(files, function(f) digest::digest(file = f), character(1)),
                      upstream = upstream_keys))
}

resolve_steps <- function(requested) {
  out <- character(0)
  add <- function(s) {
    for (d in model3_steps[[s]]$deps) add(d)
    if (!s %in% out) out <<- c(out, s)
  }
  for (s in requested) add(s)
  names(model3_steps)[names(model3_steps) %in% out]
}

run_model3 <- function(data = "test", steps = NULL, quick = FALSE, force = FALSE, plots_only = FALSE,
                     cores = NULL, make_testdata = FALSE) {
  root <- normalizePath(file.path(.model3_code_dir, ".."))
  cfg_real <- make_config("real", root = root)
  load_model3(cfg_real)
  ensure_real_derived(cfg_real)

  if (make_testdata || (data == "test" && !file.exists(file.path(cfg_real$paths$testdata, "patients.csv")))) {
    make_testdata(make_config("test", root = root))
    if (make_testdata && is.null(steps)) return(invisible(NULL))
  }

  cfg <- make_config(data, quick = quick, cores = cores, root = root)
  ensure_result_dirs(cfg)
  log <- make_logger(cfg)
  write_session_info(cfg)
  log("Model3 run | data: ", data, " | quick: ", quick, " | cores: ", cfg$runtime$cores)

  requested <- if (is.null(steps)) names(model3_steps) else steps
  if (data == "real") requested <- setdiff(requested, "V")
  todo <- resolve_steps(requested)

  data_cache <- NULL
  get_data <- function() {
    if (is.null(data_cache)) data_cache <<- read_model3_data(cfg)
    data_cache
  }
  design_cache <- NULL
  get_design <- function() {
    if (is.null(design_cache)) design_cache <<- read_model3_data(cfg_real)
    design_cache
  }
  data_hash <- data_files_hash(cfg)
  results <- list()
  keys <- list()

  for (s in todo) {
    st <- model3_steps[[s]]
    up_keys <- unlist(keys[st$deps])
    key <- step_key(s, cfg, data_hash, up_keys)
    keys[[s]] <- key

    if (plots_only) {
      results[[s]] <- load_step_result(s, cfg)
      next
    }
    # ---- guards ----
    if (s == "PR" && data == "real" && prereg_exists(cfg, "PR.sha256")) {
      log("  PR kept: the pre-registration already exists and is never regenerated for real data")
      results[[s]] <- load_step_result(s, cfg)
      keys[[s]] <- attr(results[[s]], "cache_key") %||% key
      next
    }
    if (s == "F1" && data == "real") {
      v <- verify_prereg(cfg, "PR.sha256")
      if (!v$ok) { log("  F1 refused: ", v$msg, " - archive PR first (--steps PR)"); break }
      log("  PR verified: ", v$msg)
    }
    if (!force && identical(cached_key(s, cfg), key)) {
      log("  ", s, ": cached")
      results[[s]] <- load_step_result(s, cfg)
      next
    }

    # ---- run ----
    for (dd in st$deps) if (is.null(results[[dd]])) results[[dd]] <- load_step_result(dd, cfg)
    ctx <- list(config = cfg, data = get_data, design = get_design, results = results, log = log)
    log("> ", s, " running")
    t0 <- Sys.time()
    res <- do.call(st$fun, list(ctx))
    if (s == "K1" && !isTRUE(res$all_passed)) {
      save_step_result(res, s, key, cfg)
      write_step_tables(res$tables, s, cfg)
      stop("K1 unit tests failed - see tables/K1_unit_tests.csv")
    }
    save_step_result(res, s, key, cfg)
    write_step_tables(res$tables, s, cfg)
    if (!is.null(res$prereg)) {
      write_prereg(res$prereg$files, cfg, "PR.sha256")
      log("  ", s, ": pre-registration written")
    }
    results[[s]] <- res
    log("  ", s, " done in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s")
  }

  # ---- plots ----
  for (s in names(results)) {
    res <- results[[s]]
    if (is.null(res)) next
    pf <- model3_steps[[s]]$plot
    plots <- tryCatch(do.call(pf, list(res, cfg, results = results)),
                      error = function(e) { log("  plot ", s, " failed: ", conditionMessage(e)); list() })
    save_step_plots(plots, s, cfg)
  }

  # ---- summary ----
  all_res <- lapply(stats::setNames(names(model3_steps), names(model3_steps)),
                    function(s) results[[s]] %||% load_step_result(s, cfg))
  writeLines(report_summary(all_res, cfg), file.path(cfg$paths$results, "summary.md"))
  log("Summary written: ", file.path(cfg$paths$results, "summary.md"))
  invisible(results)
}

if (sys.nframe() == 0L) {
  a <- parse_args(commandArgs(trailingOnly = TRUE))
  do.call(run_model3, a)
}
