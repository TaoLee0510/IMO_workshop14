# =============================================================================
# io_write.R
# Purpose : Everything that touches Results/: directories, cached step
#           objects, tables, figures, logs, and pre-registration files with
#           SHA-256 hashes. This layer writes; it does not compute.
# Inputs  : step results (lists with optional $tables, $prereg), ggplots.
# Outputs : Results/<source>/{objects,tables,figures,prereg,logs}/...
# =============================================================================

result_dirs <- c("objects", "tables", "figures", "prereg", "logs")

ensure_result_dirs <- function(config) {
  for (d in result_dirs) dir.create(file.path(config$paths$results, d), recursive = TRUE,
                                    showWarnings = FALSE)
  invisible(TRUE)
}

res_path <- function(config, sub, file) file.path(config$paths$results, sub, file)

# ---- logging -------------------------------------------------------------------------------------
make_logger <- function(config) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  path <- res_path(config, "logs", paste0("run_", stamp, ".log"))
  function(...) {
    msg <- paste0(format(Sys.time(), "%H:%M:%S"), "  ", paste0(...))
    message(msg)
    cat(msg, "\n", file = path, append = TRUE, sep = "")
  }
}

write_session_info <- function(config) {
  writeLines(utils::capture.output(utils::sessionInfo()),
             res_path(config, "logs", "sessionInfo.txt"))
}

# ---- cached step objects ---------------------------------------------------------------------------
save_step_result <- function(res, step, key, config) {
  attr(res, "cache_key") <- key
  attr(res, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  saveRDS(res, res_path(config, "objects", paste0(step, ".rds")))
  invisible(res)
}

load_step_result <- function(step, config) {
  p <- res_path(config, "objects", paste0(step, ".rds"))
  if (!file.exists(p)) return(NULL)
  readRDS(p)
}

cached_key <- function(step, config) {
  r <- load_step_result(step, config)
  if (is.null(r)) NA_character_ else attr(r, "cache_key")
}

# ---- tables ----------------------------------------------------------------------------------------------
write_step_tables <- function(tables, step, config) {
  tables <- Filter(Negate(is.null), tables)
  if (!length(tables)) return(character(0))
  vapply(names(tables), function(nm) {
    p <- res_path(config, "tables", paste0(step, "_", nm, ".csv"))
    data.table::fwrite(data.table::as.data.table(tables[[nm]]), p)
    p
  }, character(1))
}

# ---- figures ---------------------------------------------------------------------------------------------
# plots: named list; each element a ggplot or list(plot = , width = , height = )
save_step_plots <- function(plots, step, config) {
  if (!length(plots)) return(character(0))
  out <- character(0)
  for (nm in names(plots)) {
    p <- plots[[nm]]
    w <- config$plot$width
    h <- config$plot$height
    if (is.list(p) && !inherits(p, "gg")) {
      w <- p$width %||% w
      h <- p$height %||% h
      p <- p$plot
    }
    base <- res_path(config, "figures", paste0(step, "_", nm))
    ok <- tryCatch({
      ggplot2::ggsave(paste0(base, ".pdf"), p, width = w, height = h, device = grDevices::pdf)
      ggplot2::ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = config$plot$dpi)
      TRUE
    }, error = function(e) {
      message("  figure ", step, "_", nm, " failed: ", conditionMessage(e))
      FALSE
    })
    if (ok) out <- c(out, paste0(base, ".pdf"))
  }
  out
}


# ---- pre-registration ---------------------------------------------------------------------------------------
sha256_file <- function(path) digest::digest(file = path, algo = "sha256")

# files: named list, name = file name in prereg/, value = data.table (csv), list (json) or character (text)
write_prereg <- function(files, config, hash_name) {
  paths <- character(0)
  for (nm in names(files)) {
    p <- res_path(config, "prereg", nm)
    x <- files[[nm]]
    if (data.table::is.data.table(x) || is.data.frame(x)) data.table::fwrite(x, p)
    else if (grepl("\\.json$", nm)) jsonlite::write_json(x, p, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    else writeLines(x, p)
    paths <- c(paths, p)
  }
  hashes <- vapply(paths, sha256_file, character(1))
  writeLines(paste(hashes, basename(paths)), res_path(config, "prereg", hash_name))
  invisible(hashes)
}

verify_prereg <- function(config, hash_name) {
  hp <- res_path(config, "prereg", hash_name)
  if (!file.exists(hp)) return(list(ok = FALSE, msg = paste("missing", hash_name)))
  lines <- readLines(hp)
  for (ln in lines) {
    parts <- strsplit(ln, " ", fixed = TRUE)[[1]]
    p <- res_path(config, "prereg", parts[2])
    if (!file.exists(p)) return(list(ok = FALSE, msg = paste("missing", parts[2])))
    if (sha256_file(p) != parts[1]) return(list(ok = FALSE, msg = paste("hash mismatch:", parts[2])))
  }
  list(ok = TRUE, msg = paste("verified", length(lines), "file(s)"))
}

prereg_exists <- function(config, hash_name) file.exists(res_path(config, "prereg", hash_name))
