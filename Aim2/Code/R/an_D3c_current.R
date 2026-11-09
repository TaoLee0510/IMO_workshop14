# =============================================================================
# an_D3c_current.R
# Purpose : D3c - what the data we have now can say about the clone. D3b fixed
#           every parameter except the clone traits, so it is an upper bound.
#           Here whole cohorts at the TestData scale (350 patients, 15
#           sequenced pre-treatment, 70% training) are simulated with
#           simulate_cohort() and the D4 procedure is repeated on every
#           training split with everything estimated:
#             M0 -> M1 -> M2 -> profile likelihood over (delta, r_inf), with the
#             whole M2 population refitted at every grid point
#           Two variants share each cohort and differ only in the f0 prior:
#             indices    : f0 covariates include RDW / MCV
#             no_indices : RDW / MCV removed (CBC without red-cell indices)
# Inputs  : ctx$config (no patient data).
# Outputs : list(per_rep, summary, surface, tables, key_numbers)
#
# Metrics per cohort x variant (truth from the simulation)
#   llr_low_vs_functional : as in D3b; truth is low-output, so > 0 is correct
#   argmax_is_truth, region_size : as in D3b
#   h1 / h3 counts        : MAP mechanism attribution (as in D4) for truly
#                           persistent training patients, at the profile maximum
#   p_h3_*                : posterior P(H3) (mechanism_probability) for the same
#                           patients; sum_p_h3 is the expected number of clonal
#                           cases, to compare with n_h3
#   cor_logit_f0          : estimated vs true clone fraction, CH carriers
#
# Parallel runs: each cohort x variant is one "part" (d3c_part). Parts can be
# computed elsewhere (Code/hpc/d3c_part.R, one Slurm array task each) and saved
# under Results/<source>/objects/D3c_parts/. If that folder holds parts, the
# step only combines them; they must all exist and carry the current code key.
# =============================================================================

d3c_labels <- c(indices = "with RDW / MCV", no_indices = "without RDW / MCV")

# Training data for one variant (red-cell covariates removed for no_indices).
d3c_variant_data <- function(tr, variant) {
  if (variant != "no_indices") return(tr)
  drop <- intersect(c("log_rdw_c", "mcv_c"), names(tr$patients))
  if (length(drop)) tr$patients <- data.table::copy(tr$patients)[, (drop) := NULL]
  tr
}

d3c_metrics <- function(r, variant, tr, truth, prof, inferred, mech_prob, is_truth, s) {
  tab <- prof$table
  j <- which.max(tab$logLik)
  low <- tab$delta <= s$low_output_max
  x <- merge(inferred[, .(patient_id, mechanism, f0)],
             truth[, .(patient_id, ch, f0_true = f0, mechanism_true)], by = "patient_id")
  pers <- merge(x[mechanism_true %in% c("H1", "H3")], mech_prob, by = "patient_id", all.x = TRUE)
  carriers <- x[ch == 1]
  data.table::data.table(
    rep = r, variant = variant, n_train = nrow(x),
    n_sequenced = if (is.null(tr$vaf)) 0L else nrow(tr$vaf),
    argmax_delta = prof$best[["delta"]], argmax_r_inf = prof$best[["r_inf"]],
    argmax_is_truth = is_truth[j],
    llr_low_vs_functional = max(tab$logLik[low]) - max(tab$logLik[!low]),
    region_size = sum(tab$logLik >= max(tab$logLik) - s$region_ll),
    n_h1 = sum(pers$mechanism_true == "H1"), n_h3 = sum(pers$mechanism_true == "H3"),
    h1_correct = sum(pers$mechanism_true == "H1" & pers$mechanism == "H1"),
    h3_correct = sum(pers$mechanism_true == "H3" & pers$mechanism == "H3"),
    h1_called_h3 = sum(pers$mechanism_true == "H1" & pers$mechanism == "H3"),
    sum_p_h3_h3 = sum(pers$p_h3[pers$mechanism_true == "H3"]),
    sum_p_h3_h1 = sum(pers$p_h3[pers$mechanism_true == "H1"]),
    cor_logit_f0 = stats::cor(stats::qlogis(carriers$f0), stats::qlogis(carriers$f0_true))
  )[, `:=`(h3_rate = if (n_h3 > 0) h3_correct / n_h3 else NA_real_,
           h1_fp_rate = if (n_h1 > 0) h1_called_h3 / n_h1 else NA_real_,
           sum_p_h3 = sum_p_h3_h3 + sum_p_h3_h1)][]
}

# Settings shared by the step and by single parts.
d3c_setup <- function(cfg) {
  s <- cfg$steps$D3c
  spec <- truth_spec("main")
  spec$calib$n <- s$calib_n
  struct <- cfg$model$struct
  struct[names(spec$struct_override)] <- spec$struct_override
  grid <- do.call(data.table::CJ, cfg$fit$profile)
  is_truth <- abs(grid$delta - struct[["delta"]]) < 1e-9 & abs(grid$r_inf - struct[["r_inf"]]) < 1e-9
  list(s = s, spec = spec, struct = struct, is_truth = is_truth)
}

d3c_calibrate <- function(cfg, st) {
  set.seed(cfg$runtime$seed + 41L)
  suppressMessages(calibrate_population(cfg, st$spec, st$struct))
}

# Key of the code and settings a part was computed with (plain text hashed, so
# it is the same on every machine running the same code).
d3c_part_key <- function(cfg) {
  code <- cfg$paths$code
  files <- c(Sys.glob(file.path(code, "R", c("model_*.R", "sim_*.R"))),
             file.path(code, "R", c("phenotype.R", "likelihood.R", "fit_its.R", "io_read.R",
                                    "an_D4_popfit.R", "an_D3c_current.R")),
             file.path(code, "src", "marrow_model.cpp"))
  txt <- c(vapply(sort(files, method = "radix"), function(f) paste(readLines(f, warn = FALSE), collapse = "\n"), character(1)),
           deparse(list(cfg$model, cfg$fit, cfg$phenotype, cfg$steps$D3c, cfg$runtime$seed, cfg$quick)))
  digest::digest(paste(txt, collapse = "\n"), algo = "md5", serialize = FALSE)
}

d3c_parts_dir <- function(cfg) file.path(cfg$paths$results, "objects", "D3c_parts")
d3c_part_file <- function(cfg, r, variant) file.path(d3c_parts_dir(cfg), sprintf("rep%02d_%s.rds", r, variant))

# One cohort x variant: simulate the cohort, fit M0 -> M1 -> M2, profile with
# refits, attribute mechanisms.
d3c_part <- function(cfg, r, variant, cal, st, log_fun = message) {
  quiet <- function(...) invisible(NULL)
  t0 <- Sys.time()
  sim <- suppressMessages(simulate_cohort(cfg, "main", seed = cfg$runtime$seed + 1000L + r,
                                          cal = cal, spec = st$spec))
  tr <- analysis_set(standardise(sim$raw, cfg), "train")
  truth <- sim$truth$patient[patient_id %in% tr$patients$patient_id]
  pds <- build_pdats(tr, cfg)
  m0 <- fit_population(pds, "M0", cfg, log_fun = quiet)
  m1 <- fit_population(pds, "M1", cfg, init = m0$pop, log_fun = quiet)
  trv <- d3c_variant_data(tr, variant)
  pdv <- build_pdats(trv, cfg)
  m2 <- fit_population(pdv, "M2", cfg, init = m1$pop, log_fun = quiet)
  prof <- profile_clone_traits(pdv, m2, cfg, log_fun = log_fun)
  inferred <- d4_inferred(prof$fit_best, trv, cfg)
  pers_ids <- truth[mechanism_true %in% c("H1", "H3"), patient_id]
  mech_prob <- mechanism_probability(prof$fit_best, trv, pers_ids, cfg)
  row <- d3c_metrics(r, variant, trv, truth, prof, inferred, mech_prob, st$is_truth, st$s)
  log_fun(sprintf("  D3c: cohort %d, %-10s argmax delta=%g r_inf=%g, LLR %.1f, H3 %d/%d (sum P(H3) %.1f), H1->H3 %d/%d, cor f0 %.2f (%.1f min)",
                  r, variant, row$argmax_delta, row$argmax_r_inf, row$llr_low_vs_functional,
                  row$h3_correct, row$n_h3, row$sum_p_h3_h3, row$h1_called_h3, row$n_h1, row$cor_logit_f0,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  list(key = d3c_part_key(cfg), rep = r, variant = variant, row = row,
       surface = prof$table[, `:=`(rep = r, variant = variant, is_truth = st$is_truth)][],
       mech_prob = merge(mech_prob, truth[, .(patient_id, mechanism_true)], by = "patient_id")[
         , `:=`(rep = r, variant = variant)][])
}

an_D3c_current <- function(ctx) {
  cfg <- ctx$config
  st <- d3c_setup(cfg)
  s <- st$s
  todo <- data.table::CJ(rep = seq_len(s$n_rep), variant = s$variants, sorted = FALSE)
  files <- d3c_part_file(cfg, todo$rep, todo$variant)
  if (any(file.exists(files))) {
    if (!all(file.exists(files)))
      stop("D3c: ", sum(!file.exists(files)), " of ", length(files), " parts missing in ", d3c_parts_dir(cfg),
           " (finish the array job, or delete the folder to compute here)")
    parts <- lapply(files, readRDS)
    key <- d3c_part_key(cfg)
    if (!all(vapply(parts, `[[`, character(1), "key") == key))
      stop("D3c: parts in ", d3c_parts_dir(cfg), " were computed with other code or settings",
           " (rerun the array job, or delete the folder to compute here)")
    ctx$log("  D3c: combining ", length(parts), " parts from ", d3c_parts_dir(cfg))
  } else {
    ctx$log("  D3c: calibrating the generating population (n = ", s$calib_n, ")")
    cal <- d3c_calibrate(cfg, st)
    parts <- lapply(seq_len(nrow(todo)), function(i)
      d3c_part(cfg, todo$rep[i], todo$variant[i], cal, st, log_fun = ctx$log))
  }
  per_rep <- data.table::rbindlist(lapply(parts, `[[`, "row"))
  surface <- data.table::rbindlist(lapply(parts, `[[`, "surface"))
  mech_prob <- data.table::rbindlist(lapply(parts, `[[`, "mech_prob"))

  summ <- per_rep[, .(cohorts = .N, n_train = round(mean(n_train)), n_sequenced = round(mean(n_sequenced), 1),
                      llr_median = round(stats::median(llr_low_vs_functional), 1),
                      llr_min = round(min(llr_low_vs_functional), 1),
                      llr_max = round(max(llr_low_vs_functional), 1),
                      clone_type_correct = round(mean(llr_low_vs_functional > 0), 2),
                      llr_above_2 = round(mean(llr_low_vs_functional > 2), 2),
                      argmax_is_truth = round(mean(argmax_is_truth), 2),
                      region_size = round(mean(region_size), 1),
                      persistent_per_cohort = round(mean(n_h1 + n_h3), 1),
                      h3_per_cohort = round(mean(n_h3), 1),
                      h3_sensitivity = round(sum(h3_correct) / max(sum(n_h3), 1), 2),
                      h1_called_h3_rate = round(sum(h1_called_h3) / max(sum(n_h1), 1), 3),
                      mean_p_h3_clonal = round(sum(sum_p_h3_h3) / max(sum(n_h3), 1), 2),
                      mean_p_h3_low_reserve = round(sum(sum_p_h3_h1) / max(sum(n_h1), 1), 3),
                      expected_over_true_h3 = round(sum(sum_p_h3) / max(sum(n_h3), 1), 2),
                      cor_logit_f0 = round(mean(cor_logit_f0), 3)), by = variant]
  summ[, variant := factor(variant, s$variants)]
  data.table::setorder(summ, variant)
  summ[, variant := as.character(variant)]
  surf_mean <- surface[, .(delta_logLik = mean(delta_logLik)), by = .(variant, delta, r_inf, is_truth)]

  list(per_rep = per_rep, summary = summ, surface = surface, surface_mean = surf_mean, mech_prob = mech_prob,
       truth_traits = st$struct[c("delta", "r_inf")], settings = s,
       tables = list(summary = summ, per_rep = per_rep, profile_surface = surface,
                     mechanism_probability = mech_prob),
       key_numbers = list(summary = summ))
}
