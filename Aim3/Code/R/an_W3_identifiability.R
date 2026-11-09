# =============================================================================
# an_W3_identifiability.R
# Purpose : Step W3 - how well chemotherapy (a) and CAR T-phase (b) stem cell
#           injury can be told apart at the real sampling design. Virtual
#           cohorts (sim_cohort.R) are refitted with all five models; the
#           population fit of Mab feeds the per-patient attribution, which is
#           compared with the true mechanism class.
#           One task = (scenario, replicate, cohort size); tasks run locally or
#           as a Slurm array (hpc/w3_part.R) and are merged here.
# Inputs  : ctx (W0 derived data).
# Outputs : model-selection and attribution summaries, the reporting route.
# =============================================================================

w3_tasks <- function(config) {
  st <- config$steps$W3
  small <- data.table::CJ(scenario = st$scenarios, rep = seq_len(st$n_rep), n = 35L)
  large <- if (st$n_rep_large > 0)
    data.table::CJ(scenario = st$scenarios, rep = seq_len(st$n_rep_large), n = as.integer(st$n_large))
  out <- rbind(small, large)
  out[, task_id := .I]
  out[]
}

# Code and settings fingerprint shared by the HPC parts and the merge step.
w3_fingerprint <- function(config) {
  files <- c(list.files(file.path(config$paths$code, "R"), full.names = TRUE),
             file.path(config$paths$code, "src", "hsc_model.cpp"))
  f <- tempfile()
  writeLines(deparse(list(model = config$model, fit = config$fit, W3 = config$steps$W3,
                          phenotype = config$phenotype, data = config$data), control = "digits17"), f)
  h <- tools::md5sum(c(sort(files), f))
  g <- tempfile()
  writeLines(paste(basename(names(h)), h, collapse = "\n"), g)
  unname(tools::md5sum(g))
}

# Resample the real designs to n virtual patients (unique ids).
expand_design <- function(d, n, seed) {
  if (n == nrow(d$patients)) return(d)
  set.seed(seed)
  src <- sample(d$patients$patient_id, n, replace = TRUE)
  new_id <- sprintf("V%03d_%s", seq_len(n), src)
  map <- data.table::data.table(patient_id = src, new_id = new_id)
  ren <- function(tb) {
    x <- merge(map, tb, by = "patient_id", allow.cartesian = TRUE)
    x[, patient_id := new_id][, new_id := NULL][]
  }
  list(patients = ren(d$patients), obs = ren(d$obs), inflam = ren(d$inflam))
}

w3_run_task <- function(task, d, config, cores = config$runtime$cores, log_fun = message) {
  st <- config$model$struct
  scen <- task$scenario
  seed <- config$runtime$seed + 1000L * task$rep + 7L * match(scen, config$steps$W3$scenarios) + task$n
  real_designs <- design_context(d, config)$pdes_list
  shift <- calibrate_shift(real_designs, scen, config)
  dx <- expand_design(d, task$n, seed)
  pdes_x <- stats::setNames(lapply(dx$patients$patient_id, function(pid)
    make_pdes(dx$patients[patient_id == pid], dx$inflam[patient_id == pid], st, config)), dx$patients$patient_id)
  sim <- simulate_cohort(dx, pdes_x, scen, shift, config, seed)
  pdats <- cohort_pdats(dx, sim, config)
  quiet <- function(...) invisible(NULL)
  fits <- fit_all_models(pdats, config, cores = cores, log_fun = quiet)
  mt <- model_table(fits)
  n_att <- min(length(pdats), config$steps$W3$attrib_max)
  at <- attribute_cohort(pdats[seq_len(n_att)], fits$Mab, config, cores = cores, seed = seed)
  at <- merge(at, sim$truth, by = "patient_id")
  log_fun(sprintf("  W3 task %d (%s, rep %d, n %d): best %s, shift %.2f", task$task_id, scen, task$rep,
                  task$n, mt[which.min(BIC), model], shift))
  list(task = task, shift = shift, models = mt, attrib = at,
       pheno_mix = sim$truth[, .N, by = pheno],
       pop = lapply(fits, function(f) f$pop))
}

w3_summarise <- function(parts, config) {
  models <- data.table::rbindlist(lapply(parts, function(p) cbind(p$task, p$models)))
  best <- models[, .SD[which.min(BIC)], by = .(task_id, scenario, rep, n)]
  best[, truth := vapply(scenario, scenario_truth_model, character(1))]
  lenient <- list(a_only = "Ma", b_only = "Mb", mix = c("Mab", "Mabx"), synergy = c("Mab", "Mabx"))
  best[, correct := model == truth]
  best[, correct_lenient := mapply(function(m, s) m %in% lenient[[s]], model, scenario)]
  sel <- best[, .(n_cohorts = .N, correct = mean(correct), correct_lenient = mean(correct_lenient),
                  chosen = paste(names(table(model)), table(model), sep = ":", collapse = " ")), by = .(scenario, n)]

  att <- data.table::rbindlist(lapply(parts, function(p) cbind(p$task[, .(task_id, scenario, rep, n)], p$attrib)),
                               fill = TRUE)
  pcols <- paste0("p_", mech_classes)
  att[, p_true := mapply(function(i, m) if (is.na(m)) NA_real_ else att[[paste0("p_", m)]][i], .I, mech)]
  indiv <- att[!is.na(mech) & mech != "none",
               .(n_patients = .N, mean_p_true = mean(p_true, na.rm = TRUE),
                 map_accuracy = mean(map_class == mech, na.rm = TRUE)), by = .(scenario, n)]
  confusion <- att[!is.na(mech), .N, by = .(n, truth = mech, map = map_class)]
  calib <- data.table::rbindlist(lapply(mech_classes, function(k) {
    x <- att[!is.na(mech), .(p = get(paste0("p_", k)), y = as.numeric(mech == k))]
    x[, bin := cut(p, c(-0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.01))]
    x[, .(class = k, mean_p = mean(p), observed = mean(y), N = .N), by = bin]
  }))
  shifts <- data.table::rbindlist(lapply(parts, function(p) cbind(p$task, shift = p$shift)))
  mix <- data.table::rbindlist(lapply(parts, function(p) cbind(p$task, p$pheno_mix)))

  s35 <- sel[n == 35]
  i35 <- indiv[n == 35]
  route_cfg <- config$steps$W3$route
  cohort_ok <- nrow(s35) > 0 && min(s35$correct_lenient) >= route_cfg$cohort_correct
  indiv_ok <- nrow(i35) > 0 && mean(i35$mean_p_true) >= route_cfg$indiv_prob
  route <- if (cohort_ok && indiv_ok) 1L else if (cohort_ok) 2L else 3L
  list(selection = sel, best = best, indiv = indiv, confusion = confusion, calibration = calib,
       attrib = att, shifts = shifts, pheno_mix = mix, route = route,
       route_detail = list(cohort_min_correct = if (nrow(s35)) min(s35$correct_lenient) else NA,
                           indiv_mean_p_true = if (nrow(i35)) mean(i35$mean_p_true) else NA))
}

an_W3_identifiability <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$results$W0$derived
  parts_dir <- res_path(cfg, "objects", "W3_parts")
  # fingerprint recorded when the HPC array was submitted, so that code added
  # later (W4-W7) does not invalidate the parts; otherwise the live code
  fp_file <- res_path(cfg, "objects", "W3_expected_fingerprint.txt")
  fp <- if (file.exists(fp_file)) readLines(fp_file)[1] else w3_fingerprint(cfg)
  tasks <- w3_tasks(cfg)
  files <- if (dir.exists(parts_dir)) list.files(parts_dir, pattern = "^W3_task_\\d+\\.rds$", full.names = TRUE) else character(0)
  if (length(files)) {
    parts <- lapply(files, readRDS)
    bad <- vapply(parts, function(p) !identical(p$fingerprint, fp), logical(1))
    if (any(bad)) stop("W3 parts with a different code/settings fingerprint: ",
                       paste(basename(files[bad]), collapse = ", "), " - rerun them or delete ", parts_dir)
    have <- vapply(parts, function(p) p$task$task_id, integer(1))
    missing <- setdiff(tasks$task_id, have)
    if (length(missing)) stop("W3 parts missing for task(s) ", paste(missing, collapse = ", "))
    ctx$log("  W3: merging ", length(parts), " HPC parts")
  } else {
    ctx$log("  W3: running ", nrow(tasks), " task(s) locally")
    parts <- lapply(seq_len(nrow(tasks)), function(i) w3_run_task(tasks[i], d, cfg, log_fun = ctx$log))
  }
  s <- w3_summarise(parts, cfg)
  ctx$log("  W3: route ", s$route, " (cohort min correct ", round(s$route_detail$cohort_min_correct, 2),
          ", mean P(true class) ", round(s$route_detail$indiv_mean_p_true, 2), ")")
  c(s[c("route", "route_detail")],
    list(tables = s[c("selection", "best", "indiv", "confusion", "calibration", "attrib", "shifts", "pheno_mix")]))
}
