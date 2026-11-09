# =============================================================================
# an_K3_ident.R
# Purpose : Step K3 - what the data can resolve, on both axes, before any real
#           count is fitted. Synthetic cohorts at the real design (real35: the
#           35 patients' sampling days, ferritin, eGFR, bridging; no counts) and
#           at the planned design (n350), with known truth:
#             mechanism  a_only | b_only | mix | synergy
#             clone      none | 35% DDR-type carriers
#           Each cohort is refitted (Q2 family + the clone model) and attributed
#           (MCMC + counterfactuals); truth classes come from the same
#           counterfactuals at the true parameters.
#           Parts run as an HPC array (hpc/k3_part.R); the step combines them.
# Inputs  : ctx (config, design = the converted real extract).
# Outputs : list(tables, summary, routes, parts)
# =============================================================================

k3_tasks <- function(cfg) {
  s <- cfg$steps$K3
  real <- data.table::CJ(rep = seq_len(s$n_rep), clone = names(s$clone), mech = s$mech, sorted = FALSE)
  real[, `:=`(design = "real35", variant = "counts")]
  n350 <- data.table::CJ(rep = seq_len(s$n350$n_rep), variant = c("vaf15", "vafall"), sorted = FALSE)
  n350[, `:=`(design = "n350", mech = "mix", clone = "clone")]
  tasks <- rbind(real, n350, use.names = TRUE)
  if (cfg$quick) tasks <- tasks[(design == "real35" & mech == "mix") | (design == "n350" & variant == "vafall")][
    , .SD[1], by = .(design, clone)]
  tasks[, task := .I]
  tasks[]
}

k3_part_file <- function(cfg, task) file.path(cfg$paths$results, "objects", "K3_parts", sprintf("task%03d.rds", task))

# Truth "fit" object so that counterfactual_draw() can run at the true parameters.
truth_fit <- function(mech, truth, config) {
  spec <- truth_spec(mech)
  list(spec = spec, re = spec$re, pop = list(struct = truth_struct(mech, config),
                                             mu = c(r_s = log(stats::median(truth$r_s)))))
}

truth_eta <- function(tr, spec) {
  e <- c(r_s = log(tr$r_s), th_ap = log(tr$th_ap), th_as = log(max(tr$th_as, 1e-12)),
         th_bs = log(max(tr$th_bs, 1e-12)), f0 = stats::qlogis(min(max(tr$f0, 1e-12), 1 - 1e-9)))
  e[spec$re]
}

# Q2 model expected under each mechanism; Q1 expected under each clone setting.
k3_expected <- function(mech, clone) {
  list(q2 = switch(mech, a_only = "Ma", b_only = "Mb", mix = "Mab", synergy = "Mabx"),
       q1 = clone == "clone")
}

k3_part <- function(cfg, task_row, d_design, log_fun = message) {
  s <- cfg$steps$K3
  seed <- 1000L * task_row$task + 7L
  p_clone <- s$clone[[task_row$clone]]
  pop <- design_population(cfg)
  if (task_row$design == "real35") {
    dz <- design_real35(d_design)
    cal <- calibrate_design(dz, task_row$mech, p_clone, pop, cfg, seed = 11L)
    co <- simulate_cohort(dz, task_row$mech, p_clone, pop, cfg, seed, cal$shift,
                          ids = sprintf("K%02d", seq_along(dz)))
    models <- c(cfg$fit$q2_models, model_id(if (task_row$mech == "synergy") "Mabx" else "Mab", "ddr"))
  } else {
    kpop <- fit_kernel_population(d_design$drivers, cfg)
    dz <- design_n350(s$n350$n, d_design, kpop, cfg, seed)
    cal <- calibrate_design(design_real35(d_design), "mix", p_clone, pop, cfg, seed = 11L)
    co <- simulate_cohort(dz, "mix", p_clone, pop, cfg, seed, cal$shift,
                          vaf_n = if (task_row$variant == "vafall") Inf else s$n350$n_vaf)
    models <- c("Mab", "MabC_ddr")
  }
  log_fun("K3 task ", task_row$task, ": ", task_row$design, " ", task_row$mech, " ", task_row$clone,
          " ", task_row$variant, " | phenotype ", paste(table(factor(co$phenotype$primary, levels = cfg$phenotype$classes)), collapse = "/"))
  pdats <- build_pdats(co, cfg)
  fits <- list()
  for (m in models) {
    fits[[m]] <- fit_nested(pdats, m, fits, cfg)
    log_fun(sprintf("  %s BIC %.1f (%.0fs)", m, fits[[m]]$BIC, fits[[m]]$runtime_s))
  }
  bic <- data.table::data.table(model = names(fits), BIC = vapply(fits, `[[`, numeric(1), "BIC"),
                                logLik = vapply(fits, `[[`, numeric(1), "logLik"))
  ex <- k3_expected(task_row$mech, task_row$clone)
  q2_sel <- if (task_row$design == "real35") bic[model %in% cfg$fit$q2_models][which.min(BIC), model] else NA_character_
  base_q1 <- if (task_row$design == "real35") (if (task_row$mech == "synergy") "Mabx" else "Mab") else "Mab"
  q1_sel <- bic[model == model_id(base_q1, "ddr"), BIC] < bic[model == base_q1, BIC]

  # per-patient attribution with the clone model (both readouts) vs truth
  fm <- fits[[model_id(base_q1, "ddr")]]
  mc <- mcmc_all(pdats, fm, cfg, seed = seed)
  pdats_s <- prep_segments(pdats, fm$pop$struct, cfg)
  tf <- truth_fit(if (task_row$design == "real35") task_row$mech else "mix", co$truth, cfg)
  pdats_t <- prep_segments(pdats, tf$pop$struct, cfg)
  att <- data.table::rbindlist(lapply(names(pdats), function(pid) {
    tab <- attribute_patient(mc[[pid]]$draws, pdats_s[[pid]], fm, cfg, n_draws = min(150, cfg$fit$attrib_draws),
                             seed = seed)
    sm <- summarise_attribution(tab, cfg)
    tr <- co$truth[patient_id == pid]
    cf <- counterfactual_draw(truth_eta(tr, tf$spec), pdats_t[[pid]], tf, cfg, 0:30)
    w_true <- when_class(cf[["actual"]], cf[["no_a"]], cf[["no_bS"]], cf[["no_abS"]])
    wt_true <- when_class(cf[["actual"]], cf[["no_a"]], cf[["no_b"]], cf[["no_ab"]])
    h_true <- who_class(cf[["actual"]], cf[["no_clone"]], cf[["reserve"]])
    sm[, `:=`(true_when = w_true, true_who = h_true, true_pheno = cf[["actual"]],
              p_true_when = get(paste0("p_when_", w_true)),
              true_when_total = wt_true, p_true_when_total = get(paste0("p_whentot_", wt_true)),
              p_true_who = if (is.na(h_true)) NA_real_ else get(paste0("p_who_", h_true)),
              carrier = tr$carrier, true_f0 = tr$f0,
              eta_r_s = fm$indiv[patient_id == pid]$eta_r_s, true_r_s = tr$r_s)]
    sm
  }), fill = TRUE)
  list(task = task_row, bic = bic, q2_expected = ex$q2, q2_selected = q2_sel,
       q1_expected = ex$q1, q1_selected = q1_sel, attribution = att,
       pop = lapply(fits, function(f) list(mu = f$pop$mu, omega = f$pop$omega, struct = f$pop$struct)),
       calib = cal$history, phenotype = table(factor(co$phenotype$primary, levels = cfg$phenotype$classes)))
}

k3_combine <- function(parts, cfg) {
  s <- cfg$steps$K3
  sel <- data.table::rbindlist(lapply(parts, function(p) data.table::data.table(
    p$task, q2_expected = p$q2_expected, q2_selected = p$q2_selected,
    q1_expected = p$q1_expected, q1_selected = p$q1_selected)))
  att <- data.table::rbindlist(lapply(parts, function(p) p$attribution[, `:=`(task = p$task$task,
    design = p$task$design, mech = p$task$mech, clone = p$task$clone, variant = p$task$variant)]), fill = TRUE)
  q2 <- sel[design == "real35", .(n = .N, q2_correct = mean(q2_selected == q2_expected)), by = .(mech, clone)]
  q1 <- sel[, .(n = .N, q1_correct = mean(q1_selected == q1_expected)), by = .(design, variant, clone)]
  indiv <- att[, .(n_patients = .N, mean_p_true_when = mean(p_true_when, na.rm = TRUE),
                   mean_p_true_when_total = mean(p_true_when_total, na.rm = TRUE),
                   mean_p_true_who = mean(p_true_who, na.rm = TRUE),
                   n_nonquick = sum(true_pheno != "quick"),
                   cor_r_s = suppressWarnings(stats::cor(eta_r_s, log(true_r_s)))),
               by = .(design, variant, mech, clone)]
  wc <- c("none", "a", "b", "ab_joint", "ab_redundant", "other")
  att[, map_when := wc[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = paste0("p_when_", wc)]
  conf_when <- att[, .N, by = .(design, mech, clone, true_when, map_when)]
  # routes (Aim 4 plan section 6), on the real design
  route_of <- function(cohort, indiv_p) {
    if (is.na(cohort)) return(NA_character_)
    if (cohort >= s$route$cohort_correct && !is.na(indiv_p) && indiv_p >= s$route$indiv_prob) "1_full"
    else if (cohort >= s$route$cohort_correct) "2_cohort_only" else "3_describe_and_measure"
  }
  r35 <- indiv[design == "real35"]
  q2_cohort <- mean(sel[design == "real35"]$q2_selected == sel[design == "real35"]$q2_expected)
  q2_indiv <- mean(att[design == "real35"]$p_true_when, na.rm = TRUE)
  q1_cohort <- mean(sel[design == "real35"]$q1_selected == sel[design == "real35"]$q1_expected)
  q1_indiv <- mean(att[design == "real35" & clone == "clone"]$p_true_who, na.rm = TRUE)
  routes <- data.table::data.table(
    question = c("Q2 when (a vs b)", "Q1 who (reserve vs clone)"),
    cohort_correct = c(q2_cohort, q1_cohort), indiv_mean_p_true = c(q2_indiv, q1_indiv),
    route = c(route_of(q2_cohort, q2_indiv), route_of(q1_cohort, q1_indiv)))
  list(selection = sel, q2 = q2, q1 = q1, indiv = indiv, conf_when = conf_when, attribution = att, routes = routes)
}

an_K3_ident <- function(ctx) {
  cfg <- ctx$config
  tasks <- k3_tasks(cfg)
  dir.create(dirname(k3_part_file(cfg, 1)), recursive = TRUE, showWarnings = FALSE)
  have <- file.exists(vapply(tasks$task, function(k) k3_part_file(cfg, k), character(1)))
  if (!all(have)) {
    ctx$log("  K3: ", sum(!have), " of ", nrow(tasks), " parts missing - computing them here ",
            "(use hpc/k3_array.sbatch for the full run)")
    dd <- ctx$design()
    for (k in tasks$task[!have]) {
      part <- k3_part(cfg, tasks[task == k], dd, log_fun = function(...) ctx$log("  ", ...))
      saveRDS(part, k3_part_file(cfg, k))
    }
  }
  parts <- lapply(tasks$task, function(k) readRDS(k3_part_file(cfg, k)))
  cmb <- k3_combine(parts, cfg)
  ctx$log("  K3 routes: ", paste(sprintf("%s -> %s (cohort %.2f, indiv %.2f)", cmb$routes$question,
                                          cmb$routes$route, cmb$routes$cohort_correct,
                                          cmb$routes$indiv_mean_p_true), collapse = "; "))
  list(tables = list(q2_selection = cmb$q2, q1_selection = cmb$q1, individual = cmb$indiv,
                     confusion_when = cmb$conf_when, routes = cmb$routes, tasks = tasks),
       summary = cmb, routes = cmb$routes, n_parts = length(parts))
}
