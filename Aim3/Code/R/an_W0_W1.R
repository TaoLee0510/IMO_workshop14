# =============================================================================
# an_W0_W1.R
# Purpose : Step W0 (convert the workbook, apply the data rules, compute the
#           locked phenotype and its variants) and step W1 (unit tests and
#           typical trajectories of the model). No file I/O here.
# Inputs  : ctx (config, data loader, upstream results).
# Outputs : result lists with $tables (and $derived for W0).
# =============================================================================

an_W0_data <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  ph <- phenotype_table(d$labs, d$patients, cfg)
  d$phenotype <- merge(ph, d$patients[, .(patient_id, label_clinician, gcsf14)], by = "patient_id")
  vars <- c("primary", names(cfg$phenotype$variants))
  counts <- data.table::rbindlist(lapply(vars, function(v) {
    tb <- table(factor(d$phenotype[[v]], levels = cfg$phenotype$classes))
    data.table::data.table(rule = v, quick = tb[["quick"]], second_wave = tb[["second_wave"]],
                           delayed = tb[["delayed"]])
  }))
  locked <- unlist(counts[rule == "primary", .(quick, second_wave, delayed)])
  ctx$log("  W0: primary phenotype ", paste(names(locked), locked, sep = " ", collapse = " / "),
          if (identical(unname(locked), c(19L, 10L, 6L))) " (matches the locked 19/10/6)" else " (DIFFERS from 19/10/6)")
  list(derived = d,
       locked_ok = identical(unname(locked), c(19L, 10L, 6L)),
       tables = list(
         phenotype = d$phenotype,
         phenotype_counts = counts,
         obs_handling = d$obs[, .N, by = .(var, reason, keep)][order(var, reason)],
         label_vs_primary = d$phenotype[, .N, by = .(label_clinician, primary)][order(label_clinician, primary)],
         gcsf_by_class = d$phenotype[, .N, by = .(primary, gcsf14)][order(primary, gcsf14)],
         patients = d$patients))
}

# Design context shared by W1-W3: per-patient designs, the median driver and
# typical set points.
design_context <- function(d, config, struct = config$model$struct) {
  ids <- d$patients$patient_id
  pdes_list <- stats::setNames(lapply(ids, function(pid)
    make_pdes(d$patients[patient_id == pid], d$inflam[patient_id == pid], struct, config)), ids)
  drv0 <- lapply(ids, function(pid) {
    sig <- pdes_list[[pid]]$sig
    make_driver(sig, struct, config, t0 = config$model$t_start)
  })
  set_med <- c(anc = stats::median(d$patients$anc_set), plt = stats::median(d$patients$plt_set),
               amc = stats::median(d$patients$amc_set))
  set_low <- c(anc = unname(stats::quantile(d$patients$anc_set, 0.25)), plt = set_med[["plt"]],
               amc = set_med[["amc"]])
  list(pdes_list = pdes_list, drv = median_driver(drv0), set_med = set_med, set_low = set_low)
}

# HSC kill rate per unit of injury for a target kill fraction (eta = 0).
theta_for_loss <- function(f, pd, struct, phase = c("a", "b")) {
  phase <- match.arg(phase)
  if (f <= 0) return(0)
  integral <- if (phase == "a") sum((pd$seg[, "end"] - pd$seg[, "start"]) * pd$seg[, "dose"]) / struct[["k_e"]]
              else sum(pd$drv$v) * pd$drv$dt
  -log(1 - f) / integral
}

typical_cases <- function(pd, struct, config, loss = 0.9) {
  mu <- exp(config$model$design$mu)
  cases <- list(
    none = c(th_ap = mu[["th_ap"]], th_as = 0, th_bs = 0, r_s = mu[["r_s"]]),
    a = c(th_ap = mu[["th_ap"]], th_as = theta_for_loss(loss, pd, struct, "a"), th_bs = 0, r_s = mu[["r_s"]]),
    b = c(th_ap = mu[["th_ap"]], th_as = 0, th_bs = theta_for_loss(loss, pd, struct, "b"), r_s = mu[["r_s"]]),
    ab = c(th_ap = mu[["th_ap"]], th_as = theta_for_loss(1 - sqrt(1 - loss), pd, struct, "a"),
           th_bs = theta_for_loss(1 - sqrt(1 - loss), pd, struct, "b"), r_s = mu[["r_s"]]))
  data.table::rbindlist(lapply(names(cases), function(k) {
    par <- build_par(cases[[k]], struct, "Mab")
    tr <- sim_daily(par, pd, config, days = seq(-10, 60, by = 0.5))
    tr[, `:=`(case = k, pheno = sim_phenotype(par, pd, config))]
    tr
  }))
}

an_W1_model <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$results$W0$derived
  dc <- design_context(d, cfg)
  st <- cfg$model$struct
  pd <- design_pdes(st, cfg, dc$drv, set = dc$set_med)
  tests <- run_unit_tests(cfg, pd)
  ctx$log("  W1: ", sum(tests$passed), "/", nrow(tests), " unit tests passed")
  typical <- typical_cases(pd, st, cfg, loss = cfg$steps$W2$loss_demo)
  list(all_passed = all(tests$passed), set_med = dc$set_med, set_low = dc$set_low,
       drv = dc$drv,
       tables = list(unit_tests = tests, typical = typical,
                     driver = data.table::data.table(day = dc$drv$t0 + (seq_along(dc$drv$v) - 1) * dc$drv$dt,
                                                     D = dc$drv$v)))
}
