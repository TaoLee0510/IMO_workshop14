# =============================================================================
# an_PR_prespecify.R
# Purpose : Step PR - freeze, before any model is fitted to the real counts,
#           the phenotype labels, the models, the tests T1-T5, the attribution
#           rule, the W3 reporting routes and the landmark evaluation. The
#           files are hashed (SHA-256) by the output layer; run_all refuses to
#           regenerate them once they exist.
# Inputs  : ctx (W0 phenotype table, config).
# Outputs : $prereg = list(files, hash) and a short table.
# =============================================================================

an_PR_prespecify <- function(ctx) {
  cfg <- ctx$config
  ph <- ctx$results$W0$derived$phenotype
  code_files <- c(list.files(file.path(cfg$paths$code, "R"), full.names = TRUE),
                  file.path(cfg$paths$code, "src", "hsc_model.cpp"))
  plan <- list(
    created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    statement = paste("Written before any model is fitted to the real blood counts.",
                      "The ANC trajectories were inspected to define the phenotype.",
                      "Disclosure: while checking model adequacy (2026-11-03), the second wave was",
                      "compared with treatment flags: steroids 8/10 second wave vs 5/25 others",
                      "(Fisher OR 14, p = 0.002), G-CSF after day 14 7/10 (p = 0.015). This led to",
                      "modelling the drug-free trajectory. Ferritin-timing and dose-response relations",
                      "were not looked at."),
    endpoint = list(primary = "quick vs non-quick (second wave or delayed), locked rule, day 0-30",
                    secondary = "three classes; platelets < 50 at day 26-30",
                    interpretation = paste("The model describes the drug-free trajectory; an observed second wave is",
                                           "read as drug withdrawal unmasking a latent deficit")),
    drug_handling = list(steroid_window = cfg$data$steroid_window,
                         anc = "upper bound under steroids (days 0-21) and under G-CSF (after day 14)",
                         amc = "lower bound under steroids (days 0-21); uninformative if no differential"),
    phenotype = list(rule = cfg$phenotype$primary, variants = cfg$phenotype$variants,
                     counts = as.list(table(factor(ph$primary, levels = cfg$phenotype$classes)))),
    models = list(names = cfg$fit$models,
                  random_effects = lapply(stats::setNames(cfg$fit$models, cfg$fit$models), model_re_names),
                  struct_estimated = cfg$fit$struct_est,
                  struct_defaults = as.list(cfg$model$struct),
                  re_prior = lapply(cfg$model$re_prior, as.list),
                  method = "ITS + Laplace population fit; individual posteriors by random-walk Metropolis"),
    tests = list(
      T1 = list(role = "descriptive", what = "platelet trough day after the first recovery regressed on the patient's ferritin-peak day (platelets: no G-CSF effect)",
                expect = "slope ~0 if chemotherapy (a) dominates, ~1 if CAR T phase (b) dominates"),
      T2 = list(role = "secondary", what = "distribution of the injury time back-calculated from the AMC and PLT footprints of the fitted drug-free trajectories",
                expect = "near day -4 for a; near the ferritin peak for b"),
      T3 = list(role = "secondary", what = "dose-response: th_aS vs bridging chemotherapy; th_bS vs ferritin AUC and CRS grade",
                method = "Spearman correlation of posterior medians; Wilcoxon for binary covariates"),
      T4 = list(role = "primary", what = "model comparison M0, Ma, Mb, Mab, Mabx",
                decision = "lowest BIC; Delta BIC > 6 counted as strong; leave-one-out log predictive density as check"),
      T5 = list(role = "primary", what = "landmark prediction of non-quick recovery (primary endpoint) at days 10 and 14",
                method = "leave-one-out AUROC of the model vs CAR-HEMATOTOX; also reported in the 23 patients without G-CSF")),
    attribution = list(classes = mech_classes,
                       rule = "counterfactual removal of th_aS / th_bS; improvement = at least one phenotype level",
                       report = "P(class) per patient, expected count per class, MAP class only as a summary"),
    routes = cfg$steps$W3$route,
    landmarks = cfg$steps$W6$landmarks,
    code_md5 = as.list(stats::setNames(unname(tools::md5sum(code_files)), basename(code_files))))
  locked <- ph[, .(patient_id, primary)]
  list(tables = list(prereg_summary = data.table::data.table(
         item = c("phenotype counts", "models", "primary tests"),
         value = c(paste(names(plan$phenotype$counts), unlist(plan$phenotype$counts), collapse = " / "),
                   paste(cfg$fit$models, collapse = ", "),
                   "T4 model comparison; T5 landmark prediction of quick vs non-quick"))),
       prereg = list(files = list("PR_plan.json" = plan, "PR_phenotype_locked.csv" = locked),
                     hash = "PR_plan.sha256"))
}
