# =============================================================================
# an_S0_data.R
# Purpose : Data overview. Outcome-related summaries (phenotypes,
#           interventions, trajectories) use the training split only, so the
#           hold-out outcomes stay unseen; data-availability summaries use all.
# Inputs  : ctx$data() - standardised data.
# Outputs : list(tables, density, traj_train, key_numbers)
# =============================================================================

an_S0_data <- function(ctx) {
  cfg <- ctx$config
  d <- ctx$data()
  p <- d$patients
  w <- cfg$steps$P1$window
  mp <- cfg$steps$P1$min_points

  overview <- data.table::data.table(
    item = c("patients", "train", "test", "excluded (flags)", "ANC observations",
             "PLT observations", "CRP observations", "ferritin observations",
             "sequenced pre-treatment", "sequenced: CH positive", "post-treatment samples planned",
             "baseline RDW available", "baseline MCV available"),
    value = c(nrow(p), sum(p$split == "train"), sum(p$split == "test"), sum(p$excluded),
              d$obs[var == "anc", .N], d$obs[var == "plt", .N], d$obs[var == "crp", .N],
              d$obs[var == "ferritin", .N],
              if (is.null(d$vaf)) 0L else nrow(d$vaf),
              if (is.null(d$vaf)) 0L else sum(d$vaf$status == "positive"),
              if (is.null(d$post_planned)) 0L else nrow(d$post_planned),
              if ("rdw_base" %in% names(p)) sum(!is.na(p$rdw_base)) else 0L,
              if ("mcv_base" %in% names(p)) sum(!is.na(p$mcv_base)) else 0L))

  # late follow-up density: the number that decides how P1 can be tested
  late <- d$obs[keep == TRUE & var %in% c("anc", "plt") & day >= w[1] & day <= w[2],
                .N, by = .(patient_id, var)]
  late <- data.table::dcast(late, patient_id ~ var, value.var = "N", fill = 0)
  late <- merge(p[, .(patient_id, split)], late, by = "patient_id", all.x = TRUE)
  for (v in c("anc", "plt")) if (!v %in% names(late)) late[, (v) := 0L]
  late[is.na(anc), anc := 0]
  late[is.na(plt), plt := 0]
  density <- late[, .(n_patients = .N,
                      n_ge_min_points = sum(anc >= mp & plt >= mp),
                      pct = round(100 * mean(anc >= mp & plt >= mp), 1)), by = split]
  density[, window := paste0("day ", w[1], "-", w[2], ", >= ", mp, " ANC and PLT")]

  censoring <- d$obs[keep == FALSE, .N, by = .(var, reason)][order(var, -N)]

  # training split only from here on
  tr <- subset_split(d, "train")
  ph <- tr$pheno[excluded == FALSE]
  pheno_tab <- ph[, .N, by = phenotype][, pct := round(100 * N / sum(N), 1)][order(phenotype)]
  agree <- cross <- NULL
  if ("phenotype_clinical" %in% names(ph)) {
    a <- ph[!is.na(phenotype) & !is.na(phenotype_clinical)]
    agree <- data.table::data.table(n = nrow(a), agreement_pct = round(100 * mean(a$phenotype == a$phenotype_clinical), 1))
    cross <- data.table::dcast(ph[, .N, by = .(rule = phenotype, recorded = phenotype_clinical)],
                               rule ~ recorded, value.var = "N", fill = 0)
  }
  interventions <- tr$events[, .(n_events = .N, n_patients = data.table::uniqueN(patient_id)), by = event][order(-n_patients)]

  traj <- merge(tr$obs[keep == TRUE & var %in% c("anc", "plt")],
                ph[, .(patient_id, phenotype)], by = "patient_id")
  sampling <- d$obs[var == "anc", .N, by = day][order(day)]

  list(
    tables = list(overview = overview, late_density = density, censoring = censoring,
                  phenotype_train = pheno_tab, interventions_train = interventions,
                  checks = d$checks, clinical_agreement_train = agree, phenotype_vs_clinical = cross),
    traj_train = traj, sampling = sampling,
    key_numbers = list(late_density = density, phenotype_train = pheno_tab))
}
