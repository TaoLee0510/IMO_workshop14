# =============================================================================
# data_prep/convert_icaht.R
# Purpose : Convert the de-identified ICAHT extract (one Excel workbook) into
#           the Aim 2 data contract (TestData/README.md) in Data/.
# Usage   : Rscript Code/data_prep/convert_icaht.R [path/to/workbook.xlsx]
#           (default: the workbook in Data/)
# Inputs  : sheets daily_labs and patient_level.
# Outputs : Data/patients.csv, labs_long.csv, inflam_long.csv, events.csv,
#           split.csv, conversion_report.md. No mutation files (none exist).
#
# Rules
#   Units        : CRP mg/dL -> mg/L (x 10); counts x10^9/L as in the source.
#   Blank ANC    : with WBC < 0.5 -> left-censored, anc_upper = WBC (the ANC
#                  cannot exceed the WBC); with WBC >= 0.5 -> missing; no CBC
#                  that day -> no row. ANC = 0 is censored at the detection
#                  limit when read (io_read.R).
#   Other blanks : PLT/Hb are blank only when the whole CBC is missing; blank
#                  CRP/ferritin = not drawn.
#   Treatments   : tocilizumab, steroids and G-CSF after day 14 are yes/no
#                  only; days from approx_treatment_events() (R/events_approx.R).
#   Transfusions : product type not recorded; every transfused day is treated
#                  as a platelet transfusion (PLT dropped for 3 days after).
#   Unknown      : prior bendamustine, LD dose (standard assumed), marrow
#                  involvement, IEC-HS -> NA / 0 as noted in the report.
#   Split        : every patient "train" (evaluation by cross-validation).
# =============================================================================

suppressPackageStartupMessages(library(data.table))
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
code_dir <- normalizePath(file.path(dirname(this_file), ".."))
root <- normalizePath(file.path(code_dir, ".."))
source(file.path(code_dir, "config.R"))
source(file.path(code_dir, "R", "events_approx.R"))

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if (length(args)) args[1] else Sys.glob(file.path(root, "Data", "*.xlsx"))[1]
out_dir <- file.path(root, "Data")
cfg <- make_config("real", root = root)

num <- function(x) suppressWarnings(as.numeric(x))
yes <- function(x) as.integer(tolower(trimws(x)) %in% c("yes", "y", "1", "true"))

labs_raw <- as.data.table(readxl::read_excel(xlsx, sheet = "daily_labs", col_types = "text"))
pat_raw <- as.data.table(readxl::read_excel(xlsx, sheet = "patient_level", col_types = "text"))

# ---- labs ----------------------------------------------------------------------------------
L <- labs_raw[, .(patient_id = study_id, day = num(day), wbc = num(wbc), anc = num(anc), plt = num(plt),
                  hgb = num(hgb), alc = num(alc), mcv = num(mcv), rdw = num(rdw),
                  crp = 10 * num(crp), ferritin = num(ferritin), transfused = yes(transfused))]
L[, anc_upper := data.table::fifelse(is.na(anc) & !is.na(wbc) & wbc < 0.5, wbc, NA_real_)]
labs <- L[!is.na(wbc) | !is.na(anc) | !is.na(plt) | !is.na(hgb),
          .(patient_id, day, anc, plt, wbc, hgb, alc, mcv, rdw, anc_upper)]
inflam <- L[!is.na(crp) | !is.na(ferritin), .(patient_id, day, crp, ferritin)]

# ---- patients ------------------------------------------------------------------------------
P <- pat_raw[, .(
  patient_id = study_id, age = num(age), sex = sex, product = product, costim = "CD28",
  prior_lines = num(prior_lines), prior_asct = yes(prior_auto_sct), prior_bendamustine = NA_integer_,
  months_since_last_chemo = NA_real_,
  bridging = data.table::fifelse(yes(bridging) == 1, tolower(bridging_type), "none"),
  ld_flu_total = NA_real_, ld_cy_total = NA_real_, ld_dose_reduced = 0L,
  crs_grade = num(max_crs), crs_onset_day = NA_real_, crs_end_day = NA_real_,
  crs_max_day = num(date_max_crs_day),
  marrow_involvement = NA_integer_, iec_hs = NA_integer_, viral_infection = NA_integer_,
  relapse_day = num(progression_date_day), subsequent_tx_day = NA_real_, death_day = num(death_date_day),
  last_followup_day = num(last_contact_date_day),
  tocilizumab = yes(tocilizumab), steroid = yes(steroids), gcsf_after_d14 = yes(need_gcsf_after_d14),
  phenotype_clinical = c(quick = "quick", intermittent = "second_wave",
                         `no recovery >1500 by d30` = "persistent")[phenotype_from_anc],
  phenotype_cohort_label = cohort_phenotype_label,
  longest_anc_lt500_days = num(longest_anc_lt500_days))]

# ---- events --------------------------------------------------------------------------------
anc_for_rule <- labs[, .(patient_id, day, anc = data.table::fifelse(is.na(anc), anc_upper, anc))]
ev <- approx_treatment_events(P[, .(patient_id, tocilizumab, steroid, gcsf_after_d14, crs_max_day)],
                              anc_for_rule, cfg$data$approx_events)
tx <- L[transfused == 1, .(patient_id, event = "plt_tx", start_day = day, end_day = day, dose = NA_real_)]
events <- rbind(ev, tx)
setorder(events, patient_id, start_day, event)

split <- P[, .(patient_id, split = "train", locked_at = "none (cross-validation)", seed = cfg$runtime$seed)]

# ---- write -----------------------------------------------------------------------------------
fw <- function(x, f) fwrite(x, file.path(out_dir, f), na = "")
fw(P, "patients.csv"); fw(labs, "labs_long.csv"); fw(inflam, "inflam_long.csv")
fw(events, "events.csv"); fw(split, "split.csv")

win <- cut(labs$day, c(-Inf, -1, 10, 21, 30, 44, 75, Inf), c("<0", "0-10", "11-21", "22-30", "31-44", "45-75", ">75"))
cover <- labs[, .(anc_values = sum(!is.na(anc)), anc_censored = sum(!is.na(anc_upper)), patients = uniqueN(patient_id)),
              by = .(window = win)][order(window)]
md_tab <- function(x) c(paste("|", paste(names(x), collapse = " | "), "|"),
                        paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
                        apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
rep_txt <- c(
  "# Conversion report (generated by Code/data_prep/convert_icaht.R)", "",
  paste0("Source: `", basename(xlsx), "` - ", nrow(P), " patients, ", nrow(labs), " CBC rows, ",
         nrow(inflam), " CRP/ferritin rows."), "",
  "ANC by window (censored = blank ANC with WBC < 0.5, kept as ANC <= WBC):", "", md_tab(cover), "",
  paste0("Blank ANC with WBC >= 0.5 (treated as missing): ", L[is.na(anc) & !is.na(wbc) & wbc >= 0.5, .N],
         "; ANC = 0 (censored at the detection limit): ", L[anc == 0, .N], "."), "",
  "Approximated treatment days (only yes/no in the source):", "",
  md_tab(events[event != "plt_tx", .N, by = event]), "",
  paste0("Transfused days (type unknown, treated as platelet transfusions): ", nrow(tx), "."), "",
  "Unknown and set to NA: prior bendamustine, months since chemotherapy, LD doses (standard dose assumed),",
  "CRS onset and end, marrow involvement, IEC-HS, viral infection, subsequent therapy.")
writeLines(rep_txt, file.path(out_dir, "conversion_report.md"))
cat(rep_txt, sep = "\n")
