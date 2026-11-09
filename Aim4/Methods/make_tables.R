# =============================================================================
# make_tables.R
# Purpose : Write the parameter tables, calibration tables and number macros
#           of methods.tex from Code/config.R and the model code, so the
#           document cannot drift from the code. Stops if a structural
#           parameter in config.R has no description here.
# Usage   : Rscript Methods/make_tables.R        (from Aim4/)
# Outputs : Methods/tables/*.tex, Methods/tables/values.tex
# =============================================================================

suppressPackageStartupMessages(library(data.table))
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(here, ".."))
source(file.path(root, "Code", "config.R"))
for (f in list.files(file.path(root, "Code", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- make_config("test", root = root)
load_aim4_model(cfg)
out_dir <- file.path(here, "tables")
dir.create(out_dir, showWarnings = FALSE)
st <- cfg$model$struct
# real derived data (v1.3): set-point SDs, G-CSF windows, ALC by day window (no fits)
cfg_real <- make_config("real", root = root)
dreal <- if (file.exists(file.path(cfg_real$paths$derived, "patients.csv"))) read_aim4_data(cfg_real) else NULL

# ---- helpers ----------------------------------------------------------------------------------
num <- function(x, d = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return("--")
    if (v == round(v) && abs(v) < 1e5) return(format(v, big.mark = ""))
    if (abs(v) >= 1e3) return(sprintf("$10^{%g}$", round(log10(v), 2)))
    format(signif(v, d), drop0trailing = TRUE, scientific = FALSE)
  }, character(1))
}
code <- function(x) paste0("\\texttt{", gsub("_", "\\\\_", x), "}")
row_tex <- function(...) paste0(paste(..., sep = " & "), " \\\\")
write_tab <- function(file, caption, label, colspec, header, rows, size = "\\small", notes = NULL) {
  lines <- c("\\begin{table}[!htbp]", "\\centering", size, "\\setlength{\\tabcolsep}{3.5pt}",
             paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
             paste0("\\begin{tabular}{", colspec, "}"), "\\toprule", header, "\\midrule", rows, "\\bottomrule",
             "\\end{tabular}",
             if (!is.null(notes)) c("\\par\\smallskip", paste0("\\parbox{\\linewidth}{\\footnotesize ", notes, "}")),
             "\\end{table}")
  writeLines(lines, file.path(out_dir, file))
}

# ---- structural parameter metadata ---------------------------------------------------------------
# role: est = estimated in F1 (population); fixed; switch; book = bookkeeping only;
#       sim = simulation / what-if only; clone = set by the clone type
est <- c(cfg$fit$struct_est$base, cfg$fit$struct_est$Mabx, cfg$fit$struct_est$route_hsc,
         cfg$fit$struct_est$route_mpp)
meta <- rbindlist(list(
  # HSC, clone, self-renewal feedback, lymphoid
  list("alpha", "$\\alpha$", "Acceleration of HSC regrowth while the MPP is depleted", "--", "v1.1", "marrow"),
  list("r_a", "$r_a$", "Clone sensitivity to chemotherapy relative to normal HSC", "--", "clone scenario", "marrow"),
  list("r_b", "$r_b$", "Clone sensitivity to the CAR T phase relative to normal HSC", "--", "clone scenario", "marrow"),
  list("delta", "$\\delta$", "Output of clonal HSC relative to normal HSC", "--", "clone scenario", "marrow"),
  list("gmax", "$G_{\\max}$", "Upper bound of the self-renewal feedback factor", "--", "v1.1", "marrow"),
  list("gamma_m", "$\\gamma_M$", "Exponent of the self-renewal feedback (shared by all compartments)", "--", "v1.1", "marrow"),
  list("fb_mpp", "$\\chi_{\\mathrm{MPP}}$", "Feedback switch, MPP (progeny: ANC, AMC, PLT)", "0/1", "design", "marrow"),
  list("fb_cmp", "$\\chi_{\\mathrm{CMP}}$", "Feedback switch, CMP (progeny: ANC, AMC, PLT)", "0/1", "design", "marrow"),
  list("fb_gmp", "$\\chi_{\\mathrm{GMP}}$", "Feedback switch, GMP (progeny: ANC, AMC)", "0/1", "design", "marrow"),
  list("fb_mep", "$\\chi_{\\mathrm{MEP}}$", "Feedback switch, MEP (progeny: PLT)", "0/1", "design", "marrow"),
  list("fb_clp", "$\\chi_{\\mathrm{CLP}}$", "Feedback switch, CLP (progeny: ALC)", "0/1", "design", "marrow"),
  list("fb_ly", "$\\chi_{L}$", "Feedback switch, lymphocyte pool (ALC)", "0/1", "design", "marrow"),
  list("th_al", "$\\theta_{aL}$", "Lymphodepletion of lymphocytes per unit exposure", "d$^{-1}$", "calibrated: real ALC", "marrow"),
  list("e_sr", "$e_{\\mathrm{sr}}$", "Suppression of MPP, CMP, GMP, MEP self-renewal by $E(I)$", "--", "prior centre", "inputs"),
  # maturation chains and blood
  list("mtt_n", "$\\mathrm{MTT}_{\\mathrm{N}}$", "Neutrophil maturation transit time (3 compartments)", "d", "lineage tree (5--7 d)", "chains"),
  list("mtt_mo", "$\\mathrm{MTT}_{\\mathrm{Mo}}$", "Monocyte maturation transit time", "d", "lineage tree (1--2 d)", "chains"),
  list("mtt_pl", "$\\mathrm{MTT}_{\\mathrm{Pl}}$", "Megakaryocyte maturation time (endomitosis)", "d", "assumption", "chains"),
  list("c_n", "$c_{\\mathrm{N}}$", "Neutrophil blood turnover", "d$^{-1}$", "\\citet{dancey1976}", "chains"),
  list("c_mo", "$c_{\\mathrm{Mo}}$", "Monocyte blood turnover", "d$^{-1}$", "v1.1", "chains"),
  list("c_pl", "$c_{\\mathrm{Pl}}$", "Platelet blood turnover", "d$^{-1}$", "\\citet{harker1969}", "chains"),
  list("c_l", "$c_{L}$", "Exchange of lymphocytes between pool and blood", "d$^{-1}$", "assumption", "chains"),
  list("gam_n", "$\\gamma_{\\mathrm{N}}$", "Production feedback at the neutrophil chain entry", "--", "v1.1", "chains"),
  list("gam_mo", "$\\gamma_{\\mathrm{Mo}}$", "Production feedback at the monocyte chain entry", "--", "v1.1", "chains"),
  list("gam_pl", "$\\gamma_{\\mathrm{Pl}}$", "Production feedback at the megakaryocyte chain entry", "--", "v1.1", "chains"),
  list("d_n", "$d_{\\mathrm{N}}$", "Divisions inside the neutrophil chain (myeloblast to myelocyte)", "--", "assumption", "chains"),
  list("d_mo", "$d_{\\mathrm{Mo}}$", "Divisions inside the monocyte chain (monoblast, promonocyte)", "--", "assumption", "chains"),
  list("d_pl", "$d_{\\mathrm{Pl}}$", "Divisions inside the megakaryocyte chain (endomitosis is not a division)", "--", "lineage tree", "chains"),
  # inputs and inflammation
  list("k_e", "$k_e$", "Elimination of the chemotherapy exposure", "d$^{-1}$", "v1.1", "inputs"),
  list("beta_renal", "$\\beta$", "Exponent of the renal scaling $(90/\\mathrm{eGFR})^{\\beta}$", "--", "v1.1", "inputs"),
  list("u_br_rgo", "$u_{\\mathrm{br,RGO}}$", "Exposure rate of R-GemOx bridging", "d$^{-1}$", "v1.1", "inputs"),
  list("u_br_cy", "$u_{\\mathrm{br,Cy}}$", "Exposure rate of other (cyclophosphamide-type) bridging", "d$^{-1}$", "v1.1", "inputs"),
  list("k_i", "$k_I$", "Rate of each stage of the two-stage marrow lag", "d$^{-1}$", "Aim 2 default", "inputs"),
  list("emax", "$E_{\\max}$", "Maximal reversible suppression by marrow inflammation", "--", "Aim 2 default", "inputs"),
  list("ic50", "$\\mathrm{IC}_{50}$", "Lagged driver giving half-maximal suppression", "--", "Aim 2 plan", "inputs"),
  list("hill", "$h_E$", "Hill coefficient of the suppression", "--", "Aim 2 default", "inputs"),
  list("e_mo", "$e_{\\mathrm{Mo}}$", "Suppression of monocyte production relative to neutrophils", "--", "assumption", "inputs"),
  list("e_pl", "$e_{\\mathrm{Pl}}$", "Suppression of megakaryocyte production relative to neutrophils", "--", "assumption", "inputs"),
  list("e_l", "$e_{L}$", "Suppression of lymphoid production relative to neutrophils", "--", "assumption", "inputs"),
  list("eta_syn", "$\\eta$", "Synergy of CAR T-phase HSC injury with prior HSC loss (Mab$\\times$ only)", "--", "prior", "inputs"),
  list("gcsf_fold", "$u_{\\mathrm{G}}$", "Fold increase of neutrophil chain entry at full G-CSF effect ($E_G=1$)", "--", "v1.1 assumption (prior centre)", "inputs"),
  list("gcsf_rel", "$u_{\\mathrm{R}}$", "Acceleration of the release of mature neutrophils ($T_{\\mathrm{N}3}\\to B_{\\mathrm{N}}$) at full G-CSF effect", "--", "assumption (prior centre)", "inputs"),
  list("k_g", "$k_G$", "Rate at which the G-CSF effect $E_G$ follows the dosing indicator", "d$^{-1}$", "assumption", "inputs"),
  list("tpo_fold", "$u_{\\mathrm{T}}$", "Fold increase of megakaryocyte chain entry under TPO-RA", "--", "assumption", "inputs"),
  # CAR T - target interaction (v1.3)
  list("car_t0", "$C_{T,0}$", "CAR T cells in blood right after the infusion (day 0)", "$10^9$/L", "assumption ($\\approx1$ cell/$\\mu$L)", "cart"),
  list("k_ag", "$A_{50}$", "Antigen load giving half-maximal occupancy (expansion and killing; relative load)", "--", "assumption", "cart"),
  list("kappa", "$k_A$", "Killing of CD19$^+$ targets per CAR T cell", "(10$^9$/L)$^{-1}$ d$^{-1}$", "calibrated (K1): typical peak on day 8", "cart"),
  list("k_t", "$k_T$", "Contraction rate of CAR T cells in blood", "d$^{-1}$", "assumption (prior centre)", "cart"),
  list("beta_a", "$\\beta_A$", "Initial antigen load per unit centred log LDH, $A_i(0)=e^{\\beta_A\\,\\mathrm{ldh}_{c,i}}$", "--", "assumption (prior centre)", "cart"),
  list("k_z", "$k_Z$", "Rate of each stage of the two-stage lag of the interaction signal $Z$", "d$^{-1}$", "assumption", "cart"),
  list("z50", "$Z_{50}$", "Signal $Z$ giving half-maximal suppression of MPP self-renewal", "$10^9$/L", "calibrated (K1): typical peak of $Z$", "cart"),
  list("th_zs", "$\\theta_{zS}$", "HSC injury per unit signal $Z$ (CAR T--target route of ``b'')", "(10$^9$/L)$^{-1}$ d$^{-1}$", "centre calibrated (K1)", "cart"),
  list("e_z", "$e_Z$", "Maximal suppression of MPP self-renewal by $Z$", "--", "assumption (prior centre)", "cart")
), use.names = FALSE)
setnames(meta, c("code", "symbol", "meaning", "unit", "source", "group"))
renew <- c("mpp", "cmp", "gmp", "mep", "clp", "ly")
covered <- c(meta$code, paste0("k_", renew), paste0("sr_", renew))
missing_meta <- setdiff(names(st), covered)
extra_meta <- setdiff(covered, names(st))
if (length(missing_meta) || length(extra_meta))
  stop("parameter metadata out of date - missing: ", paste(missing_meta, collapse = ", "),
       "; unknown: ", paste(extra_meta, collapse = ", "))
role_of <- function(cd) {
  if (cd == "eta_syn") return("estimated (Mab$\\times$)")
  if (cd == "th_zs") return("estimated (b models; not \\_Zm)")
  if (cd == "e_z") return("estimated (not \\_Zh)")
  if (cd %in% est) return("estimated")
  if (startsWith(cd, "fb_")) return("switch")
  if (startsWith(cd, "d_")) return("bookkeeping")
  if (cd %in% c("tpo_fold")) return("simulation")
  if (cd %in% c("r_a", "r_b", "delta")) return("clone type")
  "fixed"
}
meta[, value := num(unname(st[code]), 4)]
meta[code == "eta_syn", value := sprintf("0 ($e^{%s}$ in Mab$\\times$)", num(cfg$fit$eta_syn_prior[["centre"]]))]
meta[code %in% c("r_a", "r_b", "delta"), value := "1 (no clone)"]
meta[, role := vapply(code, role_of, character(1))]
hdr <- "Symbol & Code & Meaning & Unit & Value & Role & Source \\\\"
spec <- "ll>{\\raggedright\\arraybackslash}p{4.0cm}>{\\raggedright\\arraybackslash}p{1.5cm}>{\\raggedright\\arraybackslash}p{1.2cm}>{\\raggedright\\arraybackslash}p{1.9cm}>{\\raggedright\\arraybackslash}p{2.5cm}"
for (g in c("marrow", "chains", "inputs", "cart")) {
  m <- meta[group == g]
  cap <- switch(g,
    marrow = "Structural parameters of the stem cell, progenitor feedback and lymphoid layers (self-renewal and exit rates are in Table~\\ref{tab:selfrenewal}). Values are prior centres for estimated parameters.",
    chains = "Structural parameters of the maturation chains and blood. $\\mathrm{MTT}_\\ell$ is converted to $\\kappa_\\ell = 3/\\mathrm{MTT}_\\ell$.",
    inputs = "Structural parameters of the inputs: chemotherapy exposure, ferritin-driven CAR T-phase inflammation and supportive care (G-CSF is an input of the fit from v1.3; TPO-RA only in what-if analyses). Values are prior centres for estimated parameters.",
    cart = "Structural parameters of the CAR T--target interaction module (v1.3). Values are prior centres for estimated parameters; the antigen load is relative (median patient $A(0)=1$).")
  write_tab(paste0("tab_struct_", g, ".tex"), cap, paste0("tab:struct_", g), spec, hdr,
            m[, row_tex(symbol, code(code), meaning, unit, value, role, source)], size = "\\footnotesize")
}

# ---- self-renewal ---------------------------------------------------------------------------------
sig <- c(mpp = "ANC, AMC, PLT", cmp = "ANC, AMC, PLT", gmp = "ANC, AMC", mep = "PLT", clp = "ALC", ly = "ALC")
sym <- c(mpp = "\\mathrm{MPP}", cmp = "\\mathrm{CMP}", gmp = "\\mathrm{GMP}", mep = "\\mathrm{MEP}", clp = "\\mathrm{CLP}", ly = "L")
lab <- c(mpp = "MPP", cmp = "CMP", gmp = "GMP", mep = "MEP", clp = "CLP", ly = "Lymphocytes ($L$)")
src <- c(mpp = "v1.1 ($k_M$, $\\phi k_M$)", cmp = "calibrated (K1)", gmp = "calibrated (K1)", mep = "calibrated (K1)",
         clp = "assumption", ly = "assumption")
rows <- c(row_tex("HSC", "$r_S$ (\\texttt{r\\_s})", "--", num(exp(cfg$model$re_prior$mu[["r_s"]])), "--", "--",
                  "niche $(1-S)$; MPP depletion", "--", "random effect"),
          vapply(renew, function(x) {
            k <- st[[paste0("k_", x)]]; s <- st[[paste0("sr_", x)]]
            row_tex(lab[[x]], sprintf("$\\sigma_{%s}$ (%s)", sym[[x]], code(paste0("sr_", x))),
                    num(k), num(s), num(s / k, 2), num(1 / (k - s)), sig[[x]],
                    if (x %in% c("mpp", "cmp", "gmp", "mep")) "yes" else "no", src[[x]])
          }, character(1)))
write_tab("tab_selfrenewal.tex",
          "The seven self-renewing cell types of the lineage tree. $k_X$: exit rate (differentiation and death, d$^{-1}$); $\\sigma_X$: self-renewal rate (d$^{-1}$, $\\sigma_X<k_X$); $\\sigma_X/k_X$: share of the compartment's production that comes from self-renewal at steady state; $1/(k_X-\\sigma_X)$: effective residence time of a cell and its self-renewed progeny (d). The feedback factor $G_X$ senses the blood counts of the compartment's own progeny; $e_{\\mathrm{sr}}$: whether inflammation suppresses the self-renewal. All rates are fixed in the primary analysis; sensitivity analyses scale every $\\sigma_X$ by 0.5 and 1.5 (capped at $0.95\\,k_X$) and switch off all feedback except the MPP's.",
          "tab:selfrenewal", "llllll>{\\raggedright\\arraybackslash}p{2.1cm}c>{\\raggedright\\arraybackslash}p{2.3cm}",
          "Cell & Self-renewal rate & $k_X$ & $\\sigma_X$ & $\\sigma_X/k_X$ & $1/(k_X-\\sigma_X)$ & Feedback signal & $e_{\\mathrm{sr}}$ & Source \\\\",
          rows, size = "\\footnotesize")

# ---- random effects -------------------------------------------------------------------------------
re <- cfg$model$re_prior
cv <- cfg$model$covariates
cvlab <- c(prior_lines_c = "prior lines", prior_asct = "prior ASCT", age_c = "age", log_rdw_c = "log RDW", mcv_c = "MCV")
re_rows <- rbindlist(list(
  list("r_s", "$r_S$", "HSC self-renewal (regrowth) rate; reserve, Q1", "log", exp(re$mu[["r_s"]])),
  list("f0", "$f_0$", "Pre-treatment clone fraction of the HSC pool, Q1", "logit", plogis(re$mu[["f0"]])),
  list("th_as", "$\\theta_{aS}$", "HSC injury per unit chemotherapy exposure (a), Q2", "log", exp(re$mu[["th_as"]])),
  list("th_ap", "$\\theta_{aP}$", "Progenitor injury per unit exposure (MPP, CMP, GMP, MEP, CLP)", "log", exp(re$mu[["th_ap"]])),
  list("th_bs", "$\\theta_{bS}$", "HSC injury per unit lagged ferritin driver (ferritin route of b), Q2", "log", exp(re$mu[["th_bs"]])),
  list("rho", "$\\rho_T$", "CAR T antigen-driven expansion rate (all parent models)", "log", exp(re$mu[["rho"]])),
  list("slope", "$\\mathrm{Slope}$", "MF: drug effect per unit exposure", "log", exp(re$mu[["slope"]])),
  list("mtt_scale", "$s_{\\mathrm{MTT}}$", "MF: individual scaling of all MTTs", "log", exp(re$mu[["mtt_scale"]]))
), use.names = FALSE)
setnames(re_rows, c("code", "symbol", "meaning", "scale", "centre"))
re_rows[, omega := num(re$omega[code])]
re_rows[, covs := vapply(code, function(k) if (length(cv[[k]])) paste(cvlab[cv[[k]]], collapse = ", ") else "--", character(1))]
write_tab("tab_re.tex",
          "Individual parameters (random effects). $\\psi_{ik}=g_k^{-1}(\\eta_{ik})$ with $\\eta_{ik}=\\mu_k+\\boldsymbol\\beta_k^\\top\\mathbf{x}_{ik}+\\omega_k z_{ik}$, $z_{ik}\\sim N(0,1)$. Centre: $g_k^{-1}(\\mu_k)$ at the start of estimation; $\\omega_k$: prior spread on the working scale (the estimate is shrunk towards it, Section~\\ref{sec:its}). The $\\theta_{aP}$ and $\\rho_T$ centres are calibrated (Section~\\ref{sec:calibration}). Each patient also has one baseline deviation per lineage (Section~\\ref{sec:baseline}).",
          "tab:re", "ll>{\\raggedright\\arraybackslash}p{5.4cm}lll>{\\raggedright\\arraybackslash}p{3.2cm}",
          "Symbol & Code & Meaning & Scale $g_k$ & Centre & $\\omega_k$ & Covariates $\\mathbf{x}$ \\\\",
          re_rows[, row_tex(symbol, code(code), meaning, scale, num(centre), omega, covs)], size = "\\footnotesize")

# ---- clone types ------------------------------------------------------------------------------------
ct <- cfg$model$clone_types
write_tab("tab_clone.tex",
          "Clone types (not estimable on the current data; each fitted as a fixed scenario).",
          "tab:clone", "lllll",
          "Clone type & Example genes & $r_a$ & $r_b$ & $\\delta$ \\\\",
          c(row_tex("DDR (chemotherapy-tolerant, low output)", "PPM1D, TP53", num(ct$ddr[["r_a"]]), num(ct$ddr[["r_b"]]), num(ct$ddr[["delta"]])),
            row_tex("DTA (inflammation-tolerant)", "DNMT3A, TET2", num(ct$dta[["r_a"]]), num(ct$dta[["r_b"]]), num(ct$dta[["delta"]]))))

# ---- MF ---------------------------------------------------------------------------------------------
ms <- cfg$model$mf_struct
mf_meta <- rbindlist(list(
  list("mtt_n", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{N}}$", "Neutrophil mean transit time ($=4/k_{\\mathrm{tr}}$)", "d"),
  list("mtt_pl", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{Pl}}$", "Platelet mean transit time", "d"),
  list("mtt_mo", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{Mo}}$", "Monocyte mean transit time", "d"),
  list("gam_mf", "$\\gamma^{\\mathrm{MF}}$", "Feedback exponent (all lineages)", "--"),
  list("slope_pl", "$\\rho_{\\mathrm{Pl}}$", "Drug effect on platelets relative to neutrophils", "--"),
  list("slope_mo", "$\\rho_{\\mathrm{Mo}}$", "Drug effect on monocytes relative to neutrophils", "--"),
  list("k_e", "$k_e$", "Exposure elimination (as parent model)", "d$^{-1}$"),
  list("beta_renal", "$\\beta$", "Renal scaling exponent (as parent model)", "--"),
  list("u_br_rgo", "$u_{\\mathrm{br,RGO}}$", "R-GemOx bridging exposure rate", "d$^{-1}$"),
  list("u_br_cy", "$u_{\\mathrm{br,Cy}}$", "Other bridging exposure rate", "d$^{-1}$")), use.names = FALSE)
setnames(mf_meta, c("code", "symbol", "meaning", "unit"))
stopifnot(setequal(mf_meta$code, names(ms)))
write_tab("tab_mf.tex", "Population parameters of the Friberg reference model MF (prior centres).", "tab:mf", "ll>{\\raggedright\\arraybackslash}p{6cm}lll",
          "Symbol & Code & Meaning & Unit & Value & Role \\\\",
          mf_meta[, row_tex(symbol, code(code), meaning, unit, num(ms[code]),
                            ifelse(code %in% cfg$fit$struct_est$MF, "estimated", "fixed"))], size = "\\footnotesize")

# ---- steady-state time and divisions since MPP entry -----------------------------------------------
ss <- steady_age_generation(st)
bk <- cfg$model$bookkeeping
ss_rows <- ss[, row_tex(compartment, num(age_days), num(generations))]
write_tab("tab_steady.tex",
          sprintf("Pre-treatment mean time (days) and mean number of divisions since MPP entry, from the closed forms of Section~\\ref{sec:moments} with the default rates. Platelets carry the generation of their megakaryocyte, which in addition undergoes about %d endomitotic rounds (%d--%d; ploidy %d--%dN) and sheds about $10^{%d}$ platelets; neither enters the normalised equations.",
                  bk$endomitosis_rounds, bk$endomitosis_range[1], bk$endomitosis_range[2], bk$ploidy_range[1], bk$ploidy_range[2],
                  round(log10(bk$platelets_per_megakaryocyte))),
          "tab:steady", "lrr", "Compartment & Time since MPP entry (d) & Divisions since MPP entry \\\\", ss_rows)

# ---- calibration and stability (model runs) --------------------------------------------------------
v11 <- v11_comparison(cfg)
write_tab("tab_v11.tex",
          "Typical-patient early nadir and recovery (no CAR T phase, no HSC injury, typical design): the v1.1 engine at its prior centre $\\theta_{aP}=1$ and the present model at the calibrated centre. Recovery: first day after the nadir with the relative count $\\geq 0.8$.",
          "tab:v11", "llrrr", "Engine & Lineage & Nadir day & Nadir (relative) & Recovery day \\\\",
          v11[, row_tex(gsub("th_aP", "$\\\\theta_{aP}$", engine), lineage, num(nadir_day), num(nadir_rel, 2), num(recovery_day))])
scr <- progenitor_screen(cfg)
write_tab("tab_screen.tex",
          "Screen of the shared self-renewal rate of CMP, GMP and MEP ($k=2$ d$^{-1}$). For each value $\\theta_{aP}$ is recalibrated; the counts are the phenotypes produced by a single chemotherapy hit over 48 combinations of $r_S\\in\\{0.01,0.02,0.04,0.08\\}$, $\\theta_{aS}\\in\\{0.4,1,1.5,2,2.5,3.5\\}$ and $\\gamma_M\\in\\{0.5,1\\}$.",
          "tab:screen", "rrrrrrr",
          "$\\sigma$ & $\\theta_{aP}$ centre & ANC nadir day & ANC recovery day & Quick & Second wave & Delayed \\\\",
          scr[, row_tex(num(sr_progenitor), num(th_ap_centre), num(anc_nadir_day), num(anc_recovery_day),
                        n_quick, n_second_wave, n_delayed)])
stab <- stability_grid(cfg)
write_tab("tab_stability.tex",
          "Linear stability of the pre-treatment steady state: largest real part of the Jacobian eigenvalues (the neutral clone-fraction direction removed) and the period of the least-damped oscillation, over the feedback exponent $\\gamma_M$ and a common scaling of every self-renewal rate (capped at $0.95\\,k_X$).",
          "tab:stability", "rrrr", "$\\gamma_M$ & $\\sigma$ scaling & $\\max\\,\\mathrm{Re}\\,\\lambda$ (d$^{-1}$) & Period (d) \\\\",
          stab[, row_tex(num(gamma_m), num(sr_scale), num(max_re, 2), num(period_days, 2))])

# ---- v1.3: physiological search box of the individual parameters ------------------------------------
lim <- cfg$model$indiv_limits
box_rows <- rbindlist(lapply(c("r_s", "th_ap", "th_as", "th_bs", "f0", "rho", "slope", "mtt_scale"), function(k) {
  cen <- if (k == "f0") plogis(re$mu[[k]]) else exp(re$mu[[k]])
  data.table(code = k, scale = if (k == "f0") "logit" else "log", lo = lim[[k]][1], hi = lim[[k]][2], centre = cen)
}))
sym_box <- c(r_s = "$r_S$", th_ap = "$\\theta_{aP}$", th_as = "$\\theta_{aS}$", th_bs = "$\\theta_{bS}$",
             f0 = "$f_0$", rho = "$\\rho_T$", slope = "Slope (MF)", mtt_scale = "$s_{\\mathrm{MTT}}$ (MF)")
bl <- cfg$model$baseline
write_tab("tab_box.tex",
          sprintf("Search box of the individual parameters (natural scale). The box is the support of each patient's parameters: the prior is truncated to it, the differential-evolution search, the importance sampling and the MCMC stay inside it, and fitted values on a bound are reported. Baseline deviations $d_{iv}$ have the box $\\pm\\max(%s,\\ %s\\,s_{iv})$, with $s_{iv}$ the prior SD of Eq.~\\ref{eq:baseprior}.",
                  num(bl$bound_min), num(bl$bound_sd)),
          "tab:box", "llrrr", "Parameter & Scale & Lower & Upper & Prior centre \\\\",
          box_rows[, row_tex(sym_box[code], scale, num(lo), num(hi), num(centre))])

# ---- v1.3: DE-ITS controls -----------------------------------------------------------------------------
de <- cfg$fit$de; ou <- cfg$fit$outer; is_ <- cfg$fit$is; mc0 <- cfg$fit$mcmc
cores_of <- function(f) {
  x <- grep("^#SBATCH -c ", readLines(file.path(root, "Code", "hpc", f)), value = TRUE)
  if (length(x)) sub("^#SBATCH -c ([0-9]+).*$", "\\1", x[1]) else "--"
}
e_cores <- cores_of("de_estep.sbatch"); m_cores <- cores_of("de_mstep.sbatch")
sci <- function(x) sprintf("$10^{%d}$", round(log10(x)))
de_rows <- c(
  row_tex("E step", "DEoptim strategy", sprintf("%d (DE/local-to-best/1/bin)", de$strategy)),
  row_tex("", "Population size NP", num(de$np)),
  row_tex("", "Differential weight $F$ / crossover CR", paste(num(de$F), "/", num(de$CR))),
  row_tex("", "Generations per chunk / convergence window", paste(num(de$chunk), "/", num(de$window))),
  row_tex("", "Best-value improvement over the window", paste0("$<$ ", sci(de$tol_best), "$\\,(|f|+1)$")),
  row_tex("", "Population objective spread (90\\% quantile $-$ min)", paste0("$<$ ", sci(de$tol_spread))),
  row_tex("", "Safety cap (generations)", num(de$max_gen)),
  row_tex("", "Polish", sprintf("L-BFGS-B, $\\le%s$ iterations", num(de$polish_maxit))),
  row_tex("", "Cores per patient (one Slurm array task each)", e_cores),
  row_tex("M step", "Structural NP", sprintf("$\\max(%s,\\ %s\\times$ dimension)", num(de$struct_np_min), num(de$struct_np_per_dim))),
  row_tex("", "Population objective spread", paste0("$<$ ", sci(de$tol_spread), " per patient")),
  row_tex("", "Polish", "L-BFGS-B, $\\le100$ iterations"),
  row_tex("", "$\\omega$ shrinkage $\\nu$ / cap", paste(num(cfg$fit$omega_nu), "/", num(cfg$fit$omega_max))),
  row_tex("", "Cores", m_cores),
  row_tex("Outer loop", "Max.\\ change of the population vector", paste0("$<$ ", num(ou$tol_par))),
  row_tex("", "Relative change of $\\sum_i\\ell_i$", paste0("$<$ ", sci(ou$tol_obj))),
  row_tex("", "Consecutive rounds required", num(ou$n_ok)),
  row_tex("", "Damping factor / minimum", paste(num(ou$damp), "/", num(ou$damp_min))),
  row_tex("", "Over-relaxation: after / growth / maximum", paste0(num(ou$relax_after), " rounds (change $>$ ",
                                                                  num(ou$relax_min), ") / ", num(ou$relax_grow), " / ", num(ou$lam_max))),
  row_tex("", "Over-relaxation off after a rise of $\\sum_i\\ell_i$ / max.\\ extrapolation", paste(num(ou$reset_rise), "/", num(ou$max_extra))),
  row_tex("", "Stationarity window / material drift", paste0(num(ou$window), " rounds / ", num(ou$tol_drift), " (and $>2$ SE)")),
  row_tex("", "Averaging phase: change of the mean", paste0("$<$ ", num(ou$tol_par), " on ", num(ou$n_ok), " rounds")),
  row_tex("", "Safety cap (rounds)", num(ou$max_iter)),
  row_tex("Marginal likelihood", "Importance draws, stage 1 / stage 2", paste(num(is_$n), "/", num(2 * is_$n))),
  row_tex("", "Proposal: multivariate $t$, df / covariance scale", paste(num(is_$df), "/", num(is_$scale))),
  row_tex("", "Minimum ESS for the IS value (else Laplace)", num(is_$ess_min)),
  row_tex("MCMC (F2)", "Chains / iterations / burn-in / thinning", paste(num(mc0$n_chains), "/", num(mc0$n_iter), "/", num(mc0$burn), "/", num(mc0$thin))),
  row_tex("", "Extension until $\\max\\hat R<$", sprintf("%s (by %s, up to %s iterations)", num(mc0$rhat_tol), num(mc0$extend), num(mc0$max_iter))))
write_tab("tab_de.tex", "Controls of the DE-ITS estimation (v1.3), from the configuration file and the cluster scripts.",
          "tab:de", "ll>{\\raggedright\\arraybackslash}p{4.6cm}", "Stage & Control & Value \\\\", de_rows, size = "\\footnotesize")

# ---- v1.3: baselines ---------------------------------------------------------------------------------------
lab_v <- c(anc = "ANC", plt = "PLT", amc = "AMC", alc = "ALC")
if (!is.null(dreal)) {
  bsd <- dreal$base_sd
  nb <- vapply(bl$vars, function(v) { x <- dreal$patients[[paste0("n_base_", v)]]
    sprintf("%s (%s--%s)", num(stats::median(x)), num(min(x)), num(max(x))) }, character(1))
  sp_rows <- vapply(bl$vars, function(v) nrow(dreal$obs[var == v & setpoint == TRUE]), numeric(1))
} else { bsd <- stats::setNames(rep(NA, 4), bl$vars); nb <- stats::setNames(rep("--", 4), bl$vars); sp_rows <- bsd }
write_tab("tab_baseline.tex",
          sprintf("Individual baselines. $s_v$: pooled within-patient SD of the log set-point values (data, fixed); $n_{iv}$: number of set-point values of a patient (median, range); these rows enter only through the prior. The set-point drift $\\tau_v$ starts at %s, is shrunk with $\\nu=%s$ towards it and is bounded below by %s.",
                  num(bl$tau_init), num(cfg$fit$omega_nu), num(bl$tau_min)),
          "tab:baseline", "lrrr", "Lineage & $s_v$ & $n_{iv}$ & Set-point rows \\\\",
          vapply(bl$vars, function(v) row_tex(lab_v[[v]], num(bsd[[v]]), nb[[v]], num(sp_rows[[v]])), character(1)))

# ---- v1.3: ALC by day window (data) vs the calibrated typical patient --------------------------------------
cart_m <- typical_cart_run(cfg)
cart_pk <- which.max(cart_m[, "cart"])
cz <- calibrate_z(cfg)
if (!is.null(dreal)) {
  a <- dreal$obs[var == "alc" & keep == TRUE][dreal$patients[, .(patient_id, alc_base)], on = "patient_id"]
  a[, rel := value / alc_base]
  brk <- c(-6, -1, 2, 5, 8, 11, 14, 21, 30, 60)
  a <- a[day > -6 & day <= 60]
  a[, w := cut(day, brk)]
  med_base <- stats::median(dreal$patients$alc_base)
  tt <- cart_m[, "t"]
  alc_rows <- a[, .(n = .N, obs = stats::median(rel)), by = w][order(w)]
  alc_rows[, c("lym", "car") := {
    lo <- brk[as.integer(w)]; hi <- brk[as.integer(w) + 1]
    sel <- tt > lo & tt <= hi
    list(mean(cart_m[sel, "alc_rel"]), mean(cart_m[sel, "cart"]) / med_base)
  }, by = w]
  write_tab("tab_alc.tex",
            sprintf("Observed ALC relative to the patient's set point (median over observations, 35 patients) and the calibrated typical patient: lymphocytes $B_L$ (relative) and CAR T cells $C_T$ divided by the median ALC set point (%s $\\times10^9$/L), averaged over the window.", num(med_base, 2)),
            "tab:alc", "lrrrrr", "Days & $n$ & Observed & Model $B_L$ & Model $C_T$ & Model total \\\\",
            alc_rows[, row_tex(gsub(",", "..", gsub("[][()]", "", as.character(w))), n, num(obs, 2), num(lym, 2), num(car, 2), num(lym + car, 2))])
}

# ---- number macros used in the text -----------------------------------------------------------------
d <- cfg$data; ph <- cfg$phenotype$primary; fit <- cfg$fit; mc <- fit$mcmc; mo <- cfg$model
mac <- c(
  valThapCentre = num(exp(re$mu[["th_ap"]])), valThal = num(st[["th_al"]]),
  valSrProg = num(st[["sr_gmp"]]), valKProg = num(st[["k_gmp"]]), valEsr = num(st[["e_sr"]]),
  valEsrSd = num(fit$e_sr_prior_sd), valStructSd = num(mo$struct_prior_sd),
  valEtaCentre = num(fit$eta_syn_prior[["centre"]]), valEtaSd = num(fit$eta_syn_prior[["sd"]]),
  valAlcTarget = num(calib_targets$alc_rel_ld), valVOneNadir = num(v11[grepl("v1.1", engine) & lineage == "ANC", nadir_rel], 3),
  valRtol = sprintf("10^{%d}", round(log10(mo$solver$rtol))), valAtol = sprintf("10^{%d}", round(log10(mo$solver$atol))),
  valHmax = num(mo$solver$hmax), valDrvDt = num(mo$drv_dt), valDrvEnd = num(mo$drv_end), valFerDecay = num(mo$fer_decay),
  valKernelH = num(mo$kernel$h), valTmin = num(mo$t_min),
  valSetFrom = num(d$setpoint_window[1]), valSetTo = num(d$setpoint_window[2]),
  valWbcNoDiff = num(d$wbc_no_diff), valManualCells = num(d$manual_diff_cells), valGcsfDay = num(d$gcsf_from_day),
  valPltBelow = num(d$plt_tx$plt_below), valPltRise = num(d$plt_tx$rise), valHgbRise = num(d$plt_tx$hgb_rise_max),
  valPltDrop = num(d$plt_tx$drop_days), valObsMax = num(d$obs_max_day),
  valBrFrom = num(d$bridging_window[1]), valBrTo = num(d$bridging_window[2]),
  valFerFrom = num(d$fer_smooth_window[1] - 3), valFerTo = num(d$fer_smooth_window[2]),
  valPhWin = num(ph$window_end), valNadFrom = num(ph$nadir_window[1]), valNadTo = num(ph$nadir_window[2]),
  valRec = num(ph$rec_thr), valDip = num(ph$dip_thr), valCensVal = num(ph$cens_value),
  valPersFrom = num(cfg$phenotype$persistent$window[1]), valPersTo = num(cfg$phenotype$persistent$window[2]),
  valPersThr = num(cfg$phenotype$persistent$thr),
  valItsMax = num(fit$its_max_iter), valItsTol = num(fit$its_tol), valStructEvery = num(fit$struct_every),
  valOmegaNu = num(fit$omega_nu), valOmegaMax = num(fit$omega_max), valMapMaxit = num(fit$map_maxit),
  valChains = num(mc$n_chains), valIter = num(mc$n_iter), valBurn = num(mc$burn), valThin = num(mc$thin),
  valScale = num(mc$scale), valAttrib = num(fit$attrib_draws),
  valLandmarks = paste(cfg$steps$F3$landmarks, collapse = ", "), valFDraws = num(cfg$steps$F3$n_draws),
  valFKernel = num(cfg$steps$F3$n_kernel),
  valGcsfFold = num(st[["gcsf_fold"]]), valTpoFold = num(st[["tpo_fold"]]),
  valNStates = num(length(a4_state_names())), valNTrack = 19, valNPar = num(length(a4_param_names())),
  valRVersion = paste(R.version$major, R.version$minor, sep = "."),
  # v1.3
  valVersion = cfg$version, valDEoptim = as.character(utils::packageVersion("DEoptim")),
  valCarTzero = num(st[["car_t0"]]), valKag = num(st[["k_ag"]]), valKappa = num(st[["kappa"]], 4),
  valKt = num(st[["k_t"]]), valBetaA = num(st[["beta_a"]]), valKz = num(st[["k_z"]]), valZfifty = num(st[["z50"]]),
  valThzs = num(st[["th_zs"]]), valEz = num(st[["e_z"]]), valKg = num(st[["k_g"]]), valGcsfRel = num(st[["gcsf_rel"]]),
  valRhoCentre = num(exp(re$mu[["rho"]])), valRhoOmega = num(re$omega[["rho"]]),
  valCartPeakDay = num(calib_targets$cart_peak_day), valCartPeak = num(calib_targets$cart_peak),
  valKernelA = num(calib_targets$kernel[["A"]]), valKernelTp = num(calib_targets$kernel[["tp"]]),
  valCartPeakDayModel = num(unname(cart_m[cart_pk, "t"])), valCartPeakModel = num(unname(cart_m[cart_pk, "cart"]), 3),
  valZfiftyModel = num(cz[["z50"]], 3), valThzsModel = num(cz[["th_zs"]], 3),
  valAucI = num(cz[["auc_i"]], 3), valAucZ = num(cz[["auc_z"]], 3),
  valRouteLoss = sprintf("%.1f", 100 * cz[["hsc_loss_route"]]),
  valThzsSd = num(cfg$model$struct_prior_sd_k[["th_zs"]]), valEzSd = num(cfg$model$struct_prior_sd_k[["e_z"]]),
  valGcsfFoldLo = num(cfg$model$struct_limits$gcsf_fold[1]), valGcsfFoldHi = num(cfg$model$struct_limits$gcsf_fold[2]),
  valTauInit = num(bl$tau_init), valTauMin = num(bl$tau_min), valBoundSd = num(bl$bound_sd), valBoundMin = num(bl$bound_min),
  valGcsfFrom = num(d$gcsf_rule$from), valGcsfThr = num(d$gcsf_rule$thr), valGcsfDays = num(d$gcsf_rule$days),
  valGcsfUntil = num(d$gcsf_rule$until),
  valGcsfWin = if (!is.null(dreal)) num(nrow(dreal$gcsf)) else "--",
  valGcsfPts = if (!is.null(dreal)) num(data.table::uniqueN(dreal$gcsf$patient_id)) else "--",
  valGcsfFlagged = if (!is.null(dreal)) num(sum(dreal$patients$gcsf_after_d14 == 1)) else "--",
  valGcsfNoWin = if (!is.null(dreal)) num(length(setdiff(dreal$patients[gcsf_after_d14 == 1]$patient_id, dreal$gcsf$patient_id))) else "--",
  valDeNp = num(de$np), valDeF = num(de$F), valDeCR = num(de$CR), valDeStrategy = num(de$strategy),
  valDeChunk = num(de$chunk), valDeWindow = num(de$window), valDeTolBest = sprintf("10^{%d}", round(log10(de$tol_best))),
  valDeTolSpread = sprintf("10^{%d}", round(log10(de$tol_spread))), valDeMaxGen = num(de$max_gen),
  valDeStructNpMin = num(de$struct_np_min), valDeStructNpDim = num(de$struct_np_per_dim), valDePolish = num(de$polish_maxit),
  valOutTolPar = num(ou$tol_par), valOutTolObj = sprintf("10^{%d}", round(log10(ou$tol_obj))), valOutNok = num(ou$n_ok),
  valOutDamp = num(ou$damp), valOutDampMin = num(ou$damp_min), valOutMax = num(ou$max_iter),
  valOutRelaxMin = num(ou$relax_min), valOutRelaxAfter = num(ou$relax_after), valOutRelaxGrow = num(ou$relax_grow),
  valOutLamMax = num(ou$lam_max), valOutWindow = num(ou$window), valOutTolDrift = num(ou$tol_drift),
  valOutResetRise = num(ou$reset_rise), valOutMaxExtra = num(ou$max_extra),
  valIsN = num(is_$n), valIsNtwo = num(2 * is_$n), valIsDf = num(is_$df), valIsScale = num(is_$scale), valIsEss = num(is_$ess_min),
  valRhatTol = num(mc0$rhat_tol), valMcmcExtend = num(mc0$extend), valMcmcMax = num(mc0$max_iter),
  valEcores = e_cores, valMcores = m_cores,
  valRcohorts = num(cfg$steps$R$n_cohorts), valRmodel = gsub("_", "\\\\_", cfg$steps$R$model),
  valRtruth = paste(sprintf("%s $\\times%s$", c(th_zs = "$\\theta_{zS}$", e_z = "$e_Z$", k_t = "$k_T$", gcsf_rel = "$u_{\\mathrm{R}}$",
                                                   beta_a = "$\\beta_A$", e_sr = "$e_{\\mathrm{sr}}$")[names(cfg$steps$R$truth_mult)],
                            num(cfg$steps$R$truth_mult)), collapse = ", "),
  valWaveOne = paste(cfg$fit$wave1, collapse = ", "))
writeLines(c("% generated by make_tables.R from Code/config.R - do not edit",
             sprintf("\\newcommand{\\%s}{%s}", names(mac), mac)), file.path(out_dir, "values.tex"))
cat("tables written to", out_dir, "\n")
