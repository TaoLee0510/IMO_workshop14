# =============================================================================
# make_tables.R
# Purpose : Write the parameter tables, calibration tables and number macros
#           of methods.tex from Code/config.R and the model code (Model3), so
#           the document cannot drift from the code. Stops if a structural
#           parameter in config.R has no description here.
# Usage   : Rscript Methods/make_tables.R        (from Model3/)
#           Reads (never writes) Code/, ../Model1/Code/config.R and ../Model2/Code/config.R
#           (Model1 and Model2 columns of the reduction tables) and the pre-registration
#           screens Results/real/prereg/N_screen.csv and P_screen.csv; the frozen Model1 and
#           Model2 engines are Code/tests/ref/model1_v10.cpp and model2_v10.cpp. The C++ code
#           is compiled into Methods/.rcpp_cache, so nothing under Code/ is written.
# Outputs : Methods/tables/*.tex, Methods/tables/values.tex
# =============================================================================

suppressPackageStartupMessages(library(data.table))
here <- normalizePath(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))))
root <- normalizePath(file.path(here, ".."))
source(file.path(root, "Code", "config.R"))
for (f in list.files(file.path(root, "Code", "R"), pattern = "\\.R$", full.names = TRUE)) source(f)
cfg <- make_config("test", root = root)
cfg$paths$rcpp_cache <- file.path(here, ".rcpp_cache")   # keep Code/ read-only
load_model3(cfg)
a4 <- load_aim4_reference(cfg)
e1 <- new.env()
Rcpp::sourceCpp(file.path(root, "Code", "tests", "ref", "model1_v10.cpp"), cacheDir = cfg$paths$rcpp_cache, env = e1)
e2 <- new.env()
Rcpp::sourceCpp(file.path(root, "Code", "tests", "ref", "model2_v10.cpp"), cacheDir = cfg$paths$rcpp_cache, env = e2)
m1env <- new.env()
sys.source(file.path(root, "..", "Model1", "Code", "config.R"), envir = m1env)
cfg_m1 <- m1env$make_config("test", root = normalizePath(file.path(root, "..", "Model1")))
st_m1 <- cfg_m1$model$struct
m1_est <- with(cfg_m1$fit$struct_est, c(base, Mabx, route_hsc, route_mpp))
m2env <- new.env()
sys.source(file.path(root, "..", "Model2", "Code", "config.R"), envir = m2env)
cfg_m2 <- m2env$make_config("test", root = normalizePath(file.path(root, "..", "Model2")))
st_m2 <- cfg_m2$model$struct
m2_est <- with(cfg_m2$fit$struct_est, c(base, Mabx, route_hsc, route_mpp))
out_dir <- file.path(here, "tables")
dir.create(out_dir, showWarnings = FALSE)
st <- cfg$model$struct
# design constants of the Model3 neutrophil chain (chosen in the pre-registration screen; not parameters)
design_k <- c("n_tr_n", "nk_n")
st_par <- st[setdiff(names(st), design_k)]
# Model2 values of the Model3 changes (the drop-one null values N, P1, P2)
m2_values <- function(s) { for (x in cfg$model$drop_null) s[names(x)] <- x; s }
# real derived data: set-point SDs, G-CSF windows, ALC by day window, phenotype counts (no fits)
cfg_real <- make_config("real", root = root)
cfg_real$paths$rcpp_cache <- cfg$paths$rcpp_cache
dreal <- if (file.exists(file.path(cfg_real$paths$derived, "patients.csv"))) read_model3_data(cfg_real) else NULL
prereg_dir <- file.path(cfg_real$paths$results, "prereg")
folds_file <- file.path(prereg_dir, "V1_folds.csv")
folds <- if (file.exists(folds_file)) fread(folds_file) else NULL

# ---- helpers ----------------------------------------------------------------------------------
num <- function(x, d = 3) {
  vapply(x, function(v) {
    if (is.na(v)) return("--")
    if (v == round(v) && abs(v) < 1e5) return(format(v, big.mark = ""))
    if (abs(v) >= 1e3) return(sprintf("$10^{%g}$", round(log10(v), 2)))
    format(signif(v, d), drop0trailing = TRUE, scientific = FALSE)
  }, character(1))
}
# signed change for the tables ("0" or "$+2$")
chg <- function(d) if (d == 0) "0" else sprintf("$%+d$", d)
dnll <- function(x) ifelse(abs(x) < 0.05, "0", sprintf("$%+.1f$", x))
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
# role: estimated in F1 (population) / fixed / clone type / simulation (what-if only) / design constant
est <- c(cfg$fit$struct_est$base, cfg$fit$struct_est$Mabx, cfg$fit$struct_est$route_hsc,
         cfg$fit$struct_est$route_mpp)
meta <- rbindlist(list(
  # stem cells, clone, MPP components, lymphocytes
  list("r_a", "$r_a$", "Clone sensitivity to chemotherapy relative to normal HSC", "--", "clone scenario", "marrow"),
  list("r_b", "$r_b$", "Clone sensitivity to the CAR T phase relative to normal HSC", "--", "clone scenario", "marrow"),
  list("delta", "$\\delta$", "Output of clonal HSC relative to normal HSC", "--", "clone scenario", "marrow"),
  list("k_p", "$k_P$", "Exit rate of every MPP component (differentiation and death)", "d$^{-1}$", "Aim4 MPP ($k_{\\mathrm{MPP}}$)", "marrow"),
  list("sr_p", "$\\sigma_P$", "Self-renewal rate of every MPP component at steady state ($\\sigma_P<k_P$)", "d$^{-1}$", "Aim4 MPP ($\\sigma_{\\mathrm{MPP}}$)", "marrow"),
  list("gmax", "$G_{\\max}$", "Upper bound of the count feedback factor", "--", "Aim4", "marrow"),
  list("gamma", "$\\gamma$", "Exponent of the count feedback (shared by the three myeloid components)", "--", "Aim4 $\\gamma_M$ (prior centre)", "marrow"),
  list("c_l", "$c_L$", "Turnover of blood lymphocytes", "d$^{-1}$", "calibrated (K1): ALC course", "marrow"),
  list("th_al", "$\\theta_{aL}$", "Lymphodepletion of blood lymphocytes per unit exposure", "d$^{-1}$", "calibrated (K1): real ALC", "marrow"),
  # maturation chains and blood
  list("mtt_n", "$\\mathrm{MTT}_{\\mathrm{N}}$", "Neutrophil maturation transit time (all $n_{\\mathrm{N}}$ compartments)", "d", "5--7 d", "chains"),
  list("n_tr_n", "$n_{\\mathrm{N}}$", "Number of compartments of the neutrophil chain (Model2: 3)", "--", "N screen", "chains"),
  list("nk_n", "$k_{\\mathrm{N}}$", "Of these, proliferating and killed by chemotherapy; the other $n_{\\mathrm{N}}-k_{\\mathrm{N}}$ form the storage pool (Model2: 2)", "--", "N screen", "chains"),
  list("mtt_mo", "$\\mathrm{MTT}_{\\mathrm{Mo}}$", "Monocyte maturation transit time", "d", "1--2 d", "chains"),
  list("mtt_pl", "$\\mathrm{MTT}_{\\mathrm{Pl}}$", "Megakaryocyte maturation time (endomitosis)", "d", "assumption", "chains"),
  list("c_n", "$c_{\\mathrm{N}}$", "Neutrophil blood turnover", "d$^{-1}$", "\\citet{dancey1976}", "chains"),
  list("c_mo", "$c_{\\mathrm{Mo}}$", "Monocyte blood turnover", "d$^{-1}$", "Aim4", "chains"),
  list("c_pl", "$c_{\\mathrm{Pl}}$", "Platelet blood turnover", "d$^{-1}$", "\\citet{harker1969}", "chains"),
  list("kappa_c", "$\\kappa_c$", "Chemotherapy kill of the proliferating chain precursors ($T_{\\mathrm{N}1}$--$T_{\\mathrm{N}k_{\\mathrm{N}}}$; $T_{\\ell1}$, $T_{\\ell2}$ of Mo and Pl) relative to the patient's MPP kill rate $\\theta_{aP}$", "--", "prior centre: as sensitive as the MPP", "chains"),
  list("th_am", "$\\theta_{aM}$", "Chemotherapy kill of circulating monocytes per unit exposure", "d$^{-1}$", "calibrated (K1): real AMC days $-1$..0", "chains"),
  list("th_cp", "$\\theta_{cP}$", "Consumption of circulating platelets per unit ferritin driver $D$ (P1; screened, not included)", "d$^{-1}$", "P1 screen: no gain", "chains"),
  # inputs and inflammation
  list("k_e", "$k_e$", "Elimination of the chemotherapy exposure", "d$^{-1}$", "Aim4", "inputs"),
  list("beta_renal", "$\\beta$", "Exponent of the renal scaling $(90/\\mathrm{eGFR})^{\\beta}$", "--", "Aim4", "inputs"),
  list("u_br_rgo", "$u_{\\mathrm{br,RGO}}$", "Exposure rate of R-GemOx bridging", "d$^{-1}$", "Aim4", "inputs"),
  list("u_br_cy", "$u_{\\mathrm{br,Cy}}$", "Exposure rate of other (cyclophosphamide-type) bridging", "d$^{-1}$", "Aim4", "inputs"),
  list("k_i", "$k_I$", "Rate of each stage of the two-stage marrow lag", "d$^{-1}$", "Aim4", "inputs"),
  list("emax", "$E_{\\max}$", "Maximal suppression of MPP self-renewal by marrow inflammation (replaces Aim4 $e_{\\mathrm{sr}}$)", "--", "screen (K1)", "inputs"),
  list("ic50", "$\\mathrm{IC}_{50}$", "Lagged driver giving half-maximal suppression", "--", "Aim4", "inputs"),
  list("hill", "$h_E$", "Hill coefficient of the suppression", "--", "Aim4", "inputs"),
  list("eta_syn", "$\\eta$", "Synergy of CAR T-phase HSC injury with prior HSC loss (Mab$\\times$ only)", "--", "prior", "inputs"),
  list("gcsf_fold", "$u_{\\mathrm{G}}$", "Fold increase of neutrophil chain entry at full G-CSF effect ($E_G=1$)", "--", "assumption (prior centre)", "inputs"),
  list("gcsf_rel", "$u_{\\mathrm{R}}$", "Acceleration of the transit through the neutrophil storage pool and the release into blood at full G-CSF effect", "--", "assumption (prior centre)", "inputs"),
  list("k_g", "$k_G$", "Rate at which the G-CSF effect $E_G$ follows the dosing indicator", "d$^{-1}$", "assumption", "inputs"),
  list("tpo_fold", "$u_{\\mathrm{T}}$", "Fold increase of megakaryocyte chain entry under TPO-RA", "--", "assumption", "inputs"),
  # CAR T - target interaction
  list("car_t0", "$C_{T,0}$", "CAR T cells in blood right after the infusion (day 0)", "$10^9$/L", "assumption ($\\approx1$ cell/$\\mu$L)", "cart"),
  list("k_ag", "$A_{50}$", "Antigen load giving half-maximal occupancy (expansion and killing; relative load)", "--", "assumption", "cart"),
  list("kappa", "$k_A$", "Killing of CD19$^+$ targets per CAR T cell", "(10$^9$/L)$^{-1}$ d$^{-1}$", "calibrated (Aim4 K1): typical peak on day 8", "cart"),
  list("k_t", "$k_T$", "Contraction rate of CAR T cells in blood", "d$^{-1}$", "assumption (prior centre)", "cart"),
  list("beta_a", "$\\beta_A$", "Initial antigen load per unit centred log LDH, $A_i(0)=e^{\\beta_A\\,\\mathrm{ldh}_{c,i}}$", "--", "assumption (prior centre)", "cart"),
  list("k_z", "$k_Z$", "Rate of each stage of the two-stage lag of the interaction signal $Z$", "d$^{-1}$", "assumption", "cart"),
  list("z50", "$Z_{50}$", "Signal $Z$ giving half-maximal suppression of MPP self-renewal", "$10^9$/L", "calibrated (Aim4 K1): typical peak of $Z$", "cart"),
  list("th_zs", "$\\theta_{zS}$", "HSC injury per unit signal $Z$ (CAR T--target route of ``b'')", "(10$^9$/L)$^{-1}$ d$^{-1}$", "centre calibrated (Aim4 K1)", "cart"),
  list("e_z", "$e_Z$", "Maximal suppression of MPP self-renewal by $Z$ (neutrophil, monocyte and lymphoid components)", "--", "assumption (prior centre)", "cart"),
  list("d_zpl", "$d_{Z\\mathrm{Pl}}$", "Logit offset of the maximal suppression of megakaryocyte MPP self-renewal by $Z$, $e_{Z,\\mathrm{Pl}}=\\mathrm{logit}^{-1}(\\mathrm{logit}\\,e_Z+d_{Z\\mathrm{Pl}})$ (P2)", "--", "prior centre: $e_{Z,\\mathrm{Pl}}=e_Z$", "cart")
), use.names = FALSE)
setnames(meta, c("code", "symbol", "meaning", "unit", "source", "group"))
missing_meta <- setdiff(names(st), meta$code)
extra_meta <- setdiff(meta$code, names(st))
if (length(missing_meta) || length(extra_meta))
  stop("parameter metadata out of date - missing: ", paste(missing_meta, collapse = ", "),
       "; unknown: ", paste(extra_meta, collapse = ", "))
role_of <- function(cd) {
  if (cd == "eta_syn") return("estimated (Mab$\\times$)")
  if (cd == "th_zs") return("estimated (b models; not \\_Zm)")
  if (cd %in% c("e_z", "d_zpl")) return("estimated (not \\_Zh)")
  if (cd %in% design_k) return("design constant")
  if (cd == "th_cp") return("fixed at 0")
  if (cd %in% est) return("estimated")
  if (cd %in% c("tpo_fold")) return("simulation")
  if (cd %in% c("r_a", "r_b", "delta")) return("clone type")
  "fixed"
}
meta[, value := num(unname(st[code]), 4)]
meta[code == "eta_syn", value := sprintf("0 ($e^{%s}$ in Mab$\\times$)", num(cfg$fit$eta_syn_prior[["centre"]]))]
meta[code %in% c("r_a", "r_b", "delta"), value := "1 (no clone)"]
meta[, role := vapply(code, role_of, character(1))]
hdr <- "Symbol & Code & Meaning & Unit & Value & Role & Source \\\\"
spec <- "ll>{\\raggedright\\arraybackslash}p{4.2cm}>{\\raggedright\\arraybackslash}p{1.5cm}>{\\raggedright\\arraybackslash}p{1.2cm}>{\\raggedright\\arraybackslash}p{1.9cm}>{\\raggedright\\arraybackslash}p{2.4cm}"
for (g in c("marrow", "chains", "inputs", "cart")) {
  m <- meta[group == g]
  cap <- switch(g,
    marrow = "Structural parameters of the stem cells (clone traits), the four MPP components and the lymphocytes. Values are prior centres for estimated parameters.",
    chains = sprintf("Structural parameters of the maturation chains and blood: the two chemotherapy kill terms of Model2, the two design constants of the Model3 neutrophil chain and the screened platelet consumption term P1 (fixed at 0). $\\mathrm{MTT}_\\ell$ is converted to $\\kappa_\\ell = n_\\ell/\\mathrm{MTT}_\\ell$ with $n_{\\mathrm{N}}=%s$ and $n_{\\mathrm{Mo}}=n_{\\mathrm{Pl}}=3$. Values are prior centres for estimated parameters.", num(st[["n_tr_n"]])),
    inputs = "Structural parameters of the inputs: chemotherapy exposure, ferritin-driven marrow inflammation and supportive care (G-CSF is an input of the fit; TPO-RA only in what-if analyses). Values are prior centres for estimated parameters.",
    cart = "Structural parameters of the CAR T--target interaction module (Aim4 v1.3) and the platelet offset $d_{Z\\mathrm{Pl}}$ of Model3 (P2, working scale: identity). Values are prior centres for estimated parameters; the antigen load is relative (median patient $A(0)=1$).")
  write_tab(paste0("tab_struct_", g, ".tex"), cap, paste0("tab:struct_", g), spec, hdr,
            m[, row_tex(symbol, code(code), meaning, unit, value, role, source)], size = "\\footnotesize")
}

# ---- self-renewing cells ----------------------------------------------------------------------------
kp <- st[["k_p"]]; sp <- st[["sr_p"]]
lab_mpp <- c(N = "MPP$_{\\mathrm{N}}$ (neutrophils)", Mo = "MPP$_{\\mathrm{Mo}}$ (monocytes)",
             Pl = "MPP$_{\\mathrm{Pl}}$ (megakaryocytes)", L = "MPP$_{L}$ (lymphocytes)")
fb_mpp <- c(N = "$B_{\\mathrm{N}}$ (ANC)", Mo = "$B_{\\mathrm{Mo}}$ (AMC)", Pl = "$B_{\\mathrm{Pl}}$ (PLT)", L = "none")
sup_mpp <- c(N = "$E(I)$; $Z$ ($e_Z$)", Mo = "$E(I)$; $Z$ ($e_Z$)", Pl = "$E(I)$; $Z$ ($e_{Z,\\mathrm{Pl}}$)",
             L = "$E(I)$; $Z$ ($e_Z$)")
rows <- c(row_tex("HSC", "$r_S$ (random effect)", "--", num(exp(cfg$model$re_prior$mu[["r_s"]])), "--", "--",
                  "niche $(1-S)$", "--"),
          vapply(names(lab_mpp), function(x)
            row_tex(lab_mpp[[x]], "$\\sigma_P$ (\\texttt{sr\\_p})", num(kp), num(sp), num(sp / kp, 2), num(1 / (kp - sp), 3),
                    fb_mpp[[x]], sup_mpp[[x]]), character(1)))
write_tab("tab_selfrenewal.tex",
          "The self-renewing cells of Model3 (as Model1 and Model2). $k_P$: exit rate (differentiation and death, d$^{-1}$); $\\sigma_P$: self-renewal rate at steady state (d$^{-1}$, $\\sigma_P<k_P$); $\\sigma_P/k_P$: share of a component's production that comes from self-renewal at steady state; $1/(k_P-\\sigma_P)$: effective residence time (d). The four MPP components share $k_P$ and $\\sigma_P$; the three myeloid components sense the blood count of their own lineage, the lymphoid component has no count feedback. All components' self-renewal is suppressed by marrow inflammation and by the CAR T--target signal; for the megakaryocyte component the maximal suppression by the signal is $e_{Z,\\mathrm{Pl}}$ (Model3 P2).",
          "tab:selfrenewal", "lllllll>{\\raggedright\\arraybackslash}p{2.0cm}>{\\raggedright\\arraybackslash}p{1.9cm}",
          "Cell & Self-renewal rate & $k_P$ & $\\sigma_P$ & $\\sigma_P/k_P$ & $1/(k_P-\\sigma_P)$ & Count feedback & Suppressed by \\\\",
          rows, size = "\\footnotesize")

# ---- parameter reduction relative to Aim4 v1.3 ------------------------------------------------------
module_of <- function(cd) {
  if (cd %in% c("alpha", "r_a", "r_b", "delta", "eta_syn", "th_zs")) return("HSC and clone ($\\alpha$, $r_a$, $r_b$, $\\delta$, $\\eta$, $\\theta_{zS}$)")
  if (grepl("^(k|sr)_(mpp|cmp|gmp|mep|clp|ly|p)$", cd)) return("Progenitor exit and self-renewal rates")
  if (cd %in% c("gmax", "gamma_m", "gamma") || startsWith(cd, "fb_")) return("Self-renewal feedback ($G_{\\max}$, exponent, switches)")
  if (cd %in% c("gam_n", "gam_mo", "gam_pl")) return("Feedback at the chain entry ($\\gamma_\\ell$)")
  if (startsWith(cd, "mtt_")) return("Maturation transit times")
  if (cd %in% c("c_n", "c_mo", "c_pl", "c_l")) return("Blood turnover")
  if (cd %in% c("k_i", "emax", "ic50", "hill", "e_mo", "e_pl", "e_l", "e_sr")) return("Marrow inflammation")
  if (cd == "th_al") return("Lymphodepletion $\\theta_{aL}$")
  if (cd %in% c("k_e", "beta_renal", "u_br_rgo", "u_br_cy")) return("Chemotherapy exposure")
  if (cd %in% c("gcsf_fold", "gcsf_rel", "k_g", "tpo_fold")) return("Supportive care (G-CSF, TPO-RA)")
  if (cd %in% c("car_t0", "k_ag", "kappa", "k_t", "beta_a", "k_z", "z50", "e_z")) return("CAR T--target module")
  if (cd %in% c("d_n", "d_mo", "d_pl")) return("Division bookkeeping (output only)")
  if (cd %in% c("kappa_c", "th_am")) return("Chemotherapy kill of chain precursors and monocytes ($\\kappa_c$, $\\theta_{aM}$)")
  if (cd == "d_zpl") return("Platelet-specific suppression by $Z$ ($d_{Z\\mathrm{Pl}}$; Model3 P2)")
  if (cd == "th_cp") return("Platelet consumption by inflammation ($\\theta_{cP}$; P1, fixed at 0)")
  stop("no module for ", cd)
}
mods <- c("HSC and clone ($\\alpha$, $r_a$, $r_b$, $\\delta$, $\\eta$, $\\theta_{zS}$)", "Progenitor exit and self-renewal rates",
          "Self-renewal feedback ($G_{\\max}$, exponent, switches)", "Feedback at the chain entry ($\\gamma_\\ell$)",
          "Maturation transit times", "Blood turnover", "Marrow inflammation", "Lymphodepletion $\\theta_{aL}$",
          "Chemotherapy exposure", "Supportive care (G-CSF, TPO-RA)", "CAR T--target module",
          "Division bookkeeping (output only)", "Chemotherapy kill of chain precursors and monocytes ($\\kappa_c$, $\\theta_{aM}$)",
          "Platelet-specific suppression by $Z$ ($d_{Z\\mathrm{Pl}}$; Model3 P2)",
          "Platelet consumption by inflammation ($\\theta_{cP}$; P1, fixed at 0)")
mod_tab <- function(nm) table(factor(vapply(nm, module_of, character(1)), levels = mods))
a4_mod <- mod_tab(names(aim4_v13_struct))
m1_mod <- mod_tab(names(st_m1))
m2_mod <- mod_tab(names(st_m2))
m3_mod <- mod_tab(names(st_par))
stopifnot(sum(a4_mod) == length(aim4_v13_struct), sum(m1_mod) == length(st_m1), sum(m2_mod) == length(st_m2),
          sum(m3_mod) == length(st_par), !any(design_k %in% names(st_m2)))
red_rows <- c(vapply(seq_along(mods), function(j)
  row_tex(mods[j], a4_mod[[j]], m1_mod[[j]], m2_mod[[j]], m3_mod[[j]], chg(m3_mod[[j]] - a4_mod[[j]])), character(1)), "\\midrule",
  row_tex("\\textbf{Total}", sum(a4_mod), sum(m1_mod), sum(m2_mod), sum(m3_mod),
          sprintf("$%+d$ ($%+d\\%%$)", sum(m3_mod) - sum(a4_mod), as.integer(round(100 * (sum(m3_mod) - sum(a4_mod)) / sum(a4_mod))))),
  row_tex("Design constants of the neutrophil chain ($n_{\\mathrm{N}}$, $k_{\\mathrm{N}}$; Model3 N)", "--", "--", "--",
          length(design_k), "--"))
write_tab("tab_reduction.tex",
          "Structural (population) parameters of Aim4 v1.3, Model1 v1.0, Model2 v1.0 and Model3 v1.0 by module. Aim4: \\texttt{Aim4/Code/config.R} of 4 November 2026; Model1, Model2 and Model3: their \\texttt{Code/config.R}. Change: Model3 minus Aim4. The two design constants of the Model3 neutrophil chain are not counted as parameters; in the earlier engines the chain length was fixed in the code (three compartments; in Model2 the first two killed by chemotherapy).",
          "tab:reduction", ">{\\raggedright\\arraybackslash}p{7.7cm}rrrrr", "Module & Aim4 v1.3 & Model1 & Model2 & Model3 & Change \\\\",
          red_rows)
# by level (largest model: Mab$\times$ with both routes of Z and a clone)
a4_est <- c("mtt_n", "mtt_pl", "mtt_mo", "ic50", "k_i", "gamma_m", "u_br_rgo", "e_sr", "k_t", "beta_a",
            "gcsf_fold", "gcsf_rel", "th_zs", "e_z", "eta_syn")
spec_big <- model_spec("MabxC_ddr")
m3_est <- struct_est_names(spec_big, cfg)
stopifnot(all(a4_est %in% names(aim4_v13_struct)), all(m1_est %in% names(st_m1)), all(m2_est %in% names(st_m2)),
          all(m3_est %in% names(st_par)))
n_re <- length(spec_big$re)
n_beta <- sum(lengths(cfg$model$covariates[spec_big$re]))
n_bl <- length(cfg$model$baseline$vars)
n_pop <- function(n_est) 2 * n_re + n_beta + n_est + 2 * n_bl
lev <- data.table(
  level = c("Structural parameters", "\\quad fixed (literature, calibration, clone type, what-if)",
            "\\quad estimated (largest model)", "Individual parameters per patient (random effects + baselines)",
            "Estimated population parameters ($\\mu$, $\\omega$, $\\beta$, structural, $\\sigma$, $\\tau$)",
            "Self-renewing cell types with their own rates", "ODE states (C++ engine)",
            "\\quad additional moment states (time and divisions since MPP entry)", "Length of the C++ parameter vector"),
  aim4 = c(length(aim4_v13_struct), length(aim4_v13_struct) - length(a4_est), length(a4_est), n_re + n_bl,
           n_pop(length(a4_est)), 7, length(a4$a4_state_names()), 38, length(a4$a4_param_names())),
  m1 = c(length(st_m1), length(st_m1) - length(m1_est), length(m1_est), n_re + n_bl, n_pop(length(m1_est)), 2,
         length(e1$m1_state_names()), 0, length(e1$m1_param_names())),
  m2 = c(length(st_m2), length(st_m2) - length(m2_est), length(m2_est), n_re + n_bl, n_pop(length(m2_est)), 2,
         length(e2$m2_state_names()), 0, length(e2$m2_param_names())),
  m3 = c(length(st_par), length(st_par) - length(m3_est), length(m3_est), n_re + n_bl, n_pop(length(m3_est)), 2,
         length(m3_state_names()), 0, length(m3_param_names())))
lev_rows <- lev[, row_tex(level, aim4, m1, m2, m3, ifelse(m3 == aim4, "0", sprintf("$%+d$", m3 - aim4)))]
lev_rows <- append(lev_rows, row_tex("\\quad design constants (neutrophil chain: $n_{\\mathrm{N}}$, $k_{\\mathrm{N}}$)",
                                     "--", "--", "--", length(design_k), "--"), after = 3)
n_tn_slots <- sum(grepl("^TN[0-9]+$", as.character(m3_state_names())))
write_tab("tab_levels.tex",
          sprintf("Model size of Aim4 v1.3, Model1 v1.0, Model2 v1.0 and Model3 v1.0 by level. ``Largest model'': Mab$\\times$ with both routes of $Z$ and a clone. Change: Model3 minus Aim4. The Model3 engine reserves %s neutrophil compartments; with $n_{\\mathrm{N}}=%s$, %s of its %s states %s constant.",
                  num(n_tn_slots), num(st[["n_tr_n"]]), num(n_tn_slots - st[["n_tr_n"]]), num(length(m3_state_names())),
                  if (n_tn_slots - st[["n_tr_n"]] == 1) "stays" else "stay"),
          "tab:levels", ">{\\raggedright\\arraybackslash}p{8.7cm}rrrrr", "Level & Aim4 v1.3 & Model1 & Model2 & Model3 & Change \\\\",
          lev_rows)

# ---- random effects -------------------------------------------------------------------------------
re <- cfg$model$re_prior
cv <- cfg$model$covariates
cvlab <- c(prior_lines_c = "prior lines", prior_asct = "prior ASCT", age_c = "age", log_rdw_c = "log RDW", mcv_c = "MCV")
re_rows <- rbindlist(list(
  list("r_s", "$r_S$", "HSC self-renewal (regrowth) rate; reserve, Q1", "log", exp(re$mu[["r_s"]])),
  list("f0", "$f_0$", "Pre-treatment clone fraction of the HSC pool, Q1", "logit", plogis(re$mu[["f0"]])),
  list("th_as", "$\\theta_{aS}$", "HSC injury per unit chemotherapy exposure (a), Q2", "log", exp(re$mu[["th_as"]])),
  list("th_ap", "$\\theta_{aP}$", "MPP injury per unit exposure (all four components; proliferating chain precursors at $\\kappa_c\\theta_{aP}$)", "log", exp(re$mu[["th_ap"]])),
  list("th_bs", "$\\theta_{bS}$", "HSC injury per unit lagged ferritin driver (ferritin route of b), Q2", "log", exp(re$mu[["th_bs"]])),
  list("rho", "$\\rho_T$", "CAR T antigen-driven expansion rate (all Model3 models)", "log", exp(re$mu[["rho"]])),
  list("slope", "$\\mathrm{Slope}$", "MF: drug effect per unit exposure", "log", exp(re$mu[["slope"]])),
  list("mtt_scale", "$s_{\\mathrm{MTT}}$", "MF: individual scaling of all MTTs", "log", exp(re$mu[["mtt_scale"]]))
), use.names = FALSE)
setnames(re_rows, c("code", "symbol", "meaning", "scale", "centre"))
re_rows[, omega := num(re$omega[code])]
re_rows[, covs := vapply(code, function(k) if (length(cv[[k]])) paste(cvlab[cv[[k]]], collapse = ", ") else "--", character(1))]
write_tab("tab_re.tex",
          "Individual parameters (random effects). $\\psi_{ik}=g_k^{-1}(\\eta_{ik})$ with $\\eta_{ik}=\\mu_k+\\boldsymbol\\beta_k^\\top\\mathbf{x}_{ik}+\\omega_k z_{ik}$, $z_{ik}\\sim N(0,1)$. Centre: the prior centre $g_k^{-1}(\\mu_k)$ (the Model3 wave-1 fits start from the Model2 population instead, Section~\\ref{sec:nested}); $\\omega_k$: prior spread on the working scale (the estimate is shrunk towards it, Section~\\ref{sec:its}). The $\\theta_{aP}$ centre is kept from Model1 and the $\\rho_T$ centre is calibrated (Section~\\ref{sec:calibration}). Each patient also has one baseline deviation per lineage (Section~\\ref{sec:baseline}).",
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
stopifnot(identical(ms, cfg_m2$model$mf_struct))   # MF is not refitted: same as Model2
mf_meta <- rbindlist(list(
  list("mtt_n", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{N}}$", "Neutrophil mean transit time ($=4/k_{\\mathrm{tr}}$)", "d"),
  list("mtt_pl", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{Pl}}$", "Platelet mean transit time", "d"),
  list("mtt_mo", "$\\mathrm{MTT}^{\\mathrm{MF}}_{\\mathrm{Mo}}$", "Monocyte mean transit time", "d"),
  list("gam_mf", "$\\gamma^{\\mathrm{MF}}$", "Feedback exponent (all lineages)", "--"),
  list("slope_pl", "$\\rho_{\\mathrm{Pl}}$", "Drug effect on platelets relative to neutrophils", "--"),
  list("slope_mo", "$\\rho_{\\mathrm{Mo}}$", "Drug effect on monocytes relative to neutrophils", "--"),
  list("k_e", "$k_e$", "Exposure elimination (as Model3)", "d$^{-1}$"),
  list("beta_renal", "$\\beta$", "Renal scaling exponent (as Model3)", "--"),
  list("u_br_rgo", "$u_{\\mathrm{br,RGO}}$", "R-GemOx bridging exposure rate", "d$^{-1}$"),
  list("u_br_cy", "$u_{\\mathrm{br,Cy}}$", "Other bridging exposure rate", "d$^{-1}$")), use.names = FALSE)
setnames(mf_meta, c("code", "symbol", "meaning", "unit"))
stopifnot(setequal(mf_meta$code, names(ms)))
write_tab("tab_mf.tex", "Population parameters of the Friberg reference model MF (prior centres; identical to Model2, whose MF fit is reused).", "tab:mf", "ll>{\\raggedright\\arraybackslash}p{6cm}lll",
          "Symbol & Code & Meaning & Unit & Value & Role \\\\",
          mf_meta[, row_tex(symbol, code(code), meaning, unit, num(ms[code]),
                            ifelse(code %in% cfg$fit$struct_est$MF, "estimated", "fixed"))], size = "\\footnotesize")

# ---- calibration and stability (model runs) --------------------------------------------------------
# Model1 and Model2 values on the Model3 engine (exact nesting, K1 tests against the frozen engines)
cfg_m2v <- cfg
cfg_m2v$model$struct <- m2_values(st)
cfg_m1v <- cfg_m2v
cfg_m1v$model$struct[c("kappa_c", "th_am")] <- 0
cmp <- rbind(
  nadir_summary(typical_aim4_run(cfg))[, engine := sprintf("Aim4 v1.3 ($\\theta_{aP}=%.2f$)", aim4_v13_th_ap_centre)],
  nadir_summary(typical_m0_run(cfg_m1v))[, engine := "Model1 values"],
  nadir_summary(typical_m0_run(cfg_m2v))[, engine := "Model2 values"],
  nadir_summary(typical_m0_run(cfg))[, engine := "Model3"])
write_tab("tab_aim4.tex",
          sprintf("Typical-patient early nadir and recovery (no CAR T phase, no HSC injury, typical design): the Aim4 v1.3 engine at its prior centre, and the Model3 engine at the prior centres ($\\theta_{aP}=%s$) with the Model1 values ($\\kappa_c=\\theta_{aM}=0$, neutrophil chain of 3 compartments), the Model2 values (3 compartments, the first 2 killed) and the Model3 values ($n_{\\mathrm{N}}=%s$, $k_{\\mathrm{N}}=%s$). Recovery: first day after the nadir with the relative count $\\geq 0.8$.",
                  num(exp(re$mu[["th_ap"]])), num(st[["n_tr_n"]]), num(st[["nk_n"]])),
          "tab:aim4", "llrrr", "Engine & Lineage & Nadir day & Nadir (relative) & Recovery day \\\\",
          cmp[, row_tex(engine, lineage, num(nadir_day), num(nadir_rel, 2), num(recovery_day))])
cal_ap <- calibrate_th_ap(cfg)
cal_ly <- calibrate_lymph(cfg)
cal_am <- calibrate_th_am(cfg)
# AMC of the typical design on days -1..0 with kappa_c = th_aM = 0 (Model1) and at the centres (Model3)
tdays <- seq(-30, 5, by = 0.25)
amc_ld <- function(s) stats::median(typical_m0_run(cfg, s, days = tdays)[tdays >= -1 & tdays <= 0, "amc_rel"])
amc_ld0 <- amc_ld(cfg_m1v$model$struct); amc_ld1 <- amc_ld(st)
# typical ANC before and after the nadir with the Model2 chain and the Model3 chain (same MTT_N)
ddays <- seq(-30, 30, by = 0.25)
anc_m2 <- typical_m0_run(cfg_m2v, days = ddays)[, "anc_rel"]
anc_m3 <- typical_m0_run(cfg, days = ddays)[, "anc_rel"]
win_mean <- function(x, a, b) mean(x[ddays >= a & ddays <= b])
scr <- phenotype_screen(cfg)
ab <- c(quick = "Q", second_wave = "SW", delayed = "D")
sb <- dcast(scr$b_only, A + th_bs ~ emax, value.var = "phenotype")
em <- setdiff(names(sb), c("A", "th_bs"))
a_cls <- unique(scr$a_only$phenotype)
a_txt <- if (identical(a_cls, "quick")) "gave quick recovery over the whole grid" else
  paste0("gave only ", paste(sort(gsub("_", " ", a_cls)), collapse = " or "), " recovery")
write_tab("tab_screen.tex",
          sprintf("Phenotypes that the CAR T phase alone can give (model Mb, typical design, ANC set point 3.5, prior centres otherwise) over the driver peak $A_D$ (peak day 9), the ferritin-route HSC injury $\\theta_{bS}$ and $E_{\\max}$. Q: quick; SW: second wave; D: delayed. A chemotherapy hit alone (model Ma, no CAR T phase; $r_S\\in\\{0.02,0.04,0.08\\}$, $\\theta_{aS}\\in\\{0.4,1,1.5,2.5\\}$) %s.", a_txt),
          "tab:screen", paste0("rr", strrep("c", length(em))),
          paste0("$A_D$ & $\\theta_{bS}$ & ", paste(sprintf("$E_{\\max}=%s$", em), collapse = " & "), " \\\\"),
          vapply(seq_len(nrow(sb)), function(i) do.call(row_tex, c(list(num(sb$A[i]), num(sb$th_bs[i])),
                                                                  as.list(ab[unlist(sb[i, em, with = FALSE])]))), character(1)))
# smallest emax of the screen giving a second wave at A_D = 2.5, th_bS = 0.06 (the rule of the emax centre)
em_sw <- scr$b_only[A == 2.5 & th_bs == 0.06 & phenotype == "second_wave", suppressWarnings(min(emax))]
stab <- stability_grid(cfg, gammas = sort(unique(c(0.25, 0.5, 0.75, cfg$model$struct_limits$gamma[2], 1.5, 2))))
stab_pos_scales <- unique(stab[max_re >= 0]$sr_scale)
stab_txt <- if (!length(stab_pos_scales)) "No value is positive." else
  sprintf("Positive values occur only with $\\sigma_P$ scaled by %s.", paste(num(sort(stab_pos_scales)), collapse = " and "))
write_tab("tab_stability.tex",
          sprintf("Linear stability of the pre-treatment steady state: largest real part of the Jacobian eigenvalues (the neutral clone-fraction and antigen directions and the unused neutrophil compartment removed) and the period of the least-damped oscillation, over the feedback exponent $\\gamma$ and a scaling of $\\sigma_P$ (capped at $0.95\\,k_P$), at the prior centres otherwise ($\\mathrm{MTT}_{\\mathrm{N}}=%s$~d). The search box of $\\gamma$ ends at %s; the larger values are shown for reference. $\\sigma_P$ is fixed in the fits. %s",
                  num(st[["mtt_n"]]), num(cfg$model$struct_limits$gamma[2]), stab_txt),
          "tab:stability", "rrrr", "$\\gamma$ & $\\sigma_P$ scaling & $\\max\\,\\mathrm{Re}\\,\\lambda$ (d$^{-1}$) & Period (d) \\\\",
          stab[, row_tex(num(gamma), num(sr_scale), num(max_re, 2), num(period_days, 2))])
# Hopf point of gamma at the configured sr_p, as a function of MTT_N, for the Model3 and the Model2 chain
hopf_of <- function(s) tryCatch(stats::uniroot(function(g) max(Re(steady_eigen(replace(s, "gamma", g)))),
                                                 c(0.3, 4), tol = 1e-6)$root, error = function(e) NA_real_)
hopf <- hopf_of(st)
hopf_m2chain <- hopf_of(m2_values(st))
# structural search box of the largest model (natural scale)
pr_big <- struct_prior_work13(m3_est, spec_big, cfg)
bx <- struct_box(pr_big, cfg)
box_lo <- work_to_struct(bx$lo, cfg); box_hi <- work_to_struct(bx$hi, cfg)
mtt_grid <- sort(unique(c(box_lo[["mtt_n"]], 2, 3.5, st[["mtt_n"]], 8.5, 12, box_hi[["mtt_n"]])))
hopf_tab <- data.table(mtt = mtt_grid)
hopf_tab[, h3 := vapply(mtt, function(m) hopf_of(replace(st, "mtt_n", m)), numeric(1))]
hopf_tab[, h2 := vapply(mtt, function(m) hopf_of(replace(m2_values(st), "mtt_n", m)), numeric(1))]
# Hopf point of the other lineages (neutrophil chain made very fast) and which lineage sets it
hopf_other <- hopf_of(replace(st, "mtt_n", 0.05))
hopf_other_nopl <- hopf_of(replace(replace(st, "mtt_n", 0.05), "mtt_pl", 0.05))
other_lin <- if (is.na(hopf_other_nopl) || hopf_other_nopl > hopf_other + 0.01) "megakaryocyte" else "monocyte"
# every corner of the MTT box at the upper gamma limit (the lineages are decoupled at the steady state)
gr <- function(k) exp(seq(log(box_lo[[k]]), log(box_hi[[k]]), length.out = 6))
corner <- CJ(mtt_n = gr("mtt_n"), mtt_pl = gr("mtt_pl"), mtt_mo = gr("mtt_mo"))
corner[, max_re := {
  s <- st
  s[c("mtt_n", "mtt_pl", "mtt_mo", "gamma")] <- c(mtt_n, mtt_pl, mtt_mo, box_hi[["gamma"]])
  max(Re(steady_eigen(s)))
}, by = .(mtt_n, mtt_pl, mtt_mo)]
box_max_re <- max(corner$max_re)
write_tab("tab_hopf.tex",
          sprintf("Value of $\\gamma$ at which the pre-treatment steady state loses stability (Hopf bifurcation), at $\\sigma_P=%s$ and the prior centres otherwise, as a function of $\\mathrm{MTT}_{\\mathrm{N}}$ (from the lower to the upper limit of its search box), for the Model3 neutrophil chain ($n_{\\mathrm{N}}=%s$, $k_{\\mathrm{N}}=%s$) and the Model2 chain (3 and 2). %s is the value set by the %s lineage, which applies when the neutrophil lineage is more stable. The search box of $\\gamma$ ends at %s. On a $6\\times6\\times6$ grid over the search boxes of $\\mathrm{MTT}_{\\mathrm{N}}$, $\\mathrm{MTT}_{\\mathrm{Pl}}$ and $\\mathrm{MTT}_{\\mathrm{Mo}}$ at $\\gamma=%s$, the largest real part of the eigenvalues is $%s$~d$^{-1}$.",
                  num(sp), num(st[["n_tr_n"]]), num(st[["nk_n"]]), sprintf("%.2f", hopf_other), other_lin,
                  num(box_hi[["gamma"]]), num(box_hi[["gamma"]]), num(box_max_re, 2)),
          "tab:hopf", "rrr", "$\\mathrm{MTT}_{\\mathrm{N}}$ (d) & Hopf $\\gamma$, Model3 chain & Hopf $\\gamma$, Model2 chain \\\\",
          hopf_tab[, row_tex(num(mtt, 3), sprintf("%.2f", h3), sprintf("%.2f", h2))])

# ---- pre-registration screens (N and P) ------------------------------------------------------------------
nsc <- fread(file.path(prereg_dir, "N_screen.csv"))
nref <- nsc[n == 3 & nk == 2]
stopifnot(nrow(nref) == 1)
nsc[, d_anc := nll_anc - nref$nll_anc]
nsel <- nsc[n == st[["n_tr_n"]] & nk == st[["nk_n"]]]
stopifnot(nrow(nsel) == 1, nsel$d_nll == min(nsc$d_nll))
bold_if <- function(x, b) ifelse(b, paste0("\\textbf{", x, "}"), x)
write_tab("tab_nscreen.tex",
          sprintf("Screen of the neutrophil chain (pre-registration, 5 November 2026; \\texttt{N\\_screen.csv}): %s variants with $n_{\\mathrm{N}}$ compartments, the first $k_{\\mathrm{N}}$ killed by chemotherapy. The individual modes of the converged Model2 Mb fit are held; $\\mathrm{MTT}_{\\mathrm{N}}$, $\\kappa_c$ and a common factor on $\\theta_{aP}$ are re-optimised (Nelder--Mead). $\\Delta$nll: change of the negative log-likelihood (counts and phenotype term, 35 patients) relative to the Model2 chain (3, 2); ANC: its ANC part. Nadir: model ANC nadir relative to baseline (median over patients) and its day (observed: 0.06 on day 7). In the screen G-CSF accelerated only the last compartment. Bold: the Model3 choice.", num(nrow(nsc))),
          "tab:nscreen", "rrrrrrrrr",
          "$n_{\\mathrm{N}}$ & $k_{\\mathrm{N}}$ & $\\Delta$nll & $\\Delta$nll (ANC) & $\\mathrm{MTT}_{\\mathrm{N}}$ (d) & $\\kappa_c$ & $\\theta_{aP}$ factor & Nadir & Nadir day \\\\",
          nsc[, row_tex(bold_if(n, n == nsel$n & nk == nsel$nk), bold_if(nk, n == nsel$n & nk == nsel$nk), dnll(d_nll),
                        dnll(d_anc), num(mtt_n, 3), num(kappa_c, 2), num(th_ap_x, 3), num(nadir, 2), num(nadir_day))])
psc <- fread(file.path(prereg_dir, "P_screen.csv"))
p_lab <- data.table(
  variant = c("neither", "P1 (th_cP x D)", "P2 (d_ZPl)", "P1 + P2", "chemo platelet loss (th x C x B_Pl)",
              "chemo kill of platelet T3", "chemo platelet loss + kill of T3", "platelet feedback exponent and cap",
              "platelet feedback + P2"),
  label = c("reference (Model2 platelet lineage)", "P1: consumption $\\theta_{cP}D\\,B_{\\mathrm{Pl}}$",
            "P2: offset $d_{Z\\mathrm{Pl}}$", "P1 + P2", "loss of blood platelets $\\theta C\\,B_{\\mathrm{Pl}}$",
            "kill of $T_{\\mathrm{Pl}3}$ (factor $k_{T3}$)", "both", "own exponent and cap of the platelet feedback",
            "platelet feedback + P2"),
  dp = c(0L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 3L))
stopifnot(all(psc$variant %in% p_lab$variant))
p_sym <- c(th_cp = "$\\theta_{cP}$", d_zpl = "$d_{Z\\mathrm{Pl}}$", th = "$\\theta$", k_T3 = "$k_{T3}$",
           gamma_pl = "$\\gamma_{\\mathrm{Pl}}$", gmax_pl = "$G_{\\max,\\mathrm{Pl}}$")
fmt_v <- function(v) {
  if (v == 0) return("0")
  if (abs(v) < 1e-3) { e <- floor(log10(abs(v))); return(sprintf("$%s\\times10^{%d}$", format(signif(v / 10^e, 3)), e)) }
  num(v, 3)
}
p_values <- function(txt) {
  kv <- strsplit(strsplit(txt, "; ", fixed = TRUE)[[1]], " ", fixed = TRUE)
  kv <- Filter(function(x) x[2] != "NA", kv)
  stopifnot(all(vapply(kv, `[`, "", 1) %in% names(p_sym)))
  paste(vapply(kv, function(x) paste0(p_sym[[x[1]]], " ", fmt_v(as.numeric(x[2]))), ""), collapse = ", ")
}
psc[, grp := match(screen, unique(screen))]
psc[, ref_nll := nll[variant == "neither"], by = screen]
psc[, ref_plt := nll_plt[variant == "neither"], by = screen]
psc <- p_lab[psc, on = "variant"]
grp_lab <- c("P1/P2" = "P1 / P2", "chemotherapy-timed" = "Chemotherapy-timed", "platelet feedback" = "Platelet feedback")
stopifnot(all(psc$screen %in% names(grp_lab)))
p_rows <- unlist(lapply(unique(psc$grp), function(g) {
  x <- psc[grp == g]
  c(if (g > 1) "\\midrule",
    x[, row_tex(ifelse(seq_len(.N) == 1, grp_lab[screen], ""), ifelse(variant == "P2 (d_ZPl)", paste0(label, "$^{\\ast}$"), label),
                vapply(value, p_values, ""), dp, dnll(nll - ref_nll), dnll(nll_plt - ref_plt))])
}))
p_mttpl <- unique(psc$mtt_pl)
write_tab("tab_pscreen.tex",
          sprintf("Screen of platelet terms (pre-registration, 5 November 2026; \\texttt{P\\_screen.csv}), with the Model3 neutrophil chain. The individual modes of the converged Model2 Mb fit are held; $\\mathrm{MTT}_{\\mathrm{N}}$, $\\kappa_c$, the $\\theta_{aP}$ factor, $\\mathrm{MTT}_{\\mathrm{Pl}}$ (floor %s~d) and the screened parameters are re-optimised. $\\Delta p$: added parameters; $\\Delta$nll: change of the negative log-likelihood relative to the reference of the same screen; PLT: its platelet part. In every variant $\\mathrm{MTT}_{\\mathrm{Pl}}$ stayed at %s~d. $^{\\ast}$: included in Model3.",
                  paste(num(p_mttpl), collapse = ", "), paste(num(p_mttpl), collapse = ", ")),
          "tab:pscreen", ">{\\raggedright\\arraybackslash}p{2.2cm}>{\\raggedright\\arraybackslash}p{4.6cm}>{\\raggedright\\arraybackslash}p{4.4cm}rrr",
          "Screen & Variant & Re-optimised values & $\\Delta p$ & $\\Delta$nll & $\\Delta$nll (PLT) \\\\",
          p_rows, size = "\\footnotesize")

# ---- physiological search box of the individual parameters ------------------------------------------
lim <- cfg$model$indiv_limits
box_rows <- rbindlist(lapply(c("r_s", "th_ap", "th_as", "th_bs", "f0", "rho", "slope", "mtt_scale"), function(k) {
  cen <- if (k == "f0") plogis(re$mu[[k]]) else exp(re$mu[[k]])
  data.table(code = k, scale = if (k == "f0") "logit" else "log", lo = lim[[k]][1], hi = lim[[k]][2], centre = cen)
}))
sym_box <- c(r_s = "$r_S$", th_ap = "$\\theta_{aP}$", th_as = "$\\theta_{aS}$", th_bs = "$\\theta_{bS}$",
             f0 = "$f_0$", rho = "$\\rho_T$", slope = "Slope (MF)", mtt_scale = "$s_{\\mathrm{MTT}}$ (MF)")
bl <- cfg$model$baseline
write_tab("tab_box.tex",
          sprintf("Search box of the individual parameters (natural scale). The box is the support of each patient's parameters: the prior is truncated to it, the differential-evolution search, the importance sampling and the MCMC stay inside it, and fitted values on a bound are reported. Round 0 (used only by fits that do not start from another fit) widens the prior but keeps this box. Baseline deviations $d_{iv}$ have the box $\\pm\\max(%s,\\ %s\\,s_{iv})$, with $s_{iv}$ the prior SD of Eq.~\\ref{eq:baseprior}.",
                  num(bl$bound_min), num(bl$bound_sd)),
          "tab:box", "llrrr", "Parameter & Scale & Lower & Upper & Prior centre \\\\",
          box_rows[, row_tex(sym_box[code], scale, num(lo), num(hi), num(centre))])

# ---- DE-ITS controls -----------------------------------------------------------------------------------
de <- cfg$fit$de; ou <- cfg$fit$outer; is_ <- cfg$fit$is; mc0 <- cfg$fit$mcmc
cores_of <- function(f) {
  x <- grep("^#SBATCH -c ", readLines(file.path(root, "Code", "hpc", f)), value = TRUE)
  if (length(x)) sub("^#SBATCH -c ([0-9]+).*$", "\\1", x[1]) else "--"
}
e_cores <- cores_of("de_estep.sbatch"); m_cores <- cores_of("de_mstep.sbatch")
sci <- function(x) sprintf("$10^{%d}$", round(log10(x)))
de_rows <- c(
  row_tex("Round 0", "Prior SD multiplier (search box unchanged; only fits without a warm start)", num(cfg$fit$round0$sd_mult)),
  row_tex("E step", "DEoptim strategy", sprintf("%d (DE/local-to-best/1/bin)", de$strategy)),
  row_tex("", "Population size NP", num(de$np)),
  row_tex("", "Differential weight $F$ / crossover CR", paste(num(de$F), "/", num(de$CR))),
  row_tex("", "Generations per chunk / convergence window", paste(num(de$chunk), "/", num(de$window))),
  row_tex("", "Best-value improvement over the window", paste0("$<$ ", sci(de$tol_best), "$\\,(|f|+1)$")),
  row_tex("", "Population objective spread (90\\% quantile $-$ min)", paste0("$<$ ", sci(de$tol_spread))),
  row_tex("", "Safety cap (generations)", sprintf("%s (Model2: %s)", num(de$max_gen), num(cfg_m2$fit$de$max_gen))),
  row_tex("", "Polish", sprintf("L-BFGS-B, $\\le%s$ iterations", num(de$polish_maxit))),
  row_tex("", "Cores per patient (one Slurm array task each)", e_cores),
  row_tex("M step", "Structural NP", sprintf("$\\max(%s,\\ %s\\times$ dimension)", num(de$struct_np_min), num(de$struct_np_per_dim))),
  row_tex("", "Population objective spread", paste0("$<$ ", sci(de$tol_spread), " per patient")),
  row_tex("", "Safety cap (generations)", sprintf("%s (Model2: %s)", num(de$max_gen), num(cfg_m2$fit$de$max_gen))),
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
write_tab("tab_de.tex", "Controls of the round-0 and DE-ITS estimation, from the configuration file and the cluster scripts.",
          "tab:de", ">{\\raggedright\\arraybackslash}p{2.2cm}>{\\raggedright\\arraybackslash}p{7.4cm}>{\\raggedright\\arraybackslash}p{4.8cm}", "Stage & Control & Value \\\\", de_rows, size = "\\footnotesize")

# ---- baselines ---------------------------------------------------------------------------------------------
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

# ---- ALC by day window (data) vs the calibrated typical patient ---------------------------------------------
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

# ---- V1 folds ------------------------------------------------------------------------------------------------
if (!is.null(folds)) {
  fl <- folds[, .(n = .N, quick = sum(phenotype == "quick"), sw = sum(phenotype == "second_wave"),
                  delayed = sum(phenotype == "delayed")), by = fold][order(fold)]
  write_tab("tab_folds.tex", "V1 folds (pre-registered): patients per fold by observed phenotype (primary rule).",
            "tab:folds", "rrrrr", "Fold & Patients & Quick & Second wave & Delayed \\\\",
            fl[, row_tex(fold, n, quick, sw, delayed)])
}

# ---- number macros used in the text -----------------------------------------------------------------
d <- cfg$data; ph <- cfg$phenotype$primary; fit <- cfg$fit; mc <- fit$mcmc; mo <- cfg$model; pt <- fit$pheno
va <- cfg$validation
nad <- function(lin, eng) cmp[engine == eng & lineage == lin]
eng_a4 <- cmp$engine[1]
ph_counts <- if (!is.null(dreal)) table(factor(dreal$phenotype$primary, levels = cfg$phenotype$classes)) else NULL
n_tests <- length(grep("^test_that\\(", readLines(file.path(root, "Code", "tests", "test_model.R"))))
psel <- psc[variant == "P2 (d_ZPl)" & screen == "P1/P2"]
mac <- c(
  valThapCentre = num(exp(re$mu[["th_ap"]])), valThapAimNadir = num(cal_ap, 3), valThapAimFour = num(aim4_v13_th_ap_centre),
  valKappaC = num(st[["kappa_c"]]), valKappaCSd = num(mo$struct_prior_sd_k[["kappa_c"]]),
  valThAm = num(st[["th_am"]]), valThAmCal = num(cal_am, 4), valThAmSd = num(mo$struct_prior_sd_k[["th_am"]]),
  valAmcTarget = num(calib_targets$amc_rel_ld), valAmcLdModelOne = num(amc_ld0, 2), valAmcLdModelThree = num(amc_ld1, 3),
  valNStructModelOne = num(length(st_m1)), valNEstModelOne = num(length(m1_est)),
  valNPopModelOne = num(n_pop(length(m1_est))),
  valNStructModelTwo = num(length(st_m2)), valNEstModelTwo = num(length(m2_est)),
  valNPopModelTwo = num(n_pop(length(m2_est))), valNFixedModelTwo = num(length(st_m2) - length(m2_est)),
  valNStatesModelTwo = num(length(e2$m2_state_names())), valNParModelTwo = num(length(e2$m2_param_names())),
  valThal = num(st[["th_al"]]), valCl = num(st[["c_l"]]), valThalCal = num(cal_ly[["th_al"]], 4),
  valClCal = num(cal_ly[["c_l"]], 4), valAlcDayThirty = num(cal_ly[["alc_d30_target"]]),
  valAncNadirAimFour = num(nad("ANC", eng_a4)$nadir_rel, 3),
  valAncNadirModel = num(nad("ANC", "Model3")$nadir_rel, 2), valAncNadirDayModel = num(nad("ANC", "Model3")$nadir_day),
  valAncNadirModelTwo = num(nad("ANC", "Model2 values")$nadir_rel, 2),
  valAncNadirDayModelTwo = num(nad("ANC", "Model2 values")$nadir_day),
  valAncNadirModelOne = num(nad("ANC", "Model1 values")$nadir_rel, 2),
  valAncNadirDayModelOne = num(nad("ANC", "Model1 values")$nadir_day),
  valAncRecAimFour = num(nad("ANC", eng_a4)$recovery_day), valAncRecModel = num(nad("ANC", "Model3")$recovery_day),
  valAncRecModelTwo = num(nad("ANC", "Model2 values")$recovery_day),
  valAmcNadirModel = num(nad("AMC", "Model3")$nadir_rel, 2), valAmcNadirDayModel = num(nad("AMC", "Model3")$nadir_day),
  valAncPreModelTwo = num(win_mean(anc_m2, -5, -1), 2), valAncPreModel = num(win_mean(anc_m3, -5, -1), 2),
  valAncLdModelTwo = num(win_mean(anc_m2, -3, 0), 2), valAncLdModel = num(win_mean(anc_m3, -3, 0), 2),
  valAncMidModelTwo = num(win_mean(anc_m2, 5, 9), 2), valAncMidModel = num(win_mean(anc_m3, 5, 9), 2),
  valKp = num(kp), valSrp = num(sp), valResP = num(1 / (kp - sp), 3), valGmax = num(st[["gmax"]]),
  valGamma = num(st[["gamma"]]),
  valGammaHi = num(box_hi[["gamma"]]), valGammaLo = num(box_lo[["gamma"]], 2),
  valGammaLimLo = num(mo$struct_limits$gamma[1]), valGammaHiModelTwo = num(cfg_m2$model$struct_limits$gamma[2]),
  valMttNLo = num(box_lo[["mtt_n"]], 2), valMttNHi = num(box_hi[["mtt_n"]]),
  valMttPlHi = num(box_hi[["mtt_pl"]], 3), valMttMoHi = num(box_hi[["mtt_mo"]], 3),
  valDzplLo = num(box_lo[["d_zpl"]]), valDzplHi = num(box_hi[["d_zpl"]]), valDzplSd = num(mo$struct_prior_sd_k[["d_zpl"]]),
  valEzPlLo = num(plogis(qlogis(st[["e_z"]]) + box_lo[["d_zpl"]]), 2),
  valEzPlHi = num(plogis(qlogis(st[["e_z"]]) + box_hi[["d_zpl"]]), 3),
  valHopf = sprintf("%.2f", hopf), valHopfModelTwoChain = sprintf("%.2f", hopf_m2chain),
  valHopfOther = sprintf("%.2f", hopf_other), valHopfOtherLineage = other_lin,
  valHopfMttHi = sprintf("%.2f", hopf_tab[mtt == max(mtt)]$h3),
  valHopfMttHiModelTwo = sprintf("%.2f", hopf_tab[mtt == max(mtt)]$h2),
  valBoxMaxRe = num(box_max_re, 2),
  valPstar = num((kp / sp)^(-1 / st[["gamma"]]), 2),
  valEmax = num(st[["emax"]]), valEmaxSd = num(mo$struct_prior_sd_k[["emax"]]),
  valEmaxScreenSW = if (is.finite(em_sw)) num(em_sw) else "--",
  valStructSd = num(mo$struct_prior_sd),
  valEtaCentre = num(fit$eta_syn_prior[["centre"]]), valEtaSd = num(fit$eta_syn_prior[["sd"]]),
  valAlcTarget = num(calib_targets$alc_rel_ld),
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
  valOmegaNu = num(fit$omega_nu), valOmegaMax = num(fit$omega_max),
  valChains = num(mc$n_chains), valIter = num(mc$n_iter), valBurn = num(mc$burn), valThin = num(mc$thin),
  valScale = num(mc$scale), valAttrib = num(fit$attrib_draws),
  valLandmarks = paste(cfg$steps$F3$landmarks, collapse = ", "), valFDraws = num(cfg$steps$F3$n_draws),
  valFKernel = num(cfg$steps$F3$n_kernel),
  valGcsfFold = num(st[["gcsf_fold"]]), valTpoFold = num(st[["tpo_fold"]]),
  valNStates = num(length(m3_state_names())), valNStatesUsed = num(length(m3_state_names()) - (n_tn_slots - st[["n_tr_n"]])),
  valNTnSlots = num(n_tn_slots),
  valNPar = num(length(m3_param_names())), valNStruct = num(length(st_par)), valNDesign = num(length(design_k)),
  valNConfigStruct = num(length(st)),
  valNStructAimFour = num(length(aim4_v13_struct)), valNEst = num(length(m3_est)),
  valNFixed = num(length(st_par) - length(m3_est)), valNFixedAimFour = num(length(aim4_v13_struct) - length(a4_est)),
  valNStatesAimFour = num(length(a4$a4_state_names())), valNPop = num(n_pop(length(m3_est))),
  valReduction = num(length(aim4_v13_struct) - length(st_par)),
  valReductionPct = num(round(100 * (length(aim4_v13_struct) - length(st_par)) / length(aim4_v13_struct))),
  # Model3 changes
  valNTrN = num(st[["n_tr_n"]]), valNkN = num(st[["nk_n"]]), valNStore = num(st[["n_tr_n"]] - st[["nk_n"]]),
  valNTrNModelTwo = num(cfg$model$drop_null$N[["n_tr_n"]]), valNkNModelTwo = num(cfg$model$drop_null$N[["nk_n"]]),
  valMttN = num(st[["mtt_n"]]), valMttNProl = num(st[["mtt_n"]] * st[["nk_n"]] / st[["n_tr_n"]], 2),
  valMttNStore = num(st[["mtt_n"]] * (st[["n_tr_n"]] - st[["nk_n"]]) / st[["n_tr_n"]], 2),
  valThcp = num(st[["th_cp"]]), valDzpl = num(st[["d_zpl"]]),
  valNScreenN = num(nrow(nsc)), valNScreenNMin = num(min(nsc$n)), valNScreenNMax = num(max(nsc$n)),
  valNScreenBest = sprintf("%.0f", nsel$d_nll), valNScreenMtt = num(nsel$mtt_n, 3), valNScreenKappa = num(nsel$kappa_c, 2),
  valNScreenNadir = num(nsel$nadir, 2), valNScreenRefNadir = num(nref$nadir, 2), valNScreenRefMtt = num(nref$mtt_n, 2),
  valPScreenDzpl = num(as.numeric(sub(".*d_zpl ([^;]+).*", "\\1", psel$value)), 3),
  valPScreenMttPl = paste(num(p_mttpl), collapse = ", "),
  valTestsN = num(n_tests),
  valMaxGen = num(de$max_gen), valMaxGenModelTwo = num(cfg_m2$fit$de$max_gen),
  valRVersion = paste(R.version$major, R.version$minor, sep = "."),
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
  valThzsSd = num(mo$struct_prior_sd_k[["th_zs"]]), valEzSd = num(mo$struct_prior_sd_k[["e_z"]]),
  valGcsfFoldLo = num(mo$struct_limits$gcsf_fold[1]), valGcsfFoldHi = num(mo$struct_limits$gcsf_fold[2]),
  valTauInit = num(bl$tau_init), valTauMin = num(bl$tau_min), valBoundSd = num(bl$bound_sd), valBoundMin = num(bl$bound_min),
  valGcsfFrom = num(d$gcsf_rule$from), valGcsfThr = num(d$gcsf_rule$thr), valGcsfDays = num(d$gcsf_rule$days),
  valGcsfUntil = num(d$gcsf_rule$until),
  valGcsfWin = if (!is.null(dreal)) num(nrow(dreal$gcsf)) else "--",
  valGcsfPts = if (!is.null(dreal)) num(data.table::uniqueN(dreal$gcsf$patient_id)) else "--",
  valGcsfFlagged = if (!is.null(dreal)) num(sum(dreal$patients$gcsf_after_d14 == 1)) else "--",
  valGcsfNoWin = if (!is.null(dreal)) num(length(setdiff(dreal$patients[gcsf_after_d14 == 1]$patient_id, dreal$gcsf$patient_id))) else "--",
  valNQuick = if (!is.null(ph_counts)) num(ph_counts[["quick"]]) else "--",
  valNSecond = if (!is.null(ph_counts)) num(ph_counts[["second_wave"]]) else "--",
  valNDelayed = if (!is.null(ph_counts)) num(ph_counts[["delayed"]]) else "--",
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
  # phenotype term, round 0, validation
  valPhW = num(pt$w), valPhS = num(pt$s), valPhEps = num(pt$eps),
  valPhSensS = paste(num(pt$sens_s), collapse = " and "), valPhSensW = paste(num(pt$sens_w), collapse = " and "),
  valPhFloor = num(-log(pt$eps), 3), valRoundZeroMult = num(fit$round0$sd_mult),
  valVk = num(va$V1$k), valVseed = num(va$V1$seed), valVper = num(35 / va$V1$k),
  valSone = num(va$S1$n_cohorts), valSoneModel = gsub("_", "\\\\_", va$S1$model),
  valSoneTruth = paste(sprintf("%s $\\times%s$", c(th_zs = "$\\theta_{zS}$", e_z = "$e_Z$", k_t = "$k_T$", gcsf_rel = "$u_{\\mathrm{R}}$",
                                                    beta_a = "$\\beta_A$", emax = "$E_{\\max}$")[names(va$S1$truth_mult)],
                             num(va$S1$truth_mult)), collapse = ", "),
  valStwoW = num(va$S2$w), valSthree = num(va$S3$n_cohorts), valSfour = num(va$S4$n_cohorts),
  valSfourModel = gsub("_", "\\\\_", va$S4$aim4_model), valSfive = num(va$S5$n_rep),
  valWaveOne = gsub("Mabx", "Mab$\\times$", paste(setdiff(cfg$fit$wave1, "MF"), collapse = ", "), fixed = TRUE),
  valDropOne = paste(cfg$fit$drop_one, collapse = " and "))
stopifnot(!anyNA(mac), all(nzchar(names(mac))))
writeLines(c("% generated by make_tables.R from Code/config.R - do not edit",
             sprintf("\\newcommand{\\%s}{%s}", names(mac), mac)), file.path(out_dir, "values.tex"))
cat("tables written to", out_dir, "\n")
