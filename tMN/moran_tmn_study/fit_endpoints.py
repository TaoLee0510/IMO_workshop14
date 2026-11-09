#!/usr/bin/env python3
"""Conditional endpoint fits for the approved nonnegative-selection scenarios.

No longitudinal observations, read-count likelihood, clinical-risk labels, or
invented detection-limit pseudocounts enter these growth fits. These are robust
descriptive fits conditional on observed positive baseline VAF, not calibrated
sequencing posteriors. The growth approximation sets mu=D0=0; transforming
alterations and existing double hits are explored separately in risk scenarios.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/tmn_fit_mpl")
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit, logit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parent
GENES = ["DNMT3A", "TP53", "PPM1D"]
COLORS = {"DNMT3A": "#3568A6", "TP53": "#D77B29", "PPM1D": "#7A5BA7"}


def kinetic_parameters(scenario):
    p = dict(scenario)
    p["alpha"] = math.log(2) / p["rise_doubling_days"]
    p["beta"] = math.log(2) / p["decay_half_life_days"]
    return p


def selection(t, s2, p):
    ld = p["s1"] * math.exp(-(t + 5) / p["tau"]) if t >= -5 else 0.0
    if t < 0:
        car = 0.0
    elif t <= p["Tmax"]:
        car = math.exp(p["alpha"] * (t - p["Tmax"]))
    else:
        car = math.exp(-p["beta"] * (t - p["Tmax"]))
    return ld + s2 * car


def propagate(z0, start, stop, s2, p, times=None):
    """Exact normalized mean-field ODE, integrated in logit coordinates."""
    z = np.atleast_1d(np.array(z0, dtype=float)).copy()
    boundaries = sorted(set([float(start), float(stop)] + [t for t in (-5.0, 0.0, float(p["Tmax"])) if start < t < stop]))
    times_arr = None if times is None else np.asarray(times, dtype=float)
    values = None if times is None else np.empty((len(z), len(times_arr)))
    if times is not None:
        values[:, times_arr == start] = z[:, None]
    for lo, hi in zip(boundaries[:-1], boundaries[1:]):
        if hi <= -5 or (p["s1"] == 0 and s2 == 0):
            if times is not None:
                use = (times_arr > lo) & (times_arr <= hi)
                values[:, use] = z[:, None]
            continue

        def rhs(t, y):
            # Right endpoints at treatment starts belong to the preceding segment
            # only during this numerical substep. The next segment starts exactly.
            tt = np.nextafter(t, lo) if t == hi and hi in (-5.0, 0.0) else t
            a = selection(tt, s2, p)
            return p["b"] * a / (1.0 + a * expit(y))

        solution = solve_ivp(rhs, (lo, hi), z, method="DOP853", rtol=2e-8, atol=1e-9, dense_output=times is not None)
        if not solution.success:
            raise RuntimeError(solution.message)
        if times is not None:
            use = (times_arr > lo) & (times_arr <= hi)
            if np.any(use):
                values[:, use] = solution.sol(times_arr[use])
        z = solution.y[:, -1]
    return z if times is None else values


def eligible(row):
    return 0 < row["first_vaf_percent"] < 50 and 0 < row["last_vaf_percent"] < 50 and row["last_day"] > row["first_day"]


def load_data():
    with (ROOT / "data/endpoints.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        for key in ("first_day", "last_day", "first_vaf_percent", "last_vaf_percent"):
            row[key] = float(row[key])
        row["source_excel_row"] = int(row["source_excel_row"])
        row["tmn_source_flag"] = int(row["tmn_source_flag"])
        row["eligible_for_conditional_fit"] = eligible(row)
    return rows


def group_rows(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["patient_id"], row["gene"])].append(row)
    return [(key, groups[key]) for key in sorted(groups, key=lambda x: (GENES.index(x[1]), int(x[0].split("_")[-1])))]


def predictions_for_rows(rows, s2, p):
    prediction = np.empty(len(rows))
    blocks = defaultdict(list)
    for index, row in enumerate(rows):
        blocks[(row["first_day"], row["last_day"])].append(index)
    for (start, stop), indices in blocks.items():
        z0 = logit(np.array([2 * rows[j]["first_vaf_percent"] / 100 for j in indices]))
        prediction[indices] = propagate(z0, start, stop, s2, p)
    return prediction


def fit_group(task):
    key, all_rows, scenario_name, scenario = task
    p = kinetic_parameters(scenario)
    rows = [r for r in all_rows if eligible(r)]
    result = {"patient_id": key[0], "gene": key[1], "scenario": scenario_name,
              "s2": None, "eligible_variants": len(rows), "total_variants": len(all_rows),
              "variant_initial_states": [{k: r[k] for k in ("variant_id", "first_day", "last_day", "first_vaf_percent", "last_vaf_percent", "source_excel_row", "quality_flags")} for r in rows],
              "status": "not_estimable_from_positive_endpoints", "boundary": False,
              "mismatch_flag": False, "flags": [], "variant_fits": []}
    if not rows:
        result["flags"] = ["no_positive_copy_neutral_compatible_endpoint_pair"]
        return result
    observed = logit(np.array([2 * r["last_vaf_percent"] / 100 for r in rows]))
    cache = {}

    def loss(s2):
        s2 = float(s2)
        if s2 not in cache:
            pred = predictions_for_rows(rows, s2, p)
            residual = pred - observed
            # Pseudo-Huber loss with delta=1 log-odds unit. A robust fitting
            # convention, not an estimate of sequencing measurement variance.
            cache[s2] = (float(np.mean(np.sqrt(1 + residual * residual) - 1)), pred)
        return cache[s2][0]

    loss_zero = loss(0.0)
    pred_zero = cache[0.0][1]
    upper = 1.0
    if np.all(pred_zero >= observed):
        best = 0.0
    else:
        loss(upper)
        while np.any(cache[upper][1] < observed) and upper < 10000:
            upper *= 10
            loss(upper)
        fitted = minimize_scalar(lambda u: loss(math.expm1(u)), bounds=(0, math.log1p(upper)), method="bounded", options={"xatol": 2e-7, "maxiter": 70})
        if not fitted.success:
            raise RuntimeError(f"Optimizer failed for {key}: {fitted.message}")
        best = math.expm1(fitted.x)
        if loss_zero <= loss(best) + 1e-11:
            best = 0.0
    objective = loss(best)
    predicted = cache[float(best)][1]
    residual = predicted - observed
    declining = any(r["last_vaf_percent"] < r["first_vaf_percent"] for r in rows)
    result.update(s2=float(best), status="conditional_descriptive_fit", boundary=best <= 1e-7,
                  mismatch_flag=declining or bool(np.max(np.abs(residual)) > 1),
                  robust_loss=objective, max_abs_logit_residual=float(np.max(np.abs(residual))),
                  rmse_logit=float(np.sqrt(np.mean(residual ** 2))), search_upper_bound=upper,
                  numerical_upper_bound_hit=best >= .999 * upper and upper >= 10000)
    if declining:
        result["flags"].append("decreasing_variant_incompatible_with_nonnegative_meanfield_selection")
    if result["boundary"]:
        result["flags"].append("s2_at_zero_boundary")
    if np.max(np.abs(residual)) > 1:
        result["flags"].append("at_least_one_absolute_logit_residual_exceeds_one")
    if result["numerical_upper_bound_hit"]:
        result["flags"].append("s2_numerical_upper_bound_hit")
    if len(rows) < len(all_rows):
        result["flags"].append("some_variants_require_detection_or_copy_number_model")
    for row, zhat, resid in zip(rows, predicted, residual):
        result["variant_fits"].append({**row, "predicted_last_vaf_percent": float(50 * expit(zhat)), "logit_residual": float(resid)})
    return result


def fit_groups(groups, scenario_name, scenario, workers):
    tasks = [(key, rows, scenario_name, scenario) for key, rows in groups]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            return list(executor.map(fit_group, tasks))
    return [fit_group(task) for task in tasks]


def save_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def flatten_fit(result):
    return {k: (";".join(v) if k == "flags" else v) for k, v in result.items() if k not in ("variant_initial_states", "variant_fits")}


def compare_genes(results, reps=2000, seed=20261105):
    estimable = [r for r in results if r["s2"] is not None]
    patients = sorted({r["patient_id"] for r in results})
    matrix = np.full((len(patients), len(GENES)), np.nan)
    for r in estimable:
        matrix[patients.index(r["patient_id"]), GENES.index(r["gene"])] = r["s2"]
    medians = np.nanmedian(matrix, axis=0)
    rng = np.random.default_rng(seed)
    draws = np.array([np.nanmedian(matrix[rng.integers(0, len(patients), len(patients))], axis=0) for _ in range(reps)])
    summaries = []
    for gi, gene in enumerate(GENES):
        sample = matrix[:, gi]
        sample = sample[np.isfinite(sample)]
        interval = np.percentile(draws[:, gi], [2.5, 97.5])
        summaries.append({"gene": gene, "n_fitted_patient_gene_groups": len(sample), "median_s2": float(medians[gi]),
                          "q25_s2": float(np.percentile(sample, 25)), "q75_s2": float(np.percentile(sample, 75)),
                          "median_cluster_bootstrap_95_low": float(interval[0]), "median_cluster_bootstrap_95_high": float(interval[1]),
                          "n_zero_boundary": sum(r["boundary"] for r in estimable if r["gene"] == gene),
                          "n_mismatch_flag": sum(r["mismatch_flag"] for r in estimable if r["gene"] == gene)})
    contrasts = []
    for left, right in (("TP53", "DNMT3A"), ("PPM1D", "DNMT3A"), ("TP53", "PPM1D")):
        a, b = GENES.index(left), GENES.index(right)
        dist = draws[:, a] - draws[:, b]
        ci = np.percentile(dist, [2.5, 97.5])
        simultaneous = np.percentile(dist, [100 * .05 / 6, 100 * (1 - .05 / 6)])
        common = np.isfinite(matrix[:, a]) & np.isfinite(matrix[:, b])
        differences = matrix[common, a] - matrix[common, b]
        paired_draws = np.array([np.median(rng.choice(differences, size=len(differences), replace=True)) for _ in range(reps)])
        paired_ci = np.percentile(paired_draws, [2.5, 97.5])
        contrasts.append({"contrast": f"{left} minus {right}", "median_difference_s2": float(medians[a] - medians[b]),
                          "cluster_bootstrap_95_low": float(ci[0]), "cluster_bootstrap_95_high": float(ci[1]),
                          "bonferroni_familywise_95_low": float(simultaneous[0]), "bonferroni_familywise_95_high": float(simultaneous[1]),
                          "paired_patients": int(np.sum(common)), "paired_median_difference": float(np.median(differences)),
                          "paired_bootstrap_95_low": float(paired_ci[0]), "paired_bootstrap_95_high": float(paired_ci[1])})
    return summaries, contrasts


def numerical_checks():
    p = kinetic_parameters({"b": 1/7, "tau": 7, "Tmax": 10, "rise_doubling_days": 2, "decay_half_life_days": 30, "s1": .1})
    comparisons = []
    for x0 in (.001, .05, .3):
        for s2 in (0, .3, 3):
            zhat = propagate([logit(x0)], -6, 180, s2, p)[0]
            x = x0
            for a, b in ((-6, -5), (-5, 0), (0, 10), (10, 180)):
                def rhs(t, y):
                    tt = np.nextafter(t, a) if t == b and b in (-5, 0) else t
                    sel = selection(tt, s2, p)
                    return p["b"] * sel * y * (1-y) / (1 + sel * y)
                sol = solve_ivp(rhs, (a, b), [x], method="RK45", rtol=1e-10, atol=1e-12)
                x = sol.y[0, -1]
            err = abs(expit(zhat) - x)
            assert err < 3e-7, err
            comparisons.append(err)
    constant = {**p, "s1": 0, "alpha": 0, "beta": 0, "Tmax": 0}
    a, t, x0 = .5, 100, .02
    z0 = logit(x0)
    rhs_value = z0 + a * np.logaddexp(0, z0) + p["b"] * a * t
    exact = brentq(lambda z: z + a * np.logaddexp(0, z) - rhs_value, z0, z0 + p["b"] * a * t)
    actual = propagate([z0], 0, t, a, constant)[0]
    assert abs(actual - exact) < 1e-7
    neutral = propagate([z0], -6, 1000, 0, {**p, "s1": 0})[0]
    assert neutral == z0
    return {"max_cell_fraction_difference_vs_independent_RK45": max(comparisons), "constant_fitness_logit_error_vs_implicit_exact_solution": abs(actual-exact), "neutral_state_preserved": True}


def plot_reference(results, groups, scenario, figure_dir):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42})
    p = kinetic_parameters(scenario)
    by_key = {(r["patient_id"], r["gene"]): r for r in results}
    with PdfPages(figure_dir / "fit_all_patient_gene_trajectories.pdf") as pdf:
        for gene in GENES:
            subset = [(key, rows) for key, rows in groups if key[1] == gene]
            for offset in range(0, len(subset), 12):
                fig, axes = plt.subplots(3, 4, figsize=(13, 8.2), squeeze=False)
                for ax, (key, rows) in zip(axes.flat, subset[offset:offset+12]):
                    fit = by_key[key]
                    for k, row in enumerate(rows):
                        color = plt.get_cmap("tab10")(k % 10)
                        valid = eligible(row)
                        ax.scatter([row["first_day"], row["last_day"]], [row["first_vaf_percent"], row["last_vaf_percent"]], s=20, color=color if valid else "#888888", marker="o" if valid else "x", zorder=4)
                        if "after" in row["quality_flags"] and "diagnos" in row["quality_flags"]:
                            ax.scatter([row["last_day"]], [row["last_vaf_percent"]], s=75, facecolors="none", edgecolors="#CC7800", linewidths=1, zorder=5)
                        if valid and fit["s2"] is not None:
                            tt = np.unique(np.r_[np.linspace(row["first_day"], min(row["last_day"], 200), 90), np.linspace(min(row["last_day"], 200), row["last_day"], 45)])
                            zz = propagate([logit(2*row["first_vaf_percent"]/100)], row["first_day"], row["last_day"], fit["s2"], p, tt)[0]
                            ax.plot(tt, 50*expit(zz), color=color, linewidth=1)
                    suffix = "not estimable" if fit["s2"] is None else f"s₂={fit['s2']:.3g}" + (" [boundary]" if fit["boundary"] else "")
                    ax.set_title(f"{key[0]} | {suffix}", fontsize=8)
                    ax.set_xlabel("Days from CAR infusion")
                    ax.set_ylabel("VAF (%)")
                    ax.set_ylim(bottom=-.02*max(1, max(max(r["first_vaf_percent"], r["last_vaf_percent"]) for r in rows)))
                    ax.axvline(0, color="#BBBBBB", linewidth=.6, linestyle=":")
                for ax in list(axes.flat)[len(subset[offset:offset+12]):]:
                    ax.axis("off")
                fig.suptitle(f"{gene}: endpoint-conditioned mean-field fits", fontsize=13)
                fig.text(.02, .015, "Colored points/lines: fitted variants. Gray ×: endpoint requires detection/copy-number model. Orange ring: last sample after a listed diagnosis.\nOnly first/last observations are used; curves between them are model interpolations. Each patient/gene shares one nonnegative s₂.", fontsize=8)
                fig.tight_layout(rect=(0, .05, 1, .95))
                fig.subplots_adjust(left=.067, right=.988, bottom=.11, top=.90, wspace=.40, hspace=.58)
                pdf.savefig(fig)
                if offset == 0:
                    fig.savefig(figure_dir / f"fit_{gene}_first_page.png", dpi=130)
                plt.close(fig)


def plot_gene_comparison(results, summaries, contrasts, out):
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1, 1.25]})
    rng = np.random.default_rng(23)
    for j, gene in enumerate(GENES):
        vals = [r["s2"] for r in results if r["gene"] == gene and r["s2"] is not None]
        s = summaries[j]
        ax.scatter(j + rng.uniform(-.14, .14, len(vals)), vals, color=COLORS[gene], s=22, alpha=.7)
        ax.errorbar(j+.24, s["median_s2"], yerr=[[s["median_s2"]-s["median_cluster_bootstrap_95_low"]], [s["median_cluster_bootstrap_95_high"]-s["median_s2"]]], fmt="D", color="black", capsize=4, markersize=5)
    ax.set_xticks(range(3), [f"{s['gene']}\nn={s['n_fitted_patient_gene_groups']}" for s in summaries])
    ax.set_ylabel("Fitted CAR selection amplitude s₂")
    ax.set_title("Patient/gene estimates and gene medians")
    bx.axvline(0, color="#999999", linewidth=.8)
    for i, c in enumerate(contrasts):
        bx.plot([c["bonferroni_familywise_95_low"], c["bonferroni_familywise_95_high"]], [i,i], color="#BBBBBB", linewidth=5)
        bx.plot([c["cluster_bootstrap_95_low"], c["cluster_bootstrap_95_high"]], [i,i], color="#333333", linewidth=2)
        bx.plot(c["median_difference_s2"], i, "o", color="#111111", markersize=5)
    bx.set_yticks(range(3), [c["contrast"].replace(" minus ", " − ") for c in contrasts])
    bx.invert_yaxis()
    bx.set_xlabel("Difference between gene medians of s₂")
    bx.set_title("Patient-cluster bootstrap contrasts")
    fig.text(.02, .01, "Conditional descriptive estimates; no clinical tMN labels enter fitting. Bars: 95% percentile intervals; pale outer bars: Bonferroni-adjusted intervals.\nBootstrap resamples patients and preserves their shared genes. Intervals exclude fixed-kinetics and baseline-VAF uncertainty.", fontsize=8)
    fig.tight_layout(rect=(0,.09,1,1))
    fig.savefig(out / "fit_gene_comparison.png", dpi=180)
    fig.savefig(out / "fit_gene_comparison.pdf")
    plt.close(fig)


def plot_sensitivity(summary_rows, out):
    names = list(dict.fromkeys(r["scenario"] for r in summary_rows))
    fig, axes = plt.subplots(1, 3, figsize=(12, 6), sharey=True, sharex=True)
    for ax, gene in zip(axes, GENES):
        rows = [r for r in summary_rows if r["gene"] == gene]
        ref = next(r["median_s2"] for r in rows if r["scenario"] == "reference")
        ax.axvline(ref, linestyle=":", color="#555555", linewidth=1)
        for row in rows:
            j = names.index(row["scenario"])
            ax.plot([ref, row["median_s2"]], [j,j], color=COLORS[gene], alpha=.4)
            ax.scatter(row["median_s2"], j, color=COLORS[gene], s=22)
        ax.set_title(gene)
        ax.set_xlabel("Median fitted s₂")
    axes[0].set_yticks(range(len(names)), names, fontsize=8)
    axes[0].invert_yaxis()
    fig.suptitle("Refitted endpoint growth under approved sensitivity scenarios", fontsize=13)
    fig.text(.02,.01,"Each row changes one shared parameter; all patient/gene amplitudes are refitted. s₁ is held fixed when b or τ changes.\nDiagnosis exclusion is a separate sensitivity because clinical reference infusions remain unresolved.",fontsize=8)
    fig.tight_layout(rect=(0,.065,1,.95))
    fig.savefig(out / "fit_sensitivity.png", dpi=170)
    fig.savefig(out / "fit_sensitivity.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--reuse-reference", action="store_true")
    args = parser.parse_args()
    results_dir, figure_dir = ROOT / "results", ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figure_dir.mkdir(exist_ok=True)
    checks = numerical_checks()
    print("Numerical checks passed.", flush=True)
    config = json.loads((ROOT / "approved_scenarios.json").read_text())
    rows = load_data()
    groups = group_rows(rows)
    reference_path = results_dir / "fit_reference.json"
    if args.reuse_reference and reference_path.exists():
        reference = json.loads(reference_path.read_text())
    else:
        reference = fit_groups(groups, "reference", config["reference"], args.workers)
    save_json(results_dir / "fit_reference.json", reference)
    write_csv(results_dir / "fit_reference.csv", [flatten_fit(r) for r in reference])
    write_csv(results_dir / "fit_variant_predictions.csv", [v for r in reference for v in r["variant_fits"]])
    write_csv(results_dir / "fit_endpoint_eligibility.csv", rows)
    print(f"Reference ready: {sum(r['s2'] is not None for r in reference)}/{len(reference)} patient/gene fits.", flush=True)
    summaries, contrasts = compare_genes(reference)
    write_csv(results_dir / "fit_gene_summary.csv", summaries)
    write_csv(results_dir / "fit_gene_contrasts.csv", contrasts)
    plot_gene_comparison(reference, summaries, contrasts, figure_dir)
    all_results = list(reference)
    sensitivity_summary = [{"scenario": "reference", **r} for r in summaries]
    for parameter, values in config["one_at_a_time_growth_sensitivity"].items():
        for value in values:
            name = f"{parameter}={value:g}"
            scenario = {**config["reference"], parameter: value}
            fits = fit_groups(groups, name, scenario, args.workers)
            all_results.extend(fits)
            simple = []
            for gene in GENES:
                valid = [r for r in fits if r["gene"] == gene and r["s2"] is not None]
                simple.append({"gene": gene, "n_fitted_patient_gene_groups": len(valid), "median_s2": float(np.median([r["s2"] for r in valid])), "n_zero_boundary": sum(r["boundary"] for r in valid), "n_mismatch_flag": sum(r["mismatch_flag"] for r in valid)})
            sensitivity_summary.extend({"scenario": name, **r} for r in simple)
            print(f"Finished sensitivity {name}.", flush=True)
    # This exclusion is only a sensitivity: the origin of clinical times is
    # unresolved and the primary objective is retrospective endpoint description.
    excluded = [r for r in rows if "after" in r["quality_flags"] and "diagnos" in r["quality_flags"]]
    retained = [r for r in rows if r not in excluded]
    diagnosis_fits = fit_groups(group_rows(retained), "exclude_after_listed_diagnosis", config["reference"], args.workers)
    all_results.extend(diagnosis_fits)
    for gene in GENES:
        valid = [r for r in diagnosis_fits if r["gene"] == gene and r["s2"] is not None]
        sensitivity_summary.append({"scenario": "exclude_after_listed_diagnosis", "gene": gene, "n_fitted_patient_gene_groups": len(valid), "median_s2": float(np.median([r["s2"] for r in valid])), "n_zero_boundary": sum(r["boundary"] for r in valid), "n_mismatch_flag": sum(r["mismatch_flag"] for r in valid)})
    # Normalize varying diagnostic fields for rectangular CSV output.
    flat = [flatten_fit(r) for r in all_results]
    fields = sorted(set().union(*(r.keys() for r in flat)))
    write_csv(results_dir / "fit_all_scenarios.csv", [{key:r.get(key) for key in fields} for r in flat])
    fields = sorted(set().union(*(r.keys() for r in sensitivity_summary)))
    write_csv(results_dir / "fit_sensitivity_summary.csv", [{key:r.get(key) for key in fields} for r in sensitivity_summary])
    plot_sensitivity(sensitivity_summary, figure_dir)
    plot_reference(reference, groups, config["reference"], figure_dir)
    report = {"method": "Exact normalized two-type Moran mean-field; pseudo-Huber endpoint logit loss, delta=1; s2 nonnegative; baseline conditioned on observed positive VAF.",
              "growth_approximation": "mu=0 and initial D=0 for conditional growth fitting only; risk scenarios independently allow mutation and pre-existing double hits.",
              "not_a_sequencing_posterior": True, "partial_pooling": False, "uncertainty": "Patient-cluster percentile bootstrap of gene summaries/contrasts; 2000 draws, seed20261105. No individual parameter confidence intervals.",
              "numerical_checks": checks, "total_variants": len(rows), "eligible_variants": sum(eligible(r) for r in rows),
              "total_patient_gene_groups": len(groups), "fitted_patient_gene_groups": sum(r["s2"] is not None for r in reference),
              "zero_boundary_groups": sum(r["boundary"] for r in reference), "mismatch_groups": sum(r["mismatch_flag"] for r in reference),
              "numeric_upper_bound_hits": sum(r.get("numerical_upper_bound_hit", False) for r in all_results),
              "after_listed_diagnosis_records_excluded_in_sensitivity": len(excluded), "gene_summary": summaries, "gene_contrasts": contrasts,
              "limitations": ["Post-diagnosis endpoints remain in retrospective primary fits and are excluded in a separate sensitivity.", "Zero and VAF>=50% endpoints are not assigned pseudocounts or capped; their patient/gene groups may be unfit.", "Variants within a patient/gene share s2 but their clonal relationships are unresolved.", "A deterministic nonnegative selection model cannot explain genuine contraction; boundary fits and residuals expose disagreement.", "Fixed kinetics, LD, baseline burden, and unresolved sample-to-progenitor mapping limit interpretation.", "No calibrated clinical tMN probability or gene-specific malignant-lineage attribution is inferred here."]}
    save_json(results_dir / "fit_validation.json", report)
    print(json.dumps({"completed": True, "eligible_variants": report["eligible_variants"], "fitted_groups": report["fitted_patient_gene_groups"], "zero_boundary_groups": report["zero_boundary_groups"], "diagnosis_exclusion_records": len(excluded), "gene_summary": summaries}), flush=True)


if __name__ == "__main__":
    main()
