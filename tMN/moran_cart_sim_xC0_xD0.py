#!/usr/bin/env python3
"""
Moran process (fixed N) of progenitor/stem cells after CAR T therapy.
One fixed gamma selection pulse; scan over initial levels X_C0 x X_D0.

Types:  W = wild type, C = single-hit mutant, D = double-hit mutant
Fitness before CAR T: f_W = f_C = f_D = 1
Fitness after  CAR T: f_W = 1, f_C = 1 + s_C(t), f_D = 1 + s_C(t) + s_D
    s_C(t) = gamma density g(t; ALPHA, THETA), with t and THETA in days,
    using area_days = 1 day, so int s_C(t) dt = 1 day, matching the fitted model.
Mutation: a dividing C cell produces a D daughter with probability mu.

Numerics: tau-leaping of Moran birth-death events (see step()).

Outputs:
    ../result/<run>.npz                  full trajectories (fractions), s_C(t), time
    ../result/<run>_summary.csv          mean / median / 5% / 95% over simulations
    ../result/params.json                all parameters
    ../result/tmn_prob_xC0_xD0.csv       1-year tMN probability per (x_C0, x_D0)
    ../figure/mutant_xC0_xD0.png         C trajectories, rows x_C0, columns x_D0
    ../figure/double_xC0_xD0.png         D trajectories, rows x_C0, columns x_D0
    ../figure/tmn_prob_xC0_xD0.png       heatmap of 1-year tMN probability
    ../figure/tmn_vs_xC0.png             P_tMN vs x_C0, one line per x_D0
    ../figure/gamma_selection_pulse.png  the (single) s_C(t) pulse with AUC_30/60/90d

1-year tMN probability: fraction of simulations whose double-mutant fraction
exceeds TMN_THRESHOLD at any time within the first year after CAR T.
"""
import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ------------------------------------------------------------------ paths
HERE = Path(__file__).resolve().parent
RESULT_DIR = (HERE / "../result").resolve()
FIG_DIR = (HERE / "../figure").resolve()

# ------------------------------------------------------------------ parameters
N = 100_000          # population size
LAM = 1 / 28         # divisions per cell per day
MU = 0.01            # C -> D probability per C division
S_D = 0.02           # extra fitness of D (on top of s_C)
T_PRE = 0            # days of neutral dynamics before CAR T
T_POST = 365         # days simulated after CAR T
DT = 1.0             # tau-leap step (days)
N_SIM = 1000         # simulations per parameter set
TMN_THRESHOLD = 0.10 # double-mutant fraction defining tMN
T_TMN = 365          # horizon (days) for the tMN probability
ALPHA = 15.848       # gamma shape (fixed)
THETA = 1.0          # gamma scale in days (fixed)
X_C0_LIST = [0.1, 0.2, 0.3]
X_D0_LIST = [0.01, 0.005, 0.0]
T_AUC_LIST = [30, 60, 90]
SEED = 2026


# ------------------------------------------------------------------ model
def step(n, f, rng):
    """One tau-leap step. n, f: arrays (R, 3) of counts and fitness [W, C, D]."""
    w = n * f
    pb = w / w.sum(axis=1, keepdims=True)
    q = pb.copy()
    q[:, 1] = pb[:, 1] * (1 - MU)
    q[:, 2] = pb[:, 2] + pb[:, 1] * MU
    pd = n / N
    rates = LAM * N * DT * q[:, :, None] * pd[:, None, :]
    k = rng.poisson(rates)
    delta = k.sum(axis=2) - k.sum(axis=1)
    n = np.maximum(n + delta, 0)
    excess = n.sum(axis=1) - N
    n[:, 0] -= excess
    bad = n[:, 0] < 0
    if bad.any():
        n[bad, 0] = 0
        n[bad] = np.round(n[bad] * N / n[bad].sum(axis=1, keepdims=True)).astype(int)
        n[bad, 0] += N - n[bad].sum(axis=1)
    return n


def s_c_of_t(t, alpha=ALPHA, theta=THETA):
    """Gamma pulse in days with area_days=1; 0 before CAR T."""
    t = np.asarray(t, float)
    days = np.maximum(t, 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        log_g = ((alpha - 1) * np.log(days) - days / theta
                 - math.lgamma(alpha) - alpha * math.log(theta))
        g = np.where(days > 0, np.exp(log_g), 1 / theta if alpha == 1 else 0.0)
    return np.where(t >= 0, g, 0.0)


def auc_s_c(t_end, dt=0.01):
    """Integrated selection over [0, t_end] days (area_days=1 times Gamma CDF)."""
    t = np.arange(0, t_end + dt, dt)
    return float(np.trapezoid(s_c_of_t(t), t))


def simulate(x_c0, x_d0, rng):
    n_c0, n_d0 = int(round(x_c0 * N)), int(round(x_d0 * N))
    n = np.tile([N - n_c0 - n_d0, n_c0, n_d0], (N_SIM, 1)).astype(np.int64)

    t = np.arange(-T_PRE, T_POST + DT, DT)
    s_c = s_c_of_t(t)
    s_mid = s_c_of_t(t[:-1] + DT / 2)
    traj = np.empty((len(t), N_SIM, 3), dtype=np.int64)
    traj[0] = n
    for i, ti in enumerate(t[:-1], start=1):
        sd = S_D if ti >= 0 else 0.0
        f = np.tile([1.0, 1.0 + s_mid[i - 1], 1.0 + s_mid[i - 1] + sd], (N_SIM, 1))
        n = step(n, f, rng)
        traj[i] = n
    x = traj / N
    return t, x[:, :, 1].T, x[:, :, 2].T, s_c


def save_run(name, t, xc, xd, s_c):
    np.savez_compressed(RESULT_DIR / f"{name}.npz", t=t, x_C=xc, x_D=xd, s_C=s_c)
    cols = [t]
    for x in (xc, xd):
        cols += [x.mean(0), np.median(x, 0), np.quantile(x, 0.05, 0), np.quantile(x, 0.95, 0)]
    header = "t_days,C_mean,C_median,C_q05,C_q95,D_mean,D_median,D_q05,D_q95"
    np.savetxt(RESULT_DIR / f"{name}_summary.csv", np.column_stack(cols),
               delimiter=",", header=header, comments="")


def tmn_probability(t, xd):
    window = (t >= 0) & (t <= T_TMN)
    return float((xd[:, window] > TMN_THRESHOLD).any(axis=1).mean())


# ------------------------------------------------------------------ plots
def pulse_label():
    return (rf"$\alpha$={ALPHA}, $\theta$={THETA} d "
            rf"(mean pulse time {ALPHA*THETA:.1f} d)")


def plot_grid(results, probs, which, color, label):
    """Rows: x_C0, columns: x_D0."""
    fig, axes = plt.subplots(len(X_C0_LIST), len(X_D0_LIST), sharex=True,
                             figsize=(4.5 * len(X_D0_LIST), 3.2 * len(X_C0_LIST)), squeeze=False)
    for (xc0, xd0), (t, xc, xd, _) in results.items():
        ax = axes[X_C0_LIST.index(xc0), X_D0_LIST.index(xd0)]
        x = xc if which == "C" else xd
        segs = np.stack([np.broadcast_to(t, x.shape), x], axis=-1)
        ax.add_collection(LineCollection(segs, colors=color, alpha=0.02, lw=0.5))
        ax.autoscale_view()
        ax.plot(t, x.mean(0), color=color, lw=2.5, label="mean")
        if which == "D":
            ax.axhline(TMN_THRESHOLD, color="k", ls="--", lw=1.0)
        if T_PRE > 0:
            ax.axvline(0, color="k", ls="--", lw=0.8)
        ax.set_title(rf"$x_C(0)$={xc0}, $x_D(0)$={xd0}, P$_{{tMN}}$={probs[(xc0, xd0)]:.2f}",
                     fontsize=10)
        ax.set_ylabel(f"{label} fraction")
        ax.legend(loc="center left", fontsize=8, frameon=False)
    for ax in axes[-1]:
        ax.set_xlabel("days after CAR T")
    fig.suptitle(f"{label} clone, {pulse_label()}, N={N:.0e}, $s_D$={S_D}, "
                 f"$\\mu$={MU:.0e}, {N_SIM} sims")
    fig.tight_layout()
    tag = "mutant" if which == "C" else "double"
    fig.savefig(FIG_DIR / f"{tag}_xC0_xD0.png", dpi=200)
    plt.close(fig)


def save_and_plot_tmn(probs):
    P = np.array([[probs[(xc0, xd0)] for xd0 in X_D0_LIST] for xc0 in X_C0_LIST])  # rows: x_C0
    rows = [[xc0] + list(P[i]) for i, xc0 in enumerate(X_C0_LIST)]
    np.savetxt(RESULT_DIR / "tmn_prob_xC0_xD0.csv", rows, delimiter=",",
               header="x_C0\\x_D0," + ",".join(str(x) for x in X_D0_LIST), comments="")

    # heatmap
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    im = ax.imshow(P, origin="lower", cmap="Reds", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(X_D0_LIST)), [str(x) for x in X_D0_LIST])
    ax.set_yticks(range(len(X_C0_LIST)), [str(x) for x in X_C0_LIST])
    ax.set_xlabel(r"$x_D(0)$")
    ax.set_ylabel(r"$x_C(0)$")
    for i in range(len(X_C0_LIST)):
        for j in range(len(X_D0_LIST)):
            ax.text(j, i, f"{P[i, j]:.2f}", ha="center", va="center", fontsize=9,
                    color="white" if P[i, j] > 0.6 else "black")
    fig.colorbar(im, ax=ax, label=f"P($x_D$ > {TMN_THRESHOLD}) within 1 year")
    ax.set_title(f"1-year tMN probability, {pulse_label()}", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "tmn_prob_xC0_xD0.png", dpi=200)
    plt.close(fig)

    # P_tMN vs x_C0, one line per x_D0
    fig, ax = plt.subplots(figsize=(6.0, 4.5))
    colors = plt.cm.Reds(np.linspace(0.4, 1.0, len(X_D0_LIST)))
    for j, xd0 in enumerate(X_D0_LIST):
        ax.plot(X_C0_LIST, P[:, j], "o-", color=colors[j], mec="k", label=rf"$x_D(0)$={xd0}")
    ax.set_xlabel(r"$x_C(0)$")
    ax.set_ylabel("1-year tMN probability")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(f"{pulse_label()}, {N_SIM} sims", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "tmn_vs_xC0.png", dpi=200)
    plt.close(fig)


def plot_gamma_pulse():
    t = np.arange(0, T_POST + DT, DT)
    s = s_c_of_t(t)
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.fill_between(t, s, color="tab:red", alpha=0.3)
    ax.plot(t, s, color="k", lw=1.2)
    for t_auc in T_AUC_LIST:
        ax.axvline(t_auc, color="k", ls=":", lw=0.8)
    ax.set_xlabel("days after CAR T")
    ax.set_ylabel(r"$s_C(t)$")
    ax.set_title(pulse_label() + "\n"
                 + ", ".join(rf"AUC$_{{{t_auc}d}}$={auc_s_c(t_auc):.2f}" for t_auc in T_AUC_LIST),
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "gamma_selection_pulse.png", dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    params = dict(N=N, LAM=LAM, MU=MU, S_D=S_D, T_PRE=T_PRE, T_POST=T_POST, DT=DT,
                  N_SIM=N_SIM, TMN_THRESHOLD=TMN_THRESHOLD, T_TMN=T_TMN,
                  ALPHA=ALPHA, THETA=THETA, X_C0_LIST=X_C0_LIST, X_D0_LIST=X_D0_LIST,
                  T_AUC_LIST=T_AUC_LIST,
                  AUC={f"{t}d": auc_s_c(t) for t in T_AUC_LIST}, SEED=SEED)
    (RESULT_DIR / "params.json").write_text(json.dumps(params, indent=2))

    rng = np.random.default_rng(SEED)
    results, probs = {}, {}
    for xc0 in X_C0_LIST:
        for xd0 in X_D0_LIST:
            t, xc, xd, s_c = simulate(xc0, xd0, rng)
            save_run(f"alpha{ALPHA}_theta{THETA}_xC0-{xc0}_xD0-{xd0}", t, xc, xd, s_c)
            results[(xc0, xd0)] = (t, xc, xd, s_c)
            probs[(xc0, xd0)] = tmn_probability(t, xd)
            print(f"x_C0={xc0} x_D0={xd0}: mean x_C(T)={xc[:, -1].mean():.3f}, "
                  f"mean x_D(T)={xd[:, -1].mean():.4f}, P_tMN(1y)={probs[(xc0, xd0)]:.3f}")

    plot_grid(results, probs, "C", "tab:blue", "Mutant")
    plot_grid(results, probs, "D", "tab:red", "Double-hit")
    save_and_plot_tmn(probs)
    plot_gamma_pulse()


if __name__ == "__main__":
    main()
