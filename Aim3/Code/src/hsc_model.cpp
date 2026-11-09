// =============================================================================
// hsc_model.cpp
// Purpose : Multi-lineage marrow model for Aim 3 (when is the stem cell pool
//           injured: chemotherapy (a) vs CAR T phase (b)). Right-hand side of
//           the ODE system plus an adaptive Dormand-Prince RK45 integrator,
//           integrated piecewise between dosing breakpoints.
// Inputs  : parameter vector (layout = hsc_param_names()), a segment matrix of
//           piecewise-constant chemotherapy dose rates, the CAR T-phase driver
//           D(t) on a uniform grid, and observation times (days).
// Outputs : matrix of model quantities at the observation times.
//
// Model (every compartment normalised to its pre-treatment value = 1)
//   Exposure   C        dC/dt = dose(t) - k_e*C   (bridging + flu-cy)
//   HSC        S        dS/dt = r_S*S*(1-S)*(1 + alpha*(1-M)+)
//                               - [th_aS*C + th_bS*D*(1 + eta*(1-S)+)]*S
//   MPP        M        dM/dt = k_M*(1-phi)*S + k_M*(phi*G - 1)*M
//                               - rho_M*[th_aP*C + q_b*D]*M      (MPP partly spared)
//                       G = min(gmax, Bbar^-gamma_M), Bbar = geometric mean of blood
//   Lineage l  P -> T1 -> T2 -> T3 -> B   (l = neutrophil, platelet, monocyte)
//              accelerated release: transit rate k_l*A_l, A_l = 1 + w_l*eps*D
//              (inflammation empties the maturation chain early: surge, then gap)
//              dP/dt  = k_l*M*B^-gam_l - k_l*A_l*P - [th_aP*C + q_b*D]*P
//              dTj/dt = k_l*A_l*(T(j-1) - Tj)
//              dB/dt  = c_l*(A_l*T3 - B)
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace hm {

// ---- parameter layout (must match hsc_param_names) ----------------------------
enum Par {
  TH_AP, TH_AS, TH_BS, R_S,                 // individual
  ALPHA, K_M, PHI, GMAX, GAMMA_M,           // stem / MPP layer
  KTR_N, KTR_PL, KTR_MO,                    // lineage transit rates (= 4 / MTT)
  C_N, C_PL, C_MO,                          // blood turnover
  GAM_N, GAM_PL, GAM_MO,                    // lineage feedback
  K_E, Q_B, ETA_SYN,                        // exposure decay, inflammatory suppression, a x b synergy
  RHO_M,                                    // MPP sensitivity relative to committed progenitors
  EPS, W_PL, W_MO,                          // accelerated release (individual) and lineage weights
  N_PAR
};

// ---- state layout -------------------------------------------------------------
enum State { S, M, CX,
             PN, TN1, TN2, TN3, BN,
             PP, TP1, TP2, TP3, BP,
             PM, TM1, TM2, TM3, BM,
             N_STATE };

struct Driver {
  double t0, dt;
  const double* v;
  int n;
  double at(double t) const {
    const double x = (t - t0) / dt;
    if (x <= 0.0) return v[0];
    if (x >= n - 1) return v[n - 1];
    const int i = static_cast<int>(x);
    const double w = x - i;
    return (1.0 - w) * v[i] + w * v[i + 1];
  }
};

inline void lineage(const double* y, double* dy, int p0, double k, double c, double gam,
                    double m, double kill, double acc) {
  const double b = std::max(y[p0 + 4], 1e-6);
  const double ka = k * acc;
  dy[p0]     = k * m * std::pow(b, -gam) - ka * y[p0] - kill * y[p0];
  dy[p0 + 1] = ka * (y[p0] - y[p0 + 1]);
  dy[p0 + 2] = ka * (y[p0 + 1] - y[p0 + 2]);
  dy[p0 + 3] = ka * (y[p0 + 2] - y[p0 + 3]);
  dy[p0 + 4] = c * (acc * y[p0 + 3] - y[p0 + 4]);
}

void rhs(double t, const double* y, double* dy, const double* p, double dose, const Driver& drv) {
  const double s = std::max(y[S], 0.0), m = std::max(y[M], 0.0), cx = std::max(y[CX], 0.0);
  const double d = std::max(drv.at(t), 0.0);

  const double bn = std::max(y[BN], 1e-6), bp = std::max(y[BP], 1e-6), bm = std::max(y[BM], 1e-6);
  const double bbar = std::cbrt(bn * bp * bm);
  const double g = std::min(p[GMAX], std::pow(bbar, -p[GAMMA_M]));

  const double stress = 1.0 + p[ALPHA] * std::max(1.0 - m, 0.0);
  const double kill_s = p[TH_AS] * cx + p[TH_BS] * d * (1.0 + p[ETA_SYN] * std::max(1.0 - s, 0.0));
  const double kill_p = p[TH_AP] * cx + p[Q_B] * d;

  dy[S]  = p[R_S] * s * (1.0 - s) * stress - kill_s * s;
  dy[M]  = p[K_M] * (1.0 - p[PHI]) * s + p[K_M] * (p[PHI] * g - 1.0) * m - p[RHO_M] * kill_p * m;
  dy[CX] = dose - p[K_E] * cx;

  const double rel = p[EPS] * d;
  lineage(y, dy, PN, p[KTR_N],  p[C_N],  p[GAM_N],  m, kill_p, 1.0 + rel);
  lineage(y, dy, PP, p[KTR_PL], p[C_PL], p[GAM_PL], m, kill_p, 1.0 + p[W_PL] * rel);
  lineage(y, dy, PM, p[KTR_MO], p[C_MO], p[GAM_MO], m, kill_p, 1.0 + p[W_MO] * rel);
}

// ---- Dormand-Prince RK45 --------------------------------------------------------
struct Ctl { double rtol, atol, hmax; long max_steps; long steps; bool failed; };

void integrate_to(double& t, double t_end, std::vector<double>& y, const double* p,
                  double dose, const Driver& drv, double& h, Ctl& ctl) {
  static const double c2 = 1.0/5, c3 = 3.0/10, c4 = 4.0/5, c5 = 8.0/9;
  static const double a21 = 1.0/5;
  static const double a31 = 3.0/40, a32 = 9.0/40;
  static const double a41 = 44.0/45, a42 = -56.0/15, a43 = 32.0/9;
  static const double a51 = 19372.0/6561, a52 = -25360.0/2187, a53 = 64448.0/6561, a54 = -212.0/729;
  static const double a61 = 9017.0/3168, a62 = -355.0/33, a63 = 46732.0/5247, a64 = 49.0/176, a65 = -5103.0/18656;
  static const double b1 = 35.0/384, b3 = 500.0/1113, b4 = 125.0/192, b5 = -2187.0/6784, b6 = 11.0/84;
  static const double e1 = 71.0/57600, e3 = -71.0/16695, e4 = 71.0/1920, e5 = -17253.0/339200,
                      e6 = 22.0/525, e7 = -1.0/40;

  const int n = N_STATE;
  double k1[N_STATE], k2[N_STATE], k3[N_STATE], k4[N_STATE], k5[N_STATE], k6[N_STATE], k7[N_STATE];
  double yt[N_STATE], yn[N_STATE];

  while (t < t_end - 1e-12) {
    if (ctl.steps++ > ctl.max_steps) { ctl.failed = true; return; }
    double hh = std::min(std::min(h, ctl.hmax), t_end - t);

    rhs(t, y.data(), k1, p, dose, drv);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * a21 * k1[i];
    rhs(t + c2 * hh, yt, k2, p, dose, drv);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
    rhs(t + c3 * hh, yt, k3, p, dose, drv);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    rhs(t + c4 * hh, yt, k4, p, dose, drv);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    rhs(t + c5 * hh, yt, k5, p, dose, drv);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    rhs(t + hh, yt, k6, p, dose, drv);
    for (int i = 0; i < n; ++i) yn[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    rhs(t + hh, yn, k7, p, dose, drv);

    double err = 0.0;
    for (int i = 0; i < n; ++i) {
      const double ei = hh * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
      const double sc = ctl.atol + ctl.rtol * std::max(std::fabs(y[i]), std::fabs(yn[i]));
      err += (ei / sc) * (ei / sc);
    }
    err = std::sqrt(err / n);

    if (!std::isfinite(err)) { ctl.failed = true; return; }
    if (err <= 1.0) {
      t += hh;
      for (int i = 0; i < n; ++i) y[i] = yn[i];
    }
    const double fac = (err == 0.0) ? 5.0 : std::min(5.0, std::max(0.2, 0.9 * std::pow(err, -0.2)));
    h = std::max(hh * fac, 1e-8);
  }
}

void initial_state(std::vector<double>& y) {
  std::fill(y.begin(), y.end(), 1.0);
  y[CX] = 0.0;
}

} // namespace hm

// [[Rcpp::export]]
CharacterVector hsc_param_names() {
  return CharacterVector::create(
    "th_ap", "th_as", "th_bs", "r_s",
    "alpha", "k_m", "phi", "gmax", "gamma_m",
    "ktr_n", "ktr_pl", "ktr_mo",
    "c_n", "c_pl", "c_mo",
    "gam_n", "gam_pl", "gam_mo",
    "k_e", "q_b", "eta_syn", "rho_m",
    "eps", "w_pl", "w_mo");
}

// [[Rcpp::export]]
CharacterVector hsc_state_names() {
  return CharacterVector::create("S", "M", "C",
                                 "PN", "TN1", "TN2", "TN3", "BN",
                                 "PP", "TP1", "TP2", "TP3", "BP",
                                 "PM", "TM1", "TM2", "TM3", "BM");
}

// Simulate one patient.
// seg   : matrix with columns start, end, dose (rows ordered in time; first start = t0)
// drv   : D(t) on the grid drv_t0 + (0:(n-1)) * drv_dt
// obs_t : ascending observation times (>= first segment start)
// value : columns t, anc_rel, plt_rel, amc_rel, S, M, C, D (all NaN if the integrator failed)
// [[Rcpp::export]]
NumericMatrix hsc_sim(NumericVector par, NumericMatrix seg, NumericVector drv, double drv_t0,
                      double drv_dt, NumericVector obs_t, double rtol = 1e-6, double atol = 1e-9,
                      double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != hm::N_PAR) stop("parameter vector has wrong length");
  if (seg.ncol() != 3) stop("segment matrix must have 3 columns");
  if (drv.size() < 2) stop("driver grid needs at least two points");
  const double* p = par.begin();
  const int nobs = obs_t.size();
  NumericMatrix out(nobs, 8);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "M", "C", "D");

  hm::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  std::vector<double> y(hm::N_STATE);
  hm::initial_state(y);
  hm::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};

  auto record = [&](int k, double t) {
    out(k, 0) = t;
    out(k, 1) = y[hm::BN];
    out(k, 2) = y[hm::BP];
    out(k, 3) = y[hm::BM];
    out(k, 4) = y[hm::S];
    out(k, 5) = y[hm::M];
    out(k, 6) = y[hm::CX];
    out(k, 7) = dr.at(t);
  };

  double t = seg(0, 0);
  double h = 0.05;
  int k = 0;
  while (k < nobs && obs_t[k] <= t + 1e-12) { record(k, obs_t[k]); ++k; }

  for (int s = 0; s < seg.nrow() && !ctl.failed; ++s) {
    const double dose = seg(s, 2);
    const double t_end = seg(s, 1);
    while (k < nobs && obs_t[k] <= t_end + 1e-12 && !ctl.failed) {
      hm::integrate_to(t, obs_t[k], y, p, dose, dr, h, ctl);
      if (!ctl.failed) { record(k, obs_t[k]); ++k; }
    }
    if (!ctl.failed) hm::integrate_to(t, t_end, y, p, dose, dr, h, ctl);
  }

  if (ctl.failed || k < nobs) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Right-hand side exposed for validation against deSolve.
// [[Rcpp::export]]
NumericVector hsc_rhs(double t, NumericVector y, NumericVector par, double dose, NumericVector drv,
                      double drv_t0, double drv_dt) {
  if (par.size() != hm::N_PAR) stop("parameter vector has wrong length");
  if (y.size() != hm::N_STATE) stop("state vector has wrong length");
  hm::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  NumericVector dy(hm::N_STATE);
  hm::rhs(t, y.begin(), dy.begin(), par.begin(), dose, dr);
  return dy;
}
