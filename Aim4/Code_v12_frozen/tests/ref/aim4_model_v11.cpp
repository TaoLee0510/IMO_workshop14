// =============================================================================
// aim4_model_v11.cpp  (frozen copy of the v1.1 engine, 2026-11-03 13:40; used only by the K1 regression test - do not edit)
// Purpose : Parent marrow model for Aim 4 (who fails to recover - reserve or
//           clone - and when the stem cell pool is hit - chemotherapy (a) or
//           CAR T phase (b)), plus the original Friberg (2002) myelosuppression
//           model as reference model MF. One adaptive Dormand-Prince RK45
//           integrator, integrated piecewise between input breakpoints.
// Inputs  : parameter vector (layout = a4_param_names() / mf_param_names()),
//           a segment matrix of piecewise-constant inputs, the CAR T-phase
//           driver D(t) on a uniform grid, and observation times (days).
// Outputs : matrix of model quantities at the observation times.
//
// Parent model (every compartment normalised to its pre-treatment value = 1)
//   Exposure  C          dC/dt = dose(t) - k_e*C        (bridging + flu-cy)
//   Marrow    I1, I      dI1/dt = k_I*(D - I1), dI/dt = k_I*(I1 - I)
//   Stem      S_W, S_X   share one niche, S = S_W + S_X
//     dS_j/dt = r_S*S_j*(1-S)*(1 + alpha*(1-M)+)
//               - [ra_j*th_aS*C + rb_j*th_bS*I*(1 + eta*(1-S)+)]*S_j
//     ra_W = rb_W = 1, ra_X = r_a, rb_X = r_b;  S_W(0) = 1 - f0, S_X(0) = f0
//   MPP       M          dM/dt = k_M*(1-phi)*F + k_M*(phi*G*(1 - E(I)) - 1)*M - th_aP*C*M
//                        F = (S_W + delta*S_X)/(1 - f0 + delta*f0)
//                        G = min(gmax, Bbar^-gamma_M), Bbar = geometric mean of blood
//   Lineage l P -> T1 -> T2 -> T3 -> B   (l = neutrophil, platelet, monocyte)
//     dP/dt = k_l*[M*B^-gam_l*(1 - e_l*E(I))*u_l(t) - P] - th_aP*C*P
//     dB/dt = c_l*(T3 - B),  E(I) = Emax*I^h/(I^h + IC50^h),  e_N = 1
//     u_N = G-CSF fold, u_Pl = TPO-RA fold while on (simulation only)
//   Segment column s_add: S_W += s_add at the segment start (CD34+ boost).
//
// MF (Friberg et al., J Clin Oncol 2002), one independent chain per lineage
//   dProl/dt = k*Prol*(1 - E)*B^-gam - k*Prol,   E = slope_l*C
//   dT_j/dt  = k*(T_{j-1} - T_j),  dB/dt = k*(T3 - B)   (kprol = kcirc = ktr = k)
//   k = ktr_l / mtt_scale
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace a4v11 {

// ---- parameter layouts (must match a4_param_names / mf_param_names) ------------
enum Par {
  R_S, F0, TH_AS, TH_AP, TH_BS,               // individual
  ALPHA, K_M, PHI, GMAX, GAMMA_M,             // stem / MPP layer
  R_A, R_B, DELTA,                            // clone traits
  KTR_N, KTR_PL, KTR_MO,                      // lineage transit rates (= 4 / MTT)
  C_N, C_PL, C_MO,                            // blood turnover
  GAM_N, GAM_PL, GAM_MO,                      // lineage feedback
  K_I, EMAX, IC50, HILL, E_PL, E_MO, ETA_SYN, // marrow inflammation, a x b synergy
  K_E, GCSF_FOLD, TPO_FOLD,                   // exposure decay, supportive care
  N_PAR
};
enum State { SW, SX, M, CX, I1, I2,
             PN, TN1, TN2, TN3, BN,
             PP, TP1, TP2, TP3, BP,
             PM, TM1, TM2, TM3, BM,
             N_STATE };

enum ParMF { MF_SLOPE, MF_MTT_SCALE, MF_KTR_N, MF_KTR_PL, MF_KTR_MO,
             MF_GAM_N, MF_GAM_PL, MF_GAM_MO, MF_SLOPE_PL, MF_SLOPE_MO, MF_K_E, MF_N_PAR };
enum StateMF { F_CX, F_PN, F_TN1, F_TN2, F_TN3, F_BN, F_PP, F_TP1, F_TP2, F_TP3, F_BP,
               F_PM, F_TM1, F_TM2, F_TM3, F_BM, MF_N_STATE };

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

// piecewise-constant inputs (one row of the segment matrix)
struct Inputs { double dose, gcsf, tpo; };

// ---- parent model -----------------------------------------------------------------
struct Parent {
  const double* p;
  const Driver* drv;
  Inputs u;
  int n_state() const { return N_STATE; }

  void lineage(const double* y, double* dy, int p0, double k, double c, double gam,
               double m, double factor, double kill) const {
    const double b = std::max(y[p0 + 4], 1e-6);
    dy[p0]     = k * (m * std::pow(b, -gam) * factor - y[p0]) - kill * y[p0];
    dy[p0 + 1] = k * (y[p0] - y[p0 + 1]);
    dy[p0 + 2] = k * (y[p0 + 1] - y[p0 + 2]);
    dy[p0 + 3] = k * (y[p0 + 2] - y[p0 + 3]);
    dy[p0 + 4] = c * (y[p0 + 3] - y[p0 + 4]);
  }

  void rhs(double t, const double* y, double* dy) const {
    const double sw = std::max(y[SW], 0.0), sx = std::max(y[SX], 0.0), s = sw + sx;
    const double m = std::max(y[M], 0.0), cx = std::max(y[CX], 0.0);
    const double i1 = std::max(y[I1], 0.0), im = std::max(y[I2], 0.0);
    const double d = std::max(drv->at(t), 0.0);

    const double bn = std::max(y[BN], 1e-6), bp = std::max(y[BP], 1e-6), bm = std::max(y[BM], 1e-6);
    const double g = std::min(p[GMAX], std::pow(std::cbrt(bn * bp * bm), -p[GAMMA_M]));

    // stem layer: shared-niche regrowth (faster when MPP is depleted), timed injuries
    const double regrow = p[R_S] * (1.0 - s) * (1.0 + p[ALPHA] * std::max(1.0 - m, 0.0));
    const double kill_a = p[TH_AS] * cx;
    const double kill_b = p[TH_BS] * im * (1.0 + p[ETA_SYN] * std::max(1.0 - s, 0.0));
    dy[SW] = regrow * sw - (kill_a + kill_b) * sw;
    dy[SX] = regrow * sx - (p[R_A] * kill_a + p[R_B] * kill_b) * sx;

    // marrow inflammation: reversible suppression of MPP self-renewal and of production
    const double ih = std::pow(im, p[HILL]);
    const double e = p[EMAX] * ih / (ih + std::pow(p[IC50], p[HILL]));

    // MPP: input from stem output, feedback-driven transient self-renewal
    const double f0n = 1.0 - p[F0] + p[DELTA] * p[F0];
    const double F = (sw + p[DELTA] * sx) / f0n;
    const double kill_p = p[TH_AP] * cx;
    dy[M] = p[K_M] * (1.0 - p[PHI]) * F + p[K_M] * (p[PHI] * g * (1.0 - e) - 1.0) * m - kill_p * m;

    // exposure and the two-stage marrow lag of the CAR T-phase driver
    dy[CX] = u.dose - p[K_E] * cx;
    dy[I1] = p[K_I] * (d - i1);
    dy[I2] = p[K_I] * (i1 - im);

    // lineages
    const double gf = (u.gcsf > 0.0) ? p[GCSF_FOLD] : 1.0;
    const double tf = (u.tpo > 0.0) ? p[TPO_FOLD] : 1.0;
    lineage(y, dy, PN, p[KTR_N],  p[C_N],  p[GAM_N],  m, (1.0 - e) * gf, kill_p);
    lineage(y, dy, PP, p[KTR_PL], p[C_PL], p[GAM_PL], m, (1.0 - p[E_PL] * e) * tf, kill_p);
    lineage(y, dy, PM, p[KTR_MO], p[C_MO], p[GAM_MO], m, (1.0 - p[E_MO] * e), kill_p);
  }

  void init(std::vector<double>& y) const {
    std::fill(y.begin(), y.end(), 1.0);
    y[SW] = 1.0 - p[F0];
    y[SX] = p[F0];
    y[CX] = 0.0;
    y[I1] = 0.0;
    y[I2] = 0.0;
  }
};

// ---- Friberg reference model --------------------------------------------------------
struct Friberg {
  const double* p;
  const Driver* drv;
  Inputs u;
  int n_state() const { return MF_N_STATE; }

  void chain(const double* y, double* dy, int p0, double k, double gam, double e) const {
    const double b = std::max(y[p0 + 4], 1e-6);
    dy[p0]     = k * y[p0] * (1.0 - e) * std::pow(b, -gam) - k * y[p0];
    dy[p0 + 1] = k * (y[p0] - y[p0 + 1]);
    dy[p0 + 2] = k * (y[p0 + 1] - y[p0 + 2]);
    dy[p0 + 3] = k * (y[p0 + 2] - y[p0 + 3]);
    dy[p0 + 4] = k * (y[p0 + 3] - y[p0 + 4]);
  }

  void rhs(double, const double* y, double* dy) const {
    const double cx = std::max(y[F_CX], 0.0);
    const double sc = p[MF_MTT_SCALE];
    dy[F_CX] = u.dose - p[MF_K_E] * cx;
    chain(y, dy, F_PN, p[MF_KTR_N] / sc,  p[MF_GAM_N],  p[MF_SLOPE] * cx);
    chain(y, dy, F_PP, p[MF_KTR_PL] / sc, p[MF_GAM_PL], p[MF_SLOPE] * p[MF_SLOPE_PL] * cx);
    chain(y, dy, F_PM, p[MF_KTR_MO] / sc, p[MF_GAM_MO], p[MF_SLOPE] * p[MF_SLOPE_MO] * cx);
  }

  void init(std::vector<double>& y) const {
    std::fill(y.begin(), y.end(), 1.0);
    y[F_CX] = 0.0;
  }
};

// ---- Dormand-Prince RK45 --------------------------------------------------------
struct Ctl { double rtol, atol, hmax; long max_steps; long steps; bool failed; };

template <class Sys>
void integrate_to(const Sys& sys, double& t, double t_end, std::vector<double>& y, double& h, Ctl& ctl) {
  static const double c2 = 1.0/5, c3 = 3.0/10, c4 = 4.0/5, c5 = 8.0/9;
  static const double a21 = 1.0/5;
  static const double a31 = 3.0/40, a32 = 9.0/40;
  static const double a41 = 44.0/45, a42 = -56.0/15, a43 = 32.0/9;
  static const double a51 = 19372.0/6561, a52 = -25360.0/2187, a53 = 64448.0/6561, a54 = -212.0/729;
  static const double a61 = 9017.0/3168, a62 = -355.0/33, a63 = 46732.0/5247, a64 = 49.0/176, a65 = -5103.0/18656;
  static const double b1 = 35.0/384, b3 = 500.0/1113, b4 = 125.0/192, b5 = -2187.0/6784, b6 = 11.0/84;
  static const double e1 = 71.0/57600, e3 = -71.0/16695, e4 = 71.0/1920, e5 = -17253.0/339200,
                      e6 = 22.0/525, e7 = -1.0/40;

  const int n = sys.n_state();
  std::vector<double> k1(n), k2(n), k3(n), k4(n), k5(n), k6(n), k7(n), yt(n), yn(n);

  while (t < t_end - 1e-12) {
    if (ctl.steps++ > ctl.max_steps) { ctl.failed = true; return; }
    double hh = std::min(std::min(h, ctl.hmax), t_end - t);

    sys.rhs(t, y.data(), k1.data());
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * a21 * k1[i];
    sys.rhs(t + c2 * hh, yt.data(), k2.data());
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
    sys.rhs(t + c3 * hh, yt.data(), k3.data());
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    sys.rhs(t + c4 * hh, yt.data(), k4.data());
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    sys.rhs(t + c5 * hh, yt.data(), k5.data());
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    sys.rhs(t + hh, yt.data(), k6.data());
    for (int i = 0; i < n; ++i) yn[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    sys.rhs(t + hh, yn.data(), k7.data());

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

// Piecewise integration over the segments, recording at obs_t.
// record(k, t, y) writes one output row.
template <class Sys, class Rec>
bool run_segments(Sys& sys, const NumericMatrix& seg, const NumericVector& obs_t, Ctl& ctl, Rec record) {
  const int nobs = obs_t.size();
  std::vector<double> y(sys.n_state());
  sys.init(y);
  double t = seg(0, 0);
  double h = 0.05;
  int k = 0;
  while (k < nobs && obs_t[k] <= t + 1e-12) { record(k, obs_t[k], y); ++k; }
  for (int s = 0; s < seg.nrow() && !ctl.failed; ++s) {
    sys.u = Inputs{seg(s, 2), seg(s, 3), seg(s, 4)};
    if (seg(s, 5) != 0.0) y[0] += seg(s, 5);  // CD34+ boost into the first (normal stem) state
    const double t_end = seg(s, 1);
    while (k < nobs && obs_t[k] <= t_end + 1e-12 && !ctl.failed) {
      integrate_to(sys, t, obs_t[k], y, h, ctl);
      if (!ctl.failed) { record(k, obs_t[k], y); ++k; }
    }
    if (!ctl.failed) integrate_to(sys, t, t_end, y, h, ctl);
  }
  return !ctl.failed && k == nobs;
}

void check_inputs(const NumericMatrix& seg, const NumericVector& drv) {
  if (seg.ncol() != 6) stop("segment matrix must have 6 columns: start, end, dose, gcsf, tpo, s_add");
  if (seg.nrow() < 1) stop("segment matrix is empty");
  if (drv.size() < 2) stop("driver grid needs at least two points");
}

} // namespace a4

// [[Rcpp::export]]
CharacterVector v11_param_names() {
  return CharacterVector::create(
    "r_s", "f0", "th_as", "th_ap", "th_bs",
    "alpha", "k_m", "phi", "gmax", "gamma_m",
    "r_a", "r_b", "delta",
    "ktr_n", "ktr_pl", "ktr_mo",
    "c_n", "c_pl", "c_mo",
    "gam_n", "gam_pl", "gam_mo",
    "k_i", "emax", "ic50", "hill", "e_pl", "e_mo", "eta_syn",
    "k_e", "gcsf_fold", "tpo_fold");
}

// [[Rcpp::export]]
CharacterVector v11_state_names() {
  return CharacterVector::create("SW", "SX", "M", "C", "I1", "I",
                                 "PN", "TN1", "TN2", "TN3", "BN",
                                 "PP", "TP1", "TP2", "TP3", "BP",
                                 "PM", "TM1", "TM2", "TM3", "BM");
}

// [[Rcpp::export]]
CharacterVector v11_mf_param_names() {
  return CharacterVector::create("slope", "mtt_scale", "ktr_n", "ktr_pl", "ktr_mo",
                                 "gam_n", "gam_pl", "gam_mo", "slope_pl", "slope_mo", "k_e");
}

// Simulate one patient with the parent model.
// seg   : matrix with columns start, end, dose, gcsf, tpo, s_add (rows ordered in time)
// drv   : D(t) on the grid drv_t0 + (0:(n-1)) * drv_dt
// obs_t : ascending observation times (>= first segment start)
// value : columns t, anc_rel, plt_rel, amc_rel, S, f_clone, M, C, D, I
//         (all NA if the integrator failed)
// [[Rcpp::export]]
NumericMatrix v11_a4_sim(NumericVector par, NumericMatrix seg, NumericVector drv, double drv_t0,
                     double drv_dt, NumericVector obs_t, double rtol = 1e-6, double atol = 1e-9,
                     double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != a4v11::N_PAR) stop("parameter vector has wrong length");
  a4v11::check_inputs(seg, drv);
  a4v11::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  a4v11::Parent sys{par.begin(), &dr, a4v11::Inputs{0, 0, 0}};
  a4v11::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  NumericMatrix out(obs_t.size(), 10);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                          "M", "C", "D", "I");
  auto record = [&](int k, double t, const std::vector<double>& y) {
    const double s = y[a4v11::SW] + y[a4v11::SX];
    out(k, 0) = t;
    out(k, 1) = y[a4v11::BN];
    out(k, 2) = y[a4v11::BP];
    out(k, 3) = y[a4v11::BM];
    out(k, 4) = s;
    out(k, 5) = (s > 1e-12) ? y[a4v11::SX] / s : 0.0;
    out(k, 6) = y[a4v11::M];
    out(k, 7) = y[a4v11::CX];
    out(k, 8) = dr.at(t);
    out(k, 9) = y[a4v11::I2];
  };
  if (!a4v11::run_segments(sys, seg, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Simulate one patient with the Friberg reference model (same output columns;
// stem-layer columns are NA).
// [[Rcpp::export]]
NumericMatrix v11_mf_sim(NumericVector par, NumericMatrix seg, NumericVector obs_t, double rtol = 1e-6,
                     double atol = 1e-9, double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != a4v11::MF_N_PAR) stop("MF parameter vector has wrong length");
  NumericVector dummy = NumericVector::create(0.0, 0.0);
  a4v11::check_inputs(seg, dummy);
  a4v11::Driver dr{0.0, 1.0, dummy.begin(), 2};
  a4v11::Friberg sys{par.begin(), &dr, a4v11::Inputs{0, 0, 0}};
  a4v11::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  NumericMatrix out(obs_t.size(), 10);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                          "M", "C", "D", "I");
  auto record = [&](int k, double t, const std::vector<double>& y) {
    out(k, 0) = t;
    out(k, 1) = y[a4v11::F_BN];
    out(k, 2) = y[a4v11::F_BP];
    out(k, 3) = y[a4v11::F_BM];
    out(k, 4) = NA_REAL;
    out(k, 5) = NA_REAL;
    out(k, 6) = y[a4v11::F_PN];
    out(k, 7) = y[a4v11::F_CX];
    out(k, 8) = 0.0;
    out(k, 9) = 0.0;
  };
  // the CD34 column has no meaning in MF
  NumericMatrix seg2 = clone(seg);
  for (int i = 0; i < seg2.nrow(); ++i) seg2(i, 5) = 0.0;
  if (!a4v11::run_segments(sys, seg2, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Right-hand side of the parent model, exposed for validation against deSolve.
// [[Rcpp::export]]
NumericVector v11_a4_rhs(double t, NumericVector y, NumericVector par, double dose, double gcsf, double tpo,
                     NumericVector drv, double drv_t0, double drv_dt) {
  if (par.size() != a4v11::N_PAR) stop("parameter vector has wrong length");
  if (y.size() != a4v11::N_STATE) stop("state vector has wrong length");
  a4v11::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  a4v11::Parent sys{par.begin(), &dr, a4v11::Inputs{dose, gcsf, tpo}};
  NumericVector dy(a4v11::N_STATE);
  sys.rhs(t, y.begin(), dy.begin());
  return dy;
}
