// =============================================================================
// marrow_model.cpp
// Purpose : Marrow digital twin for post-CAR T cytopenia (Aim 2).
//           Right-hand side of the ODE system plus an adaptive Dormand-Prince
//           RK45 integrator, integrated piecewise between input breakpoints.
// Inputs  : parameter vector (layout = marrow_param_names()), a segment matrix
//           of piecewise-constant inputs, and observation times (days).
// Outputs : matrix of model quantities at the observation times.
//
// Model (all lineage compartments normalised to their pre-treatment value = 1)
//   Stem layer   xW, xM  share one niche (capacity 1); logistic regrowth a*rho
//   LD exposure  C      dC/dt = ld - kel*C ; kills stem cells and progenitors
//   Inflammation Is (systemic, drives CRP) -> Im1 -> Im (marrow; two-stage
//                lag). Im suppresses progenitor production (reversible) and
//                causes stem attrition (selective, clone-specific).
//   Neutrophils  PN -> TN1 -> TN2 -> TN3 -> N    feedback N^-gammaN
//   Platelets    PP -> TP1 -> TP2 -> TP3 -> PL   feedback PL^-gammaP
//   Stem output  F = (xW + delta*xM) / F0, F0 = (1 - f0) + delta*f0
// =============================================================================

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

namespace mm {

// ---- parameter layout (must match marrow_param_names) -----------------------
enum Par {
  RHO, G, F0,                       // individual: reserve, inflammatory gain, clone fraction
  A, KEL_LD, KLDS, KLDP,            // stem regrowth, LD elimination, LD kill (stem, progenitor)
  R_LD, R_INF, DELTA,               // clone traits
  Q, EMAX, IC50, HILL,              // marrow inflammation: stem attrition, progenitor suppression
  K_OUT, K_M, T_PK, H_SHAPE,        // inflammation kinetics and CAR T drive kernel
  KTR_N, KCIRC_N, GAMMA_N,          // neutrophil lineage
  KTR_P, KCIRC_P, GAMMA_P, E_PLT,   // platelet lineage
  GCSF_FOLD, TPO_FOLD,              // supportive care
  TOCI_RED, STER_RED, STER_BOOST,   // anti-inflammatory treatment
  N_PAR
};

// ---- state layout -------------------------------------------------------------
enum State { XW, XM, CLD, IS, IM1, IM, PN, TN1, TN2, TN3, NN, PP, TP1, TP2, TP3, PL, N_STATE };

// ---- piecewise-constant inputs (one row of the segment matrix) ----------------
struct Inputs { double ld, gcsf, toci, ster, tpo; };

inline double car_drive(double t, const double* p) {
  if (t <= 0.0) return 0.0;
  const double s = t / p[T_PK];
  return std::pow(s, p[H_SHAPE]) * std::exp(p[H_SHAPE] * (1.0 - s));
}

inline double clone_f0_norm(const double* p) {
  return (1.0 - p[F0]) + p[DELTA] * p[F0];
}

void rhs(double t, const double* y, double* dy, const double* p, const Inputs& u, double f0n) {
  const double xw = std::max(y[XW], 0.0), xm = std::max(y[XM], 0.0);
  const double xt = xw + xm;
  const double c  = std::max(y[CLD], 0.0);
  const double is = std::max(y[IS], 0.0), im1 = std::max(y[IM1], 0.0), im = std::max(y[IM], 0.0);

  // stem layer: shared-niche logistic regrowth, LD kill, inflammatory attrition
  const double regrow = p[A] * p[RHO] * (1.0 - xt);
  dy[XW] = regrow * xw - p[KLDS] * c * xw - p[Q] * im * xw;
  dy[XM] = regrow * xm - p[R_LD] * p[KLDS] * c * xm - p[R_INF] * p[Q] * im * xm;

  // lymphodepletion exposure
  dy[CLD] = u.ld - p[KEL_LD] * c;

  // inflammation: CAR T driven systemic signal and two-stage lagged marrow signal
  const double drive = p[G] * car_drive(t, p) *
                       (1.0 - p[TOCI_RED] * u.toci) * (1.0 - p[STER_RED] * u.ster);
  dy[IS] = drive - p[K_OUT] * (1.0 + p[STER_BOOST] * u.ster) * is;
  dy[IM1] = p[K_M] * (is - im1);
  dy[IM]  = p[K_M] * (im1 - im);

  // lineages
  const double F    = (xw + p[DELTA] * xm) / f0n;
  const double imh  = std::pow(im, p[HILL]);
  const double eac  = p[EMAX] * imh / (imh + std::pow(p[IC50], p[HILL]));
  const double n    = std::max(y[NN], 1e-6), pl = std::max(y[PL], 1e-6);
  const double fbN  = std::pow(n, -p[GAMMA_N]);
  const double fbP  = std::pow(pl, -p[GAMMA_P]);
  const double gf   = (u.gcsf > 0.0) ? p[GCSF_FOLD] : 1.0;
  const double tf   = (u.tpo  > 0.0) ? p[TPO_FOLD]  : 1.0;

  dy[PN]  = p[KTR_N] * (F * fbN * (1.0 - eac) * gf - y[PN]) - p[KLDP] * c * y[PN];
  dy[TN1] = p[KTR_N] * (y[PN]  - y[TN1]);
  dy[TN2] = p[KTR_N] * (y[TN1] - y[TN2]);
  dy[TN3] = p[KTR_N] * (y[TN2] - y[TN3]);
  dy[NN]  = p[KCIRC_N] * (y[TN3] - y[NN]);

  dy[PP]  = p[KTR_P] * (F * fbP * (1.0 - p[E_PLT] * eac) * tf - y[PP]) - p[KLDP] * c * y[PP];
  dy[TP1] = p[KTR_P] * (y[PP]  - y[TP1]);
  dy[TP2] = p[KTR_P] * (y[TP1] - y[TP2]);
  dy[TP3] = p[KTR_P] * (y[TP2] - y[TP3]);
  dy[PL]  = p[KCIRC_P] * (y[TP3] - y[PL]);
}

// ---- Dormand-Prince RK45 --------------------------------------------------------
struct Ctl { double rtol, atol, hmax; long max_steps; long steps; bool failed; };

void integrate_to(double& t, double t_end, std::vector<double>& y, const double* p,
                  const Inputs& u, double f0n, double& h, Ctl& ctl) {
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

    rhs(t, y.data(), k1, p, u, f0n);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * a21 * k1[i];
    rhs(t + c2 * hh, yt, k2, p, u, f0n);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a31 * k1[i] + a32 * k2[i]);
    rhs(t + c3 * hh, yt, k3, p, u, f0n);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
    rhs(t + c4 * hh, yt, k4, p, u, f0n);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
    rhs(t + c5 * hh, yt, k5, p, u, f0n);
    for (int i = 0; i < n; ++i) yt[i] = y[i] + hh * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
    rhs(t + hh, yt, k6, p, u, f0n);
    for (int i = 0; i < n; ++i) yn[i] = y[i] + hh * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
    rhs(t + hh, yn, k7, p, u, f0n);

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

void initial_state(std::vector<double>& y, const double* p) {
  std::fill(y.begin(), y.end(), 1.0);
  y[XW] = 1.0 - p[F0];
  y[XM] = p[F0];
  y[CLD] = 0.0;
  y[IS] = 0.0;
  y[IM1] = 0.0;
  y[IM] = 0.0;
}

} // namespace mm

// [[Rcpp::export]]
CharacterVector marrow_param_names() {
  return CharacterVector::create(
    "rho", "g", "f0",
    "a", "kel_ld", "klds", "kldp",
    "r_ld", "r_inf", "delta",
    "q", "emax", "ic50", "hill",
    "k_out", "k_m", "t_pk", "h_shape",
    "ktr_n", "kcirc_n", "gamma_n",
    "ktr_p", "kcirc_p", "gamma_p", "e_plt",
    "gcsf_fold", "tpo_fold",
    "toci_red", "ster_red", "ster_boost");
}

// [[Rcpp::export]]
CharacterVector marrow_state_names() {
  return CharacterVector::create("xW", "xM", "cLD", "Is", "Im1", "Im",
                                 "PN", "TN1", "TN2", "TN3", "N",
                                 "PP", "TP1", "TP2", "TP3", "PL");
}

// Simulate one patient.
// seg   : matrix with columns start, end, ld, gcsf, toci, ster, tpo (rows ordered in time)
// obs_t : ascending observation times
// value : matrix with columns t, anc_rel, plt_rel, is, im, f_clone, x_tot, f_rel, cld
//         (all NaN if the integrator failed)
// [[Rcpp::export]]
NumericMatrix marrow_sim(NumericVector par, NumericMatrix seg, NumericVector obs_t,
                         double rtol = 1e-6, double atol = 1e-9, double hmax = 0.5,
                         int max_steps = 200000) {
  if (par.size() != mm::N_PAR) stop("parameter vector has wrong length");
  if (seg.ncol() != 7) stop("segment matrix must have 7 columns");
  const double* p = par.begin();
  const int nobs = obs_t.size();
  NumericMatrix out(nobs, 9);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "is", "im",
                                          "f_clone", "x_tot", "f_rel", "cld");

  std::vector<double> y(mm::N_STATE);
  mm::initial_state(y, p);
  const double f0n = mm::clone_f0_norm(p);
  mm::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};

  auto record = [&](int k, double t) {
    const double xw = std::max(y[mm::XW], 0.0), xm = std::max(y[mm::XM], 0.0);
    out(k, 0) = t;
    out(k, 1) = y[mm::NN];
    out(k, 2) = y[mm::PL];
    out(k, 3) = y[mm::IS];
    out(k, 4) = y[mm::IM];
    out(k, 5) = (xw + xm > 0) ? xm / (xw + xm) : NA_REAL;
    out(k, 6) = xw + xm;
    out(k, 7) = (xw + p[mm::DELTA] * xm) / f0n;
    out(k, 8) = y[mm::CLD];
  };

  double t = seg(0, 0);
  double h = 0.05;
  int k = 0;
  while (k < nobs && obs_t[k] <= t + 1e-12) { record(k, obs_t[k]); ++k; }

  // stop once every observation time is recorded (later integration cannot
  // change the output)
  for (int s = 0; s < seg.nrow() && !ctl.failed && k < nobs; ++s) {
    mm::Inputs u{seg(s, 2), seg(s, 3), seg(s, 4), seg(s, 5), seg(s, 6)};
    const double t_end = seg(s, 1);
    while (k < nobs && obs_t[k] <= t_end + 1e-12 && !ctl.failed) {
      mm::integrate_to(t, obs_t[k], y, p, u, f0n, h, ctl);
      if (!ctl.failed) { record(k, obs_t[k]); ++k; }
    }
    if (!ctl.failed) mm::integrate_to(t, t_end, y, p, u, f0n, h, ctl);
  }

  if (ctl.failed || k < nobs) {
    std::fill(out.begin(), out.end(), NA_REAL);
  }
  return out;
}

// Right-hand side exposed for validation against deSolve.
// u : ld, gcsf, toci, ster, tpo
// [[Rcpp::export]]
NumericVector marrow_rhs(double t, NumericVector y, NumericVector par, NumericVector u) {
  if (par.size() != mm::N_PAR) stop("parameter vector has wrong length");
  if (y.size() != mm::N_STATE) stop("state vector has wrong length");
  mm::Inputs in{u[0], u[1], u[2], u[3], u[4]};
  NumericVector dy(mm::N_STATE);
  mm::rhs(t, y.begin(), dy.begin(), par.begin(), in, mm::clone_f0_norm(par.begin()));
  return dy;
}

// Initial state exposed for validation against deSolve.
// [[Rcpp::export]]
NumericVector marrow_init(NumericVector par) {
  std::vector<double> y(mm::N_STATE);
  mm::initial_state(y, par.begin());
  return NumericVector(y.begin(), y.end());
}
