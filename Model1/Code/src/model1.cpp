// =============================================================================
// model1.cpp
// Purpose : Model1 marrow model (Aim4 v1.3 with a Friberg-like haematopoietic
//           layer): stem cells with a pre-existing clone, four MPP components
//           (one per blood lineage), Friberg maturation chains, CAR T - target
//           module and G-CSF input; plus the original Friberg (2002)
//           myelosuppression model as reference model MF. One adaptive
//           Dormand-Prince RK45 integrator, integrated piecewise between input
//           breakpoints.
// Inputs  : parameter vector (layout = m1_param_names() / mf_param_names()),
//           a segment matrix of piecewise-constant inputs, the CAR T-phase
//           driver D(t) on a uniform grid, and observation times (days).
// Outputs : matrix of model quantities at the observation times.
//
// Model1 (every compartment normalised to its pre-treatment value = 1)
//   Exposure  C          dC/dt = dose(t) - k_e*C        (bridging + flu-cy)
//   Marrow    I1, I      dI1/dt = k_I*(D - I1), dI/dt = k_I*(I1 - I)
//                        E(I) = Emax*I^h/(I^h + IC50^h)
//   HSC       S_W, S_X   share one niche, S = S_W + S_X; self-renewal rate r_S
//     dS_j/dt = r_S*S_j*(1-S) - [ra_j*th_aS*C + rb_j*(th_bS*I + th_zS*Z)*(1 + eta*(1-S)+)]*S_j
//     ra_W = rb_W = 1, ra_X = r_a, rb_X = r_b;  S_W(0) = 1 - f0, S_X(0) = f0
//     F = (S_W + delta*S_X)/(1 - f0 + delta*f0)                (stem cell output)
//   MPP components P_l, l in {N, Mo, Pl, L} (one per blood lineage)
//     dP_l/dt = (k_P - sr_P)*F + (sr_P*G_l - k_P)*P_l - th_aP*C*P_l
//     G_l = min(gmax, B_l^-gamma) * (1 - E(I)) * (1 - e_Z*Z/(Z + Z50)),  l in {N, Mo, Pl}
//     G_L = (1 - E(I)) * (1 - e_Z*Z/(Z + Z50))   (no count feedback: lymphodepletion keeps
//     B_L low for weeks and a feedback on P_L made ALC overshoot to ~1.3x baseline by day 30;
//     K1 2026-11-04)
//     each myeloid component senses only the blood count of its own lineage (Friberg);
//     with sr_P = k_P and F = 1, P_l is the proliferating compartment of
//     Friberg et al. (k_prol = k_P)
//   Maturation chain l in {N, Mo, Pl}
//     dT1/dt = ktr_l*[P_l*u_l(t) - T1],  dT2/dt = ktr_l*(T1 - T2)
//     dT3/dt = ktr_l*(T2 - rel_l*T3),    dB_l/dt = c_l*(rel_l*T3 - B_l)
//     u_N = 1 + (u_G - 1)*GC (G-CSF), u_Pl = TPO-RA fold while on (simulation
//     only), u_Mo = 1; rel_N = 1 + u_R*GC, rel_Mo = rel_Pl = 1
//   Lymphocytes (no maturation chain; long-lived blood pool)
//     dB_L/dt = c_L*(P_L - B_L) - th_aL*C*B_L
//   No compartment of the myeloid lineages, the stem cells, the CAR T - target
//   module or the inputs contains P_L or B_L.
//   Segment column s_add: S_W += s_add at the segment start (CD34+ boost).
//
// CAR T - target interaction (CT: CAR T cells in blood, 10^9/L, absolute;
// AG: CD19 antigen load, relative; occ = AG/(AG + K_A))
//   dCT/dt = (rho*occ - k_T)*CT,   dAG/dt = -kappa*CT*occ,   act = CT*occ
//   dZ1/dt = k_Z*(act - Z1),       dZ = k_Z*(Z1 - Z)   (two-stage lag)
//   infusion: segment column car: CT += car at the segment start (day 0);
//   AG(0) = a0 (per patient), CT = Z1 = Z = 0 before the infusion
// G-CSF effect GC: dGC/dt = k_G*(g(t) - GC), g = segment column gcsf (0/1)
// ALC observed = b_ALC * B_L + CT (assembled in R)
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

namespace m1 {

// ---- parameter layouts (must match m1_param_names / mf_param_names) ------------
enum Par {
  R_S, F0, TH_AS, TH_AP, TH_BS, RHO, A0,                // individual (A0 from LDH)
  R_A, R_B, DELTA, ETA_SYN, TH_ZS,                      // clone traits, a x b synergy, Z route to HSC
  K_P, SR_P, GMAX, GAMMA,                               // MPP components: exit, self-renewal, feedback
  KTR_N, KTR_MO, KTR_PL,                                // maturation chains (= 3 / MTT)
  C_N, C_MO, C_PL, C_L,                                 // blood turnover
  K_I, EMAX, IC50, HILL,                                // marrow inflammation
  TH_AL,                                                // lymphodepletion of blood lymphocytes
  K_E, GCSF_FOLD, TPO_FOLD, K_G, GCSF_REL,              // exposure decay, supportive care
  K_T, K_AG, KAPPA, K_Z, Z50, E_Z,                      // CAR T - target interaction
  N_PAR
};
enum State { SW, SX, PN, PMO, PPL, PL, CX, I1, I2,
             TN1, TN2, TN3, BN,
             TM1, TM2, TM3, BM,
             TP1, TP2, TP3, BP,
             BL,
             CT, AG, Z1, Z2, GC, N_STATE };

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

// ---- Model1 ---------------------------------------------------------------------------
struct Parent {
  const double* p;
  const Driver* drv;
  Inputs u;
  int n_state() const { return N_STATE; }

  // rel: release factor of the last maturation compartment into blood (G-CSF; 1 otherwise)
  static void chain(const double* y, double* dy, int t1, double k, double c, double inflow, double rel = 1.0) {
    dy[t1]     = k * (inflow - y[t1]);
    dy[t1 + 1] = k * (y[t1] - y[t1 + 1]);
    dy[t1 + 2] = k * (y[t1 + 1] - rel * y[t1 + 2]);
    dy[t1 + 3] = c * (rel * y[t1 + 2] - y[t1 + 3]);
  }

  // events at a segment start: CD34+ boost into the normal stem cells, CAR T infusion
  static void event(std::vector<double>& y, double s_add, double car) {
    y[SW] += s_add;
    y[CT] += car;
  }

  void rhs(double t, const double* y, double* dy) const {
    const double sw = std::max(y[SW], 0.0), sx = std::max(y[SX], 0.0), s = sw + sx;
    const double cx = std::max(y[CX], 0.0), i1 = std::max(y[I1], 0.0), im = std::max(y[I2], 0.0);
    const double d = std::max(drv->at(t), 0.0);
    const double bn = std::max(y[BN], 1e-6), bm = std::max(y[BM], 1e-6), bp = std::max(y[BP], 1e-6);
    const double bl = std::max(y[BL], 1e-6);
    const double ct = std::max(y[CT], 0.0), ag = std::max(y[AG], 0.0);
    const double zz = std::max(y[Z2], 0.0), gc = std::max(y[GC], 0.0);

    // CAR T - target interaction: antigen-driven expansion, killing, lagged signal Z
    const double occ = ag / (ag + p[K_AG]);
    const double act = ct * occ;
    dy[CT] = (p[RHO] * occ - p[K_T]) * ct;
    dy[AG] = -p[KAPPA] * act;
    dy[Z1] = p[K_Z] * (act - y[Z1]);
    dy[Z2] = p[K_Z] * (y[Z1] - y[Z2]);
    dy[GC] = p[K_G] * (u.gcsf - y[GC]);

    // HSC: shared-niche regrowth, timed injuries; CAR T phase (b): ferritin-driven
    // inflammation I and the CAR T - target signal Z
    const double regrow = p[R_S] * (1.0 - s);
    const double kill_a = p[TH_AS] * cx;
    const double kill_b = (p[TH_BS] * im + p[TH_ZS] * zz) * (1.0 + p[ETA_SYN] * std::max(1.0 - s, 0.0));
    dy[SW] = regrow * sw - (kill_a + kill_b) * sw;
    dy[SX] = regrow * sx - (p[R_A] * kill_a + p[R_B] * kill_b) * sx;
    const double F = (sw + p[DELTA] * sx) / (1.0 - p[F0] + p[DELTA] * p[F0]);

    // MPP components: self-renewal feedback from the own lineage's blood count, suppressed
    // by marrow inflammation and by the CAR T - target signal; chemotherapy kills
    const double ih = std::pow(im, p[HILL]);
    const double e = p[EMAX] * ih / (ih + std::pow(p[IC50], p[HILL]));
    const double sup = (1.0 - e) * (1.0 - p[E_Z] * zz / (zz + p[Z50]));
    const double kin = p[K_P] - p[SR_P], kp = p[TH_AP] * cx;
    // fb = 0: no count feedback (lymphoid component, see the header)
    auto mpp = [&](int x, double b, double fb) {
      const double g = (fb > 0.0 ? std::min(p[GMAX], std::pow(b, -p[GAMMA])) : 1.0) * sup;
      const double px = std::max(y[x], 0.0);
      dy[x] = kin * F + (p[SR_P] * g - p[K_P]) * px - kp * px;
    };
    mpp(PN, bn, 1.0);
    mpp(PMO, bm, 1.0);
    mpp(PPL, bp, 1.0);
    mpp(PL, bl, 0.0);

    // lymphocytes: long-lived blood pool, depleted directly by lymphodepletion
    const double kl = p[TH_AL] * cx;
    dy[BL] = p[C_L] * (std::max(y[PL], 0.0) - y[BL]) - kl * y[BL];

    // exposure and the two-stage marrow lag of the CAR T-phase driver
    dy[CX] = u.dose - p[K_E] * cx;
    dy[I1] = p[K_I] * (d - i1);
    dy[I2] = p[K_I] * (i1 - im);

    // maturation chains (pure delays); supportive care at entry; G-CSF also accelerates
    // the release of mature neutrophils into blood
    const double gf = 1.0 + (p[GCSF_FOLD] - 1.0) * gc;
    const double rel_n = 1.0 + p[GCSF_REL] * gc;
    const double tf = (u.tpo > 0.0) ? p[TPO_FOLD] : 1.0;
    chain(y, dy, TN1, p[KTR_N],  p[C_N],  std::max(y[PN], 0.0) * gf, rel_n);
    chain(y, dy, TM1, p[KTR_MO], p[C_MO], std::max(y[PMO], 0.0));
    chain(y, dy, TP1, p[KTR_PL], p[C_PL], std::max(y[PPL], 0.0) * tf);
  }

  // Pre-treatment steady state
  void init(std::vector<double>& y) const {
    std::fill(y.begin(), y.end(), 1.0);
    y[SW] = 1.0 - p[F0];
    y[SX] = p[F0];
    y[CX] = 0.0;
    y[I1] = 0.0;
    y[I2] = 0.0;
    y[CT] = 0.0;
    y[AG] = p[A0];
    y[Z1] = 0.0;
    y[Z2] = 0.0;
    y[GC] = 0.0;
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

  // no stem cell layer and no CAR T compartment in MF
  static void event(std::vector<double>&, double, double) {}
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
    if (seg(s, 5) != 0.0 || seg(s, 6) != 0.0) sys.event(y, seg(s, 5), seg(s, 6));  // CD34+ boost, CAR T infusion
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
  if (seg.ncol() != 7) stop("segment matrix must have 7 columns: start, end, dose, gcsf, tpo, s_add, car");
  if (seg.nrow() < 1) stop("segment matrix is empty");
  if (drv.size() < 2) stop("driver grid needs at least two points");
}

} // namespace m1

// [[Rcpp::export]]
CharacterVector m1_param_names() {
  return CharacterVector::create(
    "r_s", "f0", "th_as", "th_ap", "th_bs", "rho", "a0",
    "r_a", "r_b", "delta", "eta_syn", "th_zs",
    "k_p", "sr_p", "gmax", "gamma",
    "ktr_n", "ktr_mo", "ktr_pl",
    "c_n", "c_mo", "c_pl", "c_l",
    "k_i", "emax", "ic50", "hill",
    "th_al",
    "k_e", "gcsf_fold", "tpo_fold", "k_g", "gcsf_rel",
    "k_t", "k_ag", "kappa", "k_z", "z50", "e_z");
}

// [[Rcpp::export]]
CharacterVector m1_state_names() {
  return CharacterVector::create("SW", "SX", "PN", "PMO", "PPL", "PL", "C", "I1", "I",
                                 "TN1", "TN2", "TN3", "BN",
                                 "TM1", "TM2", "TM3", "BM",
                                 "TP1", "TP2", "TP3", "BP",
                                 "BL",
                                 "CT", "AG", "Z1", "Z", "GC");
}

// [[Rcpp::export]]
CharacterVector mf_param_names() {
  return CharacterVector::create("slope", "mtt_scale", "ktr_n", "ktr_pl", "ktr_mo",
                                 "gam_n", "gam_pl", "gam_mo", "slope_pl", "slope_mo", "k_e");
}

// Simulate one patient with Model1.
// seg   : matrix with columns start, end, dose, gcsf, tpo, s_add, car (rows ordered in time)
// drv   : D(t) on the grid drv_t0 + (0:(n-1)) * drv_dt
// obs_t : ascending observation times (>= first segment start)
// value : columns t, anc_rel, plt_rel, amc_rel, S, f_clone, M (mean of the three myeloid
//         MPP components), C, D, I, alc_rel, PN, PMO, PPL, PL (MPP components), cart (CAR T
//         cells in blood, 10^9/L), ag (antigen load), z (CAR T - target signal Z), gcsf
//         (G-CSF effect) (all NA if the integrator failed)
// [[Rcpp::export]]
NumericMatrix m1_sim(NumericVector par, NumericMatrix seg, NumericVector drv, double drv_t0,
                     double drv_dt, NumericVector obs_t, double rtol = 1e-6, double atol = 1e-9,
                     double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != m1::N_PAR) stop("parameter vector has wrong length");
  m1::check_inputs(seg, drv);
  m1::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  m1::Parent sys{par.begin(), &dr, m1::Inputs{0, 0, 0}};
  m1::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  NumericMatrix out(obs_t.size(), 19);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                          "M", "C", "D", "I", "alc_rel", "PN", "PMO", "PPL",
                                          "PL", "cart", "ag", "z", "gcsf");
  auto record = [&](int k, double t, const std::vector<double>& y) {
    const double s = y[m1::SW] + y[m1::SX];
    out(k, 0) = t;
    out(k, 1) = y[m1::BN];
    out(k, 2) = y[m1::BP];
    out(k, 3) = y[m1::BM];
    out(k, 4) = s;
    out(k, 5) = (s > 1e-12) ? y[m1::SX] / s : 0.0;
    out(k, 6) = (y[m1::PN] + y[m1::PMO] + y[m1::PPL]) / 3.0;
    out(k, 7) = y[m1::CX];
    out(k, 8) = dr.at(t);
    out(k, 9) = y[m1::I2];
    out(k, 10) = y[m1::BL];
    out(k, 11) = y[m1::PN];
    out(k, 12) = y[m1::PMO];
    out(k, 13) = y[m1::PPL];
    out(k, 14) = y[m1::PL];
    out(k, 15) = y[m1::CT];
    out(k, 16) = y[m1::AG];
    out(k, 17) = y[m1::Z2];
    out(k, 18) = y[m1::GC];
  };
  if (!m1::run_segments(sys, seg, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Simulate one patient with the Friberg reference model (same first ten output
// columns; stem-layer columns are NA).
// [[Rcpp::export]]
NumericMatrix mf_sim(NumericVector par, NumericMatrix seg, NumericVector obs_t, double rtol = 1e-6,
                     double atol = 1e-9, double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != m1::MF_N_PAR) stop("MF parameter vector has wrong length");
  NumericVector dummy = NumericVector::create(0.0, 0.0);
  m1::check_inputs(seg, dummy);
  m1::Driver dr{0.0, 1.0, dummy.begin(), 2};
  m1::Friberg sys{par.begin(), &dr, m1::Inputs{0, 0, 0}};
  m1::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  NumericMatrix out(obs_t.size(), 10);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                          "M", "C", "D", "I");
  auto record = [&](int k, double t, const std::vector<double>& y) {
    out(k, 0) = t;
    out(k, 1) = y[m1::F_BN];
    out(k, 2) = y[m1::F_BP];
    out(k, 3) = y[m1::F_BM];
    out(k, 4) = NA_REAL;
    out(k, 5) = NA_REAL;
    out(k, 6) = y[m1::F_PN];
    out(k, 7) = y[m1::F_CX];
    out(k, 8) = 0.0;
    out(k, 9) = 0.0;
  };
  // the CD34 and CAR T columns have no meaning in MF (its event() ignores them)
  if (!m1::run_segments(sys, seg, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Right-hand side of Model1 (27 states), exposed for validation against deSolve and
// for the steady-state Jacobian.
// [[Rcpp::export]]
NumericVector m1_rhs(double t, NumericVector y, NumericVector par, double dose, double gcsf, double tpo,
                     NumericVector drv, double drv_t0, double drv_dt) {
  if (par.size() != m1::N_PAR) stop("parameter vector has wrong length");
  if (y.size() != m1::N_STATE) stop("state vector has wrong length");
  m1::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  m1::Parent sys{par.begin(), &dr, m1::Inputs{dose, gcsf, tpo}};
  NumericVector dy(m1::N_STATE);
  sys.rhs(t, y.begin(), dy.begin());
  return dy;
}
