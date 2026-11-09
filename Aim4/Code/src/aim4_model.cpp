// =============================================================================
// aim4_model.cpp
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
// Parent model v1.2 (lineage tree of Aim4/image (2).png; every compartment
// normalised to its pre-treatment value = 1)
//   Exposure  C          dC/dt = dose(t) - k_e*C        (bridging + flu-cy)
//   Marrow    I1, I      dI1/dt = k_I*(D - I1), dI/dt = k_I*(I1 - I)
//                        E(I) = Emax*I^h/(I^h + IC50^h)
//   HSC       S_W, S_X   share one niche, S = S_W + S_X; self-renewal rate r_S
//     dS_j/dt = r_S*S_j*(1-S)*(1 + alpha*(1-MPP)+)
//               - [ra_j*th_aS*C + rb_j*th_bS*I*(1 + eta*(1-S)+)]*S_j
//     ra_W = rb_W = 1, ra_X = r_a, rb_X = r_b;  S_W(0) = 1 - f0, S_X(0) = f0
//     F = (S_W + delta*S_X)/(1 - f0 + delta*f0)                (stem cell output)
//   Self-renewing progenitor X in {MPP, CMP, GMP, MEP, CLP}, upstream U_X:
//     dX/dt = (k_X - sr_X)*U_X + (sr_X*G_X - k_X)*X - th_aP*C*X
//     U: MPP <- F, CMP <- MPP, CLP <- MPP, GMP <- CMP, MEP <- CMP
//     k_X exit rate, sr_X self-renewal rate (sr_X < k_X)
//     G_X = min(gmax, s_X^-gamma_M) (= 1 when switched off), s_X = geometric
//     mean of the blood counts of X's own progeny: MPP, CMP: ANC, AMC, PLT;
//     GMP: ANC, AMC; MEP: PLT; CLP, lymphocytes: ALC; for MPP, CMP, GMP and
//     MEP G_X is multiplied by (1 - e_sr*E) (0 = inflammation does not act on
//     self-renewal)
//   Maturation chain l in {N, Mo, Pl}, parent P_N = P_Mo = GMP, P_Pl = MEP
//     dT1/dt = ktr_l*[P_l*B_l^-gam_l*(1 - e_l*E)*u_l(t) - T1],  e_N = 1
//     dTj/dt = ktr_l*(T_{j-1} - Tj), j = 2, 3;   dB_l/dt = c_l*(T3 - B_l)
//     u_N = 1 + (u_G - 1)*GC (G-CSF effect, v1.3 below), u_Pl = TPO-RA fold while on (simulation only)
//   Lymphocytes (self-renewing pool L, blood B_L)
//     dL/dt   = (k_L - sr_L)*CLP*(1 - e_L*E) + (sr_L*G_L - k_L)*L - th_aL*C*L
//     dB_L/dt = c_L*(L - B_L) - th_aL*C*B_L
//   The lymphoid branch feeds nothing upstream: ANC, PLT and AMC do not depend
//   on any lymphoid parameter.
//   Segment column s_add: S_W += s_add at the segment start (CD34+ boost).
//
// v1.3 additions
//   CAR T - target interaction (CT: CAR T cells in blood, 10^9/L, absolute;
//   AG: CD19 antigen load, relative; occ = AG/(AG + K_A))
//     dCT/dt = (rho*occ - k_T)*CT,   dAG/dt = -kappa*CT*occ,   act = CT*occ
//     dZ1/dt = k_Z*(act - Z1),       dZ = k_Z*(Z1 - Z)   (two-stage lag)
//     infusion: segment column car: CT += car at the segment start (day 0);
//     AG(0) = a0 (per patient), CT = Z1 = Z = 0 before the infusion
//   Effects: HSC kill_b = (th_bS*I + th_zS*Z)*(1 + eta*(1-S)+) (clone: r_b*kill_b);
//     MPP self-renewal G_MPP *= (1 - e_Z*Z/(Z + Z50))
//   G-CSF effect GC: dGC/dt = k_G*(g(t) - GC), g = segment column gcsf (0/1)
//     neutrophil chain entry x (1 + (u_G - 1)*GC)          (u_G = gcsf_fold)
//     release T3 -> blood at rate x rel, rel = 1 + u_R*GC  (u_R = gcsf_rel):
//     dT3/dt = ktr*(T2 - rel*T3), dB_N/dt = c_N*(rel*T3 - B_N)
//   ALC observed = b_ALC * B_L + CT (assembled in R)
//
// Age and generation bookkeeping (track = true)
//   For each compartment downstream of HSC, first moments of the time since
//   entry into the MPP compartment (age mass) and of the number of divisions
//   since that entry (generation mass). Every transfer between compartments
//   is one division; a self-renewal division adds one generation to both
//   daughters; chain l adds d_l divisions over T1->T2->T3 (endomitosis of
//   megakaryocytes is not a division: d_Pl = 0). Means = mass / compartment.
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

namespace a4 {

// ---- parameter layouts (must match a4_param_names / mf_param_names) ------------
enum Par {
  R_S, F0, TH_AS, TH_AP, TH_BS,                         // individual
  ALPHA,                                                // HSC activation by MPP depletion
  K_MPP, SR_MPP, K_CMP, SR_CMP, K_GMP, SR_GMP,          // progenitor exit / self-renewal rates
  K_MEP, SR_MEP, K_CLP, SR_CLP, K_LY, SR_LY,
  GMAX, GAMMA_M, FB_MPP, FB_CMP, FB_GMP, FB_MEP, FB_CLP, FB_LY,  // self-renewal feedback
  R_A, R_B, DELTA,                                      // clone traits
  KTR_N, KTR_MO, KTR_PL,                                // maturation chains (= 3 / MTT)
  C_N, C_MO, C_PL, C_L,                                 // blood turnover
  GAM_N, GAM_MO, GAM_PL,                                // production feedback at chain entry
  K_I, EMAX, IC50, HILL, E_MO, E_PL, E_L, ETA_SYN,      // marrow inflammation, a x b synergy
  E_SR,                                                 // inflammation on MPP, CMP, GMP, MEP self-renewal
  TH_AL,                                                // lymphodepletion of lymphocytes
  K_E, GCSF_FOLD, TPO_FOLD,                             // exposure decay, supportive care
  D_N, D_MO, D_PL,                                      // divisions inside the chains (bookkeeping)
  RHO, A0,                                              // v1.3: CAR T expansion (individual), antigen load
  K_T, K_AG, KAPPA, K_Z, Z50, TH_ZS, E_Z,               // v1.3: CAR T - target interaction
  K_G, GCSF_REL,                                        // v1.3: G-CSF effect
  N_PAR
};
enum State { SW, SX, MPP, CMP, GMP, MEP, CLP, LY, CX, I1, I2,
             TN1, TN2, TN3, BN,
             TM1, TM2, TM3, BM,
             TP1, TP2, TP3, BP,
             BL,
             CT, AG, Z1, Z2, GC, N_STATE };
// compartments with age / generation moments (track = true)
enum Trk { Q_MPP, Q_CMP, Q_GMP, Q_MEP, Q_CLP, Q_LY,
           Q_TN1, Q_TN2, Q_TN3, Q_BN, Q_TM1, Q_TM2, Q_TM3, Q_BM,
           Q_TP1, Q_TP2, Q_TP3, Q_BP, Q_BL, N_TRK };

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
  bool track;
  int n_state() const { return N_STATE + (track ? 2 * N_TRK : 0); }

  // self-renewal feedback from the blood counts of the compartment's progeny
  double feedback(double signal, double on) const {
    if (on <= 0.0) return 1.0;
    return std::min(p[GMAX], std::pow(std::max(signal, 1e-6), -p[GAMMA_M]));
  }

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
    const double mpp = std::max(y[MPP], 0.0), cmp = std::max(y[CMP], 0.0), gmp = std::max(y[GMP], 0.0);
    const double mep = std::max(y[MEP], 0.0), clp = std::max(y[CLP], 0.0), ly = std::max(y[LY], 0.0);
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

    // HSC: shared-niche regrowth (faster when the MPP is depleted), timed injuries;
    // CAR T phase (b): ferritin-driven inflammation I and the CAR T - target signal Z
    const double regrow = p[R_S] * (1.0 - s) * (1.0 + p[ALPHA] * std::max(1.0 - mpp, 0.0));
    const double kill_a = p[TH_AS] * cx;
    const double kill_b = (p[TH_BS] * im + p[TH_ZS] * zz) * (1.0 + p[ETA_SYN] * std::max(1.0 - s, 0.0));
    dy[SW] = regrow * sw - (kill_a + kill_b) * sw;
    dy[SX] = regrow * sx - (p[R_A] * kill_a + p[R_B] * kill_b) * sx;
    const double F = (sw + p[DELTA] * sx) / (1.0 - p[F0] + p[DELTA] * p[F0]);

    // marrow inflammation: reversible suppression of production at each lineage entry
    // and (strength e_sr, 0 = off) of MPP and myeloid progenitor self-renewal
    const double ih = std::pow(im, p[HILL]);
    const double e = p[EMAX] * ih / (ih + std::pow(p[IC50], p[HILL]));

    // self-renewing progenitors; each senses the blood counts of its own progeny
    const double my = std::cbrt(bn * bm * bp);
    const double sr_inf = 1.0 - p[E_SR] * e;
    const double sr_z = 1.0 - p[E_Z] * zz / (zz + p[Z50]);   // CAR T - target signal on MPP self-renewal
    const double g_mpp = feedback(my, p[FB_MPP]) * sr_inf * sr_z, g_cmp = feedback(my, p[FB_CMP]) * sr_inf;
    const double g_gmp = feedback(std::sqrt(bn * bm), p[FB_GMP]) * sr_inf, g_mep = feedback(bp, p[FB_MEP]) * sr_inf;
    const double g_clp = feedback(bl, p[FB_CLP]), g_ly = feedback(bl, p[FB_LY]);
    const double kp = p[TH_AP] * cx;
    dy[MPP] = (p[K_MPP] - p[SR_MPP]) * F   + (p[SR_MPP] * g_mpp - p[K_MPP]) * mpp - kp * mpp;
    dy[CMP] = (p[K_CMP] - p[SR_CMP]) * mpp + (p[SR_CMP] * g_cmp - p[K_CMP]) * cmp - kp * cmp;
    dy[GMP] = (p[K_GMP] - p[SR_GMP]) * cmp + (p[SR_GMP] * g_gmp - p[K_GMP]) * gmp - kp * gmp;
    dy[MEP] = (p[K_MEP] - p[SR_MEP]) * cmp + (p[SR_MEP] * g_mep - p[K_MEP]) * mep - kp * mep;
    dy[CLP] = (p[K_CLP] - p[SR_CLP]) * mpp + (p[SR_CLP] * g_clp - p[K_CLP]) * clp - kp * clp;

    // lymphocytes: self-renewing pool and blood, depleted directly by lymphodepletion
    const double kl = p[TH_AL] * cx;
    const double f_ly = 1.0 - p[E_L] * e;
    dy[LY] = (p[K_LY] - p[SR_LY]) * clp * f_ly + (p[SR_LY] * g_ly - p[K_LY]) * ly - kl * ly;
    dy[BL] = p[C_L] * (ly - y[BL]) - kl * y[BL];

    // exposure and the two-stage marrow lag of the CAR T-phase driver
    dy[CX] = u.dose - p[K_E] * cx;
    dy[I1] = p[K_I] * (d - i1);
    dy[I2] = p[K_I] * (i1 - im);

    // maturation chains: production feedback, inflammation, supportive care at entry;
    // G-CSF also accelerates the release of mature neutrophils into blood
    const double gf = 1.0 + (p[GCSF_FOLD] - 1.0) * gc;
    const double rel_n = 1.0 + p[GCSF_REL] * gc;
    const double tf = (u.tpo > 0.0) ? p[TPO_FOLD] : 1.0;
    const double f_n  = std::pow(bn, -p[GAM_N])  * (1.0 - e) * gf;
    const double f_mo = std::pow(bm, -p[GAM_MO]) * (1.0 - p[E_MO] * e);
    const double f_pl = std::pow(bp, -p[GAM_PL]) * (1.0 - p[E_PL] * e) * tf;
    chain(y, dy, TN1, p[KTR_N],  p[C_N],  gmp * f_n, rel_n);
    chain(y, dy, TM1, p[KTR_MO], p[C_MO], gmp * f_mo);
    chain(y, dy, TP1, p[KTR_PL], p[C_PL], mep * f_pl);

    if (!track) return;

    // ---- age (time since MPP entry) and generation moments ----
    const double* A = y + N_STATE;
    const double* G = y + N_STATE + N_TRK;
    double* dA = dy + N_STATE;
    double* dG = dy + N_STATE + N_TRK;
    // self-renewing compartment q (state x): inflow coefficient `in` times the upstream
    // masses (aU, gU) and upstream size U; one division on entry (inc = 1, 0 for MPP)
    auto renewing = [&](int q, double x, double k, double sr, double gx, double kill,
                        double in, double aU, double gU, double U, double inc) {
      const double r = sr * gx;
      dA[q] = x + in * aU + (r - k - kill) * A[q];
      dG[q] = in * (gU + inc * U) + r * (G[q] + 2.0 * x) - (k + kill) * G[q];
    };
    renewing(Q_MPP, mpp, p[K_MPP], p[SR_MPP], g_mpp, kp, p[K_MPP] - p[SR_MPP], 0.0, 0.0, F, 0.0);
    renewing(Q_CMP, cmp, p[K_CMP], p[SR_CMP], g_cmp, kp, p[K_CMP] - p[SR_CMP], A[Q_MPP], G[Q_MPP], mpp, 1.0);
    renewing(Q_GMP, gmp, p[K_GMP], p[SR_GMP], g_gmp, kp, p[K_GMP] - p[SR_GMP], A[Q_CMP], G[Q_CMP], cmp, 1.0);
    renewing(Q_MEP, mep, p[K_MEP], p[SR_MEP], g_mep, kp, p[K_MEP] - p[SR_MEP], A[Q_CMP], G[Q_CMP], cmp, 1.0);
    renewing(Q_CLP, clp, p[K_CLP], p[SR_CLP], g_clp, kp, p[K_CLP] - p[SR_CLP], A[Q_MPP], G[Q_MPP], mpp, 1.0);
    renewing(Q_LY,  ly,  p[K_LY],  p[SR_LY],  g_ly,  kl, (p[K_LY] - p[SR_LY]) * f_ly, A[Q_CLP], G[Q_CLP], clp, 1.0);
    dA[Q_BL] = y[BL] + p[C_L] * A[Q_LY] - (p[C_L] + kl) * A[Q_BL];
    dG[Q_BL] = p[C_L] * G[Q_LY] - (p[C_L] + kl) * G[Q_BL];
    // maturation chain: entry division (+1), d_l/2 divisions on each of T1->T2, T2->T3;
    // rel scales the outflow of T3 into blood (G-CSF release, neutrophils only)
    auto mchain = [&](int q1, int t1, double k, double c, double f, int qp, double P, double dl, double rel) {
      dA[q1]     = y[t1] + k * f * A[qp] - k * A[q1];
      dG[q1]     = k * f * (G[qp] + P) - k * G[q1];
      dA[q1 + 1] = y[t1 + 1] + k * A[q1] - k * A[q1 + 1];
      dG[q1 + 1] = k * (G[q1] + 0.5 * dl * y[t1]) - k * G[q1 + 1];
      dA[q1 + 2] = y[t1 + 2] + k * A[q1 + 1] - k * rel * A[q1 + 2];
      dG[q1 + 2] = k * (G[q1 + 1] + 0.5 * dl * y[t1 + 1]) - k * rel * G[q1 + 2];
      dA[q1 + 3] = y[t1 + 3] + c * rel * A[q1 + 2] - c * A[q1 + 3];
      dG[q1 + 3] = c * rel * G[q1 + 2] - c * G[q1 + 3];
    };
    mchain(Q_TN1, TN1, p[KTR_N],  p[C_N],  f_n,  Q_GMP, gmp, p[D_N], rel_n);
    mchain(Q_TM1, TM1, p[KTR_MO], p[C_MO], f_mo, Q_GMP, gmp, p[D_MO], 1.0);
    mchain(Q_TP1, TP1, p[KTR_PL], p[C_PL], f_pl, Q_MEP, mep, p[D_PL], 1.0);
  }

  // Pre-treatment steady state; moments from the closed forms
  //   renewing compartment: age = a_in + 1/(k - sr), generation = g_in + inc + 2 sr/(k - sr)
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
    if (!track) return;
    double* A = y.data() + N_STATE;
    double* G = y.data() + N_STATE + N_TRK;
    auto renewing = [&](int q, double k, double sr, double aU, double gU, double inc) {
      A[q] = aU + 1.0 / (k - sr);
      G[q] = gU + inc + 2.0 * sr / (k - sr);
    };
    renewing(Q_MPP, p[K_MPP], p[SR_MPP], 0.0, 0.0, 0.0);
    renewing(Q_CMP, p[K_CMP], p[SR_CMP], A[Q_MPP], G[Q_MPP], 1.0);
    renewing(Q_GMP, p[K_GMP], p[SR_GMP], A[Q_CMP], G[Q_CMP], 1.0);
    renewing(Q_MEP, p[K_MEP], p[SR_MEP], A[Q_CMP], G[Q_CMP], 1.0);
    renewing(Q_CLP, p[K_CLP], p[SR_CLP], A[Q_MPP], G[Q_MPP], 1.0);
    renewing(Q_LY,  p[K_LY],  p[SR_LY],  A[Q_CLP], G[Q_CLP], 1.0);
    A[Q_BL] = A[Q_LY] + 1.0 / p[C_L];
    G[Q_BL] = G[Q_LY];
    auto mchain = [&](int q1, double k, double c, int qp, double dl) {
      A[q1] = A[qp] + 1.0 / k;          G[q1] = G[qp] + 1.0;
      A[q1 + 1] = A[q1] + 1.0 / k;      G[q1 + 1] = G[q1] + 0.5 * dl;
      A[q1 + 2] = A[q1 + 1] + 1.0 / k;  G[q1 + 2] = G[q1 + 1] + 0.5 * dl;
      A[q1 + 3] = A[q1 + 2] + 1.0 / c;  G[q1 + 3] = G[q1 + 2];
    };
    mchain(Q_TN1, p[KTR_N],  p[C_N],  Q_GMP, p[D_N]);
    mchain(Q_TM1, p[KTR_MO], p[C_MO], Q_GMP, p[D_MO]);
    mchain(Q_TP1, p[KTR_PL], p[C_PL], Q_MEP, p[D_PL]);
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

} // namespace a4

// [[Rcpp::export]]
CharacterVector a4_param_names() {
  return CharacterVector::create(
    "r_s", "f0", "th_as", "th_ap", "th_bs",
    "alpha",
    "k_mpp", "sr_mpp", "k_cmp", "sr_cmp", "k_gmp", "sr_gmp",
    "k_mep", "sr_mep", "k_clp", "sr_clp", "k_ly", "sr_ly",
    "gmax", "gamma_m", "fb_mpp", "fb_cmp", "fb_gmp", "fb_mep", "fb_clp", "fb_ly",
    "r_a", "r_b", "delta",
    "ktr_n", "ktr_mo", "ktr_pl",
    "c_n", "c_mo", "c_pl", "c_l",
    "gam_n", "gam_mo", "gam_pl",
    "k_i", "emax", "ic50", "hill", "e_mo", "e_pl", "e_l", "eta_syn",
    "e_sr",
    "th_al",
    "k_e", "gcsf_fold", "tpo_fold",
    "d_n", "d_mo", "d_pl",
    "rho", "a0",
    "k_t", "k_ag", "kappa", "k_z", "z50", "th_zs", "e_z",
    "k_g", "gcsf_rel");
}

// [[Rcpp::export]]
CharacterVector a4_state_names() {
  return CharacterVector::create("SW", "SX", "MPP", "CMP", "GMP", "MEP", "CLP", "L", "C", "I1", "I",
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

// Simulate one patient with the parent model.
// seg   : matrix with columns start, end, dose, gcsf, tpo, s_add, car (rows ordered in time)
// drv   : D(t) on the grid drv_t0 + (0:(n-1)) * drv_dt
// obs_t : ascending observation times (>= first segment start)
// track : also integrate the age / generation moments
// value : columns t, anc_rel, plt_rel, amc_rel, S, f_clone, M (= MPP), C, D, I,
//         alc_rel, CMP, GMP, MEP, CLP, L, cart (CAR T cells in blood, 10^9/L),
//         ag (antigen load), z (CAR T - target signal Z), gcsf (G-CSF effect);
//         with track also age_anc, age_plt, age_amc, age_alc (mean days since MPP
//         entry) and gen_anc, gen_plt, gen_amc, gen_alc (mean divisions since MPP
//         entry) of the blood cells (all NA if the integrator failed)
// [[Rcpp::export]]
NumericMatrix a4_sim(NumericVector par, NumericMatrix seg, NumericVector drv, double drv_t0,
                     double drv_dt, NumericVector obs_t, double rtol = 1e-6, double atol = 1e-9,
                     double hmax = 0.5, int max_steps = 200000, bool track = false) {
  if (par.size() != a4::N_PAR) stop("parameter vector has wrong length");
  a4::check_inputs(seg, drv);
  a4::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  a4::Parent sys{par.begin(), &dr, a4::Inputs{0, 0, 0}, track};
  a4::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  const int nc = track ? 28 : 20;
  NumericMatrix out(obs_t.size(), nc);
  CharacterVector cn = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                               "M", "C", "D", "I", "alc_rel", "CMP", "GMP", "MEP",
                                               "CLP", "L", "cart", "ag", "z");
  cn.push_back("gcsf");
  if (track) {
    const char* extra[] = {"age_anc", "age_plt", "age_amc", "age_alc", "gen_anc", "gen_plt", "gen_amc", "gen_alc"};
    for (int j = 0; j < 8; ++j) cn.push_back(extra[j]);
  }
  colnames(out) = cn;
  auto record = [&](int k, double t, const std::vector<double>& y) {
    const double s = y[a4::SW] + y[a4::SX];
    out(k, 0) = t;
    out(k, 1) = y[a4::BN];
    out(k, 2) = y[a4::BP];
    out(k, 3) = y[a4::BM];
    out(k, 4) = s;
    out(k, 5) = (s > 1e-12) ? y[a4::SX] / s : 0.0;
    out(k, 6) = y[a4::MPP];
    out(k, 7) = y[a4::CX];
    out(k, 8) = dr.at(t);
    out(k, 9) = y[a4::I2];
    out(k, 10) = y[a4::BL];
    out(k, 11) = y[a4::CMP];
    out(k, 12) = y[a4::GMP];
    out(k, 13) = y[a4::MEP];
    out(k, 14) = y[a4::CLP];
    out(k, 15) = y[a4::LY];
    out(k, 16) = y[a4::CT];
    out(k, 17) = y[a4::AG];
    out(k, 18) = y[a4::Z2];
    out(k, 19) = y[a4::GC];
    if (track) {
      const int A0 = a4::N_STATE, G0 = a4::N_STATE + a4::N_TRK;
      const int qb[4] = {a4::Q_BN, a4::Q_BP, a4::Q_BM, a4::Q_BL};
      const int xb[4] = {a4::BN, a4::BP, a4::BM, a4::BL};
      for (int j = 0; j < 4; ++j) {
        const double x = std::max(y[xb[j]], 1e-12);
        out(k, 20 + j) = y[A0 + qb[j]] / x;
        out(k, 24 + j) = y[G0 + qb[j]] / x;
      }
    }
  };
  if (!a4::run_segments(sys, seg, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Simulate one patient with the Friberg reference model (same first ten output
// columns; stem-layer columns are NA).
// [[Rcpp::export]]
NumericMatrix mf_sim(NumericVector par, NumericMatrix seg, NumericVector obs_t, double rtol = 1e-6,
                     double atol = 1e-9, double hmax = 0.5, int max_steps = 200000) {
  if (par.size() != a4::MF_N_PAR) stop("MF parameter vector has wrong length");
  NumericVector dummy = NumericVector::create(0.0, 0.0);
  a4::check_inputs(seg, dummy);
  a4::Driver dr{0.0, 1.0, dummy.begin(), 2};
  a4::Friberg sys{par.begin(), &dr, a4::Inputs{0, 0, 0}};
  a4::Ctl ctl{rtol, atol, hmax, max_steps, 0, false};
  NumericMatrix out(obs_t.size(), 10);
  colnames(out) = CharacterVector::create("t", "anc_rel", "plt_rel", "amc_rel", "S", "f_clone",
                                          "M", "C", "D", "I");
  auto record = [&](int k, double t, const std::vector<double>& y) {
    out(k, 0) = t;
    out(k, 1) = y[a4::F_BN];
    out(k, 2) = y[a4::F_BP];
    out(k, 3) = y[a4::F_BM];
    out(k, 4) = NA_REAL;
    out(k, 5) = NA_REAL;
    out(k, 6) = y[a4::F_PN];
    out(k, 7) = y[a4::F_CX];
    out(k, 8) = 0.0;
    out(k, 9) = 0.0;
  };
  // the CD34 and CAR T columns have no meaning in MF (its event() ignores them)
  if (!a4::run_segments(sys, seg, obs_t, ctl, record)) std::fill(out.begin(), out.end(), NA_REAL);
  return out;
}

// Right-hand side of the parent model (29 states, no moments), exposed for
// validation against deSolve and for the steady-state Jacobian.
// [[Rcpp::export]]
NumericVector a4_rhs(double t, NumericVector y, NumericVector par, double dose, double gcsf, double tpo,
                     NumericVector drv, double drv_t0, double drv_dt) {
  if (par.size() != a4::N_PAR) stop("parameter vector has wrong length");
  if (y.size() != a4::N_STATE) stop("state vector has wrong length");
  a4::Driver dr{drv_t0, drv_dt, drv.begin(), static_cast<int>(drv.size())};
  a4::Parent sys{par.begin(), &dr, a4::Inputs{dose, gcsf, tpo}, false};
  NumericVector dy(a4::N_STATE);
  sys.rhs(t, y.begin(), dy.begin());
  return dy;
}
