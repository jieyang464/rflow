#include "commutator.h"
#include <algorithm>

namespace commutator {

static scfcpp::DenseMatrix exp_taylor(const scfcpp::DenseMatrix& M, int order) {
  const int n = static_cast<int>(M.rows());
  scfcpp::DenseMatrix U = scfcpp::DenseMatrix::Identity(n, n);
  scfcpp::DenseMatrix Mk = scfcpp::DenseMatrix::Identity(n, n);
  double fact = 1.0;
  for (int k = 1; k <= order; ++k) {
    Mk = Mk * M;
    fact *= static_cast<double>(k);
    U += (1.0 / fact) * Mk;
  }
  return U;
}

static scfcpp::DenseMatrix conjugate_bch(const scfcpp::DenseMatrix& R,
                                         const scfcpp::DenseMatrix& C,
                                         int order) {
  scfcpp::DenseMatrix result = R;
  scfcpp::DenseMatrix term = R;
  double fact = 1.0;
  for (int k = 1; k <= order; ++k) {
    term = term * C - C * term;
    fact *= static_cast<double>(k);
    result += (1.0 / fact) * term;
  }
  return result;
}

static void purify_S_metric(scfcpp::DenseMatrix& P, const scfcpp::DenseMatrix& S, int steps) {
  for (int i = 0; i < steps; ++i) {
    scfcpp::DenseMatrix PSP = P * S * P;
    scfcpp::DenseMatrix PSPSP = PSP * S * P;
    P = 3.0 * PSP - 2.0 * PSPSP;
    P = 0.5 * (P + P.transpose());
  }
}

static void apply_conjugation(scfcpp::DenseMatrix& R,
                              const scfcpp::DenseMatrix& C,
                              const scfcpp::DenseMatrix& S,
                              Eigen::LLT<scfcpp::DenseMatrix>& Sllt,
                              const Options& opts) {
  if (opts.method == ConjugationMethod::Injected && opts.conjugator) {
    opts.conjugator(C, S, opts, R);
    return;
  }
  if (opts.method == ConjugationMethod::BCH) {
    R = conjugate_bch(R, C, opts.bch_order);
    return;
  }
  scfcpp::DenseMatrix U = exp_taylor(-C, opts.exp_order);
  scfcpp::DenseMatrix UtS = U.transpose() * S;
  scfcpp::DenseMatrix Uinv = Sllt.solve(UtS);
  R = U * R * Uinv;
}

void update_step_rhf(const scfcpp::DenseMatrix& F,
                     const scfcpp::DenseMatrix& S,
                     const Options& opts,
                     scfcpp::DenseMatrix& P) {
  Eigen::LLT<scfcpp::DenseMatrix> Sllt(S);
  scfcpp::DenseMatrix Ph = 0.5 * P;
  scfcpp::DenseMatrix KS = F * Ph * S - S * Ph * F;
  scfcpp::DenseMatrix A = Sllt.solve(KS);
  scfcpp::DenseMatrix C = opts.alpha * A;
  apply_conjugation(Ph, C, S, Sllt, opts);
  if (opts.do_purify) purify_S_metric(Ph, S, opts.purify_steps);
  Ph = 0.5 * (Ph + Ph.transpose());
  P = 2.0 * Ph;
}

void update_step_uhf(const scfcpp::DenseMatrix& Fa,
                     const scfcpp::DenseMatrix& Fb,
                     const scfcpp::DenseMatrix& S,
                     const Options& opts,
                     scfcpp::DenseMatrix& Pa,
                     scfcpp::DenseMatrix& Pb) {
  Eigen::LLT<scfcpp::DenseMatrix> Sllt(S);
  struct UpdateSpin {
    Eigen::LLT<scfcpp::DenseMatrix>& Sllt;
    const scfcpp::DenseMatrix& S;
    const Options& opts;
    void operator()(scfcpp::DenseMatrix& Pspin, const scfcpp::DenseMatrix& Fspin) const {
      scfcpp::DenseMatrix KS = Fspin * Pspin * S - S * Pspin * Fspin;
      scfcpp::DenseMatrix A = Sllt.solve(KS);
      scfcpp::DenseMatrix C = opts.alpha * A;
      apply_conjugation(Pspin, C, S, Sllt, opts);
      if (opts.do_purify) purify_S_metric(Pspin, S, opts.purify_steps);
      Pspin = 0.5 * (Pspin + Pspin.transpose());
    }
  } update_spin{Sllt, S, opts};
  update_spin(Pa, Fa);
  update_spin(Pb, Fb);
}

}
