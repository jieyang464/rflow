#include "scf.h"
#include "molecule.h"
#include "integral_provider.h"
#include "commutator.h"
#include "vxc_evaluator.h"
#include "vxc_libxc_grid.h"
#include <cmath>
#include <random>
#include <algorithm>

SCFEngine::SCFEngine(const IMolecule& mol,
                     const IIntegralProvider& integrals,
                     Mode mode,
                     const XCConfig& xc,
                     double conv_thresh,
                     int max_iter,
                     double init_perturb)
  : mol_(&mol), integrals_(&integrals), mode_(mode), xc_(xc),
    conv_thresh_(conv_thresh), max_iter_(max_iter), init_perturb_(init_perturb) {
  const size_t n = nbf();
  Fa_.resize(n, n);
  Fb_.resize(n, n);
  Ca_.resize(n, n);
  Cb_.resize(n, n);
  Pa_.resize(n, n);
  Pb_.resize(n, n);
  P_.resize(n, n);
  J_.resize(n, n);
  Ka_.resize(n, n);
  Kb_.resize(n, n);
  Fa_.setZero();
  Fb_.setZero();
  Ca_.setZero();
  Cb_.setZero();
  Pa_.setZero();
  Pb_.setZero();
  P_.setZero();
  J_.setZero();
  Ka_.setZero();
  Kb_.setZero();
  vxc_functor_ = make_libxc_vxc_functor();

  std::pair<int, int> ab = mol_->alpha_beta_electrons();
  int N = mol_->electron_count();
  if (mode_ == Mode::RKS && ((N % 2 != 0) || mol_->multiplicity() != 1)) {
    mode_ = Mode::UKS;
  }
  if (mode_ == Mode::UKS) {
    nelec_alpha_ = ab.first;
    nelec_beta_  = ab.second;
  } else {
    nelec_alpha_ = nelec_beta_ = N / 2;
  }
}

SCFEngine::~SCFEngine() = default;

size_t SCFEngine::nbf() const { return integrals_->nbf(*mol_); }

const IMolecule& SCFEngine::molecule() const { return *mol_; }

// Build spin density: P = sum_i C(:,i) C(:,i)^T for nocc orbitals (occupancy 1)
static void build_density_from_C(const scfcpp::DenseMatrix& C, int nocc, scfcpp::DenseMatrix& P) {
  const size_t nbf = C.rows();
  P.setZero(nbf, nbf);
  if (nocc <= 0) return;
  scfcpp::DenseMatrix Cocc = C.leftCols(nocc);
  P = Cocc * Cocc.transpose();
}

inline size_t SCFEngine::eri_index(size_t mu, size_t nu, size_t lam, size_t sig, size_t n) {
  return ((mu * n + nu) * n + lam) * n + sig;
}

void SCFEngine::build_J(const std::vector<double>& eri, const scfcpp::Matrix& P, size_t nbf,
                        scfcpp::Matrix& J) {
  scfcpp::DenseMatrix Pd = scfcpp::to_dense(P);
  scfcpp::DenseMatrix Jd(nbf, nbf);
  Jd.setZero();
  for (size_t mu = 0; mu < nbf; ++mu) {
    for (size_t nu = 0; nu < nbf; ++nu) {
      double acc = 0.0;
      for (size_t lam = 0; lam < nbf; ++lam) {
        for (size_t sig = 0; sig < nbf; ++sig) {
          acc += Pd(lam, sig) * eri[eri_index(mu, nu, lam, sig, nbf)];
        }
      }
      Jd(mu, nu) = acc;
    }
  }
  J = scfcpp::to_matrix(Jd);
}

void SCFEngine::build_K(const std::vector<double>& eri, const scfcpp::Matrix& Pspin, size_t nbf,
                        scfcpp::Matrix& K) {
  scfcpp::DenseMatrix Pd = scfcpp::to_dense(Pspin);
  scfcpp::DenseMatrix Kd(nbf, nbf);
  Kd.setZero();
  for (size_t mu = 0; mu < nbf; ++mu) {
    for (size_t nu = 0; nu < nbf; ++nu) {
      double acc = 0.0;
      for (size_t lam = 0; lam < nbf; ++lam) {
        for (size_t sig = 0; sig < nbf; ++sig) {
          acc += Pd(lam, sig) * eri[eri_index(mu, lam, nu, sig, nbf)];
        }
      }
      Kd(mu, nu) = acc;
    }
  }
  K = scfcpp::to_matrix(Kd);
}

scfcpp::DenseMatrix SCFEngine::exp_taylor(const scfcpp::DenseMatrix& M, int order) {
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

scfcpp::DenseMatrix SCFEngine::solve_FC_equals_SCe(const scfcpp::DenseMatrix& F, const scfcpp::DenseMatrix& S) {
  Eigen::GeneralizedSelfAdjointEigenSolver<scfcpp::DenseMatrix> ges(F, S);
  return ges.eigenvectors();
}

void SCFEngine::purify_S_metric(scfcpp::DenseMatrix& P, const scfcpp::DenseMatrix& S, int steps) {
  for (int i = 0; i < steps; ++i) {
    scfcpp::DenseMatrix PSP = P * S * P;
    scfcpp::DenseMatrix PSPSP = PSP * S * P;
    P = 3.0 * PSP - 2.0 * PSPSP;
    P = 0.5 * (P + P.transpose());
  }
}

double SCFEngine::compute_energy() const {
  scfcpp::DenseMatrix h = scfcpp::to_dense(const_cast<SCFEngine*>(this)->get_hcore());
  if (mode_ == Mode::UKS) {
    scfcpp::DenseMatrix Ha = h + scfcpp::to_dense(Fa_);
    scfcpp::DenseMatrix Hb = h + scfcpp::to_dense(Fb_);
    scfcpp::DenseMatrix Pa = scfcpp::to_dense(Pa_);
    scfcpp::DenseMatrix Pb = scfcpp::to_dense(Pb_);
    double Ea = 0.5 * (Pa.cwiseProduct(Ha)).sum();
    double Eb = 0.5 * (Pb.cwiseProduct(Hb)).sum();
    return Ea + Eb + (last_E_xc_ - last_tr_PVxc_);
  } else {
    scfcpp::DenseMatrix F = scfcpp::to_dense(Fa_);
    scfcpp::DenseMatrix sumHF = h + F;
    scfcpp::DenseMatrix P = scfcpp::to_dense(P_);
    return 0.5 * (P.cwiseProduct(sumHF)).sum() + (last_E_xc_ - last_tr_PVxc_);
  }
}

scfcpp::Matrix SCFEngine::get_hcore() const { return integrals_->hcore(*mol_); }
std::vector<double> SCFEngine::get_eri() const { return integrals_->eri_tensor(*mol_); }

bool SCFEngine::check_convergence() const {
  if (iter_ == 0) return false;
  double dE = std::abs(energy_ - prev_energy_);
  scfcpp::DenseMatrix Pdiff = scfcpp::to_dense(P_) - scfcpp::to_dense(prev_P_);
  double dP_rms = Pdiff.norm() / std::sqrt(static_cast<double>(Pdiff.size()));
  return (dE < conv_thresh_) && (dP_rms < conv_thresh_);
}

bool SCFEngine::iterate() {
  const size_t n = nbf();
  prev_P_ = P_;
  prev_energy_ = energy_;

  scfcpp::Matrix h = get_hcore();
  std::vector<double> eri = get_eri();
  P_ = Pa_ + Pb_;
  SCFEngine::build_J(eri, P_, n, J_);
  if (mode_ == Mode::RKS) {
    SCFEngine::build_K(eri, 0.5 * P_, n, Ka_);
    Kb_ = Ka_;
  } else {
    SCFEngine::build_K(eri, Pa_, n, Ka_);
    SCFEngine::build_K(eri, Pb_, n, Kb_);
  }

  scfcpp::Matrix Va(n, n), Vb(n, n);
  double E_xc = 0.0, tr_PVxc = 0.0;
  compute_Vxc(Pa_, Pb_, Va, Vb, E_xc, tr_PVxc);
  double ax = std::max(0.0, std::min(1.0, xc_.exact_exchange_fraction));
  scfcpp::DenseMatrix Fa_dense = scfcpp::to_dense(h) + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Ka_) + scfcpp::to_dense(Va);
  scfcpp::DenseMatrix Fb_dense = scfcpp::to_dense(h) + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Kb_) + scfcpp::to_dense(Vb);
  Fa_ = scfcpp::to_matrix(Fa_dense);
  Fb_ = scfcpp::to_matrix(Fb_dense);
  last_E_xc_ = E_xc;
  last_tr_PVxc_ = tr_PVxc;

  scfcpp::DenseMatrix S = scfcpp::to_dense(integrals_->overlap(*mol_));
  if (mode_ == Mode::RKS) {
    scfcpp::DenseMatrix F = scfcpp::to_dense(Fa_);
    Ca_ = Cb_ = SCFEngine::solve_FC_equals_SCe(F, S);
  } else {
    Ca_ = SCFEngine::solve_FC_equals_SCe(scfcpp::to_dense(Fa_), S);
    Cb_ = SCFEngine::solve_FC_equals_SCe(scfcpp::to_dense(Fb_), S);
  }

  scfcpp::DenseMatrix Pa_dense;
  scfcpp::DenseMatrix Pb_dense;
  build_density_from_C(Ca_, nelec_alpha_, Pa_dense);
  build_density_from_C(Cb_, nelec_beta_,  Pb_dense);
  Pa_ = scfcpp::to_matrix(Pa_dense);
  Pb_ = scfcpp::to_matrix(Pb_dense);
  P_ = Pa_ + Pb_;

  energy_ = compute_energy();
  ++iter_;
  converged_ = check_convergence();
  last_step_was_diag_ = true;
  return converged_;
}

void SCFEngine::run() {
  initialize_from_hcore();
  while (iter_ < max_iter_ && !converged_) {
    for (int k = 0; k < comm_steps_default_ && iter_ < max_iter_ && !converged_; ++k) {
      commutator_step(comm_opts_);
      if (check_convergence()) { converged_ = true; break; }
    }
    if (converged_ || iter_ >= max_iter_) break;
    if (iterate()) { converged_ = true; break; }
  }
  if (converged_ && !last_step_was_diag_) update_C_from_current_fock();
}

void SCFEngine::commutator_step(double alpha, int exp_order, bool do_purify, int purify_steps) {
  CommOptions opts;
  opts.alpha = alpha;
  opts.exp_order = exp_order;
  opts.do_purify = do_purify;
  opts.purify_steps = purify_steps;
  commutator_step(opts);
}

void SCFEngine::commutator_step(const CommOptions& opts) {
  const size_t n = nbf();
  prev_P_ = P_;
  prev_energy_ = energy_;

  scfcpp::Matrix h = get_hcore();
  std::vector<double> eri = get_eri();
  P_ = Pa_ + Pb_;
  SCFEngine::build_J(eri, P_, n, J_);
  if (mode_ == Mode::RKS) {
    SCFEngine::build_K(eri, 0.5 * P_, n, Ka_);
    Kb_ = Ka_;
  } else {
    SCFEngine::build_K(eri, Pa_, n, Ka_);
    SCFEngine::build_K(eri, Pb_, n, Kb_);
  }

  scfcpp::Matrix Va(n, n), Vb(n, n);
  double E_xc = 0.0, tr_PVxc = 0.0;
  compute_Vxc(Pa_, Pb_, Va, Vb, E_xc, tr_PVxc);
  double ax = std::max(0.0, std::min(1.0, xc_.exact_exchange_fraction));
  scfcpp::DenseMatrix Fa_dense = scfcpp::to_dense(h) + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Ka_) + scfcpp::to_dense(Va);
  scfcpp::DenseMatrix Fb_dense = scfcpp::to_dense(h) + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Kb_) + scfcpp::to_dense(Vb);
  Fa_ = scfcpp::to_matrix(Fa_dense);
  Fb_ = scfcpp::to_matrix(Fb_dense);
  last_E_xc_ = E_xc;
  last_tr_PVxc_ = tr_PVxc;

  scfcpp::DenseMatrix S = scfcpp::to_dense(integrals_->overlap(*mol_));
  commutator::Options copts;
  copts.alpha = opts.alpha;
  copts.exp_order = opts.exp_order;
  copts.bch_order = opts.bch_order;
  copts.do_purify = opts.do_purify;
  copts.purify_steps = opts.purify_steps;
  copts.method = opts.method;
  copts.conjugator = opts.conjugator;
  if (mode_ == Mode::RKS) {
    scfcpp::DenseMatrix Ptot = scfcpp::to_dense(P_);
    commutator::update_step_rhf(scfcpp::to_dense(Fa_), S, copts, Ptot);
    Pa_ = scfcpp::to_matrix(0.5 * Ptot);
    Pb_ = Pa_;
  } else {
    scfcpp::DenseMatrix Pa_dense = scfcpp::to_dense(Pa_);
    scfcpp::DenseMatrix Pb_dense = scfcpp::to_dense(Pb_);
    commutator::update_step_uhf(scfcpp::to_dense(Fa_), scfcpp::to_dense(Fb_), S, copts, Pa_dense, Pb_dense);
    Pa_ = scfcpp::to_matrix(Pa_dense);
    Pb_ = scfcpp::to_matrix(Pb_dense);
  }
  P_ = Pa_ + Pb_;
  energy_ = compute_energy();
  ++iter_;
  converged_ = check_convergence();
  last_step_was_diag_ = false;
}

void SCFEngine::run_mixed(int comm_steps, double alpha, int exp_order, bool do_purify) {
  CommOptions opts;
  opts.alpha = alpha;
  opts.exp_order = exp_order;
  opts.do_purify = do_purify;
  run_mixed(comm_steps, opts);
}

void SCFEngine::run_mixed(int comm_steps, const CommOptions& opts) {
  initialize_from_hcore();
  for (; iter_ < max_iter_ && !converged_; ) {
    for (int k = 0; k < comm_steps && !converged_ && iter_ < max_iter_; ++k) {
      commutator_step(opts);
      if (check_convergence()) { converged_ = true; break; }
    }
    if (converged_ || iter_ >= max_iter_) break;
    if (iterate()) { converged_ = true; break; }
  }
  if (converged_ && !last_step_was_diag_) update_C_from_current_fock();
}

void SCFEngine::run_commutator_only(const CommOptions& opts) {
  initialize_from_hcore();
  while (iter_ < max_iter_ && !converged_) {
    commutator_step(opts);
  }
  if (converged_ && !last_step_was_diag_) update_C_from_current_fock();
}

void SCFEngine::initialize_from_hcore() {
  const int n = static_cast<int>(nbf());
  std::pair<int, int> ab = mol_->alpha_beta_electrons();
  int N = mol_->electron_count();
  if (mode_ == Mode::RKS && ((N % 2 != 0) || mol_->multiplicity() != 1)) mode_ = Mode::UKS;
  if (mode_ == Mode::UKS) { nelec_alpha_ = ab.first; nelec_beta_ = ab.second; }
  else { nelec_alpha_ = nelec_beta_ = N/2; }

  scfcpp::DenseMatrix h = scfcpp::to_dense(get_hcore());
  scfcpp::DenseMatrix S = scfcpp::to_dense(integrals_->overlap(*mol_));
  scfcpp::DenseMatrix C0 = SCFEngine::solve_FC_equals_SCe(h, S);
  Ca_ = Cb_ = C0;
  if (mode_ == Mode::RKS) {
    scfcpp::DenseMatrix Phalf;
    build_density_from_C(Ca_, nelec_alpha_, Phalf);
    Pa_ = scfcpp::to_matrix(Phalf);
    Pb_ = Pa_;
  } else {
    scfcpp::DenseMatrix Pa_dense;
    scfcpp::DenseMatrix Pb_dense;
    build_density_from_C(Ca_, nelec_alpha_, Pa_dense);
    build_density_from_C(Cb_, nelec_beta_,  Pb_dense);
    Pa_ = scfcpp::to_matrix(Pa_dense);
    Pb_ = scfcpp::to_matrix(Pb_dense);
  }
  P_ = Pa_ + Pb_;

  std::vector<double> eri = get_eri();
  SCFEngine::build_J(eri, P_, n, J_);
  if (mode_ == Mode::RKS) {
    SCFEngine::build_K(eri, 0.5 * P_, n, Ka_);
    Kb_ = Ka_;
  } else {
    SCFEngine::build_K(eri, Pa_, n, Ka_);
    SCFEngine::build_K(eri, Pb_, n, Kb_);
  }

  scfcpp::Matrix Va(n, n), Vb(n, n);
  double E_xc = 0.0, tr_PVxc = 0.0;
  compute_Vxc(Pa_, Pb_, Va, Vb, E_xc, tr_PVxc);
  double ax = std::max(0.0, std::min(1.0, xc_.exact_exchange_fraction));
  scfcpp::DenseMatrix Fa_dense = h + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Ka_) + scfcpp::to_dense(Va);
  scfcpp::DenseMatrix Fb_dense = h + scfcpp::to_dense(J_) - ax * scfcpp::to_dense(Kb_) + scfcpp::to_dense(Vb);
  Fa_ = scfcpp::to_matrix(Fa_dense);
  Fb_ = scfcpp::to_matrix(Fb_dense);
  last_E_xc_ = E_xc;
  last_tr_PVxc_ = tr_PVxc;

  energy_ = compute_energy();
  prev_energy_ = energy_;
  prev_P_ = P_;
  iter_ = 0;
  converged_ = false;
}

void SCFEngine::update_C_from_current_fock() {
  scfcpp::DenseMatrix S = scfcpp::to_dense(integrals_->overlap(*mol_));
  if (mode_ == Mode::RKS) {
    Ca_ = Cb_ = SCFEngine::solve_FC_equals_SCe(scfcpp::to_dense(Fa_), S);
  } else {
    Ca_ = SCFEngine::solve_FC_equals_SCe(scfcpp::to_dense(Fa_), S);
    Cb_ = SCFEngine::solve_FC_equals_SCe(scfcpp::to_dense(Fb_), S);
  }
}

void SCFEngine::compute_Vxc(const scfcpp::Matrix& Pa,
                            const scfcpp::Matrix& Pb,
                            scfcpp::Matrix& Va,
                            scfcpp::Matrix& Vb,
                            double& E_xc,
                            double& tr_PVxc) const {
  const size_t n = nbf();
  if (xc_.exchange == "HF" && (xc_.correlation == "NONE" || xc_.correlation.empty()) && xc_.exact_exchange_fraction >= 1.0 - 1e-12) {
    Va = scfcpp::Matrix(n, n);
    Vb = scfcpp::Matrix(n, n);
    Va.setZero();
    Vb.setZero();
    E_xc = 0.0;
    tr_PVxc = 0.0;
    return;
  }
  if (!vxc_functor_) {
    Va = scfcpp::Matrix(n, n);
    Vb = scfcpp::Matrix(n, n);
    Va.setZero();
    Vb.setZero();
    E_xc = 0.0;
    tr_PVxc = 0.0;
    return;
  }
  VxcInputs in;
  in.mol = mol_;
  in.integrals = integrals_;
  in.Pa = &Pa;
  in.Pb = &Pb;
  in.xc = xc_;
  VxcResult res = vxc_functor_(in);
  Va = std::move(res.Va);
  Vb = std::move(res.Vb);
  E_xc = res.Exc;
  tr_PVxc = res.tr_PVxc;
}
