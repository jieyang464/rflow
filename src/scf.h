#pragma once

#include "commutator.h"
#include "matrix.h"
#include "xc_config.h"
#include "vxc_evaluator.h"

class IMolecule;
class IIntegralProvider;

class SCFEngine {
public:
  enum class Mode { UKS, RKS };

  struct CommOptions {
    double alpha = 0.1;
    int    exp_order = 4;
    int    bch_order = 4;
    bool   do_purify = false;
    int    purify_steps = 1;
    commutator::ConjugationMethod method = commutator::ConjugationMethod::Taylor;
    commutator::ConjugationFn conjugator;
  };

  SCFEngine(const IMolecule& mol,
            const IIntegralProvider& integrals,
            Mode mode,
            const XCConfig& xc,
            double conv_thresh = 1e-8,
            int max_iter = 50,
            double init_perturb = 1e-4);
  ~SCFEngine();

  void run();
  void run_mixed(int comm_steps,
                 double alpha = 0.1,
                 int exp_order = 4,
                 bool do_purify = false);
  void run_mixed(int comm_steps, const CommOptions& opts);
  void run_commutator_only(const CommOptions& opts);

  bool iterate();
  void commutator_step(double alpha = 0.1, int exp_order = 4, bool do_purify = false, int purify_steps = 1);
  void commutator_step(const CommOptions& opts);

  void initialize_random_density(unsigned seed = 42);
  void initialize_from_hcore();

  double compute_energy() const;
  bool check_convergence() const;

  void set_commutator_cadence(int steps_between_diags) { comm_steps_default_ = steps_between_diags; }
  void set_commutator_options(const CommOptions& opts) { comm_opts_ = opts; }
  void set_vxc_functor(VxcFunctor fn) { vxc_functor_ = std::move(fn); }

  const scfcpp::Matrix& Fa() const { return Fa_; }
  const scfcpp::Matrix& Fb() const { return Fb_; }
  const scfcpp::DenseMatrix& Ca() const { return Ca_; }
  const scfcpp::DenseMatrix& Cb() const { return Cb_; }
  const scfcpp::Matrix& Pa() const { return Pa_; }
  const scfcpp::Matrix& Pb() const { return Pb_; }
  const scfcpp::Matrix& P()  const { return P_; }
  double energy() const { return energy_; }
  bool converged() const { return converged_; }
  size_t nbf() const;
  const IMolecule& molecule() const;
  int nelec_alpha() const { return nelec_alpha_; }
  int nelec_beta()  const { return nelec_beta_;  }
  Mode mode() const { return mode_; }

private:
  const IMolecule* mol_ = nullptr;
  const IIntegralProvider* integrals_ = nullptr;
  Mode mode_ = Mode::UKS;
  XCConfig xc_{};
  VxcFunctor vxc_functor_{};

  scfcpp::Matrix Fa_, Fb_, Pa_, Pb_, P_, J_, Ka_, Kb_;
  scfcpp::DenseMatrix Ca_, Cb_;
  double energy_ = 0.0;
  double prev_energy_ = 0.0;
  bool converged_ = false;
  double conv_thresh_ = 1e-8;
  int max_iter_ = 50;
  int iter_ = 0;
  int nelec_alpha_ = 0;
  int nelec_beta_  = 0;
  double init_perturb_ = 1e-4;
  scfcpp::Matrix prev_P_;
  bool last_step_was_diag_ = false;
  double last_E_xc_ = 0.0;
  double last_tr_PVxc_ = 0.0;

  int comm_steps_default_ = 10;
  CommOptions comm_opts_{};

  scfcpp::Matrix get_hcore() const;
  std::vector<double> get_eri() const;
  static inline size_t eri_index(size_t mu, size_t nu, size_t lam, size_t sig, size_t n);
  static void build_J(const std::vector<double>& eri, const scfcpp::Matrix& P, size_t nbf, scfcpp::Matrix& J);
  static void build_K(const std::vector<double>& eri, const scfcpp::Matrix& Pspin, size_t nbf, scfcpp::Matrix& K);
  static scfcpp::DenseMatrix exp_taylor(const scfcpp::DenseMatrix& M, int order);
  static void purify_S_metric(scfcpp::DenseMatrix& P, const scfcpp::DenseMatrix& S, int steps = 1);
  static scfcpp::DenseMatrix solve_FC_equals_SCe(const scfcpp::DenseMatrix& F, const scfcpp::DenseMatrix& S);
  void update_C_from_current_fock();
  void compute_Vxc(const scfcpp::Matrix& Pa,
                   const scfcpp::Matrix& Pb,
                   scfcpp::Matrix& Va,
                   scfcpp::Matrix& Vb,
                   double& E_xc,
                   double& tr_PVxc) const;
};
