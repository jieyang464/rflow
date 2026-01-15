#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include "scf.h"
#include <iostream>
#include <iomanip>

static void print_matrix(const scfcpp::DenseMatrix& M, const std::string& name) {
  std::cout << name << " (" << M.rows() << "x" << M.cols() << ")\n";
  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(8);
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) std::cout << std::setw(14) << M(i,j);
    std::cout << '\n';
  }
}

int main() {
  Atom He{2, 0.0, 0.0, 0.0};
  Atom H {1, 0.0, 0.0, 1.4632};
  Molecule mol({He, H}, +1, "sto-3g", "");

  auto integrals = make_libint_integral_provider();
  XCConfig xc; xc.exchange = "HF"; xc.correlation = "NONE"; xc.exact_exchange_fraction = 1.0;
  SCFEngine hf_eng(mol, *integrals, SCFEngine::Mode::UKS, xc, /*conv_thresh*/1e-8, /*max_iter*/200, /*init_perturb*/0.0);

  const int N = 10;
  hf_eng.set_commutator_cadence(N);

  SCFEngine::CommOptions opts;
  opts.alpha = 0.02;
  opts.exp_order = 6;
  opts.do_purify = true;
  opts.purify_steps = 1;
  hf_eng.set_commutator_options(opts);

  hf_eng.initialize_from_hcore();

  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(12);

  int iter = 0;
  while (!hf_eng.converged() && iter < 200) {
    for (int k = 0; k < N && !hf_eng.converged() && iter < 200; ++k) {
      std::cout << "Iter " << iter+1 << " (commutator step)\n";
      std::cout << "Energy = " << hf_eng.compute_energy() << "\n";
      print_matrix(scfcpp::to_dense(hf_eng.Pa()), "Pa");
      print_matrix(scfcpp::to_dense(hf_eng.Pb()), "Pb");
      print_matrix(scfcpp::to_dense(hf_eng.Fa()), "Fa");
      print_matrix(scfcpp::to_dense(hf_eng.Fb()), "Fb");
      hf_eng.commutator_step(opts);
      ++iter;
    }
    if (hf_eng.converged() || iter >= 200) break;

    std::cout << "Iter " << iter+1 << " (diagonalization)\n";
    std::cout << "Energy = " << hf_eng.compute_energy() << "\n";
    print_matrix(scfcpp::to_dense(hf_eng.Pa()), "Pa");
    print_matrix(scfcpp::to_dense(hf_eng.Pb()), "Pb");
    print_matrix(scfcpp::to_dense(hf_eng.Fa()), "Fa");
    print_matrix(scfcpp::to_dense(hf_eng.Fb()), "Fb");
    hf_eng.iterate();
    ++iter;
  }

  std::cout << "Converged=" << std::boolalpha << hf_eng.converged() << ", total iters=" << iter
            << ", final energy=" << hf_eng.energy() << "\n";

  return 0;
}
