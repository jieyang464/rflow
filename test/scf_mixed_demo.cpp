#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include "scf.h"
#include <iostream>
#include <iomanip>

int main() {
  Atom He{2, 0.0, 0.0, 0.0};
  Atom H {1, 0.0, 0.0, 1.4632};
  Molecule mol({He, H}, +1, "sto-3g", "");

  auto integrals = make_libint_integral_provider();
  XCConfig xc; xc.exchange = "HF"; xc.correlation = "NONE"; xc.exact_exchange_fraction = 1.0;
  SCFEngine eng(mol, *integrals, SCFEngine::Mode::UKS, xc, 1e-8, 200, 1e-4);

  int comm_steps = 10;
  double alpha = 0.02;
  int exp_order = 6;
  bool purify = true;

  eng.run_mixed(comm_steps, alpha, exp_order, purify);
  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(12);
  std::cout << "Mixed SCF (HF) converged=" << std::boolalpha << eng.converged()
            << ", energy=" << eng.energy() << "\n";
  scfcpp::DenseMatrix Pa = scfcpp::to_dense(eng.Pa());
  scfcpp::DenseMatrix Pb = scfcpp::to_dense(eng.Pb());
  std::cout << "Trace Pa=" << Pa.trace() << ", Trace Pb=" << Pb.trace() << "\n";
  std::cout << "||Pa-Pb||_F=" << (Pa - Pb).norm() << "\n";
  return 0;
}
