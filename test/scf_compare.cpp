#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include "scf.h"
#include <iostream>
#include <iomanip>

static void print_matrix(const scfcpp::DenseMatrix& M, const std::string& name) {
  std::cout << name << " (" << M.rows() << "x" << M.cols() << ")\n";
  std::cout.setf(std::ios::fixed); std::cout << std::setprecision(8);
  for (int i = 0; i < M.rows(); ++i) {
    for (int j = 0; j < M.cols(); ++j) std::cout << std::setw(14) << M(i,j);
    std::cout << '\n';
  }
}

int main() {
  Atom He{2, 0.0, 0.0, 0.0};
  Atom H {1, 0.0, 0.0, 1.4632};
  Molecule mol({He, H}, +1, "sto-3g", "");

  {
    auto integrals = make_libint_integral_provider();
    XCConfig xc; xc.exchange = "HF"; xc.correlation = "NONE"; xc.exact_exchange_fraction = 1.0;
    SCFEngine eng(mol, *integrals, SCFEngine::Mode::UKS, xc, 1e-8, 100, 1e-4);
    eng.run();
    scfcpp::DenseMatrix Pa = scfcpp::to_dense(eng.Pa());
    scfcpp::DenseMatrix Pb = scfcpp::to_dense(eng.Pb());
    std::cout << "Standard HF (UHF)\n";
    std::cout << "converged=" << std::boolalpha << eng.converged() << ", energy=" << std::setprecision(12) << eng.energy() << "\n";
    std::cout << "||Pa-Pb||_F = " << (Pa - Pb).norm() << "\n";
    print_matrix(Pa, "Pa (standard)");
    print_matrix(Pb, "Pb (standard)");
  }

  return 0;
}
