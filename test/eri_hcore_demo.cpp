#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include <cmath>
#include <iostream>
#include <iomanip>

int main() {
  Atom He{2, 0.0, 0.0, 0.0};
  Atom H{1, 0.0, 0.0, 1.4632};
  Molecule mol({He, H}, +1, "sto-3g", "");

  try {
    auto integrals = make_libint_integral_provider();
    scfcpp::DenseMatrix hcore = scfcpp::to_dense(integrals->hcore(mol));
    auto eri = integrals->eri_tensor(mol);
    size_t nbf = integrals->nbf(mol);
    std::cout << "nbf=" << nbf << "\n";

    std::cout << "Hcore (row-major):\n";
    for (size_t i = 0; i < nbf; ++i) {
      for (size_t j = 0; j < nbf; ++j) std::cout << std::setw(15) << hcore(i, j);
      std::cout << '\n';
    }

    std::cout << "\nERI tensor (non-zero list)\n";
    for (size_t mu = 0; mu < nbf; ++mu)
      for (size_t nu = 0; nu < nbf; ++nu)
        for (size_t lam = 0; lam < nbf; ++lam)
          for (size_t sig = 0; sig < nbf; ++sig) {
            size_t idx = ((mu * nbf + nu) * nbf + lam) * nbf + sig;
            double v = eri[idx];
            if (std::abs(v) > 1e-12) {
              std::cout << "("<<mu<<","<<nu<<"|"<<lam<<","<<sig<<") = "<< v <<"\n";
            }
          }
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
  return 0;
}
