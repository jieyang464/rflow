#include "molecule.h"
#include "integral_provider.h"
#include "scf.h"
#include <cmath>
#include <iostream>

using ::Atom;
using ::Molecule;

int main() {
  std::vector<Atom> atoms = {
    {2, 0.0, 0.0, 0.0, ""},
    {1, 0.0, 0.0, 1.4632, ""}
  };
  Molecule mol(atoms, +1, "sto-3g");
  mol.set_multiplicity(1);

  auto integrals = make_libint_integral_provider();
  XCConfig hf_xc; hf_xc.exchange = "HF"; hf_xc.correlation = "NONE"; hf_xc.exact_exchange_fraction = 1.0;
  SCFEngine hf_eng(mol, *integrals, SCFEngine::Mode::UKS, hf_xc, 1e-8, 100, 1e-4);
  hf_eng.run();

  XCConfig dft_xc;
  dft_xc.exchange = "LDA_X";
  dft_xc.correlation = "LDA_C_VWN";
  dft_xc.exact_exchange_fraction = 0.0;
  SCFEngine dft_eng(mol, *integrals, SCFEngine::Mode::UKS, dft_xc, 1e-8, 100, 1e-4);
  dft_eng.run();

  std::cout << "HF energy        = " << hf_eng.energy() << "\n";
  std::cout << "DFT(LDA) E       = " << dft_eng.energy() << "\n";
  std::cout << "Converged(HF,DFT)= " << hf_eng.converged() << ", " << dft_eng.converged() << "\n";
  std::cout << "|E diff|         = " << std::abs(hf_eng.energy() - dft_eng.energy()) << "\n";
  return 0;
}
