#include "molecule.h"
#include <stdexcept>

void Molecule::set_multiplicity(int mult) {
  if (mult <= 0) throw std::runtime_error("Multiplicity must be a positive integer (2S+1 >= 1)");
  multiplicity_ = mult;
}

int Molecule::electron_count() const {
  int Zsum = 0;
  for (const Atom& a : atoms_) Zsum += a.Z;
  return Zsum - net_charge_;
}

std::pair<int,int> Molecule::alpha_beta_electrons() const {
  int N = electron_count();
  int M = multiplicity_; // 2S+1
  int S2 = M - 1;        // 2S
  // Parity consistency: N and S2 must have same parity for Na, Nb to be integers
  if (((N - S2) & 1) != 0) {
    throw std::runtime_error("Inconsistent charge/multiplicity: N and (M-1) parity mismatch. Adjust net charge or multiplicity.");
  }
  int Na = (N + S2) / 2;
  int Nb = N - Na; // (N - S2)/2
  if (Na < 0 || Nb < 0) {
    throw std::runtime_error("Inconsistent charge/multiplicity: negative spin electron count. Adjust settings.");
  }
  return {Na, Nb};
}

void Molecule::set_uniform_basis(const std::string& name, const std::string& file_path) {
  basis_name_ = name;
  basis_file_path_ = file_path;
  for (Atom& a : atoms_) a.basis_label = name;
}
