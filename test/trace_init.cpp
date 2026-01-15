#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include <iostream>
#include <iomanip>

// Print initial Pa and Pb traces constructed from P_total/2 by diagonalizing H_core
// In AO metric (non-orthonormal), the correct electron count is Tr(P S).
// For HeH+ in STO-3G, N=2, so Tr(Pa S)=Tr(Pb S)=1 should hold for the initial RHF-like guess.

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

  auto integrals = make_libint_integral_provider();
  scfcpp::DenseMatrix h = scfcpp::to_dense(integrals->hcore(mol));
  scfcpp::DenseMatrix S = scfcpp::to_dense(integrals->overlap(mol));

  Eigen::GeneralizedSelfAdjointEigenSolver<scfcpp::DenseMatrix> ges(h, S);
  scfcpp::DenseMatrix C = ges.eigenvectors();

  int Zsum = 0; for (const auto& a : mol.atoms()) Zsum += a.Z;
  int N = Zsum - mol.net_charge();
  int nocc = N / 2;

  // Build RHF total density P = 2 * C_occ C_occ^T, then split Pa=Pb=P/2
  scfcpp::DenseMatrix Cocc = C.leftCols(nocc);
  scfcpp::DenseMatrix Ptot = 2.0 * (Cocc * Cocc.transpose());
  scfcpp::DenseMatrix Pa = 0.5 * Ptot;
  scfcpp::DenseMatrix Pb = 0.5 * Ptot;

  double trPa = Pa.trace();
  double trPb = Pb.trace();
  double trPaS = (Pa * S).trace();
  double trPbS = (Pb * S).trace();

  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(12);
  std::cout << "Trace(Pa)    = " << trPa  << "\n";
  std::cout << "Trace(Pb)    = " << trPb  << "\n";
  std::cout << "Trace(Pa S) = " << trPaS << "\n";
  std::cout << "Trace(Pb S) = " << trPbS << "\n";

  // Compute core-Hamiltonian energy (no two-electron terms):
  // Using RHF identity E_core = 0.5 * tr(P (h + F)), but for the core guess F = h, so E_core = tr(P h).
  double E_core_via_half = 0.5 * (Ptot.cwiseProduct(h + h)).sum();
  double E_core_direct    = (Ptot.cwiseProduct(h)).sum();
  std::cout << "E_core via 0.5*tr(P*(h+h)) = " << E_core_via_half << "\n";
  std::cout << "E_core via tr(P*h)        = " << E_core_direct    << "\n";

  // Optional: print matrices for inspection
  // print_matrix(Pa, "Pa (initial from Ptot/2)");
  // print_matrix(Pb, "Pb (initial from Ptot/2)");

  return 0;
}
