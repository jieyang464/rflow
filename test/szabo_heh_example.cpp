#include "molecule.h"
#include "integral_provider.h"
#include "matrix.h"
#include "scf.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

scfcpp::DenseMatrix make_2x2(double a00, double a01, double a10, double a11) {
  scfcpp::DenseMatrix m(2, 2);
  m(0, 0) = a00;
  m(0, 1) = a01;
  m(1, 0) = a10;
  m(1, 1) = a11;
  return m;
}

double max_abs_diff(const scfcpp::DenseMatrix& a, const scfcpp::DenseMatrix& b) {
  double max_diff = 0.0;
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      double diff = std::abs(a(i, j) - b(i, j));
      if (diff > max_diff) max_diff = diff;
    }
  }
  return max_diff;
}

size_t eri_index(size_t mu, size_t nu, size_t lam, size_t sig, size_t nbf) {
  return ((mu * nbf + nu) * nbf + lam) * nbf + sig;
}

void set_eri_sym(std::vector<double>& eri, size_t nbf,
                 size_t mu, size_t nu, size_t lam, size_t sig, double val) {
  size_t a0 = mu, a1 = nu;
  size_t b0 = lam, b1 = sig;
  size_t a_vals[2] = {a0, a1};
  size_t b_vals[2] = {b0, b1};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      size_t p0 = a_vals[i];
      size_t p1 = a_vals[1 - i];
      size_t q0 = b_vals[j];
      size_t q1 = b_vals[1 - j];
      eri[eri_index(p0, p1, q0, q1, nbf)] = val;
      eri[eri_index(q0, q1, p0, p1, nbf)] = val;
    }
  }
}

bool check_eri(const std::vector<double>& eri, size_t nbf,
               size_t mu, size_t nu, size_t lam, size_t sig,
               double expected, double tol) {
  double v = eri[eri_index(mu, nu, lam, sig, nbf)];
  return std::abs(v - expected) <= tol;
}

}

int main() {
  Molecule mol({{2, 0.0, 0.0, 0.0, ""}, {1, 0.0, 0.0, 1.4632, ""}}, +1, "sto-3g");
  mol.set_multiplicity(1);

  const scfcpp::DenseMatrix S = make_2x2(1.0, 0.4508, 0.4508, 1.0);
  const scfcpp::DenseMatrix T = make_2x2(2.1643, 0.1670, 0.1670, 0.7600);
  const scfcpp::DenseMatrix V1 = make_2x2(-4.1398, -1.1029, -1.1029, -1.2652);
  const scfcpp::DenseMatrix V2 = make_2x2(-0.6772, -0.4113, -0.4113, -1.2266);
  const scfcpp::DenseMatrix Hcore = T + V1 + V2;
  const scfcpp::DenseMatrix H_expected = make_2x2(-2.6527, -1.3472, -1.3472, -1.7318);

  std::vector<double> eri(2 * 2 * 2 * 2, 0.0);
  set_eri_sym(eri, 2, 0, 0, 0, 0, 1.3072);
  set_eri_sym(eri, 2, 1, 1, 0, 0, 0.6057);
  set_eri_sym(eri, 2, 1, 0, 0, 0, 0.4373);
  set_eri_sym(eri, 2, 1, 1, 1, 0, 0.3118);
  set_eri_sym(eri, 2, 1, 0, 0, 1, 0.1773);
  set_eri_sym(eri, 2, 1, 1, 1, 1, 0.7746);

  IntegralFunctions fns;
  fns.nbf = [](const IMolecule&) { return static_cast<size_t>(2); };
  fns.overlap = [S](const IMolecule&) { return scfcpp::to_matrix(S); };
  fns.kinetic = [T](const IMolecule&) { return scfcpp::to_matrix(T); };
  fns.nuclear_attraction = [V1, V2](const IMolecule&) { return scfcpp::to_matrix(V1 + V2); };
  fns.hcore = [Hcore](const IMolecule&) { return scfcpp::to_matrix(Hcore); };
  fns.eri_tensor = [eri](const IMolecule&) { return eri; };
  FunctionIntegralProvider provider(std::move(fns));

  const double kTol = 5e-4;
  scfcpp::DenseMatrix H_from = scfcpp::to_dense(provider.hcore(mol));
  if (max_abs_diff(H_from, H_expected) > kTol) {
    std::cerr << "Hcore mismatch (max diff > " << kTol << ")\n";
    return 1;
  }

  const std::vector<double> eri_out = provider.eri_tensor(mol);
  if (!check_eri(eri_out, 2, 0, 0, 0, 0, 1.3072, kTol)) return 1;
  if (!check_eri(eri_out, 2, 1, 1, 0, 0, 0.6057, kTol)) return 1;
  if (!check_eri(eri_out, 2, 1, 0, 0, 0, 0.4373, kTol)) return 1;
  if (!check_eri(eri_out, 2, 1, 1, 1, 0, 0.3118, kTol)) return 1;
  if (!check_eri(eri_out, 2, 1, 0, 0, 1, 0.1773, kTol)) return 1;
  if (!check_eri(eri_out, 2, 1, 1, 1, 1, 0.7746, kTol)) return 1;

  XCConfig xc;
  xc.exchange = "HF";
  xc.correlation = "NONE";
  xc.exact_exchange_fraction = 1.0;
  SCFEngine eng(mol, provider, SCFEngine::Mode::RKS, xc, 1e-10, 100, 0.0);
  eng.set_commutator_cadence(0);
  eng.run();
  if (!eng.converged()) {
    std::cerr << "SCF did not converge\n";
    return 1;
  }

  const double R = 1.4632;
  const double total_energy = eng.energy() + 2.0 / R;
  const double expected_total = -2.860662;
  if (std::abs(total_energy - expected_total) > 5e-3) {
    std::cerr << "Total energy mismatch: " << total_energy << " vs " << expected_total << "\n";
    return 1;
  }

  return 0;
}
