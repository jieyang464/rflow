#include "commutator.h"
#include "integral_provider.h"
#include "matrix.h"
#include "molecule.h"
#include "scf.h"

#include <cmath>
#include <iomanip>
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

FunctionIntegralProvider make_szabo_provider() {
  const scfcpp::DenseMatrix S = make_2x2(1.0, 0.4508, 0.4508, 1.0);
  const scfcpp::DenseMatrix T = make_2x2(2.1643, 0.1670, 0.1670, 0.7600);
  const scfcpp::DenseMatrix V1 = make_2x2(-4.1398, -1.1029, -1.1029, -1.2652);
  const scfcpp::DenseMatrix V2 = make_2x2(-0.6772, -0.4113, -0.4113, -1.2266);
  const scfcpp::DenseMatrix Hcore = T + V1 + V2;

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
  fns.hcore = [Hcore](const IMolecule&) { return scfcpp::to_matrix(Hcore); };
  fns.eri_tensor = [eri](const IMolecule&) { return eri; };
  return FunctionIntegralProvider(std::move(fns));
}

int run_case(commutator::ConjugationMethod method, const char* label) {
  Molecule mol({{2, 0.0, 0.0, 0.0, ""}, {1, 0.0, 0.0, 1.4632, ""}}, +1, "sto-3g");
  mol.set_multiplicity(1);

  FunctionIntegralProvider provider = make_szabo_provider();
  XCConfig xc;
  xc.exchange = "HF";
  xc.correlation = "NONE";
  xc.exact_exchange_fraction = 1.0;

  SCFEngine eng(mol, provider, SCFEngine::Mode::RKS, xc, 1e-6, 20000, 0.0);
  SCFEngine::CommOptions opts;
  opts.alpha = 0.02;
  opts.exp_order = 6;
  opts.bch_order = 6;
  opts.do_purify = true;
  opts.purify_steps = 1;
  opts.method = method;

  eng.run_commutator_only(opts);
  if (!eng.converged()) {
    std::cerr << label << " did not converge\n";
    return 1;
  }

  const double R = 1.4632;
  const double expected_total = -2.860662;
  const double total_energy = eng.energy() + 2.0 / R;
  std::cout.setf(std::ios::fixed);
  std::cout << label << " converged=true total_energy=" << std::setprecision(12) << total_energy << "\n";
  if (std::abs(total_energy - expected_total) > 5e-3) {
    std::cerr << label << " energy mismatch: " << total_energy << " vs " << expected_total << "\n";
    return 1;
  }
  return 0;
}

}

int main() {
  int rc = 0;
  rc |= run_case(commutator::ConjugationMethod::Taylor, "Taylor");
  rc |= run_case(commutator::ConjugationMethod::BCH, "BCH");
  return rc;
}
