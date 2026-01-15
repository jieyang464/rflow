#include "integral_provider.h"
#include <cmath>
#include <stdexcept>
#include <string>

namespace {

void require_fn(bool ok, const char* name) {
  if (!ok) {
    throw std::runtime_error(std::string("Integral function not set: ") + name);
  }
}

}

size_t FunctionIntegralProvider::nbf(const IMolecule& mol) const {
  require_fn(static_cast<bool>(fns_.nbf), "nbf");
  return fns_.nbf(mol);
}

scfcpp::Matrix FunctionIntegralProvider::overlap(const IMolecule& mol) const {
  require_fn(static_cast<bool>(fns_.overlap), "overlap");
  return fns_.overlap(mol);
}

scfcpp::Matrix FunctionIntegralProvider::kinetic(const IMolecule& mol) const {
  require_fn(static_cast<bool>(fns_.kinetic), "kinetic");
  return fns_.kinetic(mol);
}

scfcpp::Matrix FunctionIntegralProvider::nuclear_attraction(const IMolecule& mol) const {
  require_fn(static_cast<bool>(fns_.nuclear_attraction), "nuclear_attraction");
  return fns_.nuclear_attraction(mol);
}

scfcpp::Matrix FunctionIntegralProvider::hcore(const IMolecule& mol) const {
  if (fns_.hcore) {
    return fns_.hcore(mol);
  }
  if (fns_.kinetic && fns_.nuclear_attraction) {
    scfcpp::DenseMatrix T = scfcpp::to_dense(fns_.kinetic(mol));
    scfcpp::DenseMatrix V = scfcpp::to_dense(fns_.nuclear_attraction(mol));
    if (T.rows() != V.rows() || T.cols() != V.cols()) {
      throw std::runtime_error("hcore fallback failed: kinetic/nuclear dimension mismatch");
    }
    return scfcpp::to_matrix(T + V);
  }
  throw std::runtime_error("Integral function not set: hcore");
}

double FunctionIntegralProvider::eri_abs_sum(const IMolecule& mol) const {
  if (fns_.eri_abs_sum) {
    return fns_.eri_abs_sum(mol);
  }
  if (fns_.eri_tensor) {
    std::vector<double> eri = fns_.eri_tensor(mol);
    double sum = 0.0;
    for (double v : eri) sum += std::abs(v);
    return sum;
  }
  throw std::runtime_error("Integral function not set: eri_abs_sum");
}

std::vector<double> FunctionIntegralProvider::eri_tensor(const IMolecule& mol) const {
  require_fn(static_cast<bool>(fns_.eri_tensor), "eri_tensor");
  return fns_.eri_tensor(mol);
}
