#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "function.h"
#include "matrix.h"
#include "molecule.h"

class IIntegralProvider {
public:
  virtual ~IIntegralProvider() = default;
  virtual size_t nbf(const IMolecule& mol) const = 0;
  virtual scfcpp::Matrix overlap(const IMolecule& mol) const = 0;
  virtual scfcpp::Matrix kinetic(const IMolecule& mol) const = 0;
  virtual scfcpp::Matrix nuclear_attraction(const IMolecule& mol) const = 0;
  virtual scfcpp::Matrix hcore(const IMolecule& mol) const = 0;
  virtual double eri_abs_sum(const IMolecule& mol) const = 0;
  virtual std::vector<double> eri_tensor(const IMolecule& mol) const = 0;
};

struct IntegralFunctions {
  scfcpp::Function<size_t(const IMolecule&)> nbf;
  scfcpp::Function<scfcpp::Matrix(const IMolecule&)> overlap;
  scfcpp::Function<scfcpp::Matrix(const IMolecule&)> kinetic;
  scfcpp::Function<scfcpp::Matrix(const IMolecule&)> nuclear_attraction;
  scfcpp::Function<scfcpp::Matrix(const IMolecule&)> hcore;
  scfcpp::Function<double(const IMolecule&)> eri_abs_sum;
  scfcpp::Function<std::vector<double>(const IMolecule&)> eri_tensor;
};

class FunctionIntegralProvider final : public IIntegralProvider {
public:
  explicit FunctionIntegralProvider(IntegralFunctions fns) : fns_(std::move(fns)) {}

  size_t nbf(const IMolecule& mol) const override;
  scfcpp::Matrix overlap(const IMolecule& mol) const override;
  scfcpp::Matrix kinetic(const IMolecule& mol) const override;
  scfcpp::Matrix nuclear_attraction(const IMolecule& mol) const override;
  scfcpp::Matrix hcore(const IMolecule& mol) const override;
  double eri_abs_sum(const IMolecule& mol) const override;
  std::vector<double> eri_tensor(const IMolecule& mol) const override;

private:
  IntegralFunctions fns_;
};

std::shared_ptr<IIntegralProvider> make_libint_integral_provider();
