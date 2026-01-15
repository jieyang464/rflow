#pragma once

#include "function.h"
#include "matrix.h"
#include "xc_config.h"

class IMolecule;
class IIntegralProvider;

struct VxcResult {
  scfcpp::Matrix Va;
  scfcpp::Matrix Vb;
  double Exc = 0.0;
  double tr_PVxc = 0.0;
};

struct VxcInputs {
  const IMolecule* mol = nullptr;
  const IIntegralProvider* integrals = nullptr;
  const scfcpp::Matrix* Pa = nullptr;
  const scfcpp::Matrix* Pb = nullptr;
  XCConfig xc{};
};

using VxcFunctor = scfcpp::Function<VxcResult(const VxcInputs&)>;

VxcFunctor make_stub_vxc_functor();
