#pragma once

#include <functional>

#include "types.h"

struct UksDensityInput {
  const T2& Da;
  const T2& Db;
  const T2& overlap;
  int xc_functional_id;
};

struct UksVxcOutput {
  T2& Vxca;
  T2& Vxcb;
  double& Exc;
  double& tr_PVxc;

  UksVxcOutput(T2& Va, T2& Vb, double& Exc_, double& tr_)
      : Vxca(Va), Vxcb(Vb), Exc(Exc_), tr_PVxc(tr_) {}
};

using UksVxcPlugin = std::function<void(const UksDensityInput&, UksVxcOutput&)>;

using VxcFunctor = UksVxcPlugin;

VxcFunctor make_stub_vxc_functor();
