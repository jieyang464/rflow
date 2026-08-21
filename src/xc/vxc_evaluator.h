// vxc_evaluator.h — Plugin contract for the exchange-correlation potential.
//
// The SCF knows nothing about grids or functionals: it hands over the two spin
// densities plus a functional selection and gets back Vxc_a, Vxc_b, E_xc and
// Tr(P Vxc).  The last quantity is what the SCF needs to correct the
// 0.5*Tr[D(H+F)] energy expression into the Kohn-Sham energy.
#pragma once

#include <functional>

#include "types.h"
#include "xc/xc_config.h"

struct UksDensityInput {
  const T2& Da;
  const T2& Db;
  const xc::FunctionalSpec& functional;
};

struct UksVxcOutput {
  T2& Vxca;
  T2& Vxcb;
  double& Exc;
  double& tr_PVxc;  // Tr(Da Vxc_a) + Tr(Db Vxc_b)

  UksVxcOutput(T2& Va, T2& Vb, double& Exc_, double& tr_)
      : Vxca(Va), Vxcb(Vb), Exc(Exc_), tr_PVxc(tr_) {}
};

using UksVxcPlugin = std::function<void(const UksDensityInput&, UksVxcOutput&)>;
using VxcFunctor = UksVxcPlugin;

// A functor that contributes nothing; used for pure Hartree-Fock.
VxcFunctor make_stub_vxc_functor();
