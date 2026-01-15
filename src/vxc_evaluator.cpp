#include "vxc_evaluator.h"
#include "integral_provider.h"

VxcFunctor make_stub_vxc_functor() {
  return [](const VxcInputs& in) {
    const size_t n = (in.mol && in.integrals) ? in.integrals->nbf(*in.mol) : 0;
    VxcResult res;
    res.Va = scfcpp::Matrix(n, n);
    res.Va.setZero();
    res.Vb = scfcpp::Matrix(n, n);
    res.Vb.setZero();
    res.Exc = 0.0;
    res.tr_PVxc = 0.0;
    return res;
  };
}
