#include "xc/vxc_evaluator.h"

VxcFunctor make_stub_vxc_functor() {
  return [](const UksDensityInput&, UksVxcOutput&) {
    // Contract: the caller zeroes the output buffers before invoking the plugin.
  };
}
