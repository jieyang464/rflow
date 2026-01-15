#pragma once

#include <functional>

namespace scfcpp {

template <typename Signature>
using Function = std::function<Signature>;

}
