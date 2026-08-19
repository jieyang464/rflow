#pragma once

#include <vector>

#include "types.h"

struct IntegralDerivatives {
    std::vector<T2> d_overlap;
    std::vector<T2> d_hcore;
    std::vector<T4> d_eri;
};

struct IIntegralDerivativeProvider {
    virtual ~IIntegralDerivativeProvider() = default;
    virtual IntegralDerivatives ComputeFirstDerivatives() const = 0;
};
