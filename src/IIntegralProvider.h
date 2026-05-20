#pragma once

#include <functional>
#include <vector>
#include <stdexcept>
#include "types.h"

struct IIntegralProvider {
    virtual ~IIntegralProvider() = default;

    virtual T2 ComputeHcore() const = 0;
    virtual T2 ComputeOverlap() const = 0;
    
    // Computes and returns the full 4D ER tensor
    virtual T4 ComputeERI() const = 0;

    // Feed ERI tensor element (i, j | k, l) on the fly
    virtual double ComputeERI(int i, int j, int k, int l) const = 0;

    virtual double ComputeNuclearRepulsionEnergy() const = 0;
};


 