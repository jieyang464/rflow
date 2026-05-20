#pragma once

#include "IIntegralDerivativeProvider.h"
#include "IIntegralProvider.h"

class SzaboHeHIntegralProvider final : public IIntegralProvider,
                                       public IIntegralDerivativeProvider {
public:
    T2 ComputeHcore() const override;
    T2 ComputeOverlap() const override;
    T4 ComputeERI() const override;
    double ComputeERI(int i, int j, int k, int l) const override;
    IntegralDerivatives ComputeFirstDerivatives() const override;
    double ComputeNuclearRepulsionEnergy() const override {
        return 1.3669;  // from Szabo's book
    }
};
