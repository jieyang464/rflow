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

    int NumBasisFunctions() const override { return 2; }

    // Fixed HeH+ geometry used by the Szabo table: He at the origin, H along z,
    // bond length 2.0/1.3669 bohr (the value test_libint and the other tests use).
    std::vector<Atom> GetAtoms() const override;
    BasisShells GetBasisShells() const override;
};
