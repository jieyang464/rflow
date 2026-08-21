// gradient.h — Assemble a nuclear gradient from integral derivatives and
// effective densities.
//
//   dE/dR_A = Tr[D_eff dHcore/dR_A]
//           - Tr[W      dS/dR_A]
//           + 1/2 sum_{munulamsig} Gamma_{munulamsig} d(mu nu|lam sig)/dR_A
//           + dVnn/dR_A
//
// Nothing here is method-specific: Hartree-Fock, MP2 and anything else differ
// only in what IEffectiveDensities supplies.
#pragma once

#include <vector>

#include "IIntegralDerivativeProvider.h"
#include "basis.h"
#include "effective_densities.h"
#include "types.h"

// One Vec3 per atom, in Hartree/bohr.
std::vector<Vec3> AssembleGradient(const IIntegralDerivativeProvider& derivatives,
                                   const IEffectiveDensities& densities,
                                   const std::vector<Atom>& atoms);

// The individual terms, for debugging a gradient that is close but not right --
// a sign error on the Pulay term or a missing degeneracy factor looks like a
// small uniform discrepancy in the total and is much easier to see split out.
struct GradientTerms {
    std::vector<Vec3> hcore;
    std::vector<Vec3> overlap;   // already carries its minus sign
    std::vector<Vec3> two_electron;
    std::vector<Vec3> nuclear;
    std::vector<Vec3> total;
};

GradientTerms AssembleGradientTerms(const IIntegralDerivativeProvider& derivatives,
                                    const IEffectiveDensities& densities,
                                    const std::vector<Atom>& atoms);
