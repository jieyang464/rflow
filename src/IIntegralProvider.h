#pragma once

#include <vector>

#include "basis.h"
#include "types.h"

// Contract for anything that can hand the SCF the one- and two-electron
// integrals over a fixed geometry and basis.
//
// Providers also publish the *geometry* and the *basis* itself, because the DFT
// part needs the value of each basis function in real space to build the
// exchange-correlation quadrature.  Keeping that on this interface is what lets
// the grid work unchanged with either the hard-coded Szabo HeH+ table or libint2.
struct IIntegralProvider {
    virtual ~IIntegralProvider() = default;

    // Number of basis functions; must equal NumBasisFunctions(GetBasisShells()).
    virtual int NumBasisFunctions() const = 0;

    virtual T2 ComputeHcore() const = 0;
    virtual T2 ComputeOverlap() const = 0;

    // Computes and returns the full 4D ERI tensor in chemists' notation (ij|kl).
    virtual T4 ComputeERI() const = 0;

    // Single ERI element (i j | k l).  Implementations are expected to make
    // repeated calls cheap (e.g. by caching), since direct J/K builds call this
    // O(N^4) times per SCF iteration.
    virtual double ComputeERI(int i, int j, int k, int l) const = 0;

    virtual double ComputeNuclearRepulsionEnergy() const = 0;

    // Atomic centres, in bohr, in the order the basis functions are blocked.
    virtual std::vector<Atom> GetAtoms() const = 0;

    // Contracted Gaussian shells, in the same order as the basis functions.
    virtual BasisShells GetBasisShells() const = 0;
};
