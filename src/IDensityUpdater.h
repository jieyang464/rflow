#pragma once

#include "types.h"

class IDensityUpdater {
public:
    virtual ~IDensityUpdater() = default;
    virtual void UpdateDensity(const T2& F, const T2& S, T2& D) const = 0;
};

class CommutatorDensityUpdater : public IDensityUpdater {
public:
    CommutatorDensityUpdater(int bch_order = 4, double step = 1.0);

    void UpdateDensity(const T2& F, const T2& S, T2& D) const override;

private:
    int bch_order_{4};
    double step_{1.0};
};

// The textbook update: diagonalize F in the S metric and rebuild D from the
// lowest occupied eigenvectors.  This is the baseline the commutator updater is
// measured against, and the thing it exists to avoid -- the diagonalization is
// O(N^3) and cannot be made to scale linearly.
//
// `damping` mixes the new density with the old:
//
//     D <- damping * D_new + (1 - damping) * D_old
//
// 1.0 is pure Roothaan iteration, which is what the method is usually written
// as.  It converges for Hartree-Fock on small closed-shell systems but is not
// unconditionally stable: with a local XC potential it can settle into a
// two-cycle limit cycle instead, where the energy stalls while the density
// alternates between two configurations.  Around 0.5 fixes that.  Damping slows
// convergence when it was not needed, so the default stays undamped.
class FockDiagonalizationDensityUpdater : public IDensityUpdater {
public:
    explicit FockDiagonalizationDensityUpdater(double damping = 1.0);

    void UpdateDensity(const T2& F, const T2& S, T2& D) const override;

private:
    double damping_{1.0};
};
