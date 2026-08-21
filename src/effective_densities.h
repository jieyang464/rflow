// effective_densities.h — The density-like objects a nuclear gradient contracts
// the integral derivatives against.
//
//   dE/dx = Tr[D_eff dHcore/dx] - Tr[W dS/dx]
//         + 1/2 sum Gamma d(mu nu|lam sig)/dx + dVnn/dx
//
// This interface is the *only* method-specific part of the gradient: swapping
// Hartree-Fock for MP2 means supplying different D_eff, W and Gamma, while the
// derivative provider and the assembler are unchanged.
#pragma once

#include "IIntegralDerivativeProvider.h"
#include "types.h"

class IEffectiveDensities {
public:
    virtual ~IEffectiveDensities() = default;

    // One-particle effective density; to be contracts with dHcore/dx.
    virtual const T2& OneParticle() const = 0; //D_eff

    // Energy-weighted density; to be contracts with -dS/dx. (the Pulay term).
    virtual const T2& EnergyWeighted() const = 0; //W

    // Two-particle effective density, one shell-quartet block at a time.  This
    // is the N^4 object, so it is produced on demand rather than stored.
    virtual void TwoParticleBlock(const ShellQuartet& quartet, double* out) const = 0; 
};



// Everything here is native to the AO basis -- no MO transform is involved.

class HartreeFockEffectiveDensities final : public IEffectiveDensities {
public:
    // Da, Db, Fa, Fb from a converged SCF.  Copied, since they are only N^2 and
    // the gradient outlives no particular SCF iteration.
    HartreeFockEffectiveDensities(const T2& Da, const T2& Db, const T2& Fa, const T2& Fb);

    const T2& OneParticle() const override { return one_particle_; }
    const T2& EnergyWeighted() const override { return energy_weighted_; }
    void TwoParticleBlock(const ShellQuartet& quartet, double* out) const override;

private:
    T2 Da_, Db_;
    T2 one_particle_;     // Da + Db
    T2 energy_weighted_;  // Da Fa Da + Db Fb Db
};
