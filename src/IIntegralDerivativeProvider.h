// IIntegralDerivativeProvider.h — First derivatives of the integrals with
// respect to nuclear coordinates.
//
// The asymmetry in this interface is deliberate and follows from size:
//
//   dS/dx, dHcore/dx   are N^2 per coordinate.  3*natoms of them is a few MB,
//                      so they are returned whole.
//   dERI/dx            is N^4 per coordinate -- 3*natoms times the size of the
//                      ERI tensor itself.  It is never materialized.
//
// The two-electron term is the only place this matters, and there both N^4
// objects in the contraction
//
//     1/2 sum_{munulamsig} Gamma_{munulamsig} d(mu nu|lam sig)/dx
//
// are needed *together and only together*.  So the provider drives the shell
// quartets and hands each derivative block to a sink, which builds that
// quartet's Gamma block from the (N^2) density and contracts on the spot.
// Neither N^4 object is ever stored.  This is also the natural place for Schwarz
// screening to live, since screening is a decision about a shell quartet.
#pragma once

#include <array>
#include <functional>
#include <vector>

#include "types.h"

// One shell quartet, with everything a sink needs to place its contribution.
struct ShellQuartet {
    std::array<int, 4> shell{{0, 0, 0, 0}};  // shell indices within the basis
    std::array<int, 4> bf{{0, 0, 0, 0}};     // first basis function of each shell
    std::array<int, 4> n{{0, 0, 0, 0}};      // number of functions in each shell
    std::array<int, 4> atom{{0, 0, 0, 0}};   // atom each shell is centred on

    // Permutational multiplicity of this canonical quartet: 1, 2, 4 or 8.
    // Quartets are enumerated canonically (s0 >= s1, s2 >= s3, (s0 s1) >= (s2 s3))
    // and the sink multiplies by this instead of visiting every permutation.
    // Correct only because Gamma carries the full 8-fold symmetry of (mu nu|lam sig).
    int degeneracy{1};

    int size() const { return n[0] * n[1] * n[2] * n[3]; }
};

// The twelve derivative buffers libint2 produces for one quartet: three
// Cartesian directions for each of the four centres.  A buffer may be null when
// that derivative is identically zero.  Each non-null buffer holds
// quartet.size() doubles in (f0, f1, f2, f3) row-major order -- the same layout
// as IEffectiveDensities::TwoParticleBlock writes.
struct EriDerivativeBlock {
    std::array<const double*, 12> buffer{};

    static constexpr int Index(int centre, int axis) { return 3 * centre + axis; }
    const double* At(int centre, int axis) const { return buffer[Index(centre, axis)]; }
};

class IIntegralDerivativeProvider {
public:
    virtual ~IIntegralDerivativeProvider() = default;

    virtual int NumBasisFunctions() const = 0;
    virtual int NumAtoms() const = 0;

    // Indexed [3 * atom + axis], axis 0/1/2 = x/y/z.
    virtual std::vector<T2> ComputeOverlapDerivatives() const = 0;
    virtual std::vector<T2> ComputeHcoreDerivatives() const = 0;

    // Stream the ERI derivatives one shell quartet at a time.
    using EriDerivativeSink =
        std::function<void(const ShellQuartet&, const EriDerivativeBlock&)>;
    virtual void ForEachEriDerivative(const EriDerivativeSink& sink) const = 0;
};
