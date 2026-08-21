#pragma once

#include <vector>

#include "IIntegralProvider.h"
#include "basis.h"
#include "geometry.h"

// The HeH+ worked example from Szabo & Ostlund, "Modern Quantum Chemistry",
// section 3.5.2: STO-3G with zeta_He = 2.0925 and zeta_H = 1.24, at R = 1.4632
// bohr.  The integrals are the values printed in the book, so this provider
// needs no integral library at all and serves as the fixed reference that the
// libint2 provider is checked against.
//
// Note the basis is *not* the tabulated STO-3G: the library value for helium is
// zeta = 1.690.  GetBasisShells() therefore returns Szabo's zeta values, which
// is what makes a like-for-like comparison with libint2 possible.
namespace szabo_heh {

// Bond length in bohr; chosen so that Vnn = Z_He * Z_H / R = 1.3669 Hartree.
constexpr double kBondLengthBohr = 1.4632;
constexpr double kZetaHe = 2.0925;
constexpr double kZetaH = 1.24;
constexpr double kNuclearRepulsion = 1.3669;

// STO-3G contraction for a unit-zeta 1s Slater function (Hehre, Stewart &
// Pople, JCP 51, 2657 (1969)).  Exponents scale as zeta^2.
inline const std::vector<double>& UnitZetaExponents() {
    static const std::vector<double> v{2.22766, 0.405771, 0.109818};
    return v;
}
inline const std::vector<double>& UnitZetaCoefficients() {
    static const std::vector<double> v{0.154329, 0.535328, 0.444635};
    return v;
}

// The two 1s shells, in the same order as the basis functions of the tables
// below: index 0 is helium, index 1 is hydrogen.
BasisShells MakeBasis();
std::vector<Atom> MakeAtoms();

}  // namespace szabo_heh

// Note this provider does *not* implement IIntegralDerivativeProvider.  It used
// to, returning zeros -- which would silently produce a zero gradient rather
// than an error.  A provider that cannot differentiate should say so by not
// offering the interface.
class SzaboHeHIntegralProvider final : public IIntegralProvider {
public:
    int NumBasisFunctions() const override { return 2; }
    T2 ComputeHcore() const override;
    T2 ComputeOverlap() const override;
    T4 ComputeERI() const override;
    double ComputeERI(int i, int j, int k, int l) const override;
    // Computed from the geometry rather than returning the book's printed
    // 1.3669: that rounding was worth 3.3e-5 Hartree of spurious disagreement
    // with libint2.  See geometry.h for why this is a free function.
    double ComputeNuclearRepulsionEnergy() const override {
        return NuclearRepulsionEnergy(GetAtoms());
    }
    std::vector<Atom> GetAtoms() const override { return szabo_heh::MakeAtoms(); }
    BasisShells GetBasisShells() const override { return szabo_heh::MakeBasis(); }
};
