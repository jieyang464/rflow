#include "SzaboHeHIntegral.h"

#include <stdexcept>

namespace szabo_heh {
namespace {

BasisShells BuildBasis() {
    auto scaled = [](double zeta) {
        std::vector<double> exponents = UnitZetaExponents();
        for (double& a : exponents) a *= zeta * zeta;
        return exponents;
    };

    BasisShells shells;
    shells.push_back(MakeShell(/*l=*/0, /*pure=*/false, Vec3{{0.0, 0.0, 0.0}},
                               scaled(kZetaHe), UnitZetaCoefficients()));
    shells.push_back(MakeShell(/*l=*/0, /*pure=*/false, Vec3{{0.0, 0.0, kBondLengthBohr}},
                               scaled(kZetaH), UnitZetaCoefficients()));
    return shells;
}

}  // namespace

BasisShells MakeBasis() { return BuildBasis(); }

std::vector<Atom> MakeAtoms() {
    return {
        Atom{"He", 2, 0.0, 0.0, 0.0},
        Atom{"H", 1, 0.0, 0.0, kBondLengthBohr},
    };
}

}  // namespace szabo_heh

T2 SzaboHeHIntegralProvider::ComputeHcore() const {
    T2 hcore(2, 2);
    hcore.setZero();
    hcore(0, 0) = -2.652744703;
    hcore(0, 1) = -1.347205024;
    hcore(1, 0) = -1.347205024;
    hcore(1, 1) = -1.731828436;
    return hcore;
}

T2 SzaboHeHIntegralProvider::ComputeOverlap() const {
    T2 overlap(2, 2);
    overlap.setZero();
    overlap(0, 0) = 1.0;
    overlap(0, 1) = 0.4508;
    overlap(1, 0) = 0.4508;
    overlap(1, 1) = 1.0;
    return overlap;
}

T4 SzaboHeHIntegralProvider::ComputeERI() const {
    T4 eri(2, 2, 2, 2);
    eri.setZero();
    eri(0, 0, 0, 0) = 1.307152;
    eri(0, 0, 0, 1) = 0.437279;
    eri(0, 0, 1, 0) = 0.437279;
    eri(0, 0, 1, 1) = 0.605703;
    eri(0, 1, 0, 0) = 0.437279;
    eri(0, 1, 0, 1) = 0.177267;
    eri(0, 1, 1, 0) = 0.177267;
    eri(0, 1, 1, 1) = 0.311795;
    eri(1, 0, 0, 0) = 0.437279;
    eri(1, 0, 0, 1) = 0.177267;
    eri(1, 0, 1, 0) = 0.177267;
    eri(1, 0, 1, 1) = 0.311795;
    eri(1, 1, 0, 0) = 0.605703;
    eri(1, 1, 0, 1) = 0.311795;
    eri(1, 1, 1, 0) = 0.311795;
    eri(1, 1, 1, 1) = 0.774608;
    return eri;
}

double SzaboHeHIntegralProvider::ComputeERI(int i, int j, int k, int l) const {
    if (i < 0 || j < 0 || k < 0 || l < 0 || i > 1 || j > 1 || k > 1 || l > 1) {
        throw std::out_of_range("SzaboHeHIntegralProvider: ERI index out of range.");
    }
    // Rebuilding a 2x2x2x2 tensor is cheap enough that caching would only add
    // state; the libint2 provider, where it matters, does cache.
    return ComputeERI()(i, j, k, l);
}

