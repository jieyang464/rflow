#include "Libint2IntegralProvider.h"
#include "SzaboHeHIntegral.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

Molecule MakeHeHPlusMolecule() {
    Molecule molecule;
    const double bond_length_bohr = 2.0 / 1.3669;

    molecule.atoms.push_back(
        AtomWithBasis{Atom{"He", 2, 0.0, 0.0, 0.0}, "STO-3G"});
    molecule.atoms.push_back(
        AtomWithBasis{Atom{"H", 1, 0.0, 0.0, bond_length_bohr}, "STO-3G"});

    return molecule;
}

double MaxAbsDiff(const T2& a, const T2& b) {
    double max_diff = 0.0;
    for (Eigen::Index i = 0; i < a.dimension(0); ++i) {
        for (Eigen::Index j = 0; j < a.dimension(1); ++j) {
            max_diff = std::max(max_diff, std::abs(a(i, j) - b(i, j)));
        }
    }
    return max_diff;
}

double MaxAbsDiff(const T4& a, const T4& b) {
    double max_diff = 0.0;
    for (Eigen::Index i = 0; i < a.dimension(0); ++i) {
        for (Eigen::Index j = 0; j < a.dimension(1); ++j) {
            for (Eigen::Index k = 0; k < a.dimension(2); ++k) {
                for (Eigen::Index l = 0; l < a.dimension(3); ++l) {
                    max_diff = std::max(max_diff, std::abs(a(i, j, k, l) - b(i, j, k, l)));
                }
            }
        }
    }
    return max_diff;
}

}  // namespace

int main() {
    const SzaboHeHIntegralProvider reference;
    const Libint2IntegralProvider libint_provider(MakeHeHPlusMolecule());

    const T2 overlap_ref = reference.ComputeOverlap();
    const T2 overlap_libint = libint_provider.ComputeOverlap();
    const T2 hcore_ref = reference.ComputeHcore();
    const T2 hcore_libint = libint_provider.ComputeHcore();
    const T4 eri_ref = reference.ComputeERI();
    const T4 eri_libint = libint_provider.ComputeERI();

    const double vnn_diff = std::abs(reference.ComputeNuclearRepulsionEnergy() -
                                     libint_provider.ComputeNuclearRepulsionEnergy());
    const double overlap_diff = MaxAbsDiff(overlap_ref, overlap_libint);
    const double hcore_diff = MaxAbsDiff(hcore_ref, hcore_libint);
    const double eri_diff = MaxAbsDiff(eri_ref, eri_libint);

    constexpr double tol = 5.0e-4;
    const bool passed = vnn_diff < tol &&
                        overlap_diff < tol &&
                        hcore_diff < tol &&
                        eri_diff < tol;

    std::cout << "HeH+ Libint vs Szabo integral comparison\n";
    std::cout << "  max |Vnn - Vnn_ref| = " << vnn_diff << "\n";
    std::cout << "  max |S - S_ref|     = " << overlap_diff << "\n";
    std::cout << "  max |H - H_ref|     = " << hcore_diff << "\n";
    std::cout << "  max |ERI - ERI_ref| = " << eri_diff << "\n";

    if (!passed) {
        std::cerr << "FAIL: Libint integrals do not match the Szabo HeH+ table within "
                  << tol << ".\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
