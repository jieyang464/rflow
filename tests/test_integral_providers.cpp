// Cross-check the two integral providers against each other on HeH+.
//
// The Szabo provider is a table of numbers printed in the book; the libint2
// provider computes them.  Agreement between them validates the libint2 shell
// loops, the symmetric-block filling, and the basis conversion in both
// directions -- and it is the only test in the suite that has an *external*
// reference rather than checking the code against itself.
//
// One subtlety drives the whole test: Szabo's HeH+ example uses zeta_He = 2.0925,
// whereas the tabulated STO-3G helium basis is zeta = 1.690.  Asking libint2 for
// "STO-3G" therefore compares two different calculations.  The test does both:
// it demands agreement on Szabo's actual basis, and separately reports the
// library-STO-3G numbers so the difference is visible rather than mysterious.

#include <libint2.hpp>

#include <array>
#include <cstdio>
#include <memory>
#include <utility>

#include "IDensityUpdater.h"
#include "Libint2IntegralProvider.h"
#include "SzaboHeHIntegral.h"
#include "fock_builders.h"
#include "jk_builder.h"
#include "scf.h"
#include "test_util.h"

namespace {

Molecule MakeHeHPlusMolecule(const std::string& basis_name) {
    Molecule molecule;
    molecule.atoms.push_back(
        AtomWithBasis{Atom{"He", 2, 0.0, 0.0, 0.0}, basis_name});
    molecule.atoms.push_back(
        AtomWithBasis{Atom{"H", 1, 0.0, 0.0, szabo_heh::kBondLengthBohr}, basis_name});
    return molecule;
}

// Converge HeH+ with a fixed updater, so the only thing varying is the provider.
SCFResults RunHeHPlusSCF(const IIntegralProvider& provider, const T2& hcore, const T2& overlap) {
    SCFSettings settings;
    settings.Na = 1;
    settings.Nb = 1;
    settings.energy_tol = 1e-11;
    settings.gradient_tol = 1e-9;

    SCFResults results;
    results.nuclear_repulsion = provider.ComputeNuclearRepulsionEnergy();
    GenerateInitialGuess(hcore, overlap, settings.Na, settings.Nb, results);

    auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
    UHFBuilder fock_builder(std::move(jk), hcore);
    const FockDiagonalizationDensityUpdater updater;

    SCFLoop(settings, fock_builder, updater, overlap, results);
    return results;
}

// Reference one-electron matrices built with no permutational symmetry at all:
// every ordered shell pair (s1, s2) computed and assigned.  Slower by 2x, but it
// cannot get the symmetric-block bookkeeping wrong, which is precisely what the
// provider's s2 <= s1 loop has to get right.
struct ReferenceOneElectron {
    T2 overlap;
    T2 hcore;
};

ReferenceOneElectron BuildReferenceOneElectron(const std::vector<libint2::Atom>& atoms,
                                               const std::string& basis_name) {
    const libint2::BasisSet bs(basis_name, atoms, true);
    const auto n = static_cast<Eigen::Index>(bs.nbf());
    const auto shell2bf = bs.shell2bf();

    std::vector<std::pair<double, std::array<double, 3>>> charges;
    for (const auto& a : atoms) {
        charges.push_back({static_cast<double>(a.atomic_number), {a.x, a.y, a.z}});
    }

    libint2::Engine s_eng(libint2::Operator::overlap, bs.max_nprim(), bs.max_l(), 0);
    libint2::Engine t_eng(libint2::Operator::kinetic, bs.max_nprim(), bs.max_l(), 0);
    libint2::Engine v_eng(libint2::Operator::nuclear, bs.max_nprim(), bs.max_l(), 0);
    v_eng.set_params(charges);

    ReferenceOneElectron out;
    out.overlap = T2(n, n);
    out.hcore = T2(n, n);
    out.overlap.setZero();
    out.hcore.setZero();

    for (std::size_t s1 = 0; s1 < bs.size(); ++s1) {
        for (std::size_t s2 = 0; s2 < bs.size(); ++s2) {  // full loop, both orders
            const std::size_t n1 = bs[s1].size(), n2 = bs[s2].size();
            const std::size_t i0 = shell2bf[s1], j0 = shell2bf[s2];

            s_eng.compute(bs[s1], bs[s2]);
            t_eng.compute(bs[s1], bs[s2]);
            v_eng.compute(bs[s1], bs[s2]);

            for (std::size_t p = 0; p < n1; ++p) {
                for (std::size_t q = 0; q < n2; ++q) {
                    const auto i = static_cast<Eigen::Index>(i0 + p);
                    const auto j = static_cast<Eigen::Index>(j0 + q);
                    if (s_eng.results()[0]) out.overlap(i, j) = s_eng.results()[0][p * n2 + q];
                    double h = 0.0;
                    if (t_eng.results()[0]) h += t_eng.results()[0][p * n2 + q];
                    if (v_eng.results()[0]) h += v_eng.results()[0][p * n2 + q];
                    out.hcore(i, j) = h;
                }
            }
        }
    }
    return out;
}

// A geometry with no symmetry whatsoever.  This matters more than it looks: with
// a symmetric molecule the off-diagonal elements of a p shell's own block vanish,
// so a fill that double-counts them still gives the right answer.  Symmetric test
// cases are exactly how that class of bug survives.
Molecule MakeAsymmetricWater(const std::string& basis_name) {
    Molecule m;
    m.atoms.push_back(AtomWithBasis{Atom{"O", 8, 0.11, 0.23, 0.07}, basis_name});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.90, 1.51, 1.13}, basis_name});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, -1.37, -0.42, 1.02}, basis_name});
    return m;
}

std::vector<libint2::Atom> ToLibintAtoms(const Molecule& m) {
    std::vector<libint2::Atom> out;
    for (const auto& a : m.atoms) {
        out.push_back(libint2::Atom{a.atom.atomic_number, a.atom.x, a.atom.y, a.atom.z});
    }
    return out;
}

}  // namespace

int main() {
    testing::TestRunner t("integral providers: Szabo HeH+ table vs libint2");

    const SzaboHeHIntegralProvider szabo;

    // ---------------------------------------------------------------------
    t.Section("libint2 on Szabo's own basis (zeta_He = 2.0925, zeta_H = 1.24)");

    const Libint2IntegralProvider libint_exact(MakeHeHPlusMolecule("STO-3G"),
                                               szabo.GetBasisShells());

    t.Check(libint_exact.NumBasisFunctions() == szabo.NumBasisFunctions(),
            "basis dimensions agree (nbf = 2)");

    // Szabo prints the overlap to 4 decimals and the rest to 6, so the tolerance
    // is set by the book's precision, not by the integrals.
    t.Near(libint_exact.ComputeNuclearRepulsionEnergy(),
           szabo.ComputeNuclearRepulsionEnergy(), 1e-4, "nuclear repulsion Vnn");

    const T2 s_ref = szabo.ComputeOverlap();
    const T2 s_libint = libint_exact.ComputeOverlap();
    t.Near(testing::MaxAbsDiff(s_ref, s_libint), 0.0, 1e-4, "max |S - S_szabo|");

    const T2 h_ref = szabo.ComputeHcore();
    const T2 h_libint = libint_exact.ComputeHcore();
    t.Near(testing::MaxAbsDiff(h_ref, h_libint), 0.0, 1e-5, "max |Hcore - Hcore_szabo|");

    const T4 eri_ref = szabo.ComputeERI();
    const T4 eri_libint = libint_exact.ComputeERI();
    t.Near(testing::MaxAbsDiff(eri_ref, eri_libint), 0.0, 1e-5, "max |ERI - ERI_szabo|");

    // The element accessor must agree with the full tensor everywhere: the
    // libint2 side reaches it through a lazily built cache, which is exactly the
    // sort of thing that can silently go stale.
    double worst_element = 0.0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k)
                for (int l = 0; l < 2; ++l)
                    worst_element = std::max(
                        worst_element,
                        std::fabs(libint_exact.ComputeERI(i, j, k, l) - eri_libint(i, j, k, l)));
    t.Near(worst_element, 0.0, 0.0, "ComputeERI(i,j,k,l) matches ComputeERI() exactly");

    // ---------------------------------------------------------------------
    t.Section("the two providers converge to the same SCF solution");

    const SCFResults szabo_scf = RunHeHPlusSCF(szabo, h_ref, s_ref);
    const SCFResults libint_scf = RunHeHPlusSCF(libint_exact, h_libint, s_libint);

    t.Check(szabo_scf.converged, "Szabo provider SCF converged");
    t.Check(libint_scf.converged, "libint2 provider SCF converged");
    t.Report("E(HeH+) from Szabo table   [Hartree]", szabo_scf.energy);
    t.Report("E(HeH+) from libint2       [Hartree]", libint_scf.energy);

    // Szabo & Ostlund quote -2.8606 Hartree for this system.
    t.Near(szabo_scf.energy, -2.8606, 5e-4, "Szabo total energy vs the book's -2.8606");
    t.Near(libint_scf.energy, szabo_scf.energy, 2e-4,
           "libint2 total energy matches the Szabo table");
    t.Near(testing::MaxAbsDiff(libint_scf.densityMatrices.Da, szabo_scf.densityMatrices.Da),
           0.0, 2e-4, "converged density matrices agree");

    // ---------------------------------------------------------------------
    t.Section("library STO-3G is a different basis (reported, not asserted)");

    const Libint2IntegralProvider libint_library(MakeHeHPlusMolecule("STO-3G"));
    const T2 h_library = libint_library.ComputeHcore();
    const T2 s_library = libint_library.ComputeOverlap();
    const SCFResults library_scf = RunHeHPlusSCF(libint_library, h_library, s_library);

    t.Report("max |S_library - S_szabo|", testing::MaxAbsDiff(s_library, s_ref));
    t.Report("E(HeH+) with library STO-3G [Hartree]", library_scf.energy);
    t.Report("  ... difference from Szabo's basis", library_scf.energy - szabo_scf.energy);
    std::printf(
        "   note  the library helium 1s is zeta = 1.690, not Szabo's 2.0925, so this\n"
        "         is a different calculation -- it is shown to make that explicit.\n");

    // ---------------------------------------------------------------------
    t.Section("shell-pair symmetry shortcuts agree with a full unsymmetrized loop");

    // p and d shells have size > 1, so their diagonal shell block is a square
    // that libint2 already returns in full -- mirroring it again would double the
    // off-diagonal elements.  Only a low-symmetry molecule exposes that.
    for (const char* basis : {"6-31G", "cc-pVDZ"}) {
        const Molecule water = MakeAsymmetricWater(basis);
        const Libint2IntegralProvider provider(water);
        const ReferenceOneElectron reference =
            BuildReferenceOneElectron(ToLibintAtoms(water), basis);

        t.Check(provider.NumBasisFunctions() == reference.overlap.dimension(0),
                std::string(basis) + ": nbf matches the reference build");
        t.Near(testing::MaxAbsDiff(provider.ComputeOverlap(), reference.overlap), 0.0, 1e-12,
               std::string(basis) + ": S matches unsymmetrized reference");
        t.Near(testing::MaxAbsDiff(provider.ComputeHcore(), reference.hcore), 0.0, 1e-12,
               std::string(basis) + ": Hcore matches unsymmetrized reference");
    }

    // ---------------------------------------------------------------------
    t.Section("ERI permutational symmetry");

    {
        const Molecule water = MakeAsymmetricWater("6-31G");
        const Libint2IntegralProvider provider(water);
        const T4 eri = provider.ComputeERI();
        const int n = provider.NumBasisFunctions();

        double worst_symmetry = 0.0;
        double worst_diagonal = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                for (int k = 0; k < n; ++k)
                    for (int l = 0; l < n; ++l) {
                        const double v = eri(i, j, k, l);
                        // All eight permutations of (ij|kl) for real orbitals.
                        for (const double other : {eri(j, i, k, l), eri(i, j, l, k),
                                                   eri(j, i, l, k), eri(k, l, i, j),
                                                   eri(l, k, i, j), eri(k, l, j, i),
                                                   eri(l, k, j, i)}) {
                            worst_symmetry = std::max(worst_symmetry, std::fabs(v - other));
                        }
                        if (i == k && j == l) worst_diagonal = std::min(worst_diagonal, v);
                    }
        t.Near(worst_symmetry, 0.0, 1e-14, "6-31G: (ij|kl) 8-fold permutational symmetry");
        // (ij|ij) is the self-repulsion of a charge distribution: never negative.
        t.Check(worst_diagonal >= 0.0, "6-31G: all (ij|ij) are non-negative");
    }

    return t.Summary();
}
