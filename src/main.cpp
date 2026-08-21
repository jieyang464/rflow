// Demonstration driver: converge the same system twice, once by diagonalizing
// the Fock matrix and once by rotating the density with e^{-XS} D e^{SX}, and
// print the two answers side by side.

#include <cstdio>
#include <memory>
#include <string>

#include "IDensityUpdater.h"
#include "Libint2IntegralProvider.h"
#include "SzaboHeHIntegral.h"
#include "dft/dft_helper.h"
#include "fock_builders.h"
#include "geometry.h"
#include "jk_builder.h"
#include "scf.h"

namespace {

double MaxAbsDiff(const T2& a, const T2& b) {
    double worst = 0.0;
    for (Eigen::Index i = 0; i < a.dimension(0); ++i)
        for (Eigen::Index j = 0; j < a.dimension(1); ++j)
            worst = std::max(worst, std::fabs(a(i, j) - b(i, j)));
    return worst;
}

Molecule MakeWater() {
    Molecule m;
    m.atoms.push_back(AtomWithBasis{Atom{"O", 8, 0.0000, 0.0000, 0.0000}, "STO-3G"});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, 1.4310, 1.1070}, "STO-3G"});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, -1.4310, 1.1070}, "STO-3G"});
    return m;
}

void PrintMatrix(const char* label, const T2& m) {
    std::printf("  %s\n", label);
    for (Eigen::Index i = 0; i < m.dimension(0); ++i) {
        std::printf("    ");
        for (Eigen::Index j = 0; j < m.dimension(1); ++j) std::printf("% .6f  ", m(i, j));
        std::printf("\n");
    }
}

}  // namespace

int main() {
    // ---------------------------------------------------------------------
    std::printf("\n=== HeH+ / STO-3G, integrals from Szabo & Ostlund section 3.5.2 ===\n\n");
    {
        const SzaboHeHIntegralProvider provider;
        const T2 hcore = provider.ComputeHcore();
        const T2 overlap = provider.ComputeOverlap();

        SCFSettings settings;
        settings.Na = 1;
        settings.Nb = 1;
        settings.energy_tol = 1e-12;
        settings.gradient_tol = 1e-10;

        const auto converge = [&](const IDensityUpdater& updater, int max_iter) {
            SCFSettings s = settings;
            s.max_iter = max_iter;
            SCFResults results;
            results.nuclear_repulsion = NuclearRepulsionEnergy(provider.GetAtoms());
            GenerateInitialGuess(hcore, overlap, s.Na, s.Nb, results);
            auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
            UHFBuilder builder(std::move(jk), hcore);
            SCFLoop(s, builder, updater, overlap, results);
            return results;
        };

        const SCFResults diagonalized = converge(FockDiagonalizationDensityUpdater(), 200);
        const SCFResults rotated = converge(CommutatorDensityUpdater(8, 1.0), 2000);

        std::printf("  Fock diagonalization : E = %.12f Ha   (%d cycles)\n",
                    diagonalized.energy, diagonalized.iteration);
        std::printf("  Commutator rotation  : E = %.12f Ha   (%d cycles)\n", rotated.energy,
                    rotated.iteration);
        std::printf("  difference           : %.2e Ha,  max |dD| = %.2e\n\n",
                    rotated.energy - diagonalized.energy,
                    MaxAbsDiff(rotated.densityMatrices.Da, diagonalized.densityMatrices.Da));
        PrintMatrix("converged alpha density (commutator path):", rotated.densityMatrices.Da);
    }

    // ---------------------------------------------------------------------
    std::printf("\n=== H2O / STO-3G, integrals from libint2 ===\n\n");
    {
        const Molecule water = MakeWater();
        const Libint2IntegralProvider provider(water);
        const T2 hcore = provider.ComputeHcore();
        const T2 overlap = provider.ComputeOverlap();
        const T4 eri = provider.ComputeERI();
        const double vnn = NuclearRepulsionEnergy(provider.GetAtoms());

        const dft::DftGridContext grid = dft::build_dft_grid_context(
            provider.GetBasisShells(), provider.GetAtoms(), GridSettings::Medium());
        std::printf("  Becke grid: %zu points\n\n", grid.grid.size());

        const auto converge = [&](bool use_dft, const IDensityUpdater& updater, int max_iter) {
            SCFSettings settings;
            settings.Na = 5;
            settings.Nb = 5;
            settings.max_iter = max_iter;
            settings.energy_tol = 1e-12;
            settings.gradient_tol = 1e-8;

            SCFResults results;
            results.nuclear_repulsion = vnn;
            GenerateInitialGuess(hcore, overlap, 5, 5, results);
            auto jk = std::make_unique<InCoreJKBuilder>(eri);
            std::unique_ptr<IFockBuilder> builder;
            if (use_dft) {
                builder = std::make_unique<UKSBuilder>(std::move(jk), hcore, grid,
                                                       xc::LsdaSpec());
            } else {
                builder = std::make_unique<UHFBuilder>(std::move(jk), hcore);
            }
            SCFLoop(settings, *builder, updater, overlap, results);
            return results;
        };

        std::printf("  %-22s %-20s %8s  %10s\n", "method", "energy / Ha", "cycles", "max |dD|");
        std::printf("  %s\n", std::string(64, '-').c_str());

        for (const bool use_dft : {false, true}) {
            const char* label = use_dft ? "UKS / SVWN5" : "UHF";
            const SCFResults diagonalized =
                converge(use_dft, FockDiagonalizationDensityUpdater(use_dft ? 0.7 : 1.0), 300);
            const SCFResults rotated = converge(use_dft, CommutatorDensityUpdater(8, 0.1), 3000);

            std::printf("  %-22s %20.12f %8d\n", (std::string(label) + ", diagonalized").c_str(),
                        diagonalized.energy, diagonalized.iteration);
            std::printf("  %-22s %20.12f %8d  %10.2e\n",
                        (std::string(label) + ", commutator").c_str(), rotated.energy,
                        rotated.iteration,
                        MaxAbsDiff(rotated.densityMatrices.Da, diagonalized.densityMatrices.Da));
        }
        std::printf(
            "\n  The commutator path never diagonalizes a Fock matrix while iterating.\n"
            "  It needs many more cycles: the update is steepest descent in the\n"
            "  exponential parameterization, with no preconditioning.\n");
    }

    std::printf("\n");
    return 0;
}
