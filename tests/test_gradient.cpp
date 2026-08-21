// Hartree-Fock nuclear gradient, checked against finite differences of the
// converged total energy.
//
// This is the only test here that means much.  Analytic gradients fail in ways
// that look entirely plausible -- a sign on the Pulay term, a missing
// permutational degeneracy factor, a shell-centre mapped to the wrong atom --
// and every one of those produces a gradient that is smooth, symmetric, and
// wrong.  Comparing against [E(x+h) - E(x-h)] / 2h catches all of them.

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "IDensityUpdater.h"
#include "Libint2DerivativeProvider.h"
#include "Libint2IntegralProvider.h"
#include "effective_densities.h"
#include "fock_builders.h"
#include "geometry.h"
#include "gradient.h"
#include "jk_builder.h"
#include "scf.h"
#include "test_util.h"

namespace {

// Deliberately asymmetric: a symmetric geometry makes whole components of the
// gradient vanish, which would hide errors rather than expose them.
Molecule MakeAsymmetricWater(const std::string& basis) {
    Molecule m;
    m.atoms.push_back(AtomWithBasis{Atom{"O", 8, 0.10, 0.20, 0.05}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.85, 1.48, 1.09}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, -1.32, -0.39, 0.98}, basis});
    return m;
}

Molecule Displace(const Molecule& m, std::size_t atom, int axis, double delta) {
    Molecule moved = m;
    switch (axis) {
        case 0: moved.atoms[atom].atom.x += delta; break;
        case 1: moved.atoms[atom].atom.y += delta; break;
        default: moved.atoms[atom].atom.z += delta; break;
    }
    return moved;
}

struct ConvergedSCF {
    SCFResults results;
    T2 hcore, overlap;
};

ConvergedSCF RunSCF(const Molecule& molecule, int na, int nb) {
    const Libint2IntegralProvider provider(molecule);

    ConvergedSCF out;
    out.hcore = provider.ComputeHcore();
    out.overlap = provider.ComputeOverlap();

    SCFSettings settings;
    settings.Na = na;
    settings.Nb = nb;
    settings.max_iter = 200;
    settings.energy_tol = 1e-12;
    settings.gradient_tol = 1e-10;

    out.results.nuclear_repulsion = NuclearRepulsionEnergy(provider.GetAtoms());
    GenerateInitialGuess(out.hcore, out.overlap, na, nb, out.results);

    auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
    UHFBuilder fock_builder(std::move(jk), out.hcore);
    const FockDiagonalizationDensityUpdater updater;
    SCFLoop(settings, fock_builder, updater, out.overlap, out.results);
    return out;
}

double TotalEnergy(const Molecule& molecule, int na, int nb) {
    return RunSCF(molecule, na, nb).results.energy;
}

}  // namespace

int main() {
    testing::TestRunner t("Hartree-Fock nuclear gradient");

    // ---------------------------------------------------------------------
    t.Section("nuclear repulsion gradient vs finite differences");
    {
        const std::vector<Atom> atoms = {
            Atom{"O", 8, 0.10, 0.20, 0.05},
            Atom{"H", 1, 0.85, 1.48, 1.09},
            Atom{"H", 1, -1.32, -0.39, 0.98},
        };
        const std::vector<Vec3> analytic = NuclearRepulsionGradient(atoms);
        const double h = 1e-5;
        double worst = 0.0;
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            for (int axis = 0; axis < 3; ++axis) {
                const double plus = NuclearRepulsionEnergy(DisplaceAtom(atoms, a, axis, +h));
                const double minus = NuclearRepulsionEnergy(DisplaceAtom(atoms, a, axis, -h));
                worst = std::max(worst, std::fabs(analytic[a][axis] - (plus - minus) / (2 * h)));
            }
        }
        t.Near(worst, 0.0, 1e-8, "max |dVnn/dx - finite difference|");

        // Translating the whole molecule cannot change the energy, so the
        // gradient must sum to zero across atoms.
        double worst_sum = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
            double sum = 0.0;
            for (const auto& g : analytic) sum += g[axis];
            worst_sum = std::max(worst_sum, std::fabs(sum));
        }
        t.Near(worst_sum, 0.0, 1e-12, "dVnn/dx satisfies translational invariance");
    }

    // ---------------------------------------------------------------------
    t.Section("libint2 build capability");
    t.Report("analytic one-body derivatives available (1 = yes)",
             Libint2HasAnalyticOneBodyDerivatives() ? 1.0 : 0.0);
    std::printf(
        "   note  with INCLUDE_ONEBODY=0 the one-electron derivatives come from\n"
        "         central differences of the integrals; the ERI derivatives are\n"
        "         analytic either way.\n");

    // ---------------------------------------------------------------------
    for (const char* basis : {"STO-3G", "6-31G"}) {
        t.Section(std::string("UHF gradient for asymmetric H2O/") + basis);

        const Molecule water = MakeAsymmetricWater(basis);
        const int na = 5, nb = 5;  // 10 electrons

        const ConvergedSCF scf = RunSCF(water, na, nb);
        t.Check(scf.results.converged, std::string(basis) + ": reference SCF converged");
        t.Report(std::string(basis) + ": E(HF) [Hartree]", scf.results.energy);

        const Libint2DerivativeProvider derivatives(water);
        const HartreeFockEffectiveDensities densities(
            scf.results.densityMatrices.Da, scf.results.densityMatrices.Db,
            scf.results.fockMatrices.Fa, scf.results.fockMatrices.Fb);

        const std::vector<Atom> atoms = AtomsOf(water);
        const GradientTerms terms = AssembleGradientTerms(derivatives, densities, atoms);

        // Finite differences of the *total converged energy*: the whole SCF is
        // rerun at each displaced geometry.
        const double h = 5e-4;
        double worst = 0.0;
        std::printf("   %-6s %-4s %14s %14s %12s\n", "atom", "axis", "analytic", "finite-diff",
                    "difference");
        for (std::size_t a = 0; a < water.atoms.size(); ++a) {
            for (int axis = 0; axis < 3; ++axis) {
                const double plus = TotalEnergy(Displace(water, a, axis, +h), na, nb);
                const double minus = TotalEnergy(Displace(water, a, axis, -h), na, nb);
                const double fd = (plus - minus) / (2 * h);
                const double analytic = terms.total[a][axis];
                worst = std::max(worst, std::fabs(analytic - fd));
                std::printf("   %-6zu %-4c % 14.9f % 14.9f % 12.2e\n", a, "xyz"[axis], analytic,
                            fd, analytic - fd);
            }
        }
        t.Near(worst, 0.0, 2e-5,
               std::string(basis) + ": max |analytic - finite difference|");

        // Translational invariance again -- but now of the assembled gradient,
        // which is a much sharper test: it fails if any shell is attributed to
        // the wrong atom.
        double worst_sum = 0.0;
        for (int axis = 0; axis < 3; ++axis) {
            double sum = 0.0;
            for (const auto& g : terms.total) sum += g[axis];
            worst_sum = std::max(worst_sum, std::fabs(sum));
        }
        t.Near(worst_sum, 0.0, 1e-6,
               std::string(basis) + ": gradient sums to zero over atoms");
    }

    return t.Summary();
}
