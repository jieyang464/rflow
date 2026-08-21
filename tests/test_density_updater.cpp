// The claim this project exists to demonstrate: the SCF can be converged by
// rotating the density matrix with e^{-XS} D e^{SX}, never diagonalizing a Fock
// matrix, and it reaches the *same density matrix* the diagonalizing SCF does.
//
// The tests are arranged so that only the updater varies.  Same provider, same
// Fock builder, same initial guess, same convergence thresholds -- swap
// FockDiagonalizationDensityUpdater for CommutatorDensityUpdater and compare.

#include <cmath>
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
#include "linalg.h"
#include "scf.h"
#include "test_util.h"

namespace {

Molecule MakeWater(const std::string& basis) {
    Molecule m;
    m.atoms.push_back(AtomWithBasis{Atom{"O", 8, 0.0000, 0.0000, 0.0000}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, 1.4310, 1.1070}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, -1.4310, 1.1070}, basis});
    return m;
}

double TraceProduct(const T2& A, const T2& B) {
    double trace = 0.0;
    for (Eigen::Index i = 0; i < A.dimension(0); ++i)
        for (Eigen::Index j = 0; j < A.dimension(1); ++j) trace += A(i, j) * B(j, i);
    return trace;
}

}  // namespace

int main() {
    testing::TestRunner t("density updater: commutator rotation vs Fock diagonalization");

    // =====================================================================
    t.Section("one step at the fixed point leaves the density alone");
    // At convergence FDS - SDF = 0, so X = 0 and e^{-XS} D e^{SX} = D.  If the
    // commutator update moved a converged density it could not share a fixed
    // point with the diagonalizing update, whatever else it did.
    {
        const SzaboHeHIntegralProvider szabo;
        const T2 hcore = szabo.ComputeHcore();
        const T2 overlap = szabo.ComputeOverlap();

        SCFSettings settings;
        settings.Na = 1;
        settings.Nb = 1;
        settings.energy_tol = 1e-13;
        settings.gradient_tol = 1e-11;

        SCFResults converged;
        converged.nuclear_repulsion = szabo.ComputeNuclearRepulsionEnergy();
        GenerateInitialGuess(hcore, overlap, 1, 1, converged);
        auto jk = std::make_unique<InCoreJKBuilder>(szabo.ComputeERI());
        UHFBuilder builder(std::move(jk), hcore);
        SCFLoop(settings, builder, FockDiagonalizationDensityUpdater(), overlap, converged);
        t.Check(converged.converged, "reference SCF converged");

        const T2 D_reference = converged.densityMatrices.Da;

        T2 D_commutator = D_reference;
        CommutatorDensityUpdater(8, 1.0).UpdateDensity(converged.fockMatrices.Fa, overlap,
                                                       D_commutator);
        t.Near(testing::MaxAbsDiff(D_commutator, D_reference), 0.0, 1e-9,
               "commutator update is a fixed point of the converged density");

        T2 D_diagonal = D_reference;
        FockDiagonalizationDensityUpdater().UpdateDensity(converged.fockMatrices.Fa, overlap,
                                                          D_diagonal);
        t.Near(testing::MaxAbsDiff(D_diagonal, D_reference), 0.0, 1e-9,
               "diagonalizing update is a fixed point of the same density");
    }

    // =====================================================================
    t.Section("the rotation preserves idempotency and electron count");
    // e^{-XS} D e^{SX} is a similarity transform in the S metric, so Tr(DS) and
    // D S D = D are preserved exactly; truncating the BCH series breaks both
    // slightly, which is what the McWeeny purification is there to repair.
    {
        const SzaboHeHIntegralProvider szabo;
        const T2 hcore = szabo.ComputeHcore();
        const T2 overlap = szabo.ComputeOverlap();

        SCFResults state;
        state.nuclear_repulsion = szabo.ComputeNuclearRepulsionEnergy();
        GenerateInitialGuess(hcore, overlap, 1, 1, state);

        auto jk = std::make_unique<InCoreJKBuilder>(szabo.ComputeERI());
        UHFBuilder builder(std::move(jk), hcore);
        FockBuildResult fock;
        builder.build({state.densityMatrices.Da, state.densityMatrices.Db}, fock);

        T2 D = state.densityMatrices.Da;
        const CommutatorDensityUpdater updater(8, 0.5);

        double worst_trace = 0.0;
        double worst_idempotency = 0.0;
        for (int step = 0; step < 20; ++step) {
            updater.UpdateDensity(fock.Fa, overlap, D);
            worst_trace = std::max(worst_trace, std::fabs(TraceProduct(D, overlap) - 1.0));
            const T2 DSD = MatMul(MatMul(D, overlap), D);
            worst_idempotency = std::max(worst_idempotency, MaxAbsElement(DSD - D));
        }
        t.Near(worst_trace, 0.0, 1e-10, "Tr(D S) stays at 1 electron over 20 rotations");
        t.Near(worst_idempotency, 0.0, 1e-10, "D S D = D holds over 20 rotations");
    }

    // =====================================================================
    // The main event: same everything, only the updater swapped.
    struct Case {
        const char* label;
        bool use_dft;
        int na, nb;
    };

    const Libint2IntegralProvider water_provider(MakeWater("STO-3G"));
    const T2 water_hcore = water_provider.ComputeHcore();
    const T2 water_overlap = water_provider.ComputeOverlap();
    const T4 water_eri = water_provider.ComputeERI();
    const dft::DftGridContext water_grid = dft::build_dft_grid_context(
        water_provider.GetBasisShells(), water_provider.GetAtoms(), GridSettings::Coarse());

    for (const Case& c : {Case{"H2O/STO-3G UHF", false, 5, 5},
                          Case{"H2O/STO-3G UKS(SVWN5)", true, 5, 5}}) {
        t.Section(std::string(c.label) + ": swapping the density updater");

        const auto run = [&](const IDensityUpdater& updater, int max_iter) {
            SCFSettings settings;
            settings.Na = c.na;
            settings.Nb = c.nb;
            settings.max_iter = max_iter;
            settings.energy_tol = 1e-12;
            settings.gradient_tol = 1e-8;

            SCFResults results;
            results.nuclear_repulsion = NuclearRepulsionEnergy(water_provider.GetAtoms());
            GenerateInitialGuess(water_hcore, water_overlap, c.na, c.nb, results);

            auto jk = std::make_unique<InCoreJKBuilder>(water_eri);
            std::unique_ptr<IFockBuilder> builder;
            if (c.use_dft) {
                builder = std::make_unique<UKSBuilder>(std::move(jk), water_hcore, water_grid,
                                                       xc::LsdaSpec());
            } else {
                builder = std::make_unique<UHFBuilder>(std::move(jk), water_hcore);
            }
            SCFLoop(settings, *builder, updater, water_overlap, results);
            return results;
        };

        const SCFResults diagonalized =
            run(FockDiagonalizationDensityUpdater(c.use_dft ? 0.7 : 1.0), 300);
        const SCFResults rotated = run(CommutatorDensityUpdater(8, 0.1), 3000);

        t.Check(diagonalized.converged, "diagonalizing SCF converged");
        t.Check(rotated.converged, "commutator SCF converged (no Fock diagonalization)");
        t.Report("E, diagonalized [Hartree]", diagonalized.energy);
        t.Report("E, commutator   [Hartree]", rotated.energy);
        t.Report("cycles, diagonalized", diagonalized.iteration);
        t.Report("cycles, commutator", rotated.iteration);

        t.Near(rotated.energy, diagonalized.energy, 1e-10, "same converged energy");
        t.Near(testing::MaxAbsDiff(rotated.densityMatrices.Da, diagonalized.densityMatrices.Da),
               0.0, 1e-6, "same converged alpha density matrix");
        t.Near(testing::MaxAbsDiff(rotated.densityMatrices.Db, diagonalized.densityMatrices.Db),
               0.0, 1e-6, "same converged beta density matrix");
        t.Near(TraceProduct(rotated.densityMatrices.Da, water_overlap), 5.0, 1e-8,
               "commutator path kept 5 alpha electrons");
    }

    // =====================================================================
    t.Section("known limitation: the fixed step is not unconditionally stable");
    // Recorded as a test so it cannot quietly change.  The update is steepest
    // descent in the exponential parameterization with no line search and no
    // trust region, so above some step length the truncated BCH series and the
    // single purification can no longer keep D on the idempotent manifold.
    {
        SCFSettings settings;
        settings.Na = 5;
        settings.Nb = 5;
        settings.max_iter = 1200;
        settings.energy_tol = 1e-12;
        settings.gradient_tol = 1e-8;

        const auto energy_for_step = [&](double step) {
            SCFResults results;
            results.nuclear_repulsion = NuclearRepulsionEnergy(water_provider.GetAtoms());
            GenerateInitialGuess(water_hcore, water_overlap, 5, 5, results);
            auto jk = std::make_unique<InCoreJKBuilder>(water_eri);
            UHFBuilder builder(std::move(jk), water_hcore);
            SCFLoop(settings, builder, CommutatorDensityUpdater(8, step), water_overlap, results);
            return results;
        };

        const SCFResults small = energy_for_step(0.1);
        const SCFResults large = energy_for_step(0.5);
        t.Check(small.converged, "step = 0.1 converges");
        t.Check(!large.converged || !std::isfinite(large.energy),
                "step = 0.5 does not converge (documented limitation, not a regression)");
        t.Report("step = 0.5 final energy", large.energy);
        std::printf(
            "   note  a line search and a trust region on the rotation would remove\n"
            "         this; the fixed-step form is what the proof of concept covers.\n");
    }

    return t.Summary();
}
