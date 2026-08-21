#include "energy_gradient.h"

#include <memory>
#include <stdexcept>
#include <utility>

#include "IDensityUpdater.h"
#include "Libint2DerivativeProvider.h"
#include "Libint2IntegralProvider.h"
#include "effective_densities.h"
#include "fock_builders.h"
#include "geometry.h"
#include "gradient.h"
#include "jk_builder.h"
#include "scf.h"

namespace {

Molecule ToMolecule(const std::vector<Atom>& atoms, const std::string& basis) {
    Molecule molecule;
    molecule.atoms.reserve(atoms.size());
    for (const auto& atom : atoms) molecule.atoms.push_back(AtomWithBasis{atom, basis});
    return molecule;
}

}  // namespace

EnergyGradientFunction MakeHartreeFockEnergyGradient(HartreeFockJobSettings settings) {
    if (settings.Na < 0 || settings.Nb < 0) {
        throw std::invalid_argument("MakeHartreeFockEnergyGradient: negative electron count.");
    }

    // Carried between calls so consecutive geometries can warm-start.  Held by
    // shared_ptr because the returned std::function is copyable.
    struct Cache {
        bool has_guess{false};
        T2 Da, Db;
    };
    auto cache = std::make_shared<Cache>();

    return [settings, cache](const std::vector<Atom>& atoms) -> EnergyAndGradient {
        const Molecule molecule = ToMolecule(atoms, settings.basis);
        const Libint2IntegralProvider provider(molecule);

        const T2 hcore = provider.ComputeHcore();
        const T2 overlap = provider.ComputeOverlap();

        SCFSettings scf_settings;
        scf_settings.Na = settings.Na;
        scf_settings.Nb = settings.Nb;
        scf_settings.max_iter = settings.max_scf_iterations;
        scf_settings.energy_tol = settings.scf_energy_tol;
        scf_settings.gradient_tol = settings.scf_gradient_tol;

        SCFResults results;
        results.nuclear_repulsion = NuclearRepulsionEnergy(atoms);
        GenerateInitialGuess(hcore, overlap, settings.Na, settings.Nb, results);

        // Warm start, but only when the basis dimension still matches -- it will
        // not if the caller changed basis between calls.
        if (settings.reuse_density_guess && cache->has_guess &&
            cache->Da.dimension(0) == hcore.dimension(0)) {
            results.densityMatrices.Da = cache->Da;
            results.densityMatrices.Db = cache->Db;
        }

        auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
        UHFBuilder fock_builder(std::move(jk), hcore);
        const FockDiagonalizationDensityUpdater updater;
        SCFLoop(scf_settings, fock_builder, updater, overlap, results);

        cache->has_guess = results.converged;
        cache->Da = results.densityMatrices.Da;
        cache->Db = results.densityMatrices.Db;

        const Libint2DerivativeProvider derivatives(molecule);
        const HartreeFockEffectiveDensities densities(
            results.densityMatrices.Da, results.densityMatrices.Db,
            results.fockMatrices.Fa, results.fockMatrices.Fb);

        EnergyAndGradient out;
        out.energy = results.energy;
        out.gradient = AssembleGradient(derivatives, densities, atoms);
        out.converged = results.converged;
        out.scf_iterations = results.iteration;
        return out;
    };
}
