#include <iostream>

#include "SzaboHeHIntegral.h"
#include "scf.h"
#include "fock_builders.h"
#include "IDensityUpdater.h"
#include "ci.h"

int main() {
    SzaboHeHIntegralProvider szabo_integral_provider;
    const IIntegralProvider& integral_provider = szabo_integral_provider;
    const double nuclear_repulsion = integral_provider.ComputeNuclearRepulsionEnergy();
    
    // Abstract out the T2s
    T2 hcore = integral_provider.ComputeHcore();
    T2 overlap = integral_provider.ComputeOverlap();
    T4 eri = integral_provider.ComputeERI();

    SCFSettings settings;
    settings.Na = 1;
    settings.Nb = 1;

    SCFResults scfResults;
    scfResults.nuclear_repulsion = nuclear_repulsion;

    GenerateInitialGuess(hcore, overlap, settings.Na, settings.Nb, scfResults);

    // Create the JK builder (using InCore since we just fetched the T4)
    auto jk_builder = std::make_unique<InCoreJKBuilder>(std::move(eri));
    
    UHFBuilder fock_builder(std::move(jk_builder));
    FockDiagonalizationDensityUpdater updater;

    SCFLoop(settings, fock_builder, updater, overlap, hcore, scfResults);

    std::cout << "SCF converged energy = " << scfResults.energy << "\n";
    std::cout << "Nuclear repulsion energy (Vnn) = " << nuclear_repulsion << "\n";
    std::cout << "SCF orbitals (C_alpha):\n";
    for (int i = 0; i < scfResults.moCoefficients.Ca.dimension(0); ++i) {
        for (int j = 0; j < scfResults.moCoefficients.Ca.dimension(1); ++j) {
            std::cout << scfResults.moCoefficients.Ca(i, j) << " ";
        }
        std::cout << "\n";
    }

    // Since we moved eri, we will need it for TransformAO2MO.
    // So let's re-fetch it from provider for the CI part.
    T4 eri_ao = integral_provider.ComputeERI();
    T2 h_mo;
    T4 eri_mo;
    TransformAO2MO(scfResults.moCoefficients.Ca,
                   hcore,
                   eri_ao,
                   h_mo,
                   eri_mo);

    CIResult ci_result;
    ci_result.space.Build(2, 1, 1); // 2 MOs, 1 alpha, 1 beta
    ci_result.hcore_mo = h_mo;
    ci_result.eri_mo = eri_mo;

    BuildCIHamiltonian(ci_result);
    DiagonalizeCI(ci_result);

    int nstates = static_cast<int>(ci_result.space.dets.size());
    std::cout << "CI dimension = " << nstates << "\n";

    std::cout << "CI basis states (alpha,beta occupations):\n";
    for (int i = 0; i < nstates; ++i) {
        const auto& det = ci_result.space.dets[i];
        std::cout << "  state " << i << ": alpha=[";
        for (int p = 0; p < (int)det.alpha.size(); ++p) {
            std::cout << det.alpha[p] << (p+1<(int)det.alpha.size() ? "," : "");
        }
        std::cout << "] beta=[";
        for (int p = 0; p < (int)det.beta.size(); ++p) {
            std::cout << det.beta[p] << (p+1<(int)det.beta.size() ? "," : "");
        }
        std::cout << "]\n";
    }

    std::cout << "CI Hamiltonian matrix:\n";
    for (int i = 0; i < nstates; ++i) {
        for (int j = 0; j < nstates; ++j) {
            std::cout << ci_result.H(i, j) << (j+1<nstates ? "\t" : "");
        }
        std::cout << "\n";
    }

    std::cout << "CI eigenvalues (electronic states):\n";
    for (int i = 0; i < nstates; ++i) {
        std::cout << "  " << i << ": " << ci_result.eigenvalues(i, i) << "\n";
    }

    std::cout << "CI total energies (electronic + Vnn):\n";
    for (int i = 0; i < nstates; ++i) {
        std::cout << "  " << i << ": "
                  << ci_result.eigenvalues(i, i) + nuclear_repulsion << "\n";
    }

    std::cout << "CI eigenvectors / CI wavefunctions (columns in determinant basis):\n";
    for (int state = 0; state < nstates; ++state) {
        std::cout << "  state " << state << ":\n";
        for (int basis = 0; basis < nstates; ++basis) {
            const auto& det = ci_result.space.dets[basis];
            std::cout << "    coeff[" << basis << "] = "
                      << ci_result.cicoeff(basis, state)
                      << " for alpha=[";
            for (int p = 0; p < static_cast<int>(det.alpha.size()); ++p) {
                std::cout << det.alpha[p] << (p + 1 < static_cast<int>(det.alpha.size()) ? "," : "");
            }
            std::cout << "] beta=[";
            for (int p = 0; p < static_cast<int>(det.beta.size()); ++p) {
                std::cout << det.beta[p] << (p + 1 < static_cast<int>(det.beta.size()) ? "," : "");
            }
            std::cout << "]\n";
        }
    }

    return 0;
}
