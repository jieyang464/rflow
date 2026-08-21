// energy_gradient.h — What a geometry optimizer needs from an electronic
// structure method: one scalar and 3N derivatives, for a given set of nuclei.
//
// This is deliberately a callable rather than an interface.  There is exactly
// one operation here, and everything that genuinely varies between methods --
// the Fock build, the effective densities -- already varies below this line
// through IFockBuilder and IEffectiveDensities.  By the time the optimizer sees
// a geometry, all of that has collapsed into a number and a gradient.
//
// A closure also gets one thing an interface would need extra methods for: it
// can carry the previous converged density forward as the initial guess for the
// next geometry, which is where most of the saving in an optimization comes from.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "basis.h"
#include "types.h"

struct EnergyAndGradient {
    double energy{0.0};              // total energy, Hartree
    std::vector<Vec3> gradient;      // dE/dR_A, Hartree/bohr, one per atom
    bool converged{false};           // did the underlying SCF converge
    int scf_iterations{0};
};

using EnergyGradientFunction =
    std::function<EnergyAndGradient(const std::vector<Atom>& atoms)>;

struct HartreeFockJobSettings {
    std::string basis{"STO-3G"};
    int Na{0};
    int Nb{0};
    int max_scf_iterations{200};
    double scf_energy_tol{1e-11};
    double scf_gradient_tol{1e-9};
    // Reuse the previous geometry's converged density as the next initial guess.
    // Turned off, every step restarts from the core Hamiltonian.
    bool reuse_density_guess{true};
};

// Build an unrestricted Hartree-Fock energy+gradient callable.
EnergyGradientFunction MakeHartreeFockEnergyGradient(HartreeFockJobSettings settings);
