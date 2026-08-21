// geometry.h — Quantities that depend on the nuclei alone.
//
// Nuclear repulsion is not an integral over basis functions, so it does not
// belong on IIntegralProvider: every provider would have to reimplement the same
// loop, and one of them would eventually disagree with the others.  (Exactly
// that happened here -- the Szabo provider used to return the book's rounded
// 1.3669 rather than 2/1.4632, a 3.3e-5 discrepancy visible in the test output.)
//
// It is also not something to contract with a density.  Both the energy and the
// gradient are simply *added* by whoever assembles a total, which is why they
// live here rather than on the integral or derivative providers.
#pragma once

#include <vector>

#include "basis.h"
#include "types.h"

// The bare nuclei of a molecule, dropping the per-atom basis assignment.
std::vector<Atom> AtomsOf(const Molecule& molecule);

// V_nn = sum_{A<B} Z_A Z_B / R_AB, in Hartree, with coordinates in bohr.
double NuclearRepulsionEnergy(const std::vector<Atom>& atoms);

// dV_nn/dR_A, one Vec3 per atom.
std::vector<Vec3> NuclearRepulsionGradient(const std::vector<Atom>& atoms);

// Copy of `atoms` with atom `index` displaced by `delta` along `axis` (0=x,1=y,
// 2=z).  Used by finite-difference derivatives and by gradient tests.
std::vector<Atom> DisplaceAtom(const std::vector<Atom>& atoms, std::size_t index,
                               int axis, double delta);
