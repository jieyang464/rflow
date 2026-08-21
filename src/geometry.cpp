#include "geometry.h"

#include <cmath>
#include <stdexcept>

namespace {

double Distance(const Atom& a, const Atom& b) {
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

}  // namespace

std::vector<Atom> AtomsOf(const Molecule& molecule) {
    std::vector<Atom> atoms;
    atoms.reserve(molecule.atoms.size());
    for (const auto& a : molecule.atoms) atoms.push_back(a.atom);
    return atoms;
}

double NuclearRepulsionEnergy(const std::vector<Atom>& atoms) {
    double energy = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            const double r = Distance(atoms[i], atoms[j]);
            if (r == 0.0) {
                throw std::runtime_error("NuclearRepulsionEnergy: coincident atoms.");
            }
            energy += static_cast<double>(atoms[i].atomic_number) *
                      static_cast<double>(atoms[j].atomic_number) / r;
        }
    }
    return energy;
}

std::vector<Vec3> NuclearRepulsionGradient(const std::vector<Atom>& atoms) {
    std::vector<Vec3> gradient(atoms.size(), Vec3{{0.0, 0.0, 0.0}});

    // d/dR_A sum_{B<A} Z_A Z_B / |R_A - R_B|  =  -sum_{B != A} Z_A Z_B (R_A - R_B)/R^3
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            const double r = Distance(atoms[i], atoms[j]);
            if (r == 0.0) {
                throw std::runtime_error("NuclearRepulsionGradient: coincident atoms.");
            }
            const double zz = static_cast<double>(atoms[i].atomic_number) *
                              static_cast<double>(atoms[j].atomic_number);
            const double scale = -zz / (r * r * r);
            const double d[3] = {atoms[i].x - atoms[j].x, atoms[i].y - atoms[j].y,
                                 atoms[i].z - atoms[j].z};
            for (int c = 0; c < 3; ++c) {
                gradient[i][c] += scale * d[c];
                gradient[j][c] -= scale * d[c];  // equal and opposite
            }
        }
    }
    return gradient;
}

std::vector<Atom> DisplaceAtom(const std::vector<Atom>& atoms, std::size_t index,
                               int axis, double delta) {
    if (index >= atoms.size()) throw std::out_of_range("DisplaceAtom: atom index out of range.");
    if (axis < 0 || axis > 2) throw std::out_of_range("DisplaceAtom: axis must be 0, 1 or 2.");

    std::vector<Atom> moved = atoms;
    switch (axis) {
        case 0: moved[index].x += delta; break;
        case 1: moved[index].y += delta; break;
        default: moved[index].z += delta; break;
    }
    return moved;
}
