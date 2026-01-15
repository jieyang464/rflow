#include "grid.h"
#include "molecule.h"
#include <stdexcept>

// Placeholder for radial and angular grid generation.
// A real implementation would use Gauss-Chebyshev for radial and Lebedev for angular.
static std::vector<GridPoint> generate_atomic_grid(const Atom& atom, int radial_points, int angular_points) {
    std::vector<GridPoint> grid;
    if (radial_points <= 0 || angular_points <= 0) return grid;

    // Simplified grid for demonstration: concentric shells
    for (int i = 1; i <= radial_points; ++i) {
        double r = static_cast<double>(i); // Simplified radius
        // Simplified angular part: points on axes
        grid.push_back({{atom.x + r, atom.y, atom.z}, 1.0});
        grid.push_back({{atom.x - r, atom.y, atom.z}, 1.0});
        grid.push_back({{atom.x, atom.y + r, atom.z}, 1.0});
        grid.push_back({{atom.x, atom.y - r, atom.z}, 1.0});
        grid.push_back({{atom.x, atom.y, atom.z + r}, 1.0});
        grid.push_back({{atom.x, atom.y, atom.z - r}, 1.0});
    }
    // Normalize weights (crude approximation)
    double total_weight = grid.size();
    for (GridPoint& p : grid) {
        p.weight /= total_weight;
    }
    return grid;
}

// Becke partitioning function (simplified)
static double becke_partition(const scfcpp::Vector3& point, const Atom& atom_i, const std::vector<Atom>& all_atoms) {
    double p_i = 1.0;
    for (const Atom& atom_j : all_atoms) {
        if (&atom_i == &atom_j) continue;
        scfcpp::Vector3 r_i = point - scfcpp::Vector3(atom_i.x, atom_i.y, atom_i.z);
        scfcpp::Vector3 r_j = point - scfcpp::Vector3(atom_j.x, atom_j.y, atom_j.z);
        double mu = (r_i.norm() - r_j.norm()) / (r_i - r_j).norm();
        // A simple polynomial mapping mu -> [0,1]
        struct F {
            double operator()(double x) const { return 0.5 * (1.0 - x); }
        } f;
        p_i *= f(mu);
    }
    return p_i;
}


std::vector<GridPoint> generate_molecular_grid(const IMolecule& mol,
                                               int radial_points,
                                               int angular_points) {
    if (radial_points <= 0 || angular_points <= 0) {
        return {};
    }
    std::vector<GridPoint> molecular_grid;
    const std::vector<Atom>& atoms = mol.atoms();

    // 1. Generate raw atomic grids
    std::vector<std::vector<GridPoint>> atomic_grids;
    for (const Atom& atom : atoms) {
        atomic_grids.push_back(generate_atomic_grid(atom, radial_points, angular_points));
    }

    // 2. Combine grids and apply Becke partitioning
    double total_weight_sum = 0;
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (const GridPoint& point : atomic_grids[i]) {
            double total_p = 0;
            std::vector<double> p_k;
            for (size_t k = 0; k < atoms.size(); ++k) {
                double p = becke_partition(point.pos, atoms[k], atoms);
                p_k.push_back(p);
                total_p += p;
            }

            if (total_p > 1e-12) {
                double w_i = p_k[i] / total_p;
                molecular_grid.push_back({point.pos, w_i * point.weight});
                total_weight_sum += w_i * point.weight;
            }
        }
    }
    
    // Renormalize to match total number of electrons (a common practice)
    double N_elec = mol.electron_count();
    if (total_weight_sum > 1e-12) {
        double scale = N_elec / total_weight_sum;
        for (GridPoint& p : molecular_grid) {
            p.weight *= scale;
        }
    }

    return molecular_grid;
}
