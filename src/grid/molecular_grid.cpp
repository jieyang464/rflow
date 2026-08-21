#include "grid/molecular_grid.h"

#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kAngstromToBohr = 1.8897261254578281;

// Bragg-Slater radii in Angstrom (Slater, JCP 41, 3199 (1964)); H is given
// Becke's value of 0.35 A rather than Slater's 0.25 A.
const double kBraggSlaterAngstrom[] = {
    /* 0  */ 0.35,
    /* H  */ 0.35, /* He */ 0.35,
    /* Li */ 1.45, /* Be */ 1.05, /* B  */ 0.85, /* C  */ 0.70, /* N  */ 0.65,
    /* O  */ 0.60, /* F  */ 0.50, /* Ne */ 0.45,
    /* Na */ 1.80, /* Mg */ 1.50, /* Al */ 1.25, /* Si */ 1.10, /* P  */ 1.00,
    /* S  */ 1.00, /* Cl */ 1.00, /* Ar */ 1.00,
    /* K  */ 2.20, /* Ca */ 1.80, /* Sc */ 1.60, /* Ti */ 1.40, /* V  */ 1.35,
    /* Cr */ 1.40, /* Mn */ 1.40, /* Fe */ 1.40, /* Co */ 1.35, /* Ni */ 1.35,
    /* Cu */ 1.35, /* Zn */ 1.35, /* Ga */ 1.30, /* Ge */ 1.25, /* As */ 1.15,
    /* Se */ 1.15, /* Br */ 1.15, /* Kr */ 1.10,
};
constexpr int kBraggSlaterCount = sizeof(kBraggSlaterAngstrom) / sizeof(double);

// Becke's iterated cutoff polynomial, eqns (19)-(21).
double BeckeStepFunction(double mu, int iterations) {
    double f = mu;
    for (int i = 0; i < iterations; ++i) f = 1.5 * f - 0.5 * f * f * f;
    return 0.5 * (1.0 - f);
}

// Gauss-Legendre nodes/weights on [-1,1] via Newton iteration on P_n.
void GaussLegendre(int n, std::vector<double>& nodes, std::vector<double>& weights) {
    nodes.assign(n, 0.0);
    weights.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        // Initial guess (Abramowitz & Stegun 25.4.30).
        double x = std::cos(kPi * (i + 0.75) / (n + 0.5));
        double dp = 0.0;
        for (int iter = 0; iter < 100; ++iter) {
            double p0 = 1.0, p1 = 0.0;
            for (int k = 0; k < n; ++k) {
                const double p2 = p1;
                p1 = p0;
                p0 = ((2.0 * k + 1.0) * x * p1 - k * p2) / (k + 1.0);
            }
            dp = n * (x * p0 - p1) / (x * x - 1.0);
            const double dx = -p0 / dp;
            x += dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        nodes[i] = x;
        weights[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
}

}  // namespace

double BraggSlaterRadius(int atomic_number) {
    if (atomic_number < 0) {
        throw std::invalid_argument("BraggSlaterRadius: negative atomic number.");
    }
    const int z = atomic_number < kBraggSlaterCount ? atomic_number : kBraggSlaterCount - 1;
    return kBraggSlaterAngstrom[z] * kAngstromToBohr;
}

GridSettings GridSettings::Coarse() { return GridSettings{32, 12, 24, 3, 1e-14}; }
GridSettings GridSettings::Medium() { return GridSettings{}; }
GridSettings GridSettings::Fine() { return GridSettings{96, 32, 64, 3, 1e-16}; }

std::vector<GridPoint> BuildBeckeGrid(const std::vector<Atom>& atoms,
                                      const GridSettings& settings) {
    if (atoms.empty()) return {};
    if (settings.n_radial < 1 || settings.n_theta < 1 || settings.n_phi < 1) {
        throw std::invalid_argument("BuildBeckeGrid: grid dimensions must be positive.");
    }

    const std::size_t natoms = atoms.size();

    // --- Becke atomic size adjustment, eqns (A2)-(A6) -------------------------
    // a_ij = u_ij / (u_ij^2 - 1), clipped to |a| <= 0.5, with u = (X-1)/(X+1)
    // and X the ratio of Bragg-Slater radii.
    std::vector<double> radii(natoms);
    for (std::size_t i = 0; i < natoms; ++i) radii[i] = BraggSlaterRadius(atoms[i].atomic_number);

    std::vector<double> a_ij(natoms * natoms, 0.0);
    for (std::size_t i = 0; i < natoms; ++i) {
        for (std::size_t j = 0; j < natoms; ++j) {
            if (i == j) continue;
            const double chi = radii[i] / radii[j];
            const double u = (chi - 1.0) / (chi + 1.0);
            double a = u / (u * u - 1.0);
            if (a > 0.5) a = 0.5;
            if (a < -0.5) a = -0.5;
            a_ij[i * natoms + j] = a;
        }
    }

    // Internuclear distances.
    std::vector<double> r_ij(natoms * natoms, 0.0);
    for (std::size_t i = 0; i < natoms; ++i) {
        for (std::size_t j = 0; j < natoms; ++j) {
            const double dx = atoms[i].x - atoms[j].x;
            const double dy = atoms[i].y - atoms[j].y;
            const double dz = atoms[i].z - atoms[j].z;
            const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (i != j && d < 1e-10) {
                throw std::runtime_error("BuildBeckeGrid: coincident atoms.");
            }
            r_ij[i * natoms + j] = d;
        }
    }

    // --- Angular quadrature (shared by all atoms and all radial shells) ------
    std::vector<double> cos_theta, w_theta;
    GaussLegendre(settings.n_theta, cos_theta, w_theta);

    struct AngularPoint {
        double ux, uy, uz, w;
    };
    std::vector<AngularPoint> angular;
    angular.reserve(static_cast<std::size_t>(settings.n_theta) * settings.n_phi);
    const double dphi = 2.0 * kPi / settings.n_phi;
    for (int t = 0; t < settings.n_theta; ++t) {
        const double ct = cos_theta[t];
        const double st = std::sqrt(std::max(0.0, 1.0 - ct * ct));
        for (int p = 0; p < settings.n_phi; ++p) {
            const double phi = p * dphi;
            angular.push_back({st * std::cos(phi), st * std::sin(phi), ct, w_theta[t] * dphi});
        }
    }

    std::vector<GridPoint> grid;
    grid.reserve(natoms * settings.n_radial * angular.size());

    std::vector<double> cell(natoms, 0.0);

    for (std::size_t a = 0; a < natoms; ++a) {
        // Becke's radial map: r = R_m (1+x)/(1-x) with R_m = half the
        // Bragg-Slater radius, except for hydrogen where R_m = R_BS.
        const double r_m =
            (atoms[a].atomic_number == 1) ? radii[a] : 0.5 * radii[a];

        for (int i = 1; i <= settings.n_radial; ++i) {
            // Gauss-Chebyshev (2nd kind) node/weight on [-1,1].
            const double angle = i * kPi / (settings.n_radial + 1);
            const double x = std::cos(angle);
            const double sin_angle = std::sin(angle);
            const double one_minus_x = 1.0 - x;
            if (one_minus_x < 1e-14) continue;

            const double r = r_m * (1.0 + x) / one_minus_x;
            // w_GC / sqrt(1-x^2) = pi/(n+1) * sin(angle); times r^2 dr/dx.
            const double dr_dx = 2.0 * r_m / (one_minus_x * one_minus_x);
            const double w_radial =
                (kPi / (settings.n_radial + 1)) * sin_angle * r * r * dr_dx;
            if (std::fabs(w_radial) < settings.weight_cutoff) continue;

            for (const auto& ang : angular) {
                const double px = atoms[a].x + r * ang.ux;
                const double py = atoms[a].y + r * ang.uy;
                const double pz = atoms[a].z + r * ang.uz;

                // --- Becke nuclear weight w_a(r) ---------------------------
                double total = 0.0;
                for (std::size_t i_at = 0; i_at < natoms; ++i_at) {
                    double p_i = 1.0;
                    const double dxi = px - atoms[i_at].x;
                    const double dyi = py - atoms[i_at].y;
                    const double dzi = pz - atoms[i_at].z;
                    const double ri = std::sqrt(dxi * dxi + dyi * dyi + dzi * dzi);
                    for (std::size_t j_at = 0; j_at < natoms; ++j_at) {
                        if (i_at == j_at) continue;
                        const double dxj = px - atoms[j_at].x;
                        const double dyj = py - atoms[j_at].y;
                        const double dzj = pz - atoms[j_at].z;
                        const double rj = std::sqrt(dxj * dxj + dyj * dyj + dzj * dzj);
                        double mu = (ri - rj) / r_ij[i_at * natoms + j_at];
                        // Atomic size adjustment: nu = mu + a(1 - mu^2).
                        mu += a_ij[i_at * natoms + j_at] * (1.0 - mu * mu);
                        if (mu > 1.0) mu = 1.0;
                        if (mu < -1.0) mu = -1.0;
                        p_i *= BeckeStepFunction(mu, settings.becke_iterations);
                        if (p_i == 0.0) break;
                    }
                    cell[i_at] = p_i;
                    total += p_i;
                }
                if (total <= 0.0) continue;

                const double weight = w_radial * ang.w * (cell[a] / total);
                if (std::fabs(weight) < settings.weight_cutoff) continue;

                grid.push_back(GridPoint{px, py, pz, weight});
            }
        }
    }

    return grid;
}
