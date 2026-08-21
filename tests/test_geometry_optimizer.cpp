// Geometry optimization, validated against a reference that does not use the
// gradient at all.
//
// For H2 the equilibrium bond length is located by scanning the energy on a fine
// grid and fitting a parabola through the minimum -- energies only, no
// derivatives.  If the analytic gradient had a systematic error the optimizer
// would settle somewhere else, and the two numbers would disagree.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "energy_gradient.h"
#include "geometry.h"
#include "geometry_optimizer.h"
#include "test_util.h"

namespace {

std::vector<Atom> MakeH2(double bond_length_bohr) {
    return {Atom{"H", 1, 0.0, 0.0, 0.0}, Atom{"H", 1, 0.0, 0.0, bond_length_bohr}};
}

double Distance(const Atom& a, const Atom& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Minimum of a fine energy scan, refined by a parabola through the three points
// bracketing the lowest sample.  Uses energies only.
double ScanForMinimum(const EnergyGradientFunction& job, double from, double to, double spacing) {
    std::vector<double> radii, energies;
    for (double r = from; r <= to + 1e-12; r += spacing) {
        radii.push_back(r);
        energies.push_back(job(MakeH2(r)).energy);
    }

    std::size_t best = 0;
    for (std::size_t i = 1; i < energies.size(); ++i) {
        if (energies[i] < energies[best]) best = i;
    }
    if (best == 0 || best + 1 >= energies.size()) return radii[best];

    // Vertex of the parabola through (r-h, r, r+h).
    const double e0 = energies[best - 1], e1 = energies[best], e2 = energies[best + 1];
    const double denominator = e0 - 2.0 * e1 + e2;
    if (std::fabs(denominator) < 1e-14) return radii[best];
    return radii[best] - 0.5 * spacing * (e2 - e0) / denominator;
}

}  // namespace

int main() {
    testing::TestRunner t("geometry optimization (steepest descent, no Hessian)");

    // ---------------------------------------------------------------------
    t.Section("H2 / STO-3G: optimizer vs an energy-only scan");
    {
        HartreeFockJobSettings job_settings;
        job_settings.basis = "STO-3G";
        job_settings.Na = 1;
        job_settings.Nb = 1;
        const EnergyGradientFunction job = MakeHartreeFockEnergyGradient(job_settings);

        const double scan_minimum = ScanForMinimum(job, 1.20, 1.60, 0.005);
        t.Report("bond length from energy scan  [bohr]", scan_minimum);

        GeometryOptimizerSettings opt;
        opt.max_force_tol = 1e-5;   // tighter than default, to compare positions
        opt.rms_force_tol = 1e-5;
        const GeometryOptimizationResult result =
            OptimizeGeometry(job, MakeH2(1.80), opt);  // start well stretched

        const double optimized = Distance(result.atoms[0], result.atoms[1]);
        t.Check(result.converged, "optimizer reported convergence");
        t.Report("bond length from optimizer     [bohr]", optimized);
        t.Report("optimized energy               [Hartree]", result.energy);
        t.Report("steps / energy evaluations", result.energy_evaluations);

        t.Near(optimized, scan_minimum, 2e-3,
               "optimized bond length matches the scan minimum");
        t.Near(result.max_force, 0.0, opt.max_force_tol, "final |F|max below threshold");

        // The optimizer must never accept an uphill step.
        bool monotonic = true;
        for (std::size_t i = 1; i < result.energy_history.size(); ++i) {
            if (result.energy_history[i] > result.energy_history[i - 1]) monotonic = false;
        }
        t.Check(monotonic, "accepted energies decrease monotonically");
        t.Check(result.energy < job(MakeH2(1.80)).energy,
                "optimized energy is below the starting energy");
    }

    // ---------------------------------------------------------------------
    t.Section("H2O / STO-3G from a distorted start");
    {
        HartreeFockJobSettings job_settings;
        job_settings.basis = "STO-3G";
        job_settings.Na = 5;
        job_settings.Nb = 5;
        const EnergyGradientFunction job = MakeHartreeFockEnergyGradient(job_settings);

        const std::vector<Atom> start = {
            Atom{"O", 8, 0.00, 0.00, 0.00},
            Atom{"H", 1, 0.00, 1.70, 1.30},
            Atom{"H", 1, 0.00, -1.70, 1.30},
        };
        const double start_energy = job(start).energy;

        GeometryOptimizerSettings opt;
        opt.max_steps = 60;
        const GeometryOptimizationResult result = OptimizeGeometry(job, start, opt);

        const double r1 = Distance(result.atoms[0], result.atoms[1]);
        const double r2 = Distance(result.atoms[0], result.atoms[2]);
        // Bond angle at the oxygen, in degrees.
        const double angle = [&] {
            const auto& o = result.atoms[0];
            const double a[3] = {result.atoms[1].x - o.x, result.atoms[1].y - o.y,
                                 result.atoms[1].z - o.z};
            const double b[3] = {result.atoms[2].x - o.x, result.atoms[2].y - o.y,
                                 result.atoms[2].z - o.z};
            const double dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            return std::acos(dot / (r1 * r2)) * 180.0 / 3.14159265358979323846;
        }();

        t.Check(result.converged, "optimizer reported convergence");
        t.Report("starting energy   [Hartree]", start_energy);
        t.Report("optimized energy  [Hartree]", result.energy);
        t.Report("O-H bond length   [bohr]", r1);
        t.Report("O-H bond length   [bohr]", r2);
        t.Report("H-O-H angle       [degrees]", angle);
        t.Report("energy evaluations", result.energy_evaluations);

        t.Check(result.energy < start_energy, "energy decreased from the distorted start");
        t.Near(result.max_force, 0.0, opt.max_force_tol, "final |F|max below threshold");
        t.Near(r1, r2, 1e-3, "the two O-H bonds relax to the same length");

        // Sanity, not a literature comparison: HF/STO-3G water is near 1.8 bohr
        // and 100 degrees, so anything far from that means the gradient is wrong.
        t.Check(r1 > 1.5 && r1 < 2.1, "O-H bond length is physically sensible");
        t.Check(angle > 90.0 && angle < 115.0, "H-O-H angle is physically sensible");
    }

    return t.Summary();
}
