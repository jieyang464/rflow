#include "geometry_optimizer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace {

std::vector<Atom> StepAlong(const std::vector<Atom>& atoms, const std::vector<Vec3>& gradient,
                            double step, double max_displacement) {
    // Cap the largest single-atom displacement so an early large gradient cannot
    // throw the geometry somewhere the SCF will not converge.
    double largest = 0.0;
    for (const auto& g : gradient) {
        for (int c = 0; c < 3; ++c) largest = std::max(largest, std::fabs(step * g[c]));
    }
    double scale = 1.0;
    if (largest > max_displacement && largest > 0.0) scale = max_displacement / largest;

    std::vector<Atom> moved = atoms;
    for (std::size_t a = 0; a < moved.size(); ++a) {
        moved[a].x -= scale * step * gradient[a][0];
        moved[a].y -= scale * step * gradient[a][1];
        moved[a].z -= scale * step * gradient[a][2];
    }
    return moved;
}

}  // namespace

double MaxForce(const std::vector<Vec3>& gradient) {
    double worst = 0.0;
    for (const auto& g : gradient) {
        for (int c = 0; c < 3; ++c) worst = std::max(worst, std::fabs(g[c]));
    }
    return worst;
}

double RmsForce(const std::vector<Vec3>& gradient) {
    if (gradient.empty()) return 0.0;
    double sum = 0.0;
    for (const auto& g : gradient) {
        for (int c = 0; c < 3; ++c) sum += g[c] * g[c];
    }
    return std::sqrt(sum / (3.0 * gradient.size()));
}

GeometryOptimizationResult OptimizeGeometry(const EnergyGradientFunction& energy_gradient,
                                            std::vector<Atom> atoms,
                                            const GeometryOptimizerSettings& settings) {
    if (!energy_gradient) {
        throw std::invalid_argument("OptimizeGeometry: empty energy/gradient function.");
    }
    if (atoms.empty()) throw std::invalid_argument("OptimizeGeometry: no atoms.");

    GeometryOptimizationResult result;
    EnergyAndGradient current = energy_gradient(atoms);
    result.energy_evaluations = 1;
    result.energy_history.push_back(current.energy);

    double step = settings.initial_step;

    for (int iteration = 0; iteration < settings.max_steps; ++iteration) {
        result.steps = iteration;
        const double max_force = MaxForce(current.gradient);
        const double rms_force = RmsForce(current.gradient);

        if (settings.verbose) {
            std::printf("  opt %2d  E = %.10f  |F|max = %.3e  |F|rms = %.3e  step = %.3f\n",
                        iteration, current.energy, max_force, rms_force, step);
        }

        if (max_force < settings.max_force_tol && rms_force < settings.rms_force_tol) {
            result.converged = true;
            break;
        }

        // Try a step; shrink and retry while the energy goes up.  Each retry is a
        // full SCF, which is why the step is grown on success rather than a line
        // search being run from scratch every iteration.
        bool accepted = false;
        for (int rejection = 0; rejection < settings.max_rejections; ++rejection) {
            const std::vector<Atom> trial =
                StepAlong(atoms, current.gradient, step, settings.max_displacement);
            const EnergyAndGradient probe = energy_gradient(trial);
            ++result.energy_evaluations;

            if (probe.converged && probe.energy < current.energy) {
                const double energy_change = current.energy - probe.energy;
                atoms = trial;
                current = probe;
                result.energy_history.push_back(current.energy);
                step = std::min(step * settings.step_grow, settings.max_step);
                accepted = true;

                // Converged if the energy has also stopped moving and the forces
                // are already small; checked at the top of the next iteration.
                if (energy_change < settings.energy_tol &&
                    MaxForce(current.gradient) < settings.max_force_tol) {
                    result.converged = true;
                }
                break;
            }
            step *= settings.step_shrink;
        }

        if (!accepted) break;   // could not find a downhill step
        if (result.converged) break;
    }

    result.atoms = std::move(atoms);
    result.energy = current.energy;
    result.gradient = current.gradient;
    result.max_force = MaxForce(current.gradient);
    result.rms_force = RmsForce(current.gradient);
    return result;
}
