// geometry_optimizer.h — Steepest descent on the nuclear coordinates.
//
// The simplest thing that works: step downhill along -g, grow the step while the
// energy keeps falling, halve it and retry when it rises.  No Hessian is formed
// and none is approximated.
//
// That is genuinely a limitation rather than a simplification.  Steepest descent
// converges slowly on molecules because stretches and bends differ in curvature
// by an order of magnitude, so a step that suits one is wrong for the other; the
// standard fix is a quasi-Newton update (BFGS) building an approximate inverse
// Hessian from the gradient history, which needs no second derivatives either.
// This is the honest baseline that such an optimizer would be measured against.
//
// No translation/rotation projection is applied: the gradient already sums to
// zero over atoms (translational invariance), so the centre of mass does not
// drift, and rotations are harmless for an energy minimization.
#pragma once

#include <vector>

#include "basis.h"
#include "energy_gradient.h"
#include "types.h"

struct GeometryOptimizerSettings {
    int max_steps{50};
    // Convergence thresholds on the force, in Hartree/bohr.  The defaults match
    // the values commonly used for a "normal" optimization.
    double max_force_tol{4.5e-4};
    double rms_force_tol{3.0e-4};
    // Additionally require the energy to have stopped changing.
    double energy_tol{1e-6};

    double initial_step{0.5};       // bohr per (Hartree/bohr) of gradient
    double max_step{2.0};
    double max_displacement{0.3};   // cap on how far any one atom may move per step
    double step_grow{1.2};          // on an accepted step
    double step_shrink{0.5};        // on a rejected one
    int max_rejections{8};          // consecutive shrinks before giving up

    bool verbose{false};
};

struct GeometryOptimizationResult {
    std::vector<Atom> atoms;         // optimized geometry
    double energy{0.0};
    std::vector<Vec3> gradient;
    double max_force{0.0};
    double rms_force{0.0};
    int steps{0};
    int energy_evaluations{0};
    bool converged{false};
    std::vector<double> energy_history;
};

GeometryOptimizationResult OptimizeGeometry(
    const EnergyGradientFunction& energy_gradient,
    std::vector<Atom> atoms,
    const GeometryOptimizerSettings& settings = GeometryOptimizerSettings());

// Largest |component| of a gradient, and its root-mean-square over all 3N.
double MaxForce(const std::vector<Vec3>& gradient);
double RmsForce(const std::vector<Vec3>& gradient);
