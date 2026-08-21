// vxc_grid.h — Exchange-correlation potential by numerical quadrature.
//
//   rho_sigma(r)  = sum_{mu,nu} D^sigma_{mu,nu} phi_mu(r) phi_nu(r)
//   E_xc          = sum_g w_g e_xc(rho_a(g), rho_b(g))
//   V^sigma_{mu,nu} = sum_g w_g (d e_xc / d rho_sigma)(g) phi_mu(g) phi_nu(g)
//
// The evaluator caches the AO values on the grid, so an SCF iteration costs one
// pass of O(n_points * nbf^2).
#pragma once

#include <cstddef>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "IIntegralProvider.h"
#include "basis.h"
#include "grid/molecular_grid.h"
#include "xc/vxc_evaluator.h"

class GridVxcEvaluator {
public:
    GridVxcEvaluator(BasisShells shells, std::vector<GridPoint> points);

    // Convenience: take geometry and basis straight from an integral provider.
    static std::shared_ptr<const GridVxcEvaluator> FromProvider(
        const IIntegralProvider& provider,
        const GridSettings& settings = GridSettings());

    void operator()(const UksDensityInput& in, UksVxcOutput& out) const;

    int nbf() const { return nbf_; }
    std::size_t num_points() const { return points_.size(); }

    // --- Diagnostics, used by the unit tests -------------------------------
    // Number of electrons recovered by integrating a density matrix: Tr over
    // the grid of sum_{mu,nu} D_{mu,nu} phi_mu phi_nu.  Should equal Tr(D S).
    double IntegrateDensity(const T2& D) const;
    // Overlap matrix reproduced by the quadrature; compare against the
    // analytic S from the integral provider to validate the grid.
    T2 GridOverlap() const;

private:
    BasisShells shells_;
    std::vector<GridPoint> points_;
    Eigen::MatrixXd ao_values_;  // (n_points x nbf)
    int nbf_{0};
};

// Wrap an evaluator as the VxcFunctor the Fock builder expects.
VxcFunctor MakeGridVxcFunctor(std::shared_ptr<const GridVxcEvaluator> evaluator);
