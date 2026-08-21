// molecular_grid.h — Atom-centred numerical integration grid for DFT.
//
// Becke's fuzzy-Voronoi scheme [JCP 88, 2547 (1988)]: the molecular integral
//   \int f(r) dr  =  sum_A \int w_A(r) f(r) dr
// is split into atomic contributions by a smooth nuclear weight function w_A,
// and each atomic piece is done in spherical polar coordinates around A.
//
//   radial   : Gauss-Chebyshev (2nd kind) on [-1,1] mapped to [0,inf) by
//              Becke's r = R_m (1+x)/(1-x)
//   angular  : Gauss-Legendre in cos(theta) x uniform trapezoid in phi
//              (exact for spherical harmonics up to l < min(2*n_theta, n_phi/2))
//
// A Lebedev grid would need far fewer points for the same accuracy; the product
// grid is used here because it is a dozen lines of code and needs no tables.
#pragma once

#include <vector>

#include "basis.h"
#include "types.h"

struct GridPoint {
    double x{0.0}, y{0.0}, z{0.0};
    double weight{0.0};
};

struct GridSettings {
    int n_radial{64};
    int n_theta{20};   // Gauss-Legendre nodes in cos(theta)
    int n_phi{40};     // uniform nodes in phi
    int becke_iterations{3};
    double weight_cutoff{1e-14};  // drop points that cannot contribute

    static GridSettings Coarse();
    static GridSettings Medium();  // the default above
    static GridSettings Fine();
};

// Build the full molecular grid.  Coordinates are in bohr, matching Atom.
std::vector<GridPoint> BuildBeckeGrid(const std::vector<Atom>& atoms,
                                      const GridSettings& settings = GridSettings());

// Bragg-Slater atomic radius in bohr (used for Becke's atomic size adjustment).
double BraggSlaterRadius(int atomic_number);
