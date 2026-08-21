// dft_helper.h — The DFT pipeline, as a sequence of plain functions.
//
//   geometry  ->  grid
//   grid      ->  AO values         (built once, reused every SCF iteration)
//   (Da, Db)  ->  spin densities on the grid
//   densities ->  XC potential and energy density on the grid
//   grid Vxc  ->  AO matrix
//
// Each step is a free function with explicit inputs, so each can be tested on
// its own: the grid against an analytic integral, the AO values against the
// overlap matrix, the functional against finite differences.  UKSBuilder is the
// only caller that runs the whole chain, and it is the only place in the project
// that knows what a functional is.
//
// Note on where the cost is
// -------------------------
// Evaluating phi_mu(r) at every grid point dominates everything else here -- it
// is n_points x nbf exponentials, against n_points x nbf^2 cheap arithmetic for
// the rest.  So the AO values are a first-class object (GridAoValues) that flows
// through the pipeline rather than being recomputed inside each step, and the
// grid plus AO values are bundled into a DftGridContext that UKSBuilder builds
// once and reuses.  Rebuilding them per SCF iteration would dominate the run.
#pragma once

#include <memory>
#include <vector>

#include "basis.h"
#include "grid/molecular_grid.h"
#include "types.h"
#include "xc/xc_config.h"

namespace dft {

// Values of every basis function at every grid point, row-major (point, mu).
struct GridAoValues {
    int nbf{0};
    std::size_t num_points{0};
    std::vector<double> values;  // num_points * nbf

    const double* at(std::size_t point) const { return values.data() + point * nbf; }
};

// Spin densities sampled on the grid.
//
// Kept as (rho_a, rho_b) rather than (rho, rho_a - rho_b): local functionals are
// defined on the two spin densities, and the spin-polarization form would only be
// converted straight back at every point.
struct GridDensity {
    std::vector<double> rho_a;
    std::vector<double> rho_b;
};

// XC quantities sampled on the grid.
struct XcOnGrid {
    std::vector<double> v_rho_a;        // d(e_xc)/d(rho_a) at each point
    std::vector<double> v_rho_b;
    std::vector<double> energy_density; // e_xc per unit *volume*, i.e. rho * eps_xc
    double exc{0.0};                    // sum_g w_g energy_density[g]
    double tr_p_vxc{0.0};               // sum_g w_g (rho_a v_a + rho_b v_b)[g]
                                        //   == Tr(Da Vxc_a) + Tr(Db Vxc_b)
};

// Everything that depends only on geometry and basis, so it survives every SCF
// iteration unchanged.
struct DftGridContext {
    std::vector<GridPoint> grid;
    GridAoValues ao;

    bool empty() const { return grid.empty(); }
};

struct DftBuildResult {
    T2 Vxc_a;
    T2 Vxc_b;
    double Exc{0.0};
    double tr_PVxc{0.0};
};

// --- 1) geometry -> grid ---------------------------------------------------
std::vector<GridPoint> build_molecular_grid(const std::vector<Atom>& atoms,
                                            const GridSettings& settings = GridSettings());

// --- 2) grid -> AO values --------------------------------------------------
GridAoValues evaluate_aos_on_grid(const BasisShells& shells,
                                  const std::vector<GridPoint>& grid);

// Steps 1 and 2 together; this is what UKSBuilder holds.
DftGridContext build_dft_grid_context(const BasisShells& shells,
                                      const std::vector<Atom>& atoms,
                                      const GridSettings& settings = GridSettings());

// --- 3) (Da, Db) -> spin densities on the grid -----------------------------
//   rho_sigma(g) = sum_{mu,nu} D^sigma_{mu,nu} phi_mu(g) phi_nu(g)
GridDensity build_spin_density_on_grid(const GridAoValues& ao, const T2& Da, const T2& Db);

// --- 4) densities -> XC on the grid ----------------------------------------
XcOnGrid evaluate_xc_on_grid(const GridDensity& density, const std::vector<GridPoint>& grid,
                             const xc::FunctionalSpec& functional);

// --- 5) grid Vxc -> AO matrix ----------------------------------------------
//   V_{mu,nu} = sum_g w_g v(g) phi_mu(g) phi_nu(g)
T2 assemble_vxc_matrix(const GridAoValues& ao, const std::vector<GridPoint>& grid,
                       const std::vector<double>& v_xc);

// --- the whole chain, for UKSBuilder ---------------------------------------
DftBuildResult build_uks_vxc(const DftGridContext& context, const T2& Da, const T2& Db,
                             const xc::FunctionalSpec& functional);

// Convenience form that builds the context too.  Fine for a one-off or a test;
// calling it from an SCF loop rebuilds the grid and the AO values every
// iteration, which is the expensive mistake this header is arranged to avoid.
DftBuildResult build_uks_vxc(const BasisShells& shells, const std::vector<Atom>& atoms,
                             const T2& Da, const T2& Db,
                             const xc::FunctionalSpec& functional,
                             const GridSettings& settings = GridSettings());

// --- diagnostics used by the tests -----------------------------------------
// Overlap matrix reproduced by the quadrature; compare with the analytic S from
// the integral provider to validate grid and AO evaluation together.
T2 grid_overlap_matrix(const GridAoValues& ao, const std::vector<GridPoint>& grid);

// Number of electrons recovered by integrating a density matrix on the grid;
// should equal Tr(D S).
double integrate_density_on_grid(const GridAoValues& ao, const std::vector<GridPoint>& grid,
                                 const T2& D);

}  // namespace dft
