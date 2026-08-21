#include "dft/dft_helper.h"

#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>

#include "xc/lsda.h"

#ifdef USE_LIBXC
#include "xc/libxc_wrapper.h"
#endif

namespace dft {
namespace {

using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

// Points are processed in blocks so the AO values can be handed to Eigen as one
// gemm instead of a rank-1 update per point, without holding an n_points x nbf
// temporary for the whole grid.
constexpr std::size_t kBlockSize = 2048;

Eigen::Map<const EigenMatrix> AsMatrix(const T2& t) {
    return Eigen::Map<const EigenMatrix>(t.data(), t.dimension(0), t.dimension(1));
}

T2 ToTensor(const EigenMatrix& m) {
    T2 out(m.rows(), m.cols());
    std::copy(m.data(), m.data() + m.size(), out.data());
    return out;
}

// AO values of one block of points as an (n_points x nbf) Eigen matrix.
EigenMatrix BlockAoMatrix(const GridAoValues& ao, std::size_t start, std::size_t count) {
    EigenMatrix block(static_cast<Eigen::Index>(count), ao.nbf);
    for (std::size_t g = 0; g < count; ++g) {
        const double* row = ao.at(start + g);
        for (int mu = 0; mu < ao.nbf; ++mu) {
            block(static_cast<Eigen::Index>(g), mu) = row[mu];
        }
    }
    return block;
}

XcPointResult EvaluatePoint(const xc::FunctionalSpec& functional, double rho_a, double rho_b) {
    XcPointResult result = EvaluateBuiltinFunctional(functional.builtin, rho_a, rho_b);

#ifdef USE_LIBXC
    for (const int id : functional.libxc_ids) {
        const SpinXCOutput xc = evaluate_spin_lda_point(id, SpinXCInput{rho_a, rho_b});
        result.energy_density += (rho_a + rho_b) * xc.eps_xc;
        result.v_rho_a += xc.v_rho_a;
        result.v_rho_b += xc.v_rho_b;
    }
#else
    if (!functional.libxc_ids.empty()) {
        throw std::runtime_error(
            "FunctionalSpec requests libxc functionals but this build has no libxc "
            "(rebuild with USE_LIBXC=1).");
    }
#endif

    return result;
}

void CheckSquare(const T2& D, int nbf, const char* name) {
    if (D.dimension(0) != nbf || D.dimension(1) != nbf) {
        throw std::invalid_argument(std::string("dft: ") + name +
                                    " does not match the number of basis functions.");
    }
}

}  // namespace

std::vector<GridPoint> build_molecular_grid(const std::vector<Atom>& atoms,
                                            const GridSettings& settings) {
    return BuildBeckeGrid(atoms, settings);
}

GridAoValues evaluate_aos_on_grid(const BasisShells& shells,
                                  const std::vector<GridPoint>& grid) {
    GridAoValues ao;
    ao.nbf = NumBasisFunctions(shells);
    if (ao.nbf <= 0) throw std::invalid_argument("evaluate_aos_on_grid: empty basis.");

    ao.num_points = grid.size();
    ao.values.assign(ao.num_points * static_cast<std::size_t>(ao.nbf), 0.0);

    for (std::size_t g = 0; g < grid.size(); ++g) {
        EvaluateAOs(shells, Vec3{{grid[g].x, grid[g].y, grid[g].z}},
                    ao.values.data() + g * ao.nbf);
    }
    return ao;
}

DftGridContext build_dft_grid_context(const BasisShells& shells, const std::vector<Atom>& atoms,
                                      const GridSettings& settings) {
    DftGridContext context;
    context.grid = build_molecular_grid(atoms, settings);
    context.ao = evaluate_aos_on_grid(shells, context.grid);
    return context;
}

GridDensity build_spin_density_on_grid(const GridAoValues& ao, const T2& Da, const T2& Db) {
    CheckSquare(Da, ao.nbf, "Da");
    CheckSquare(Db, ao.nbf, "Db");

    const auto Da_m = AsMatrix(Da);
    const auto Db_m = AsMatrix(Db);

    GridDensity density;
    density.rho_a.assign(ao.num_points, 0.0);
    density.rho_b.assign(ao.num_points, 0.0);

    for (std::size_t start = 0; start < ao.num_points; start += kBlockSize) {
        const std::size_t count = std::min(kBlockSize, ao.num_points - start);
        const EigenMatrix Phi = BlockAoMatrix(ao, start, count);

        // rho_sigma(g) = sum_mu phi_mu(g) (D phi)_mu(g)
        const EigenMatrix Ta = Phi * Da_m;
        const EigenMatrix Tb = Phi * Db_m;
        const Eigen::VectorXd ra = (Phi.array() * Ta.array()).rowwise().sum();
        const Eigen::VectorXd rb = (Phi.array() * Tb.array()).rowwise().sum();

        for (std::size_t g = 0; g < count; ++g) {
            // Quadrature noise can push a near-vacuum density slightly negative;
            // the functionals are undefined there.
            density.rho_a[start + g] = std::max(0.0, ra(static_cast<Eigen::Index>(g)));
            density.rho_b[start + g] = std::max(0.0, rb(static_cast<Eigen::Index>(g)));
        }
    }
    return density;
}

XcOnGrid evaluate_xc_on_grid(const GridDensity& density, const std::vector<GridPoint>& grid,
                             const xc::FunctionalSpec& functional) {
    if (density.rho_a.size() != grid.size() || density.rho_b.size() != grid.size()) {
        throw std::invalid_argument("evaluate_xc_on_grid: density and grid sizes disagree.");
    }

    XcOnGrid out;
    out.v_rho_a.assign(grid.size(), 0.0);
    out.v_rho_b.assign(grid.size(), 0.0);
    out.energy_density.assign(grid.size(), 0.0);

    for (std::size_t g = 0; g < grid.size(); ++g) {
        const double rho_a = density.rho_a[g];
        const double rho_b = density.rho_b[g];
        if (rho_a + rho_b < kDensityThreshold) continue;

        const XcPointResult xc = EvaluatePoint(functional, rho_a, rho_b);
        out.v_rho_a[g] = xc.v_rho_a;
        out.v_rho_b[g] = xc.v_rho_b;
        out.energy_density[g] = xc.energy_density;

        const double w = grid[g].weight;
        out.exc += w * xc.energy_density;
        // Same number as Tr(Da Vxc_a) + Tr(Db Vxc_b), but available here without
        // needing the assembled matrices.
        out.tr_p_vxc += w * (rho_a * xc.v_rho_a + rho_b * xc.v_rho_b);
    }
    return out;
}

T2 assemble_vxc_matrix(const GridAoValues& ao, const std::vector<GridPoint>& grid,
                       const std::vector<double>& v_xc) {
    if (v_xc.size() != grid.size() || ao.num_points != grid.size()) {
        throw std::invalid_argument("assemble_vxc_matrix: grid, AO values and v_xc disagree.");
    }

    EigenMatrix V = EigenMatrix::Zero(ao.nbf, ao.nbf);
    for (std::size_t start = 0; start < ao.num_points; start += kBlockSize) {
        const std::size_t count = std::min(kBlockSize, ao.num_points - start);
        const EigenMatrix Phi = BlockAoMatrix(ao, start, count);

        EigenMatrix Weighted(static_cast<Eigen::Index>(count), ao.nbf);
        for (std::size_t g = 0; g < count; ++g) {
            Weighted.row(static_cast<Eigen::Index>(g)) =
                grid[start + g].weight * v_xc[start + g] * Phi.row(static_cast<Eigen::Index>(g));
        }
        V.noalias() += Phi.transpose() * Weighted;
    }

    // Symmetric up to rounding by construction; make it exactly so.
    V = (0.5 * (V + V.transpose())).eval();
    return ToTensor(V);
}

DftBuildResult build_uks_vxc(const DftGridContext& context, const T2& Da, const T2& Db,
                             const xc::FunctionalSpec& functional) {
    DftBuildResult result;

    const Eigen::Index nbf = Da.dimension(0);
    result.Vxc_a = T2(nbf, nbf);
    result.Vxc_b = T2(nbf, nbf);
    result.Vxc_a.setZero();
    result.Vxc_b.setZero();

    if (!functional.NeedsGrid()) return result;
    if (context.empty()) {
        throw std::invalid_argument(
            "build_uks_vxc: the functional needs a grid but the context is empty.");
    }

    const GridDensity density = build_spin_density_on_grid(context.ao, Da, Db);
    const XcOnGrid xc = evaluate_xc_on_grid(density, context.grid, functional);

    result.Vxc_a = assemble_vxc_matrix(context.ao, context.grid, xc.v_rho_a);
    result.Vxc_b = assemble_vxc_matrix(context.ao, context.grid, xc.v_rho_b);
    result.Exc = xc.exc;
    result.tr_PVxc = xc.tr_p_vxc;
    return result;
}

DftBuildResult build_uks_vxc(const BasisShells& shells, const std::vector<Atom>& atoms,
                             const T2& Da, const T2& Db, const xc::FunctionalSpec& functional,
                             const GridSettings& settings) {
    return build_uks_vxc(build_dft_grid_context(shells, atoms, settings), Da, Db, functional);
}

T2 grid_overlap_matrix(const GridAoValues& ao, const std::vector<GridPoint>& grid) {
    return assemble_vxc_matrix(ao, grid, std::vector<double>(grid.size(), 1.0));
}

double integrate_density_on_grid(const GridAoValues& ao, const std::vector<GridPoint>& grid,
                                 const T2& D) {
    const GridDensity density = build_spin_density_on_grid(ao, D, D);
    double total = 0.0;
    for (std::size_t g = 0; g < grid.size(); ++g) total += grid[g].weight * density.rho_a[g];
    return total;
}

}  // namespace dft
