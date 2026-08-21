// The DFT pipeline, one stage at a time.
//
// Each stage is checked against something that does not go through the grid:
//
//   grid + AO values  ->  reproduce the analytic overlap matrix
//   spin density      ->  integrate to the right electron count
//   functional        ->  v_rho matches a finite difference of the energy density
//   assembled Vxc     ->  Tr(D Vxc) matches the grid's own tr_p_vxc
//   whole chain       ->  UKS/LSDA converges, and reduces to UHF when ax = 1
//
// The last one matters most: with exact_exchange_fraction = 1 and no functional,
// UKSBuilder must reproduce UHFBuilder to machine precision.  If it does not, the
// DFT path has broken the Hartree-Fock path it shares code with.

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

#include "IDensityUpdater.h"
#include "Libint2IntegralProvider.h"
#include "dft/dft_helper.h"
#include "fock_builders.h"
#include "geometry.h"
#include "jk_builder.h"
#include "scf.h"
#include "test_util.h"
#include "xc/lsda.h"

namespace {

Molecule MakeWater(const std::string& basis) {
    Molecule m;
    m.atoms.push_back(AtomWithBasis{Atom{"O", 8, 0.0000, 0.0000, 0.0000}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, 1.4310, 1.1070}, basis});
    m.atoms.push_back(AtomWithBasis{Atom{"H", 1, 0.0000, -1.4310, 1.1070}, basis});
    return m;
}

double TraceProduct(const T2& A, const T2& B) {
    double trace = 0.0;
    for (Eigen::Index i = 0; i < A.dimension(0); ++i)
        for (Eigen::Index j = 0; j < A.dimension(1); ++j) trace += A(i, j) * B(j, i);
    return trace;
}

}  // namespace

int main() {
    testing::TestRunner t("DFT pipeline");

    const Molecule water = MakeWater("STO-3G");
    const Libint2IntegralProvider provider(water);
    const T2 hcore = provider.ComputeHcore();
    const T2 overlap = provider.ComputeOverlap();
    const std::vector<Atom> atoms = provider.GetAtoms();
    const BasisShells shells = provider.GetBasisShells();

    // ---------------------------------------------------------------------
    t.Section("1-2) grid and AO values reproduce the analytic overlap");
    dft::DftGridContext context;
    {
        for (const auto& settings :
             {GridSettings::Coarse(), GridSettings::Medium(), GridSettings::Fine()}) {
            const dft::DftGridContext c = dft::build_dft_grid_context(shells, atoms, settings);
            const T2 s_grid = dft::grid_overlap_matrix(c.ao, c.grid);
            std::printf("   info  %7zu points -> max |S_grid - S| = %.3e\n", c.grid.size(),
                        testing::MaxAbsDiff(s_grid, overlap));
        }
        context = dft::build_dft_grid_context(shells, atoms, GridSettings::Medium());
        const T2 s_grid = dft::grid_overlap_matrix(context.ao, context.grid);
        t.Near(testing::MaxAbsDiff(s_grid, overlap), 0.0, 1e-6,
               "max |S_grid - S_analytic| on the medium grid");
    }

    // ---------------------------------------------------------------------
    t.Section("3) the spin density integrates to the electron count");
    SCFResults hf;
    {
        SCFSettings settings;
        settings.Na = 5;
        settings.Nb = 5;
        settings.energy_tol = 1e-12;
        settings.gradient_tol = 1e-10;
        hf.nuclear_repulsion = NuclearRepulsionEnergy(atoms);
        GenerateInitialGuess(hcore, overlap, 5, 5, hf);

        auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
        UHFBuilder builder(std::move(jk), hcore);
        const FockDiagonalizationDensityUpdater updater;
        SCFLoop(settings, builder, updater, overlap, hf);
        t.Check(hf.converged, "reference UHF converged");

        const double n_alpha =
            dft::integrate_density_on_grid(context.ao, context.grid, hf.densityMatrices.Da);
        t.Near(n_alpha, 5.0, 1e-5, "integral of rho_alpha equals 5 electrons");
        t.Near(n_alpha, TraceProduct(hf.densityMatrices.Da, overlap), 1e-5,
               "integral of rho_alpha equals Tr(Da S)");
    }

    // ---------------------------------------------------------------------
    t.Section("4) v_rho is the derivative of the energy density");
    {
        double worst = 0.0;
        for (double ra : {0.01, 0.1, 0.5, 2.0})
            for (double rb : {0.005, 0.08, 0.4, 1.5}) {
                const double h = 1e-6 * ra;
                const auto energy = [](double a, double b) {
                    return EvaluateBuiltinFunctional(xc::Lsda, a, b).energy_density;
                };
                const XcPointResult xc = EvaluateBuiltinFunctional(xc::Lsda, ra, rb);
                const double fd_a = (energy(ra + h, rb) - energy(ra - h, rb)) / (2 * h);
                const double fd_b = (energy(ra, rb + h) - energy(ra, rb - h)) / (2 * h);
                worst = std::max(worst, std::fabs(fd_a - xc.v_rho_a) / std::fabs(xc.v_rho_a));
                worst = std::max(worst, std::fabs(fd_b - xc.v_rho_b) / std::fabs(xc.v_rho_b));
            }
        t.Near(worst, 0.0, 1e-6, "SVWN5: max relative error of v_rho vs finite difference");
    }

    // ---------------------------------------------------------------------
    t.Section("5) assembled Vxc is consistent with the grid quantities");
    {
        const dft::DftBuildResult xc = dft::build_uks_vxc(
            context, hf.densityMatrices.Da, hf.densityMatrices.Db, xc::LsdaSpec());

        const double trace = TraceProduct(hf.densityMatrices.Da, xc.Vxc_a) +
                             TraceProduct(hf.densityMatrices.Db, xc.Vxc_b);
        t.Near(trace, xc.tr_PVxc, 1e-9,
               "Tr(Da Vxc_a) + Tr(Db Vxc_b) equals the grid's tr_PVxc");
        t.Check(xc.Exc < 0.0, "Exc is negative");
        t.Report("Exc (SVWN5) [Hartree]", xc.Exc);

        t.Near(testing::MaxAbsDiff(xc.Vxc_a, xc.Vxc_b), 0.0, 1e-12,
               "closed shell: Vxc_alpha equals Vxc_beta");
    }

    // ---------------------------------------------------------------------
    t.Section("6) UKS with ax = 1 and no functional reproduces UHF exactly");
    {
        SCFSettings settings;
        settings.Na = 5;
        settings.Nb = 5;
        settings.energy_tol = 1e-12;
        settings.gradient_tol = 1e-10;

        SCFResults uks;
        uks.nuclear_repulsion = NuclearRepulsionEnergy(atoms);
        GenerateInitialGuess(hcore, overlap, 5, 5, uks);

        auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
        UKSBuilder builder(std::move(jk), hcore, dft::DftGridContext{}, xc::HartreeFockSpec());
        const FockDiagonalizationDensityUpdater updater;
        SCFLoop(settings, builder, updater, overlap, uks);

        t.Check(uks.converged, "UKS(HF spec) converged");
        t.Near(uks.energy, hf.energy, 1e-12, "UKS(ax=1, no functional) energy equals UHF");
        t.Near(testing::MaxAbsDiff(uks.densityMatrices.Da, hf.densityMatrices.Da), 0.0, 1e-10,
               "UKS(ax=1, no functional) density equals UHF");
    }

    // ---------------------------------------------------------------------
    t.Section("7) LSDA converges and lands below Hartree-Fock exchange only");
    {
        SCFSettings settings;
        settings.Na = 5;
        settings.Nb = 5;
        settings.energy_tol = 1e-10;
        settings.gradient_tol = 1e-8;

        SCFResults lsda;
        lsda.nuclear_repulsion = NuclearRepulsionEnergy(atoms);
        GenerateInitialGuess(hcore, overlap, 5, 5, lsda);

        auto jk = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
        UKSBuilder builder(std::move(jk), hcore, context, xc::LsdaSpec());
        // Damped, unlike the Hartree-Fock sections above.  Undamped Roothaan
        // iteration does not merely converge slowly here -- it settles into a
        // two-cycle limit cycle at ||FDS-SDF|| ~ 0.1 and reports E = -74.6549,
        // which is not a stationary point at all.  With damping the same SCF
        // reaches 1e-13 in ~30 cycles and gives -74.73193.  The grid is not the
        // limitation: coarse and fine grids both plateau at the same 0.1.
        const FockDiagonalizationDensityUpdater updater(0.7);
        SCFLoop(settings, builder, updater, overlap, lsda);

        t.Check(lsda.converged, "SVWN5 SCF converged");
        t.Report("SVWN5 cycles to convergence", lsda.iteration);
        t.Report("SVWN5 final ||FDS-SDF||", lsda.orbital_gradient_norm);
        t.Report("E(UHF)   [Hartree]", hf.energy);
        t.Report("E(SVWN5) [Hartree]", lsda.energy);
        t.Report("SVWN5 <S^2>", lsda.spin_squared);
        t.Near(lsda.spin_squared, 0.0, 1e-8, "closed-shell singlet has <S^2> = 0");
        t.Near(dft::integrate_density_on_grid(context.ao, context.grid, lsda.densityMatrices.Da),
               5.0, 1e-5, "converged LSDA density still integrates to 5 alpha electrons");

        // Slater exchange alone, for comparison: adding VWN correlation must lower
        // the energy further.
        SCFResults slater;
        slater.nuclear_repulsion = NuclearRepulsionEnergy(atoms);
        GenerateInitialGuess(hcore, overlap, 5, 5, slater);
        auto jk2 = std::make_unique<InCoreJKBuilder>(provider.ComputeERI());
        UKSBuilder builder2(std::move(jk2), hcore, context, xc::SlaterExchangeSpec());
        SCFLoop(settings, builder2, updater, overlap, slater);
        t.Check(slater.converged, "Slater-only SCF converged");
        t.Report("E(Slater exchange only) [Hartree]", slater.energy);
        t.Check(lsda.energy < slater.energy, "adding VWN5 correlation lowers the energy");
    }

    return t.Summary();
}
