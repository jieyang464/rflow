#include "scf.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "IDensityUpdater.h"
#include "fock_builders.h"
#include "linalg.h"

namespace {

double TraceProduct(const T2& A, const T2& B) {
    const Eigen::Index nbf = A.dimension(0);
    double trace = 0.0;
    for (Eigen::Index i = 0; i < nbf; ++i) {
        for (Eigen::Index j = 0; j < nbf; ++j) {
            trace += A(i, j) * B(j, i);
        }
    }
    return trace;
}

}  // namespace

void GenerateInitialGuess(const T2& Hcore, const T2& S, int Na, int Nb, SCFResults& scfResults) {
    const DiagResult diag = DiagonalizeInSMetric(Hcore, S);
    scfResults.moCoefficients.Ca = diag.C;
    scfResults.moCoefficients.Cb = diag.C;
    scfResults.fockEigenvalues.eps_a = diag.eps;
    scfResults.fockEigenvalues.eps_b = diag.eps;

    const Eigen::Index nbf = Hcore.dimension(0);
    const int occ_a = std::max<int>(0, std::min<int>(Na, static_cast<int>(nbf)));
    const int occ_b = std::max<int>(0, std::min<int>(Nb, static_cast<int>(nbf)));

    scfResults.densityMatrices.Da = BuildDensityMatrix(scfResults.moCoefficients.Ca, occ_a);
    scfResults.densityMatrices.Db = BuildDensityMatrix(scfResults.moCoefficients.Cb, occ_b);
}

void SCFLoop(const SCFSettings& settings,
             IFockBuilder& fock_builder,
             const IDensityUpdater& updater,
             const T2& S,
             SCFResults& scfResults) {
    double prev_energy = 0.0;
    scfResults.converged = false;

    // Allocated once and reused every cycle: the builder writes straight into it
    // and Eigen::Tensor does not reallocate when the shape is unchanged.  This is
    // the SCF's scratch, which is why it is a local and not part of SCFResults.
    FockBuildResult fock;

    for (int iter = 0; iter < settings.max_iter; ++iter) {
        scfResults.iteration = iter;

        const FockBuildInput build_in{scfResults.densityMatrices.Da,
                                      scfResults.densityMatrices.Db};
        fock_builder.build(build_in, fock);

        scfResults.electronic_energy = fock.electronic_energy;
        scfResults.energy = fock.electronic_energy + scfResults.nuclear_repulsion;

        // The gradient of the energy with respect to an orbital rotation is the
        // S-metric commutator FDS - SDF; its vanishing is the real convergence
        // test, and it means the same thing for either density updater.
        const double err_a =
            MaxAbsElement(ComputeDIISError(fock.Fa, scfResults.densityMatrices.Da, S));
        const double err_b =
            MaxAbsElement(ComputeDIISError(fock.Fb, scfResults.densityMatrices.Db, S));
        scfResults.orbital_gradient_norm = std::max(err_a, err_b);

        const double delta_energy = std::abs(scfResults.energy - prev_energy);

        if (settings.verbose) {
            std::cout << "  SCF " << iter << "  E = " << scfResults.energy
                      << "  dE = " << (iter > 0 ? delta_energy : 0.0)
                      << "  ||FDS-SDF|| = " << scfResults.orbital_gradient_norm << "\n";
        }

        if (iter > 0 && delta_energy < settings.energy_tol &&
            scfResults.orbital_gradient_norm < settings.gradient_tol) {
            scfResults.converged = true;
            break;
        }
        prev_energy = scfResults.energy;

        updater.UpdateDensity(fock.Fa, S, scfResults.densityMatrices.Da);
        updater.UpdateDensity(fock.Fb, S, scfResults.densityMatrices.Db);
    }

    // One copy, at exit: the converged Fock matrices are an output (gradients
    // need them), unlike the per-iteration buffer above.
    scfResults.fockMatrices.Fa = fock.Fa;
    scfResults.fockMatrices.Fb = fock.Fb;

    scfResults.spin_squared = ComputeSpinSquared(scfResults.densityMatrices.Da,
                                                 scfResults.densityMatrices.Db, S);

    if (settings.finalize_orbitals) {
        // Post-processing only, deliberately outside the iteration: the
        // commutator path never diagonalizes while converging.
        const DiagResult diag_a = DiagonalizeInSMetric(fock.Fa, S);
        const DiagResult diag_b = DiagonalizeInSMetric(fock.Fb, S);
        scfResults.moCoefficients.Ca = diag_a.C;
        scfResults.moCoefficients.Cb = diag_b.C;
        scfResults.fockEigenvalues.eps_a = diag_a.eps;
        scfResults.fockEigenvalues.eps_b = diag_b.eps;
    }
}

double ComputeSpinSquared(const T2& Da, const T2& Db, const T2& S) {
    const double n_alpha = TraceProduct(Da, S);
    const double n_beta = TraceProduct(Db, S);
    const double sz = 0.5 * (n_alpha - n_beta);

    // Tr(Da S Db S) counts the overlap between the occupied alpha and beta spaces.
    const double overlap_term = TraceProduct(MatMul(MatMul(Da, S), Db), S);
    return sz * (sz + 1.0) + n_beta - overlap_term;
}

void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao, T2& h_mo, T4& eri_mo) {
    const Eigen::Index nbf = C.dimension(0);

    // One-electron part: h_mo = C^T h_ao C.
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_second = {Eigen::IndexPair<int>(1, 0)};
    const T2 h_half = h_ao.contract(C, contract_second);
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_first = {Eigen::IndexPair<int>(0, 0)};
    h_mo = C.contract(h_half, contract_first);

    // Two-electron part: four quarter transforms, O(N^5) each rather than the
    // O(N^8) of a direct four-index contraction.
    T4 t1(nbf, nbf, nbf, nbf);
    T4 t2(nbf, nbf, nbf, nbf);
    T4 t3(nbf, nbf, nbf, nbf);
    eri_mo.resize(nbf, nbf, nbf, nbf);
    t1.setZero();
    t2.setZero();
    t3.setZero();
    eri_mo.setZero();

    for (Eigen::Index p = 0; p < nbf; ++p)
        for (Eigen::Index nu = 0; nu < nbf; ++nu)
            for (Eigen::Index lam = 0; lam < nbf; ++lam)
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index mu = 0; mu < nbf; ++mu)
                        value += C(mu, p) * eri_ao(mu, nu, lam, sig);
                    t1(p, nu, lam, sig) = value;
                }

    for (Eigen::Index p = 0; p < nbf; ++p)
        for (Eigen::Index q = 0; q < nbf; ++q)
            for (Eigen::Index lam = 0; lam < nbf; ++lam)
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index nu = 0; nu < nbf; ++nu)
                        value += C(nu, q) * t1(p, nu, lam, sig);
                    t2(p, q, lam, sig) = value;
                }

    for (Eigen::Index p = 0; p < nbf; ++p)
        for (Eigen::Index q = 0; q < nbf; ++q)
            for (Eigen::Index r = 0; r < nbf; ++r)
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index lam = 0; lam < nbf; ++lam)
                        value += C(lam, r) * t2(p, q, lam, sig);
                    t3(p, q, r, sig) = value;
                }

    for (Eigen::Index p = 0; p < nbf; ++p)
        for (Eigen::Index q = 0; q < nbf; ++q)
            for (Eigen::Index r = 0; r < nbf; ++r)
                for (Eigen::Index s = 0; s < nbf; ++s) {
                    double value = 0.0;
                    for (Eigen::Index sig = 0; sig < nbf; ++sig)
                        value += C(sig, s) * t3(p, q, r, sig);
                    eri_mo(p, q, r, s) = value;
                }
}
