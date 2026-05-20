#include "scf.h"
#include "fock_builders.h"
#include "IDensityUpdater.h"
#include "linalg.h"
#include <algorithm>
#include <cmath>
#include <iostream>

inline double TraceProduct(const T2& A, const T2& B) {
    const Eigen::Index nbf = A.dimension(0);
    double trace = 0.0;
    for (Eigen::Index i = 0; i < nbf; ++i) {
        for (Eigen::Index j = 0; j < nbf; ++j) {
            trace += A(i, j) * B(j, i);
        }
    }
    return trace;
}

void GenerateInitialGuess(const T2& Hcore, const T2& S, int Na, int Nb, SCFResults& scfResults) {
    auto diag = DiagonalizeInSMetric(Hcore, S);
    scfResults.moCoefficients.Ca = diag.C;
    scfResults.moCoefficients.Cb = diag.C;
    scfResults.fockEigenvalues.eps_a = diag.eps;
    scfResults.fockEigenvalues.eps_b = diag.eps;

    const Eigen::Index nbf = Hcore.dimension(0);
    int occ_a = std::max<int>(0, std::min<int>(Na, nbf));
    int occ_b = std::max<int>(0, std::min<int>(Nb, nbf));

    scfResults.densityMatrices.Da = BuildDensityMatrix(scfResults.moCoefficients.Ca, occ_a);
    scfResults.densityMatrices.Db = BuildDensityMatrix(scfResults.moCoefficients.Cb, occ_b);
}

void SCFLoop(const SCFSettings& settings,
             IFockBuilder& fock_builder,
             const IDensityUpdater& updater,
             const T2& S, const T2& Hcore,
             SCFResults& scfResults) {

    double prev_energy = 0.0;
    const int max_iter = settings.max_iter;
    const double energy_tol = settings.energy_tol;

    for (int iter = 0; iter < max_iter; ++iter) {
        scfResults.iteration = iter;

        // Build Fock matrix
        FockBuildInput build_in{Hcore, S, scfResults.densityMatrices.Da, scfResults.densityMatrices.Db};
        auto build_out = fock_builder.build(build_in);

        scfResults.fockMatrices.Fa = build_out.Fa;
        scfResults.fockMatrices.Fb = build_out.Fb;

        // Compute electronic energy
        const double e_elec = 0.5 * (
            TraceProduct(scfResults.densityMatrices.Da, Hcore + scfResults.fockMatrices.Fa) +
            TraceProduct(scfResults.densityMatrices.Db, Hcore + scfResults.fockMatrices.Fb)
        );
        scfResults.energy = e_elec + scfResults.nuclear_repulsion + build_out.energy_correction;

        std::cout << "SCF cycle " << iter << " energy=" << scfResults.energy << "\n";

        if (iter > 0) {
            if (std::abs(scfResults.energy - prev_energy) < energy_tol) {
                break;
            }
        }
        prev_energy = scfResults.energy;

        // Update density matrices using the provided IDensityUpdater
        updater.UpdateDensity(scfResults.fockMatrices.Fa, S, scfResults.densityMatrices.Da);
        updater.UpdateDensity(scfResults.fockMatrices.Fb, S, scfResults.densityMatrices.Db);
    }
}

void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao, T2& h_mo, T4& eri_mo) {
    const Eigen::Index nbf = C.dimension(0);
    h_mo.resize(nbf, nbf);
    eri_mo.resize(nbf, nbf, nbf, nbf);

    const Eigen::array<Eigen::IndexPair<int>, 1> h1 = {
        Eigen::IndexPair<int>(1, 0)};
    T2 t = h_ao.contract(C, h1); // (mu,q)

    const Eigen::array<Eigen::IndexPair<int>, 1> h2 = {
        Eigen::IndexPair<int>(0, 0)};
    h_mo = C.contract(t, h2); // (p,q)

    T4 t1(nbf, nbf, nbf, nbf);  // (p, nu, lam, sig)
    T4 t2(nbf, nbf, nbf, nbf);  // (p, q, lam, sig)
    T4 t3(nbf, nbf, nbf, nbf);  // (p, q, r, sig)
    t1.setZero();
    t2.setZero();
    t3.setZero();
    eri_mo.setZero();

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index nu = 0; nu < nbf; ++nu) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index mu = 0; mu < nbf; ++mu) {
                        value += C(mu, p) * eri_ao(mu, nu, lam, sig);
                    }
                    t1(p, nu, lam, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index nu = 0; nu < nbf; ++nu) {
                        value += C(nu, q) * t1(p, nu, lam, sig);
                    }
                    t2(p, q, lam, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index r = 0; r < nbf; ++r) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double value = 0.0;
                    for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                        value += C(lam, r) * t2(p, q, lam, sig);
                    }
                    t3(p, q, r, sig) = value;
                }
            }
        }
    }

    for (Eigen::Index p = 0; p < nbf; ++p) {
        for (Eigen::Index q = 0; q < nbf; ++q) {
            for (Eigen::Index r = 0; r < nbf; ++r) {
                for (Eigen::Index s = 0; s < nbf; ++s) {
                    double value = 0.0;
                    for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                        value += C(sig, s) * t3(p, q, r, sig);
                    }
                    eri_mo(p, q, r, s) = value;
                }
            }
        }
    }
}
