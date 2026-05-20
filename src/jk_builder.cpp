#include "jk_builder.h"

void DirectJKBuilder::build_JK(const T2& Da, const T2& Db, T2& J, T2& Ka, T2& Kb) const {
    const Eigen::Index nbf = Da.dimension(0);

    for (Eigen::Index mu = 0; mu < nbf; ++mu) {
        for (Eigen::Index nu = 0; nu < nbf; ++nu) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double eri_val = provider_.ComputeERI(mu, nu, lam, sig);
                    double eri_exch = provider_.ComputeERI(mu, lam, nu, sig);

                    J(mu, nu) += (Da(lam, sig) + Db(lam, sig)) * eri_val;
                    Ka(mu, nu) += Da(lam, sig) * eri_exch;
                    Kb(mu, nu) += Db(lam, sig) * eri_exch;
                }
            }
        }
    }
}

void InCoreJKBuilder::build_JK(const T2& Da, const T2& Db, T2& J, T2& Ka, T2& Kb) const {
    const Eigen::Index nbf = Da.dimension(0);

    // Can use Eigen contraction for more performance later, but this loops over in-memory T4
    for (Eigen::Index mu = 0; mu < nbf; ++mu) {
        for (Eigen::Index nu = 0; nu < nbf; ++nu) {
            for (Eigen::Index lam = 0; lam < nbf; ++lam) {
                for (Eigen::Index sig = 0; sig < nbf; ++sig) {
                    double eri_val = eri_(mu, nu, lam, sig);
                    double eri_exch = eri_(mu, lam, nu, sig);

                    J(mu, nu) += (Da(lam, sig) + Db(lam, sig)) * eri_val;
                    Ka(mu, nu) += Da(lam, sig) * eri_exch;
                    Kb(mu, nu) += Db(lam, sig) * eri_exch;
                }
            }
        }
    }
}
