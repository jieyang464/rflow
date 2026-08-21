#include "effective_densities.h"

#include <stdexcept>

#include "linalg.h"

HartreeFockEffectiveDensities::HartreeFockEffectiveDensities(const T2& Da, const T2& Db,
                                                             const T2& Fa, const T2& Fb)
    : Da_(Da), Db_(Db) {
    const Eigen::Index nbf = Da.dimension(0);
    if (Db.dimension(0) != nbf || Fa.dimension(0) != nbf || Fb.dimension(0) != nbf) {
        throw std::invalid_argument(
            "HartreeFockEffectiveDensities: Da, Db, Fa and Fb must have the same shape.");
    }

    one_particle_ = Da + Db;

    // W_sigma = D_sigma F_sigma D_sigma.  Equivalent to sum_i eps_i C_i C_i^T,
    // but expressed without ever forming the orbitals.
    energy_weighted_ = MatMul(MatMul(Da, Fa), Da) + MatMul(MatMul(Db, Fb), Db);
}

void HartreeFockEffectiveDensities::TwoParticleBlock(const ShellQuartet& quartet,
                                                     double* out) const {
    // Two-electron energy in chemists' notation, with K^sigma_{mu nu} =
    // sum_{lam sig} D^sigma_{lam sig} (mu lam|nu sig):
    //
    //   E_2 = 1/2 sum_{abcd} (ab|cd) [ D_ab D_cd - Da_ac Da_bd - Db_ac Db_bd ]
    //
    // The exchange products are written symmetrized,
    //   -1/2 [ D_ac D_bd + D_ad D_bc ],
    // which leaves the energy unchanged -- relabelling c <-> d maps one onto the
    // other -- but gives Gamma the full 8-fold symmetry of (ab|cd).  That is
    // what makes ShellQuartet::degeneracy valid instead of an approximation.
    const T2& Da = Da_;
    const T2& Db = Db_;
    const T2& D = one_particle_;

    std::size_t index = 0;
    for (int f0 = 0; f0 < quartet.n[0]; ++f0) {
        const Eigen::Index a = quartet.bf[0] + f0;
        for (int f1 = 0; f1 < quartet.n[1]; ++f1) {
            const Eigen::Index b = quartet.bf[1] + f1;
            for (int f2 = 0; f2 < quartet.n[2]; ++f2) {
                const Eigen::Index c = quartet.bf[2] + f2;
                for (int f3 = 0; f3 < quartet.n[3]; ++f3) {
                    const Eigen::Index d = quartet.bf[3] + f3;

                    const double coulomb = D(a, b) * D(c, d);
                    const double exchange_a = 0.5 * (Da(a, c) * Da(b, d) + Da(a, d) * Da(b, c));
                    const double exchange_b = 0.5 * (Db(a, c) * Db(b, d) + Db(a, d) * Db(b, c));

                    out[index++] = coulomb - exchange_a - exchange_b;
                }
            }
        }
    }
}
