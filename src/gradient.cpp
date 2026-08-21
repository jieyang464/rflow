#include "gradient.h"

#include <stdexcept>
#include <vector>

#include "geometry.h"

namespace {

double TraceProduct(const T2& A, const T2& B) {
    const Eigen::Index n = A.dimension(0);
    double trace = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        for (Eigen::Index j = 0; j < n; ++j) trace += A(i, j) * B(j, i);
    }
    return trace;
}

std::vector<Vec3> ZeroGradient(std::size_t natoms) {
    return std::vector<Vec3>(natoms, Vec3{{0.0, 0.0, 0.0}});
}

}  // namespace

GradientTerms AssembleGradientTerms(const IIntegralDerivativeProvider& derivatives,
                                    const IEffectiveDensities& densities,
                                    const std::vector<Atom>& atoms) {
    const std::size_t natoms = atoms.size();
    if (static_cast<int>(natoms) != derivatives.NumAtoms()) {
        throw std::invalid_argument("AssembleGradient: atom count does not match the provider.");
    }

    GradientTerms terms;
    terms.hcore = ZeroGradient(natoms);
    terms.overlap = ZeroGradient(natoms);
    terms.two_electron = ZeroGradient(natoms);

    // --- one-electron terms -------------------------------------------------
    {
        const std::vector<T2> d_hcore = derivatives.ComputeHcoreDerivatives();
        const std::vector<T2> d_overlap = derivatives.ComputeOverlapDerivatives();
        if (d_hcore.size() != 3 * natoms || d_overlap.size() != 3 * natoms) {
            throw std::runtime_error(
                "AssembleGradient: expected 3*natoms one-electron derivative matrices.");
        }

        const T2& D = densities.OneParticle();
        const T2& W = densities.EnergyWeighted();
        for (std::size_t a = 0; a < natoms; ++a) {
            for (int axis = 0; axis < 3; ++axis) {
                const std::size_t k = 3 * a + axis;
                terms.hcore[a][axis] = TraceProduct(D, d_hcore[k]);
                terms.overlap[a][axis] = -TraceProduct(W, d_overlap[k]);
            }
        }
    }

    // --- two-electron term --------------------------------------------------
    // Gamma and dERI/dx are both N^4 and neither is stored: for each shell
    // quartet the derivative buffers arrive, the matching Gamma block is built
    // from the N^2 densities, and the contraction happens immediately.
    {
        std::vector<double> gamma;
        derivatives.ForEachEriDerivative(
            [&](const ShellQuartet& quartet, const EriDerivativeBlock& block) {
                const int block_size = quartet.size();
                gamma.resize(static_cast<std::size_t>(block_size));
                densities.TwoParticleBlock(quartet, gamma.data());

                const double prefactor = 0.5 * quartet.degeneracy;
                for (int centre = 0; centre < 4; ++centre) {
                    const int atom = quartet.atom[centre];
                    for (int axis = 0; axis < 3; ++axis) {
                        const double* buf = block.At(centre, axis);
                        if (buf == nullptr) continue;

                        double accumulated = 0.0;
                        for (int i = 0; i < block_size; ++i) accumulated += gamma[i] * buf[i];
                        terms.two_electron[atom][axis] += prefactor * accumulated;
                    }
                }
            });
    }

    // --- nuclear repulsion --------------------------------------------------
    terms.nuclear = NuclearRepulsionGradient(atoms);

    terms.total = ZeroGradient(natoms);
    for (std::size_t a = 0; a < natoms; ++a) {
        for (int axis = 0; axis < 3; ++axis) {
            terms.total[a][axis] = terms.hcore[a][axis] + terms.overlap[a][axis] +
                                   terms.two_electron[a][axis] + terms.nuclear[a][axis];
        }
    }
    return terms;
}

std::vector<Vec3> AssembleGradient(const IIntegralDerivativeProvider& derivatives,
                                   const IEffectiveDensities& densities,
                                   const std::vector<Atom>& atoms) {
    return AssembleGradientTerms(derivatives, densities, atoms).total;
}
