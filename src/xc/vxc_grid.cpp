#include "xc/vxc_grid.h"

#include <algorithm>
#include <stdexcept>

#include "xc/lsda.h"

#ifdef USE_LIBXC
#include "xc/libxc_wrapper.h"
#endif

namespace {

constexpr Eigen::Index kChunkSize = 2048;

using EigenMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

Eigen::Map<const EigenMatrix> AsMatrix(const T2& t) {
    return Eigen::Map<const EigenMatrix>(t.data(), t.dimension(0), t.dimension(1));
}

T2 ToTensor(const EigenMatrix& m) {
    T2 out(m.rows(), m.cols());
    std::copy(m.data(), m.data() + m.size(), out.data());
    return out;
}

XcPointResult EvaluateFunctional(const xc::FunctionalSpec& spec, double rho_a, double rho_b) {
    XcPointResult result = EvaluateBuiltinFunctional(spec.builtin, rho_a, rho_b);

#ifdef USE_LIBXC
    for (const int id : spec.libxc_ids) {
        const SpinXCOutput xc = evaluate_spin_lda_point(id, SpinXCInput{rho_a, rho_b});
        result.energy_density += (rho_a + rho_b) * xc.eps_xc;
        result.v_rho_a += xc.v_rho_a;
        result.v_rho_b += xc.v_rho_b;
    }
#else
    if (!spec.libxc_ids.empty()) {
        throw std::runtime_error(
            "FunctionalSpec requests libxc functionals but this build has no libxc "
            "(rebuild with USE_LIBXC=1).");
    }
#endif

    return result;
}

}  // namespace

GridVxcEvaluator::GridVxcEvaluator(BasisShells shells, std::vector<GridPoint> points)
    : shells_(std::move(shells)), points_(std::move(points)) {
    nbf_ = NumBasisFunctions(shells_);
    if (nbf_ <= 0) throw std::invalid_argument("GridVxcEvaluator: empty basis.");

    ao_values_.resize(static_cast<Eigen::Index>(points_.size()), nbf_);
    std::vector<double> buffer(nbf_);
    for (std::size_t g = 0; g < points_.size(); ++g) {
        EvaluateAOs(shells_, Vec3{{points_[g].x, points_[g].y, points_[g].z}}, buffer.data());
        for (int mu = 0; mu < nbf_; ++mu) {
            ao_values_(static_cast<Eigen::Index>(g), mu) = buffer[mu];
        }
    }
}

std::shared_ptr<const GridVxcEvaluator> GridVxcEvaluator::FromProvider(
    const IIntegralProvider& provider, const GridSettings& settings) {
    return std::make_shared<const GridVxcEvaluator>(provider.GetBasisShells(),
                                                    BuildBeckeGrid(provider.GetAtoms(), settings));
}

void GridVxcEvaluator::operator()(const UksDensityInput& in, UksVxcOutput& out) const {
    if (in.Da.dimension(0) != nbf_ || in.Da.dimension(1) != nbf_ ||
        in.Db.dimension(0) != nbf_ || in.Db.dimension(1) != nbf_) {
        throw std::invalid_argument("GridVxcEvaluator: density matrix shape does not match basis.");
    }
    if (!in.functional.NeedsGrid()) return;

    const auto Da = AsMatrix(in.Da);
    const auto Db = AsMatrix(in.Db);

    EigenMatrix Va = EigenMatrix::Zero(nbf_, nbf_);
    EigenMatrix Vb = EigenMatrix::Zero(nbf_, nbf_);
    double Exc = 0.0;

    const Eigen::Index npts = ao_values_.rows();
    for (Eigen::Index start = 0; start < npts; start += kChunkSize) {
        const Eigen::Index n = std::min(kChunkSize, npts - start);
        const auto Phi = ao_values_.middleRows(start, n);          // n x nbf
        const EigenMatrix Ta = Phi * Da;                            // n x nbf
        const EigenMatrix Tb = Phi * Db;

        // rho_sigma(g) = sum_mu phi_mu(g) * (D phi)_mu(g)
        const Eigen::VectorXd rho_a = (Phi.array() * Ta.array()).rowwise().sum();
        const Eigen::VectorXd rho_b = (Phi.array() * Tb.array()).rowwise().sum();

        EigenMatrix Ba(n, nbf_);
        EigenMatrix Bb(n, nbf_);
        for (Eigen::Index i = 0; i < n; ++i) {
            const double w = points_[static_cast<std::size_t>(start + i)].weight;
            const double ra = std::max(0.0, rho_a(i));
            const double rb = std::max(0.0, rho_b(i));

            if (ra + rb < kDensityThreshold) {
                Ba.row(i).setZero();
                Bb.row(i).setZero();
                continue;
            }

            const XcPointResult xc = EvaluateFunctional(in.functional, ra, rb);
            Exc += w * xc.energy_density;
            Ba.row(i) = (w * xc.v_rho_a) * Phi.row(i);
            Bb.row(i) = (w * xc.v_rho_b) * Phi.row(i);
        }

        Va.noalias() += Phi.transpose() * Ba;
        Vb.noalias() += Phi.transpose() * Bb;
    }

    // Symmetrize: the quadrature makes Va symmetric up to rounding.
    Va = 0.5 * (Va + Va.transpose()).eval();
    Vb = 0.5 * (Vb + Vb.transpose()).eval();

    out.Vxca += ToTensor(Va);
    out.Vxcb += ToTensor(Vb);
    out.Exc += Exc;
    out.tr_PVxc += (Da.array() * Va.array()).sum() + (Db.array() * Vb.array()).sum();
}

double GridVxcEvaluator::IntegrateDensity(const T2& D) const {
    const auto Dm = AsMatrix(D);
    double total = 0.0;
    const Eigen::Index npts = ao_values_.rows();
    for (Eigen::Index start = 0; start < npts; start += kChunkSize) {
        const Eigen::Index n = std::min(kChunkSize, npts - start);
        const auto Phi = ao_values_.middleRows(start, n);
        const EigenMatrix T = Phi * Dm;
        const Eigen::VectorXd rho = (Phi.array() * T.array()).rowwise().sum();
        for (Eigen::Index i = 0; i < n; ++i) {
            total += points_[static_cast<std::size_t>(start + i)].weight * rho(i);
        }
    }
    return total;
}

T2 GridVxcEvaluator::GridOverlap() const {
    EigenMatrix S = EigenMatrix::Zero(nbf_, nbf_);
    const Eigen::Index npts = ao_values_.rows();
    for (Eigen::Index start = 0; start < npts; start += kChunkSize) {
        const Eigen::Index n = std::min(kChunkSize, npts - start);
        const auto Phi = ao_values_.middleRows(start, n);
        EigenMatrix B(n, nbf_);
        for (Eigen::Index i = 0; i < n; ++i) {
            B.row(i) = points_[static_cast<std::size_t>(start + i)].weight * Phi.row(i);
        }
        S.noalias() += Phi.transpose() * B;
    }
    return ToTensor(S);
}

VxcFunctor MakeGridVxcFunctor(std::shared_ptr<const GridVxcEvaluator> evaluator) {
    if (!evaluator) throw std::invalid_argument("MakeGridVxcFunctor: null evaluator.");
    return [evaluator](const UksDensityInput& in, UksVxcOutput& out) {
        (*evaluator)(in, out);
    };
}
