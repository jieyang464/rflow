#include "fock_builders.h"

#include <stdexcept>
#include <string>
#include <utility>

#include "linalg.h"

namespace {

// Tr(A B) for symmetric A and B.
double TraceProduct(const T2& A, const T2& B) {
    const Eigen::Index n = A.dimension(0);
    double trace = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) {
        for (Eigen::Index j = 0; j < n; ++j) {
            trace += A(i, j) * B(j, i);
        }
    }
    return trace;
}

// Eigen::Tensor::resize does not reallocate when the total size is unchanged,
// so this reuses the existing buffers across SCF iterations.
void ResizeAndZero(T2& m, Eigen::Index n) {
    m.resize(n, n);
    m.setZero();
}

void CheckDensityShape(const T2& D, Eigen::Index nbf, const char* name) {
    if (D.dimension(0) != nbf || D.dimension(1) != nbf) {
        throw std::invalid_argument(std::string("Fock build: ") + name +
                                    " does not match the shape of Hcore.");
    }
}

}  // namespace

UHFBuilder::UHFBuilder(std::unique_ptr<IJKBuilder> jk_builder, const T2& hcore)
    : jk_builder_(std::move(jk_builder)), hcore_(hcore) {
    if (!jk_builder_) throw std::invalid_argument("UHFBuilder: null JK builder.");
}

void UHFBuilder::build(const FockBuildInput& in, FockBuildResult& out) {
    const Eigen::Index nbf = hcore_.dimension(0);
    CheckDensityShape(in.Da, nbf, "Da");
    CheckDensityShape(in.Db, nbf, "Db");

    ResizeAndZero(scratch_.J, nbf);
    ResizeAndZero(scratch_.Ka, nbf);
    ResizeAndZero(scratch_.Kb, nbf);

    jk_builder_->build_JK(in.Da, in.Db, scratch_.J, scratch_.Ka, scratch_.Kb);

    out.Fa = hcore_ + scratch_.J - scratch_.Ka;
    out.Fb = hcore_ + scratch_.J - scratch_.Kb;

    const T2 D_total = in.Da + in.Db;
    out.electronic_energy = TraceProduct(D_total, hcore_) +
                            0.5 * TraceProduct(D_total, scratch_.J) -
                            0.5 * (TraceProduct(in.Da, scratch_.Ka) +
                                   TraceProduct(in.Db, scratch_.Kb));
}

UKSBuilder::UKSBuilder(std::unique_ptr<IJKBuilder> jk_builder,
                       const T2& hcore,
                       dft::DftGridContext grid,
                       xc::FunctionalSpec functional)
    : jk_builder_(std::move(jk_builder)),
      hcore_(hcore),
      grid_(std::move(grid)),
      functional_(std::move(functional)) {
    if (!jk_builder_) throw std::invalid_argument("UKSBuilder: null JK builder.");
    if (functional_.NeedsGrid() && grid_.empty()) {
        throw std::invalid_argument(
            "UKSBuilder: the selected functional needs a grid but the context is empty.");
    }
}

void UKSBuilder::build(const FockBuildInput& in, FockBuildResult& out) {
    const Eigen::Index nbf = hcore_.dimension(0);
    CheckDensityShape(in.Da, nbf, "Da");
    CheckDensityShape(in.Db, nbf, "Db");

    ResizeAndZero(scratch_.J, nbf);
    ResizeAndZero(scratch_.Ka, nbf);
    ResizeAndZero(scratch_.Kb, nbf);

    jk_builder_->build_JK(in.Da, in.Db, scratch_.J, scratch_.Ka, scratch_.Kb);

    const double ax = functional_.exact_exchange_fraction;

    out.Fa = hcore_ + scratch_.J - ax * scratch_.Ka;
    out.Fb = hcore_ + scratch_.J - ax * scratch_.Kb;

    ResizeAndZero(scratch_.Vxc_a, nbf);
    ResizeAndZero(scratch_.Vxc_b, nbf);
    scratch_.Exc = 0.0;
    scratch_.tr_PVxc = 0.0;

    if (functional_.NeedsGrid()) {
        // geometry -> grid -> AO values happened once, at construction.  Per
        // iteration only: density on the grid -> functional -> AO matrix.
        const dft::DftBuildResult xc =
            dft::build_uks_vxc(grid_, in.Da, in.Db, functional_);

        scratch_.Vxc_a = xc.Vxc_a;
        scratch_.Vxc_b = xc.Vxc_b;
        scratch_.Exc = xc.Exc;
        scratch_.tr_PVxc = xc.tr_PVxc;

        out.Fa += scratch_.Vxc_a;
        out.Fb += scratch_.Vxc_b;
    }

    const T2 D_total = in.Da + in.Db;
    out.electronic_energy = TraceProduct(D_total, hcore_) +
                            0.5 * TraceProduct(D_total, scratch_.J) -
                            0.5 * ax * (TraceProduct(in.Da, scratch_.Ka) +
                                        TraceProduct(in.Db, scratch_.Kb)) +
                            scratch_.Exc;
}
