#include "fock_builders.h"
#include "scf.h"
#include "xc/vxc_evaluator.h"

#include <algorithm>

UHFBuilder::UHFBuilder(std::unique_ptr<IJKBuilder> jk_builder)
    : jk_builder_(std::move(jk_builder)) {}

FockBuildResult UHFBuilder::build(const FockBuildInput& in) {
    const T2& H = in.Hcore;
    const Eigen::Index nbf = H.dimension(0);

    scratch_.J.resize(nbf, nbf);
    scratch_.Ka.resize(nbf, nbf);
    scratch_.Kb.resize(nbf, nbf);
    scratch_.J.setZero();
    scratch_.Ka.setZero();
    scratch_.Kb.setZero();

    jk_builder_->build_JK(in.Da, in.Db, scratch_.J, scratch_.Ka, scratch_.Kb);

    FockBuildResult result;
    result.Fa = H + scratch_.J - scratch_.Ka;
    result.Fb = H + scratch_.J - scratch_.Kb;
    result.energy_correction = 0.0;
    return result;
}

UKSBuilder::UKSBuilder(std::unique_ptr<IJKBuilder> jk_builder,
                       VxcFunctor vxc_functor,
                       int xc_functional_id,
                       double exact_exchange_fraction)
    : jk_builder_(std::move(jk_builder)),
      vxc_functor_(std::move(vxc_functor)),
      xc_functional_id_(xc_functional_id),
      exact_exchange_fraction_(exact_exchange_fraction) {}

FockBuildResult UKSBuilder::build(const FockBuildInput& in) {
    const T2& H = in.Hcore;
    const T2& S = in.S;
    const Eigen::Index nbf = H.dimension(0);

    scratch_.J.resize(nbf, nbf);
    scratch_.Ka.resize(nbf, nbf);
    scratch_.Kb.resize(nbf, nbf);
    scratch_.J.setZero();
    scratch_.Ka.setZero();
    scratch_.Kb.setZero();

    jk_builder_->build_JK(in.Da, in.Db, scratch_.J, scratch_.Ka, scratch_.Kb);

    const double ax = exact_exchange_fraction_;
    FockBuildResult result;
    result.Fa = H + scratch_.J - ax * scratch_.Ka;
    result.Fb = H + scratch_.J - ax * scratch_.Kb;

    scratch_.Vxc_a.resize(nbf, nbf);
    scratch_.Vxc_b.resize(nbf, nbf);
    scratch_.Vxc_a.setZero();
    scratch_.Vxc_b.setZero();
    scratch_.Exc = 0.0;
    scratch_.tr_PVxc = 0.0;

    UksDensityInput din{ in.Da, in.Db, S, xc_functional_id_ };
    UksVxcOutput dout{scratch_.Vxc_a, scratch_.Vxc_b, scratch_.Exc, scratch_.tr_PVxc};

    vxc_functor_(din, dout);

    result.Fa += scratch_.Vxc_a;
    result.Fb += scratch_.Vxc_b;
    result.energy_correction = scratch_.Exc - 0.5 * scratch_.tr_PVxc;
    return result;
}
