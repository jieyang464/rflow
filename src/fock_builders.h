#pragma once

#include "jk_builder.h"
#include "types.h"
#include "xc/vxc_evaluator.h"
#include <memory>

struct FockBuildInput {
    const T2& Hcore;
    const T2& S;
    const T2& Da;
    const T2& Db;
};

struct FockBuildResult {
    T2 Fa;
    T2 Fb;
    double energy_correction{0.0};
};

class IFockBuilder {
public:
    virtual ~IFockBuilder() = default;
    virtual FockBuildResult build(const FockBuildInput& in) = 0;
};

class UHFBuilder : public IFockBuilder {
public:
    explicit UHFBuilder(std::unique_ptr<IJKBuilder> jk_builder);
    FockBuildResult build(const FockBuildInput& in) override;

private:
    std::unique_ptr<IJKBuilder> jk_builder_;

    struct Scratch {
        T2 J;
        T2 Ka;
        T2 Kb;
    } scratch_;
};

class UKSBuilder : public IFockBuilder {
public:
    UKSBuilder(std::unique_ptr<IJKBuilder> jk_builder,
               VxcFunctor vxc_functor,
               int xc_functional_id,
               double exact_exchange_fraction = 0.0);

    FockBuildResult build(const FockBuildInput& in) override;

private:
    std::unique_ptr<IJKBuilder> jk_builder_;
    VxcFunctor vxc_functor_;
    int xc_functional_id_{0};
    double exact_exchange_fraction_{0.0};

    struct Scratch {
        T2 J;
        T2 Ka;
        T2 Kb;
        T2 Vxc_a;
        T2 Vxc_b;
        double Exc{0.0};
        double tr_PVxc{0.0};
    } scratch_;
};
