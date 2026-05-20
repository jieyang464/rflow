#pragma once

#include "types.h"
#include "IIntegralProvider.h"
#include <utility>

class IJKBuilder {
public:
    virtual ~IJKBuilder() = default;

    // Fills J, Ka, and Kb matrices based on the densities Da and Db
    virtual void build_JK(const T2& Da, const T2& Db, T2& J, T2& Ka, T2& Kb) const = 0;
};

// Direct (on-the-fly) JK builder using a stateless provider
class DirectJKBuilder : public IJKBuilder {
public:
    explicit DirectJKBuilder(const IIntegralProvider& provider) : provider_(provider) {}
    void build_JK(const T2& Da, const T2& Db, T2& J, T2& Ka, T2& Kb) const override;
private:
    const IIntegralProvider& provider_;
};

// In-Core (precomputed) JK builder taking ownership of the massive T4 tensor
class InCoreJKBuilder : public IJKBuilder {
public:
    explicit InCoreJKBuilder(T4&& eriTensor) : eri_(std::move(eriTensor)) {}
    // Or allow copying if absolutely necessary, but prefer moving exactly as you requested
    explicit InCoreJKBuilder(const T4& eriTensor) : eri_(eriTensor) {}

    void build_JK(const T2& Da, const T2& Db, T2& J, T2& Ka, T2& Kb) const override;
private:
    T4 eri_;
};
