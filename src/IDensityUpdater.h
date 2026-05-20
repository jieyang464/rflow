#pragma once

#include "types.h"

class IDensityUpdater {
public:
    virtual ~IDensityUpdater() = default;
    virtual void UpdateDensity(const T2& F, const T2& S, T2& D) const = 0;
};

class CommutatorDensityUpdater : public IDensityUpdater {
public:
    CommutatorDensityUpdater(int bch_order = 4, double step = 1.0);

    void UpdateDensity(const T2& F, const T2& S, T2& D) const override;

private:
    int bch_order_{4};
    double step_{1.0};
};

class FockDiagonalizationDensityUpdater : public IDensityUpdater {
public:
    void UpdateDensity(const T2& F, const T2& S, T2& D) const override;
};
