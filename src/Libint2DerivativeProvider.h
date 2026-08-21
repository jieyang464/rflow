#pragma once

#include <memory>
#include <string>
#include <vector>

#include "IIntegralDerivativeProvider.h"
#include "basis.h"
#include "types.h"

struct DerivativeBuildOptions {
    // Central-difference step, in bohr, for the one-electron derivatives.
    //
    // Why finite differences at all: libint2 must be compiled with one-body
    // derivative support, and the common packaged builds are not -- this one
    // reports INCLUDE_ONEBODY 0 and aborts on Engine(overlap, ..., deriv=1).
    // dS/dx and dHcore/dx are only N^2 and the one-electron integrals are cheap,
    // so 6*natoms extra builds costs far less than the analytic ERI derivative
    // pass that follows.  Accurate to roughly 1e-9 at the default step.
    //
    // The ERI derivatives are always analytic: libint2 does provide those.
    double finite_difference_step = 1.0e-4;

    // Set true if libint2 was built with one-body derivatives; then dS/dx and
    // dHcore/dx are taken analytically instead.
    bool use_analytic_onebody = false;
};

// Reports whether this libint2 build can differentiate the one-electron
// operators, so callers can pick a strategy instead of hitting an assertion.
bool Libint2HasAnalyticOneBodyDerivatives();

class Libint2DerivativeProvider final : public IIntegralDerivativeProvider {
public:
    Libint2DerivativeProvider(const Molecule& molecule,
                              std::vector<std::string> basis_by_atom = {},
                              DerivativeBuildOptions options = {});
    // Explicit-basis form, matching Libint2IntegralProvider.  The shells are
    // translated with the atom they sit on when the geometry is displaced.
    Libint2DerivativeProvider(const Molecule& molecule,
                              BasisShells explicit_basis,
                              DerivativeBuildOptions options = {});
    ~Libint2DerivativeProvider();

    Libint2DerivativeProvider(Libint2DerivativeProvider&&) noexcept;
    Libint2DerivativeProvider& operator=(Libint2DerivativeProvider&&) noexcept;

    int NumBasisFunctions() const override;
    int NumAtoms() const override;

    std::vector<T2> ComputeOverlapDerivatives() const override;
    std::vector<T2> ComputeHcoreDerivatives() const override;
    void ForEachEriDerivative(const EriDerivativeSink& sink) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
