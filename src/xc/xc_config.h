// xc_config.h — Selection of the exchange-correlation treatment.
//
// Hartree-Fock is just the special case builtin == None, exact_exchange = 1.
#pragma once

#include <string>
#include <vector>

namespace xc {

// Functionals implemented directly in this project, so that DFT works with no
// external dependency.  Values are small integers deliberately disjoint from
// the libxc id space used by FunctionalSpec::libxc_ids.
enum BuiltinFunctionalId {
    None = 0,
    SlaterExchange = 1,     // Dirac/Slater LSDA exchange
    Vwn5Correlation = 2,    // Vosko-Wilk-Nusair fit V correlation
    Lsda = 3,               // SlaterExchange + Vwn5Correlation ("SVWN5")
};

struct FunctionalSpec {
    BuiltinFunctionalId builtin{None};

    // Optional extra local functionals evaluated through libxc; their
    // contributions are added to the builtin ones.  Ignored (with a build-time
    // warning path) unless compiled with -DUSE_LIBXC.
    std::vector<int> libxc_ids{};

    // Fraction of exact (Hartree-Fock) exchange in the Fock matrix.
    double exact_exchange_fraction{0.0};

    bool NeedsGrid() const { return builtin != None || !libxc_ids.empty(); }
};

// Pure Hartree-Fock: no XC functional, full exact exchange.
inline FunctionalSpec HartreeFockSpec() {
    FunctionalSpec spec;
    spec.builtin = None;
    spec.exact_exchange_fraction = 1.0;
    return spec;
}

// SVWN5 local spin-density approximation: no exact exchange.
inline FunctionalSpec LsdaSpec() {
    FunctionalSpec spec;
    spec.builtin = Lsda;
    spec.exact_exchange_fraction = 0.0;
    return spec;
}

// Slater exchange only (the "X-alpha"-like functional with alpha = 2/3).
inline FunctionalSpec SlaterExchangeSpec() {
    FunctionalSpec spec;
    spec.builtin = SlaterExchange;
    spec.exact_exchange_fraction = 0.0;
    return spec;
}

}  // namespace xc
