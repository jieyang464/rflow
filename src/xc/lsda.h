// lsda.h — Local spin-density functionals evaluated at a single grid point.
#pragma once

#include "xc/xc_config.h"

struct XcPointResult {
    // Exchange-correlation energy *density* (energy per unit volume), i.e.
    // rho * eps_xc, so that E_xc = sum_g w_g * energy_density(g).
    double energy_density{0.0};
    // Functional derivatives d(rho*eps_xc)/d(rho_sigma).
    double v_rho_a{0.0};
    double v_rho_b{0.0};

    XcPointResult& operator+=(const XcPointResult& other) {
        energy_density += other.energy_density;
        v_rho_a += other.v_rho_a;
        v_rho_b += other.v_rho_b;
        return *this;
    }
};

// Dirac/Slater local spin-density exchange:
//   e_x = -(3/4) (6/pi)^{1/3} (rho_a^{4/3} + rho_b^{4/3})
XcPointResult SlaterExchange(double rho_a, double rho_b);

// Vosko-Wilk-Nusair correlation, parameterization "V" (Can. J. Phys. 58, 1200
// (1980), eqns 4.4 and Table 5) -- i.e. VWN5, the fit used by "SVWN5"/"LSDA".
XcPointResult Vwn5Correlation(double rho_a, double rho_b);

// Dispatch on a builtin functional id.
XcPointResult EvaluateBuiltinFunctional(xc::BuiltinFunctionalId id,
                                        double rho_a, double rho_b);

// Densities below this are treated as vacuum (all outputs zero).
constexpr double kDensityThreshold = 1e-14;
