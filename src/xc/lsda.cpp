#include "xc/lsda.h"

#include <cmath>
#include <stdexcept>

namespace {

constexpr double kPi = 3.14159265358979323846;

// f(zeta) = [(1+z)^{4/3} + (1-z)^{4/3} - 2] / (2^{4/3} - 2)
double SpinInterpolation(double zeta) {
    const double denominator = std::pow(2.0, 4.0 / 3.0) - 2.0;
    return (std::pow(1.0 + zeta, 4.0 / 3.0) + std::pow(1.0 - zeta, 4.0 / 3.0) - 2.0) / denominator;
}

double SpinInterpolationDerivative(double zeta) {
    const double denominator = std::pow(2.0, 4.0 / 3.0) - 2.0;
    return (4.0 / 3.0) * (std::cbrt(1.0 + zeta) - std::cbrt(1.0 - zeta)) / denominator;
}

// One of the three VWN Pade fits; returns eps and d(eps)/d(rs).
struct VwnFit {
    double A, x0, b, c;
};

void EvaluateVwnFit(const VwnFit& fit, double rs, double& eps, double& deps_drs) {
    const double x = std::sqrt(rs);
    const double X = x * x + fit.b * x + fit.c;
    const double X0 = fit.x0 * fit.x0 + fit.b * fit.x0 + fit.c;
    const double Q = std::sqrt(4.0 * fit.c - fit.b * fit.b);
    const double atan_arg = Q / (2.0 * x + fit.b);
    const double atan_term = std::atan(atan_arg);
    const double bx0_over_X0 = fit.b * fit.x0 / X0;

    eps = fit.A * (std::log(x * x / X) + (2.0 * fit.b / Q) * atan_term -
                   bx0_over_X0 * (std::log((x - fit.x0) * (x - fit.x0) / X) +
                                  (2.0 * (fit.b + 2.0 * fit.x0) / Q) * atan_term));

    // d/dx of the bracket above.
    const double denom_atan = (2.0 * x + fit.b) * (2.0 * x + fit.b) + Q * Q;
    const double deps_dx =
        fit.A * (2.0 / x - (2.0 * x + fit.b) / X - 4.0 * fit.b / denom_atan -
                 bx0_over_X0 * (2.0 / (x - fit.x0) - (2.0 * x + fit.b) / X -
                                4.0 * (fit.b + 2.0 * fit.x0) / denom_atan));

    deps_drs = deps_dx / (2.0 * x);  // dx/drs = 1/(2 sqrt(rs))
}

}  // namespace

XcPointResult SlaterExchange(double rho_a, double rho_b) {
    XcPointResult out;
    if (rho_a < 0.0) rho_a = 0.0;
    if (rho_b < 0.0) rho_b = 0.0;
    if (rho_a + rho_b < kDensityThreshold) return out;

    // E_x[rho_a, rho_b] = (1/2)(E_x[2 rho_a] + E_x[2 rho_b]) with the
    // unpolarized Dirac functional, which collapses to the constant below.
    const double c_x = -0.75 * std::cbrt(6.0 / kPi);

    const auto pow43 = [](double r) { return r > 0.0 ? r * std::cbrt(r) : 0.0; };
    const auto pow13 = [](double r) { return r > 0.0 ? std::cbrt(r) : 0.0; };

    out.energy_density = c_x * (pow43(rho_a) + pow43(rho_b));
    out.v_rho_a = (4.0 / 3.0) * c_x * pow13(rho_a);
    out.v_rho_b = (4.0 / 3.0) * c_x * pow13(rho_b);
    return out;
}

XcPointResult Vwn5Correlation(double rho_a, double rho_b) {
    XcPointResult out;
    if (rho_a < 0.0) rho_a = 0.0;
    if (rho_b < 0.0) rho_b = 0.0;

    const double rho = rho_a + rho_b;
    if (rho < kDensityThreshold) return out;

    // VWN Table 5 ("fit V"), converted from Rydberg to Hartree (A -> A/2).
    static const VwnFit kParamagnetic{0.0310907, -0.10498, 3.72744, 12.9352};
    static const VwnFit kFerromagnetic{0.01554535, -0.32500, 7.06042, 18.0578};
    static const VwnFit kSpinStiffness{-1.0 / (6.0 * kPi * kPi), -0.0047584, 1.13107, 13.0045};

    const double rs = std::cbrt(3.0 / (4.0 * kPi * rho));
    double zeta = (rho_a - rho_b) / rho;
    if (zeta > 1.0) zeta = 1.0;
    if (zeta < -1.0) zeta = -1.0;

    double eps_p = 0.0, deps_p = 0.0;
    double eps_f = 0.0, deps_f = 0.0;
    double alpha_c = 0.0, dalpha_c = 0.0;
    EvaluateVwnFit(kParamagnetic, rs, eps_p, deps_p);
    EvaluateVwnFit(kFerromagnetic, rs, eps_f, deps_f);
    EvaluateVwnFit(kSpinStiffness, rs, alpha_c, dalpha_c);

    // f''(0) = 4 / (9 (2^{1/3} - 1))
    const double fpp0 = 4.0 / (9.0 * (std::cbrt(2.0) - 1.0));
    const double f = SpinInterpolation(zeta);
    const double df = SpinInterpolationDerivative(zeta);
    const double z3 = zeta * zeta * zeta;
    const double z4 = z3 * zeta;

    // VWN eqn (4.4): the spin stiffness carries the (1 - zeta^4) piece.
    const double eps_c = eps_p + alpha_c * (f / fpp0) * (1.0 - z4) + (eps_f - eps_p) * f * z4;

    const double deps_c_drs =
        deps_p + dalpha_c * (f / fpp0) * (1.0 - z4) + (deps_f - deps_p) * f * z4;

    const double deps_c_dzeta = alpha_c / fpp0 * (df * (1.0 - z4) - 4.0 * z3 * f) +
                                (eps_f - eps_p) * (df * z4 + 4.0 * z3 * f);

    // v_sigma = d(rho eps_c)/d(rho_sigma)
    //         = eps_c - (rs/3) d(eps_c)/d(rs) + (delta_{sigma,a} - zeta) d(eps_c)/d(zeta)
    out.energy_density = rho * eps_c;
    const double common = eps_c - (rs / 3.0) * deps_c_drs - zeta * deps_c_dzeta;
    out.v_rho_a = common + deps_c_dzeta;
    out.v_rho_b = common - deps_c_dzeta;
    return out;
}

XcPointResult EvaluateBuiltinFunctional(xc::BuiltinFunctionalId id,
                                        double rho_a, double rho_b) {
    switch (id) {
        case xc::None:
            return XcPointResult{};
        case xc::SlaterExchange:
            return SlaterExchange(rho_a, rho_b);
        case xc::Vwn5Correlation:
            return Vwn5Correlation(rho_a, rho_b);
        case xc::Lsda: {
            XcPointResult out = SlaterExchange(rho_a, rho_b);
            out += Vwn5Correlation(rho_a, rho_b);
            return out;
        }
    }
    throw std::invalid_argument("EvaluateBuiltinFunctional: unknown builtin functional id.");
}
