#include "basis.h"

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace {

constexpr int kMaxAngularMomentum = 6;  // i functions; plenty for this project

// n!  (double avoids overflow bookkeeping; exact for n <= 22 in IEEE doubles)
double Factorial(int n) {
    static double table[21] = {0.0};
    static bool initialized = false;
    if (!initialized) {
        table[0] = 1.0;
        for (int i = 1; i <= 20; ++i) table[i] = table[i - 1] * i;
        initialized = true;
    }
    if (n < 0 || n > 20) {
        throw std::out_of_range("Factorial argument out of range: " + std::to_string(n));
    }
    return table[n];
}

// (k-1)!!  with the convention (-1)!! = 0!! = 1, matching libint2::Shell::df_Kminus1.
double DoubleFactorialKMinus1(int k) {
    static double table[31] = {0.0};
    static bool initialized = false;
    if (!initialized) {
        table[0] = 1.0;
        table[1] = 1.0;
        for (int i = 2; i <= 30; ++i) table[i] = table[i - 2] * (i - 1);
        initialized = true;
    }
    if (k < 0 || k > 30) {
        throw std::out_of_range("Double factorial argument out of range: " + std::to_string(k));
    }
    return table[k];
}

double Binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

// (-1)^i, replicating libint2's parity() including its behaviour for i < 0.
int Parity(int i) { return (i % 2) ? -1 : 1; }

// Cartesian components of a shell in libint2's "standard" (CCA) ordering:
//   for lx = l..0 { for ly = l-lx..0 { lz = l-lx-ly } }
struct CartesianComponent {
    int lx, ly, lz;
};

std::vector<CartesianComponent> CartesianComponents(int l) {
    std::vector<CartesianComponent> components;
    components.reserve((l + 1) * (l + 2) / 2);
    for (int lx = l; lx >= 0; --lx) {
        for (int ly = l - lx; ly >= 0; --ly) {
            components.push_back({lx, ly, l - lx - ly});
        }
    }
    return components;
}

// Normalization constant of a single primitive x^l exp(-alpha r^2) so that the
// (l,0,0) Cartesian component has unit self-overlap.
double PrimitiveNormalization(double alpha, int l) {
    const double sqrt_pi_cubed = 5.56832799683170784528481798212;  // pi^{3/2}
    const double two_alpha = 2.0 * alpha;
    const double two_alpha_to_am32 = std::pow(two_alpha, l + 1) * std::sqrt(two_alpha);
    return std::sqrt(std::pow(2.0, l) * two_alpha_to_am32 /
                     (sqrt_pi_cubed * DoubleFactorialKMinus1(2 * l)));
}

double IntegerPow(double base, int exponent) {
    double result = 1.0;
    for (int i = 0; i < exponent; ++i) result *= base;
    return result;
}

}  // namespace

int NumBasisFunctions(const BasisShells& shells) {
    int nbf = 0;
    for (const auto& shell : shells) nbf += shell.size();
    return nbf;
}

std::vector<int> ShellToBasisFunction(const BasisShells& shells) {
    std::vector<int> shell2bf;
    shell2bf.reserve(shells.size());
    int offset = 0;
    for (const auto& shell : shells) {
        shell2bf.push_back(offset);
        offset += shell.size();
    }
    return shell2bf;
}

GaussianShell MakeShell(int l, bool pure, const Vec3& origin,
                        std::vector<double> exponents,
                        std::vector<double> tabulated_coefficients) {
    if (l < 0 || l > kMaxAngularMomentum) {
        throw std::invalid_argument("MakeShell: unsupported angular momentum " + std::to_string(l));
    }
    if (exponents.empty() || exponents.size() != tabulated_coefficients.size()) {
        throw std::invalid_argument(
            "MakeShell: exponents and coefficients must be non-empty and equal in length.");
    }

    GaussianShell shell;
    shell.l = l;
    shell.pure = pure;
    shell.origin = origin;
    shell.exponents = std::move(exponents);
    shell.coefficients = std::move(tabulated_coefficients);

    // Step 1: fold the primitive normalization into the coefficients.
    for (std::size_t p = 0; p < shell.nprim(); ++p) {
        if (shell.exponents[p] <= 0.0) {
            throw std::invalid_argument("MakeShell: primitive exponents must be positive.");
        }
        shell.coefficients[p] *= PrimitiveNormalization(shell.exponents[p], l);
    }

    // Step 2: rescale the whole contraction to unit self-overlap.
    const double sqrt_pi_cubed = 5.56832799683170784528481798212;
    double norm = 0.0;
    for (std::size_t p = 0; p < shell.nprim(); ++p) {
        for (std::size_t q = 0; q <= p; ++q) {
            const double gamma = shell.exponents[p] + shell.exponents[q];
            norm += (p == q ? 1.0 : 2.0) * DoubleFactorialKMinus1(2 * l) * sqrt_pi_cubed *
                    shell.coefficients[p] * shell.coefficients[q] /
                    (std::pow(2.0, l) * std::pow(gamma, l + 1) * std::sqrt(gamma));
        }
    }
    if (!(norm > 0.0)) {
        throw std::invalid_argument("MakeShell: contraction has non-positive self-overlap.");
    }
    const double scale = 1.0 / std::sqrt(norm);
    for (auto& c : shell.coefficients) c *= scale;

    return shell;
}

double SolidHarmonicCoefficient(int l, int m, int lx, int ly, int lz) {
    if (l < 0 || l > kMaxAngularMomentum) {
        throw std::out_of_range("SolidHarmonicCoefficient: unsupported angular momentum.");
    }
    if (lx + ly + lz != l || std::abs(m) > l) return 0.0;

    const int abs_m = std::abs(m);
    if ((lx + ly - abs_m) % 2) return 0.0;

    const int j = (lx + ly - abs_m) / 2;
    if (j < 0) return 0.0;

    const int comp = (m >= 0) ? 1 : -1;
    const int i_shift = abs_m - lx;
    if (comp != Parity(std::abs(i_shift))) return 0.0;

    double pfac = std::sqrt(
        ((Factorial(2 * lx) * Factorial(2 * ly) * Factorial(2 * lz)) / Factorial(2 * l)) *
        (Factorial(l - abs_m) / Factorial(l)) * (1.0 / Factorial(l + abs_m)) *
        (1.0 / (Factorial(lx) * Factorial(ly) * Factorial(lz))));
    pfac /= static_cast<double>(1L << l);
    pfac *= (m < 0) ? Parity((i_shift - 1) / 2) : Parity(i_shift / 2);

    const int i_min = j;
    const int i_max = (l - abs_m) / 2;
    double sum = 0.0;
    for (int i = i_min; i <= i_max; ++i) {
        double pfac1 = Binomial(l, i) * Binomial(i, j);
        pfac1 *= (Parity(i) * Factorial(2 * (l - i))) / Factorial(l - abs_m - 2 * i);
        double sum1 = 0.0;
        const int k_min = std::max((lx - abs_m) / 2, 0);
        const int k_max = std::min(j, lx / 2);
        for (int k = k_min; k <= k_max; ++k) {
            if (lx - 2 * k <= abs_m) {
                sum1 += Binomial(j, k) * Binomial(abs_m, lx - 2 * k) * Parity(k);
            }
        }
        sum += pfac1 * sum1;
    }
    sum *= std::sqrt(DoubleFactorialKMinus1(2 * l) /
                     (DoubleFactorialKMinus1(2 * lx) * DoubleFactorialKMinus1(2 * ly) *
                      DoubleFactorialKMinus1(2 * lz)));

    return (m == 0) ? pfac * sum : M_SQRT2 * pfac * sum;
}

void EvaluateAOs(const BasisShells& shells, const Vec3& r, double* out) {
    int offset = 0;
    std::vector<double> cartesian;

    for (const auto& shell : shells) {
        const double dx = r[0] - shell.origin[0];
        const double dy = r[1] - shell.origin[1];
        const double dz = r[2] - shell.origin[2];
        const double r2 = dx * dx + dy * dy + dz * dz;

        double radial = 0.0;
        for (std::size_t p = 0; p < shell.nprim(); ++p) {
            const double arg = shell.exponents[p] * r2;
            // exp() underflows to 0 well before this, but skipping saves the call.
            if (arg < 700.0) radial += shell.coefficients[p] * std::exp(-arg);
        }

        const auto components = CartesianComponents(shell.l);
        cartesian.assign(components.size(), 0.0);
        for (std::size_t c = 0; c < components.size(); ++c) {
            const auto& comp = components[c];
            cartesian[c] = radial * IntegerPow(dx, comp.lx) * IntegerPow(dy, comp.ly) *
                           IntegerPow(dz, comp.lz);
        }

        if (!shell.pure) {
            for (std::size_t c = 0; c < components.size(); ++c) out[offset + c] = cartesian[c];
            offset += static_cast<int>(components.size());
        } else {
            for (int m = -shell.l; m <= shell.l; ++m) {
                double value = 0.0;
                for (std::size_t c = 0; c < components.size(); ++c) {
                    const auto& comp = components[c];
                    const double coefficient =
                        SolidHarmonicCoefficient(shell.l, m, comp.lx, comp.ly, comp.lz);
                    if (coefficient != 0.0) value += coefficient * cartesian[c];
                }
                out[offset++] = value;
            }
        }
    }
}
