// basis.h — Provider-agnostic description of a contracted Gaussian basis.
//
// Both integral providers (the hard-coded Szabo HeH+ table and libint2) publish
// their basis through this type so that anything that needs the *value* of a
// basis function in real space — most importantly the DFT quadrature grid — can
// work without knowing which provider produced it.
//
// Normalization convention (identical to libint2::Shell after renorm()):
//   the stored `coefficients` multiply *normalization-free* primitives, i.e.
//
//     chi_{lx,ly,lz}(r) = x^lx y^ly z^lz * sum_p coefficients[p] * exp(-exponents[p] r^2)
//
//   with x,y,z measured from `origin`, and the coefficients scaled so that the
//   (l,0,0) Cartesian component has unit self-overlap.  Solid-harmonic ("pure")
//   shells are linear combinations of those Cartesians with the coefficients
//   returned by SolidHarmonicCoefficient(), which reproduces libint2's
//   convention exactly (IJQC 54, 83 (1995), eqn 15).
#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "types.h"

using Vec3 = std::array<double, 3>;

struct GaussianShell {
    int l{0};
    bool pure{false};
    Vec3 origin{{0.0, 0.0, 0.0}};
    std::vector<double> exponents;
    std::vector<double> coefficients;  // normalization-free convention, see above

    int cartesian_size() const { return (l + 1) * (l + 2) / 2; }
    int size() const { return pure ? 2 * l + 1 : cartesian_size(); }
    std::size_t nprim() const { return exponents.size(); }
};

using BasisShells = std::vector<GaussianShell>;

// Total number of basis functions spanned by `shells`.
int NumBasisFunctions(const BasisShells& shells);

// Index of the first basis function of each shell (length shells.size()).
std::vector<int> ShellToBasisFunction(const BasisShells& shells);

// Build a shell from contraction coefficients as they are *tabulated* in basis
// set libraries (i.e. referred to unit-normalized primitives) and convert them
// to the normalization-free convention documented above.
GaussianShell MakeShell(int l, bool pure, const Vec3& origin,
                        std::vector<double> exponents,
                        std::vector<double> tabulated_coefficients);

// Coefficient of the Cartesian monomial (lx,ly,lz) in the real solid harmonic
// Y_{l,m}. `m` runs -l..l (libint2 "standard" solid-harmonic ordering).
double SolidHarmonicCoefficient(int l, int m, int lx, int ly, int lz);

// Evaluate every basis function of `shells` at `r`.
// `out` must point to at least NumBasisFunctions(shells) doubles.
void EvaluateAOs(const BasisShells& shells, const Vec3& r, double* out);
