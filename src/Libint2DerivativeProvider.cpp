#include "Libint2DerivativeProvider.h"

#include <libint2.hpp>

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "Libint2IntegralProvider.h"
#include "geometry.h"

namespace {

// Build a provider at a displaced geometry, moving explicit shells with their
// atom so that the finite difference includes the basis-function motion (the
// Pulay part) and not only the operator's explicit dependence.
std::unique_ptr<Libint2IntegralProvider> MakeDisplacedProvider(
    const Molecule& molecule, const std::vector<std::string>& basis_by_atom,
    const BasisShells& explicit_basis, std::size_t atom_index, int axis, double delta) {
    Molecule moved = molecule;
    switch (axis) {
        case 0: moved.atoms[atom_index].atom.x += delta; break;
        case 1: moved.atoms[atom_index].atom.y += delta; break;
        default: moved.atoms[atom_index].atom.z += delta; break;
    }

    if (explicit_basis.empty()) {
        return std::make_unique<Libint2IntegralProvider>(moved, basis_by_atom);
    }

    BasisShells shells = explicit_basis;
    const Atom& original = molecule.atoms[atom_index].atom;
    for (auto& shell : shells) {
        const double dx = shell.origin[0] - original.x;
        const double dy = shell.origin[1] - original.y;
        const double dz = shell.origin[2] - original.z;
        if (dx * dx + dy * dy + dz * dz < 1e-16) shell.origin[axis] += delta;
    }
    return std::make_unique<Libint2IntegralProvider>(moved, shells);
}

}  // namespace

bool Libint2HasAnalyticOneBodyDerivatives() {
#ifdef INCLUDE_ONEBODY
    return INCLUDE_ONEBODY > 0;
#else
    return false;
#endif
}

struct Libint2DerivativeProvider::Impl {
    Molecule molecule;
    std::vector<std::string> basis_by_atom;
    BasisShells explicit_basis;  // empty unless constructed from explicit shells
    DerivativeBuildOptions options;

    std::unique_ptr<Libint2IntegralProvider> reference;  // at the undisplaced geometry
    std::vector<libint2::Atom> atoms;
    std::vector<libint2::Shell> shells;
    std::vector<int> shell2atom;
    std::vector<int> shell2bf;
    int nbf{0};
    std::size_t max_nprim{0};
    int max_l{0};

    void Build();

    // Central difference of a one-electron matrix with respect to one coordinate.
    T2 OneElectronDerivative(std::size_t atom_index, int axis,
                             T2 (IIntegralProvider::*compute)() const) const {
        const double h = options.finite_difference_step;
        const auto plus = MakeDisplacedProvider(molecule, basis_by_atom, explicit_basis,
                                                atom_index, axis, +h);
        const auto minus = MakeDisplacedProvider(molecule, basis_by_atom, explicit_basis,
                                                 atom_index, axis, -h);
        const T2 forward = ((*plus).*compute)();
        const T2 backward = ((*minus).*compute)();
        return (forward - backward) * (0.5 / h);
    }
};

void Libint2DerivativeProvider::Impl::Build() {
    if (options.use_analytic_onebody && !Libint2HasAnalyticOneBodyDerivatives()) {
        throw std::runtime_error(
            "Libint2DerivativeProvider: analytic one-body derivatives were requested but this "
            "libint2 build has INCLUDE_ONEBODY=0. Rebuild libint2 with one-body derivative "
            "support, or leave use_analytic_onebody false to use central differences.");
    }
    if (options.finite_difference_step <= 0.0) {
        throw std::invalid_argument("Libint2DerivativeProvider: step must be positive.");
    }

    reference = explicit_basis.empty()
                    ? std::make_unique<Libint2IntegralProvider>(molecule, basis_by_atom)
                    : std::make_unique<Libint2IntegralProvider>(molecule, explicit_basis);

    // Rebuild the libint2 shell list here so the ERI derivative loop can drive
    // shells directly; the integral provider keeps its own copy behind its pimpl.
    const BasisShells basis = reference->GetBasisShells();
    const std::vector<Atom> geometry = reference->GetAtoms();

    atoms.clear();
    for (const auto& a : geometry) {
        atoms.push_back(libint2::Atom{a.atomic_number, a.x, a.y, a.z});
    }

    shells.clear();
    shell2atom.clear();
    shell2bf.clear();
    nbf = 0;
    max_nprim = 0;
    max_l = 0;

    for (const auto& gs : basis) {
        libint2::svector<double> alpha(gs.exponents.begin(), gs.exponents.end());
        libint2::Shell shell{alpha,
                             {{gs.l, gs.pure, libint2::svector<double>(alpha.size(), 1.0)}},
                             {{gs.origin[0], gs.origin[1], gs.origin[2]}}};
        shell.contr[0].coeff.assign(gs.coefficients.begin(), gs.coefficients.end());
        shell.max_ln_coeff.resize(alpha.size());
        for (std::size_t p = 0; p < alpha.size(); ++p) {
            shell.max_ln_coeff[p] = std::log(std::abs(gs.coefficients[p]));
        }

        int owner = -1;
        for (std::size_t a = 0; a < geometry.size(); ++a) {
            const double dx = gs.origin[0] - geometry[a].x;
            const double dy = gs.origin[1] - geometry[a].y;
            const double dz = gs.origin[2] - geometry[a].z;
            if (dx * dx + dy * dy + dz * dz < 1e-16) { owner = static_cast<int>(a); break; }
        }
        if (owner < 0) {
            throw std::runtime_error(
                "Libint2DerivativeProvider: a shell is not centred on any atom.");
        }

        shell2atom.push_back(owner);
        shell2bf.push_back(nbf);
        nbf += static_cast<int>(shell.size());
        max_nprim = std::max(max_nprim, shell.nprim());
        max_l = std::max(max_l, gs.l);
        shells.push_back(std::move(shell));
    }
}

Libint2DerivativeProvider::Libint2DerivativeProvider(const Molecule& molecule,
                                                     std::vector<std::string> basis_by_atom,
                                                     DerivativeBuildOptions options)
    : impl_(std::make_unique<Impl>()) {
    impl_->molecule = molecule;
    impl_->basis_by_atom = std::move(basis_by_atom);
    impl_->options = options;
    impl_->Build();
}

Libint2DerivativeProvider::Libint2DerivativeProvider(const Molecule& molecule,
                                                     BasisShells explicit_basis,
                                                     DerivativeBuildOptions options)
    : impl_(std::make_unique<Impl>()) {
    impl_->molecule = molecule;
    impl_->explicit_basis = std::move(explicit_basis);
    impl_->options = options;
    impl_->Build();
}

Libint2DerivativeProvider::~Libint2DerivativeProvider() = default;
Libint2DerivativeProvider::Libint2DerivativeProvider(Libint2DerivativeProvider&&) noexcept =
    default;
Libint2DerivativeProvider& Libint2DerivativeProvider::operator=(
    Libint2DerivativeProvider&&) noexcept = default;

int Libint2DerivativeProvider::NumBasisFunctions() const { return impl_->nbf; }
int Libint2DerivativeProvider::NumAtoms() const {
    return static_cast<int>(impl_->atoms.size());
}

std::vector<T2> Libint2DerivativeProvider::ComputeOverlapDerivatives() const {
    std::vector<T2> derivatives;
    derivatives.reserve(impl_->atoms.size() * 3);
    for (std::size_t a = 0; a < impl_->atoms.size(); ++a) {
        for (int axis = 0; axis < 3; ++axis) {
            derivatives.push_back(
                impl_->OneElectronDerivative(a, axis, &IIntegralProvider::ComputeOverlap));
        }
    }
    return derivatives;
}

std::vector<T2> Libint2DerivativeProvider::ComputeHcoreDerivatives() const {
    std::vector<T2> derivatives;
    derivatives.reserve(impl_->atoms.size() * 3);
    for (std::size_t a = 0; a < impl_->atoms.size(); ++a) {
        for (int axis = 0; axis < 3; ++axis) {
            // Displacing an atom moves both its basis functions and its nucleus,
            // so this captures the Hellmann-Feynman and Pulay parts together.
            derivatives.push_back(
                impl_->OneElectronDerivative(a, axis, &IIntegralProvider::ComputeHcore));
        }
    }
    return derivatives;
}

void Libint2DerivativeProvider::ForEachEriDerivative(const EriDerivativeSink& sink) const {
    if (!sink) throw std::invalid_argument("ForEachEriDerivative: empty sink.");

    const auto& shells = impl_->shells;
    libint2::Engine engine(libint2::Operator::coulomb, impl_->max_nprim, impl_->max_l,
                           /*deriv_order=*/1);
    const auto& results = engine.results();

    for (std::size_t s0 = 0; s0 < shells.size(); ++s0) {
        for (std::size_t s1 = 0; s1 <= s0; ++s1) {
            for (std::size_t s2 = 0; s2 <= s0; ++s2) {
                const std::size_t s3_max = (s0 == s2) ? s1 : s2;
                for (std::size_t s3 = 0; s3 <= s3_max; ++s3) {
                    engine.compute(shells[s0], shells[s1], shells[s2], shells[s3]);

                    ShellQuartet quartet;
                    const std::size_t s[4] = {s0, s1, s2, s3};
                    for (int i = 0; i < 4; ++i) {
                        quartet.shell[i] = static_cast<int>(s[i]);
                        quartet.bf[i] = impl_->shell2bf[s[i]];
                        quartet.n[i] = static_cast<int>(shells[s[i]].size());
                        quartet.atom[i] = impl_->shell2atom[s[i]];
                    }
                    // Multiplicity of the permutations this canonical quartet
                    // stands in for: mu<->nu, lam<->sig, and (mu nu)<->(lam sig).
                    quartet.degeneracy = (s0 == s1 ? 1 : 2) * (s2 == s3 ? 1 : 2) *
                                         ((s0 == s2 && s1 == s3) ? 1 : 2);

                    EriDerivativeBlock block;
                    bool any = false;
                    for (int i = 0; i < 12; ++i) {
                        block.buffer[i] = (i < static_cast<int>(results.size())) ? results[i]
                                                                                 : nullptr;
                        if (block.buffer[i]) any = true;
                    }
                    if (!any) continue;  // screened out entirely

                    sink(quartet, block);
                }
            }
        }
    }
}
