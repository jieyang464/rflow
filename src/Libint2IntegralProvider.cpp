#include "Libint2IntegralProvider.h"

#include <libint2.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

// --- libint2 global state -------------------------------------------------
// libint2::initialize() must be called before any engine is used and
// libint2::finalize() after the last one is destroyed.  Reference counting lets
// several providers coexist without either double-initializing or tearing the
// library down while another provider is still alive.
class LibintSession {
public:
    LibintSession() {
        std::lock_guard<std::mutex> lock(Mutex());
        if (RefCount()++ == 0) {
            ConfigureDataPath();
            libint2::initialize();
        }
    }
    ~LibintSession() {
        std::lock_guard<std::mutex> lock(Mutex());
        if (--RefCount() == 0) {
            libint2::finalize();
        }
    }
    LibintSession(const LibintSession&) : LibintSession() {}
    LibintSession& operator=(const LibintSession&) { return *this; }

private:
    static std::mutex& Mutex() {
        static std::mutex m;
        return m;
    }
    static int& RefCount() {
        static int count = 0;
        return count;
    }
    // libint2 looks up basis-set files through LIBINT_DATA_PATH.  Point it at
    // the copy vendored by install_basis_sets.sh unless the user set it.
    static void ConfigureDataPath() {
#ifdef SCFCXX_BASIS_PATH
        if (std::getenv("LIBINT_DATA_PATH") == nullptr) {
            setenv("LIBINT_DATA_PATH", SCFCXX_BASIS_PATH, /*overwrite=*/0);
        }
#endif
    }
};

int AtomicNumberFromSymbol(const std::string& symbol) {
    static const std::unordered_map<std::string, int> table = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},
        {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17}, {"Ar", 18},
        {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
        {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
    };

    const auto it = table.find(symbol);
    if (it == table.end()) {
        throw std::invalid_argument("Unknown atomic symbol for libint2 conversion: " + symbol);
    }
    return it->second;
}

int ResolveAtomicNumber(const Atom& atom) {
    if (atom.atomic_number > 0) return atom.atomic_number;
    if (atom.symbol.empty()) {
        throw std::invalid_argument(
            "Atom must provide either atomic_number or symbol for libint2 conversion.");
    }
    return AtomicNumberFromSymbol(atom.symbol);
}

std::vector<libint2::Atom> ToLibintAtoms(const Molecule& molecule) {
    std::vector<libint2::Atom> atoms;
    atoms.reserve(molecule.atoms.size());
    for (const auto& atom_with_basis : molecule.atoms) {
        const auto& atom = atom_with_basis.atom;
        libint2::Atom libint_atom;
        libint_atom.atomic_number = ResolveAtomicNumber(atom);
        libint_atom.x = atom.x;
        libint_atom.y = atom.y;
        libint_atom.z = atom.z;
        atoms.push_back(libint_atom);
    }
    return atoms;
}

std::vector<std::pair<double, std::array<double, 3>>> BuildNuclearCharges(
    const std::vector<libint2::Atom>& atoms) {
    std::vector<std::pair<double, std::array<double, 3>>> nuclear_charges;
    nuclear_charges.reserve(atoms.size());
    for (const auto& atom : atoms) {
        nuclear_charges.push_back(
            {static_cast<double>(atom.atomic_number), {atom.x, atom.y, atom.z}});
    }
    return nuclear_charges;
}

std::vector<std::string> ResolveBasisNames(const Molecule& molecule,
                                           const std::vector<std::string>& basis_by_atom) {
    if (!basis_by_atom.empty() && basis_by_atom.size() != 1 &&
        basis_by_atom.size() != molecule.atoms.size()) {
        throw std::invalid_argument(
            "basis_by_atom must be empty, contain one uniform basis name, or match the atom count.");
    }

    std::vector<std::string> basis_names;
    basis_names.reserve(molecule.atoms.size());
    for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
        std::string basis_name;
        if (basis_by_atom.size() == 1) {
            basis_name = basis_by_atom.front();
        } else if (basis_by_atom.size() == molecule.atoms.size()) {
            basis_name = basis_by_atom[i];
        }
        if (basis_name.empty()) basis_name = molecule.atoms[i].basis_set;
        if (basis_name.empty()) basis_name = "STO-3G";
        basis_names.push_back(std::move(basis_name));
    }
    return basis_names;
}

// Convert this project's provider-agnostic shell description into libint2's.
// GaussianShell already stores coefficients in libint2's normalization-free
// convention, so the numbers pass through untouched; the only reason to bypass
// libint2::Shell's constructor (which would renormalize them a second time) is
// to set the coefficients directly.
libint2::Shell ToLibintShell(const GaussianShell& shell) {
    libint2::svector<double> alpha(shell.exponents.begin(), shell.exponents.end());
    libint2::svector<double> coeff(shell.coefficients.begin(), shell.coefficients.end());

    // Build with unit coefficients, then overwrite: libint2::Shell's constructor
    // calls renorm(), which must not be applied to already-normalized values.
    libint2::Shell out{alpha,
                       {{shell.l, shell.pure, libint2::svector<double>(alpha.size(), 1.0)}},
                       {{shell.origin[0], shell.origin[1], shell.origin[2]}}};
    out.contr[0].coeff = coeff;

    // renorm() also fills max_ln_coeff, which the screening machinery reads.
    out.max_ln_coeff.resize(alpha.size());
    for (std::size_t p = 0; p < alpha.size(); ++p) {
        out.max_ln_coeff[p] = std::log(std::abs(coeff[p]));
    }
    return out;
}

// Write one shell-pair block of a symmetric matrix.  `mirror` must be true only
// when the two shells differ: for a diagonal block (s1 == s2) libint2 already
// returns the full square, so mirroring it would double the off-diagonal
// elements inside that block.
void AccumulateBlock(T2& mat, std::size_t i0, std::size_t j0, int n1, int n2,
                     const double* buf, bool mirror) {
    for (int p = 0; p < n1; ++p) {
        for (int q = 0; q < n2; ++q) {
            const double value = buf[p * n2 + q];
            mat(static_cast<Eigen::Index>(i0 + p), static_cast<Eigen::Index>(j0 + q)) += value;
            if (mirror) {
                mat(static_cast<Eigen::Index>(j0 + q), static_cast<Eigen::Index>(i0 + p)) += value;
            }
        }
    }
}

}  // namespace

struct Libint2IntegralProvider::Impl {
    LibintSession session;  // must outlive every engine created below

    Molecule molecule;
    std::vector<std::string> basis_by_atom;
    IntegralBuildOptions options;

    std::vector<libint2::Atom> atoms;
    std::vector<std::string> basis_names;      // resolved, one per atom
    std::vector<libint2::Shell> shells;
    std::vector<std::size_t> shell2atom;       // per-shell index into `atoms`
    std::vector<std::size_t> shell2bf;         // per-shell first basis function
    std::vector<std::pair<double, std::array<double, 3>>> nuclear_charges;
    int nbf{0};
    std::size_t max_nprim{0};
    int max_l{0};

    mutable std::unique_ptr<T4> eri_cache;

    // Construct on an explicitly supplied basis.  shell2atom is resolved by
    // matching each shell's origin against the atomic centres.
    Impl(const Molecule& molecule_, const BasisShells& explicit_basis,
         IntegralBuildOptions options_)
        : molecule(molecule_), options(options_) {
        atoms = ToLibintAtoms(molecule);
        nuclear_charges = BuildNuclearCharges(atoms);
        basis_names.assign(molecule.atoms.size(), "<explicit>");

        for (const auto& shell : explicit_basis) {
            shells.push_back(ToLibintShell(shell));

            std::size_t owner = atoms.size();
            for (std::size_t a = 0; a < atoms.size(); ++a) {
                const double dx = shell.origin[0] - atoms[a].x;
                const double dy = shell.origin[1] - atoms[a].y;
                const double dz = shell.origin[2] - atoms[a].z;
                if (dx * dx + dy * dy + dz * dz < 1e-16) { owner = a; break; }
            }
            if (owner == atoms.size()) {
                throw std::invalid_argument(
                    "Libint2IntegralProvider: an explicit shell is not centred on any atom.");
            }
            shell2atom.push_back(owner);
        }

        Finalize();
    }

    Impl(const Molecule& molecule_, std::vector<std::string> basis_by_atom_,
         IntegralBuildOptions options_)
        : molecule(molecule_),
          basis_by_atom(std::move(basis_by_atom_)),
          options(options_) {
        atoms = ToLibintAtoms(molecule);
        nuclear_charges = BuildNuclearCharges(atoms);
        basis_names = ResolveBasisNames(molecule, basis_by_atom);

        // Build the basis one atom at a time so that different atoms may carry
        // different basis sets.
        for (std::size_t a = 0; a < atoms.size(); ++a) {
            const std::vector<libint2::Atom> single_atom{atoms[a]};
            libint2::BasisSet atom_basis(basis_names[a], single_atom, /*throw_if_no_match=*/true);
            for (auto& shell : atom_basis) {
                shells.push_back(std::move(shell));
                shell2atom.push_back(a);
            }
        }

        Finalize();
    }

    // Derive the per-shell offsets and engine sizing shared by both constructors.
    void Finalize() {
        for (const auto& shell : shells) {
            shell2bf.push_back(static_cast<std::size_t>(nbf));
            nbf += static_cast<int>(shell.size());
            max_nprim = std::max(max_nprim, shell.nprim());
            for (const auto& contraction : shell.contr) {
                max_l = std::max(max_l, contraction.l);
            }
        }

        if (nbf == 0) {
            throw std::runtime_error(
                "Libint2IntegralProvider: the resolved basis contains no functions.");
        }
    }

    // Copying must not carry the cache's ownership; recompute lazily instead.
    Impl(const Impl& other)
        : session(other.session),
          molecule(other.molecule),
          basis_by_atom(other.basis_by_atom),
          options(other.options),
          atoms(other.atoms),
          basis_names(other.basis_names),
          shells(other.shells),
          shell2atom(other.shell2atom),
          shell2bf(other.shell2bf),
          nuclear_charges(other.nuclear_charges),
          nbf(other.nbf),
          max_nprim(other.max_nprim),
          max_l(other.max_l) {}

    Impl& operator=(const Impl& other) {
        if (this != &other) {
            Impl copy(other);
            std::swap(*this, copy);
        }
        return *this;
    }
    Impl(Impl&&) = default;
    Impl& operator=(Impl&&) = default;

    const T4& Eri() const {
        if (!eri_cache) eri_cache = std::make_unique<T4>(BuildEri());
        return *eri_cache;
    }

    T4 BuildEri() const;
};

Libint2IntegralProvider::Libint2IntegralProvider(const Molecule& molecule,
                                                 std::vector<std::string> basis_by_atom,
                                                 IntegralBuildOptions options)
    : impl_(std::make_unique<Impl>(molecule, std::move(basis_by_atom), options)) {}

Libint2IntegralProvider::Libint2IntegralProvider(const Molecule& molecule,
                                                 BasisShells explicit_basis,
                                                 IntegralBuildOptions options)
    : impl_(std::make_unique<Impl>(molecule, explicit_basis, options)) {}

Libint2IntegralProvider::~Libint2IntegralProvider() = default;

Libint2IntegralProvider::Libint2IntegralProvider(const Libint2IntegralProvider& other)
    : impl_(std::make_unique<Impl>(*other.impl_)) {}

Libint2IntegralProvider& Libint2IntegralProvider::operator=(const Libint2IntegralProvider& other) {
    if (this != &other) impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}

Libint2IntegralProvider::Libint2IntegralProvider(Libint2IntegralProvider&&) noexcept = default;
Libint2IntegralProvider& Libint2IntegralProvider::operator=(Libint2IntegralProvider&&) noexcept =
    default;

void Libint2IntegralProvider::SetMolecule(const Molecule& molecule) {
    // Build the replacement from local copies: `impl_` is about to be replaced.
    auto basis_by_atom = impl_->basis_by_atom;
    const auto options = impl_->options;
    impl_ = std::make_unique<Impl>(molecule, std::move(basis_by_atom), options);
}

void Libint2IntegralProvider::SetBasisByAtom(std::vector<std::string> basis_by_atom) {
    const Molecule molecule = impl_->molecule;
    const auto options = impl_->options;
    impl_ = std::make_unique<Impl>(molecule, std::move(basis_by_atom), options);
}

void Libint2IntegralProvider::SetOptions(IntegralBuildOptions options) {
    const Molecule molecule = impl_->molecule;
    auto basis_by_atom = impl_->basis_by_atom;
    impl_ = std::make_unique<Impl>(molecule, std::move(basis_by_atom), options);
}

int Libint2IntegralProvider::NumBasisFunctions() const { return impl_->nbf; }

std::vector<std::string> Libint2IntegralProvider::BasisNames() const { return impl_->basis_names; }

std::vector<Atom> Libint2IntegralProvider::GetAtoms() const {
    std::vector<Atom> atoms;
    atoms.reserve(impl_->molecule.atoms.size());
    for (std::size_t a = 0; a < impl_->molecule.atoms.size(); ++a) {
        Atom atom = impl_->molecule.atoms[a].atom;
        atom.atomic_number = impl_->atoms[a].atomic_number;
        atoms.push_back(std::move(atom));
    }
    return atoms;
}

BasisShells Libint2IntegralProvider::GetBasisShells() const {
    BasisShells out;
    out.reserve(impl_->shells.size());
    for (const auto& shell : impl_->shells) {
        // A libint2 shell may hold several contractions (a general contraction);
        // its basis functions are laid out contraction by contraction, so emit
        // one GaussianShell each, in the same order.
        for (const auto& contraction : shell.contr) {
            GaussianShell gs;
            gs.l = contraction.l;
            gs.pure = contraction.pure;
            gs.origin = {{shell.O[0], shell.O[1], shell.O[2]}};
            gs.exponents.assign(shell.alpha.begin(), shell.alpha.end());
            // libint2 has already folded the normalization into these, which is
            // exactly the convention GaussianShell documents.
            gs.coefficients.assign(contraction.coeff.begin(), contraction.coeff.end());
            out.push_back(std::move(gs));
        }
    }
    return out;
}

double Libint2IntegralProvider::ComputeNuclearRepulsionEnergy() const {
    double energy = 0.0;
    const auto& atoms = impl_->atoms;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            const double dx = atoms[i].x - atoms[j].x;
            const double dy = atoms[i].y - atoms[j].y;
            const double dz = atoms[i].z - atoms[j].z;
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r == 0.0) {
                throw std::runtime_error("Nuclear repulsion is undefined for coincident atoms.");
            }
            energy += static_cast<double>(atoms[i].atomic_number) *
                      static_cast<double>(atoms[j].atomic_number) / r;
        }
    }
    return energy;
}

T2 Libint2IntegralProvider::ComputeOverlap() const {
    const auto& shells = impl_->shells;
    libint2::Engine engine(libint2::Operator::overlap, impl_->max_nprim, impl_->max_l, 0);

    T2 overlap(impl_->nbf, impl_->nbf);
    overlap.setZero();

    const auto& buf = engine.results();
    for (std::size_t s1 = 0; s1 < shells.size(); ++s1) {
        for (std::size_t s2 = 0; s2 <= s1; ++s2) {
            engine.compute(shells[s1], shells[s2]);
            if (buf[0] == nullptr) continue;
            AccumulateBlock(overlap, impl_->shell2bf[s1], impl_->shell2bf[s2],
                            static_cast<int>(shells[s1].size()),
                            static_cast<int>(shells[s2].size()), buf[0], s1 != s2);
        }
    }
    return overlap;
}

T2 Libint2IntegralProvider::ComputeHcore() const {
    const auto& shells = impl_->shells;
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, impl_->max_nprim, impl_->max_l, 0);
    libint2::Engine nuclear_engine(libint2::Operator::nuclear, impl_->max_nprim, impl_->max_l, 0);
    nuclear_engine.set_params(impl_->nuclear_charges);

    T2 hcore(impl_->nbf, impl_->nbf);
    hcore.setZero();

    const auto& kinetic_buf = kinetic_engine.results();
    const auto& nuclear_buf = nuclear_engine.results();
    for (std::size_t s1 = 0; s1 < shells.size(); ++s1) {
        for (std::size_t s2 = 0; s2 <= s1; ++s2) {
            const int n1 = static_cast<int>(shells[s1].size());
            const int n2 = static_cast<int>(shells[s2].size());
            const bool mirror = (s1 != s2);

            kinetic_engine.compute(shells[s1], shells[s2]);
            if (kinetic_buf[0] != nullptr) {
                AccumulateBlock(hcore, impl_->shell2bf[s1], impl_->shell2bf[s2], n1, n2,
                                kinetic_buf[0], mirror);
            }
            nuclear_engine.compute(shells[s1], shells[s2]);
            if (nuclear_buf[0] != nullptr) {
                AccumulateBlock(hcore, impl_->shell2bf[s1], impl_->shell2bf[s2], n1, n2,
                                nuclear_buf[0], mirror);
            }
        }
    }
    return hcore;
}

T4 Libint2IntegralProvider::Impl::BuildEri() const {
    T4 eri(nbf, nbf, nbf, nbf);
    eri.setZero();

    libint2::Engine engine(libint2::Operator::coulomb, max_nprim, max_l, 0);
    const auto& buf = engine.results();

    // Canonical shell-quartet loop: s1 >= s2, s3 >= s4, (s1,s2) >= (s3,s4).
    for (std::size_t s1 = 0; s1 < shells.size(); ++s1) {
        for (std::size_t s2 = 0; s2 <= s1; ++s2) {
            for (std::size_t s3 = 0; s3 <= s1; ++s3) {
                const std::size_t s4_max = (s1 == s3) ? s2 : s3;
                for (std::size_t s4 = 0; s4 <= s4_max; ++s4) {
                    engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);
                    if (buf[0] == nullptr) continue;  // screened out by libint2

                    const std::size_t i0 = shell2bf[s1], n1 = shells[s1].size();
                    const std::size_t j0 = shell2bf[s2], n2 = shells[s2].size();
                    const std::size_t k0 = shell2bf[s3], n3 = shells[s3].size();
                    const std::size_t l0 = shell2bf[s4], n4 = shells[s4].size();

                    const double* ptr = buf[0];
                    for (std::size_t f1 = 0; f1 < n1; ++f1) {
                        const auto bf1 = static_cast<Eigen::Index>(i0 + f1);
                        for (std::size_t f2 = 0; f2 < n2; ++f2) {
                            const auto bf2 = static_cast<Eigen::Index>(j0 + f2);
                            for (std::size_t f3 = 0; f3 < n3; ++f3) {
                                const auto bf3 = static_cast<Eigen::Index>(k0 + f3);
                                for (std::size_t f4 = 0; f4 < n4; ++f4) {
                                    const auto bf4 = static_cast<Eigen::Index>(l0 + f4);
                                    const double value = *ptr++;
                                    // All eight permutations of (ij|kl) for real
                                    // orbitals; assigning (rather than adding)
                                    // keeps repeats among them harmless.
                                    eri(bf1, bf2, bf3, bf4) = value;
                                    eri(bf2, bf1, bf3, bf4) = value;
                                    eri(bf1, bf2, bf4, bf3) = value;
                                    eri(bf2, bf1, bf4, bf3) = value;
                                    eri(bf3, bf4, bf1, bf2) = value;
                                    eri(bf4, bf3, bf1, bf2) = value;
                                    eri(bf3, bf4, bf2, bf1) = value;
                                    eri(bf4, bf3, bf2, bf1) = value;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return eri;
}

T4 Libint2IntegralProvider::ComputeERI() const {
    if (!impl_->options.cache_eri) return impl_->BuildEri();
    return impl_->Eri();
}

double Libint2IntegralProvider::ComputeERI(int i, int j, int k, int l) const {
    const int nbf = impl_->nbf;
    if (i < 0 || j < 0 || k < 0 || l < 0 || i >= nbf || j >= nbf || k >= nbf || l >= nbf) {
        throw std::out_of_range("ERI index out of range.");
    }
    if (!impl_->options.cache_eri) return impl_->BuildEri()(i, j, k, l);
    return impl_->Eri()(i, j, k, l);
}
