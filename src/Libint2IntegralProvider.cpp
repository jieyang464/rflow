#include "Libint2IntegralProvider.h"

#include <libint2.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

int AtomicNumberFromSymbol(const std::string& symbol) {
    static const std::unordered_map<std::string, int> table = {
        {"H", 1},   {"He", 2},  {"Li", 3}, {"Be", 4}, {"B", 5},  {"C", 6},
        {"N", 7},   {"O", 8},   {"F", 9},  {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18},
        {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23}, {"Cr", 24},
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
    if (atom.atomic_number > 0) {
        return atom.atomic_number;
    }
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
        nuclear_charges.push_back({
            static_cast<double>(atom.atomic_number),
            {atom.x, atom.y, atom.z},
        });
    }

    return nuclear_charges;
}

std::vector<std::string> ResolveBasisNames(
    const Molecule& molecule,
    const std::vector<std::string>& basis_by_atom) {
    std::vector<std::string> basis_names;
    basis_names.reserve(molecule.atoms.size());

    if (!basis_by_atom.empty() &&
        basis_by_atom.size() != 1 &&
        basis_by_atom.size() != molecule.atoms.size()) {
        throw std::invalid_argument(
            "basis_by_atom must be empty, contain one uniform basis name, or match the atom count.");
    }

    for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
        std::string basis_name;
        if (basis_by_atom.size() == 1) {
            basis_name = basis_by_atom.front();
        } else if (basis_by_atom.size() == molecule.atoms.size()) {
            basis_name = basis_by_atom[i];
        }

        if (basis_name.empty()) {
            basis_name = molecule.atoms[i].basis_set;
        }
        if (basis_name.empty()) {
            basis_name = "STO-3G";
        }

        basis_names.push_back(std::move(basis_name));
    }

    return basis_names;
}

std::string RequireUniformBasisName(const std::vector<std::string>& basis_names) {
    if (basis_names.empty()) {
        return "STO-3G";
    }

    const std::string& basis_name = basis_names.front();
    if (!std::all_of(basis_names.begin(), basis_names.end(),
                     [&](const std::string& name) { return name == basis_name; })) {
        throw std::invalid_argument(
            "Libint2IntegralProvider currently supports one uniform basis set name per molecule.");
    }

    return basis_name;
}

}  // namespace

struct Libint2IntegralProvider::Impl {
    Molecule molecule;
    std::vector<std::string> basis_by_atom;
    IntegralBuildOptions options;

    std::vector<libint2::Atom> atoms;
    libint2::BasisSet basis;
    std::vector<std::pair<double, std::array<double, 3>>> nuclear_charges;
    std::vector<long> shell2atom;  // per-shell mapping to molecule atom index

    Impl(const Molecule& molecule_,
         std::vector<std::string> basis_by_atom_,
         IntegralBuildOptions options_)
        : molecule(molecule_),
          basis_by_atom(std::move(basis_by_atom_)),
          options(options_) {
        libint2::initialize();

        atoms = ToLibintAtoms(molecule);
        nuclear_charges = BuildNuclearCharges(atoms);

        const auto basis_names = ResolveBasisNames(molecule, basis_by_atom);
        const auto basis_name = RequireUniformBasisName(basis_names);
        basis = libint2::BasisSet(basis_name, atoms, true);
        shell2atom = basis.shell2atom(atoms);
    }
};

Libint2IntegralProvider::Libint2IntegralProvider(
    const Molecule& molecule, std::vector<std::string> basis_by_atom,
    IntegralBuildOptions options)
    : impl_(std::make_unique<Impl>(molecule, std::move(basis_by_atom),
                                   options)) {}

Libint2IntegralProvider::~Libint2IntegralProvider() = default;

Libint2IntegralProvider::Libint2IntegralProvider(const Libint2IntegralProvider& other)
    : impl_(std::make_unique<Impl>(*other.impl_)) {}

Libint2IntegralProvider& Libint2IntegralProvider::operator=(
    const Libint2IntegralProvider& other) {
    impl_ = std::make_unique<Impl>(*other.impl_);
    return *this;
}

Libint2IntegralProvider::Libint2IntegralProvider(Libint2IntegralProvider&&) noexcept = default;
Libint2IntegralProvider& Libint2IntegralProvider::operator=(
    Libint2IntegralProvider&&) noexcept = default;

void Libint2IntegralProvider::SetMolecule(const Molecule& molecule) {
    impl_->molecule = molecule;
    // Rebuild basis and nuclear charges to keep in sync.
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

void Libint2IntegralProvider::SetBasisByAtom(std::vector<std::string> basis_by_atom) {
    impl_->basis_by_atom = std::move(basis_by_atom);
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

void Libint2IntegralProvider::SetOptions(IntegralBuildOptions options) {
    impl_->options = options;
    *impl_ = Impl(impl_->molecule, impl_->basis_by_atom, impl_->options);
}

static void FillSymmetricMatrix(Eigen::Tensor<double, 2>& mat,
                                size_t i0,
                                size_t j0,
                                int n1,
                                int n2,
                                const double* buf) {
    for (int p = 0; p < n1; ++p) {
        for (int q = 0; q < n2; ++q) {
            const double val = buf[p * n2 + q];
            mat(i0 + p, j0 + q) = val;
            if (i0 + p != j0 + q) {
                mat(j0 + q, i0 + p) = val;
            }
        }
    }
}

static void AccumulateSymmetricMatrix(Eigen::Tensor<double, 2>& mat,
                                      size_t i0,
                                      size_t j0,
                                      int n1,
                                      int n2,
                                      const double* buf) {
    for (int p = 0; p < n1; ++p) {
        for (int q = 0; q < n2; ++q) {
            const double val = buf[p * n2 + q];
            mat(i0 + p, j0 + q) += val;
            if (i0 + p != j0 + q) {
                mat(j0 + q, i0 + p) += val;
            }
        }
    }
}

double Libint2IntegralProvider::ComputeERI(int i, int j, int k, int l) const {
    const auto nbf = static_cast<int>(impl_->basis.nbf());
    if (i < 0 || j < 0 || k < 0 || l < 0 ||
        i >= nbf || j >= nbf || k >= nbf || l >= nbf) {
        throw std::out_of_range("ERI index out of range.");
    }

    const T4 eri = ComputeERI();
    return eri(i, j, k, l);
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
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const auto max_nprim = basis.max_nprim();
    const auto max_l_i = static_cast<int>(basis.max_l());
    libint2::Engine overlap_engine(libint2::Operator::overlap, max_nprim, max_l_i, 0);

    T2 overlap(nbf, nbf);
    overlap.setZero();
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            auto ov_results = overlap_engine.compute(shell1, shell2);
            if (ov_results[0]) {
                FillSymmetricMatrix(overlap, i0, j0, shell1.size(), shell2.size(), ov_results[0]);
            }
        }
    }
    return overlap;
}

T2 Libint2IntegralProvider::ComputeHcore() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const auto max_nprim = basis.max_nprim();
    const auto max_l_i = static_cast<int>(basis.max_l());
    
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, max_nprim, max_l_i, 0);
    libint2::Engine nuclear_engine(libint2::Operator::nuclear, max_nprim, max_l_i, 0);
    nuclear_engine.set_params(impl_->nuclear_charges);

    T2 hcore(nbf, nbf);
    hcore.setZero();
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            
            auto kin_results = kinetic_engine.compute(shell1, shell2);
            if (kin_results[0]) {
                AccumulateSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), kin_results[0]);
            }
            auto nuc_results = nuclear_engine.compute(shell1, shell2);
            if (nuc_results[0]) {
                AccumulateSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), nuc_results[0]);
            }
        }
    }
    return hcore;
}

T4 Libint2IntegralProvider::ComputeERI() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    T4 eri(nbf, nbf, nbf, nbf);
    eri.setZero();

    libint2::Engine coulomb_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l(), 0);
    const auto& shell2bf = basis.shell2bf();
    const auto nshell = basis.size();

    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            for (size_t s3 = 0; s3 <= s1; ++s3) {
                const auto s4_max = (s1 == s3) ? s2 : s3;
                for (size_t s4 = 0; s4 <= s4_max; ++s4) {
                    auto buf = coulomb_engine.compute(basis[s1], basis[s2], basis[s3], basis[s4]);
                    if (!buf.empty() && buf[0] != nullptr) {
                        const size_t i0 = shell2bf[s1];
                        const size_t j0 = shell2bf[s2];
                        const size_t k0 = shell2bf[s3];
                        const size_t l0 = shell2bf[s4];
                        const size_t n1 = basis[s1].size();
                        const size_t n2 = basis[s2].size();
                        const size_t n3 = basis[s3].size();
                        const size_t n4 = basis[s4].size();

                        const double* ptr = buf[0];
                        for (size_t f1 = 0; f1 < n1; ++f1) {
                            const Eigen::Index bf1 = static_cast<Eigen::Index>(i0 + f1);
                            for (size_t f2 = 0; f2 < n2; ++f2) {
                                const Eigen::Index bf2 = static_cast<Eigen::Index>(j0 + f2);
                                for (size_t f3 = 0; f3 < n3; ++f3) {
                                    const Eigen::Index bf3 = static_cast<Eigen::Index>(k0 + f3);
                                    for (size_t f4 = 0; f4 < n4; ++f4) {
                                        const Eigen::Index bf4 = static_cast<Eigen::Index>(l0 + f4);
                                        double val = *ptr++;
                                        eri(bf1, bf2, bf3, bf4) = val;
                                        eri(bf2, bf1, bf3, bf4) = val;
                                        eri(bf1, bf2, bf4, bf3) = val;
                                        eri(bf2, bf1, bf4, bf3) = val;
                                        eri(bf3, bf4, bf1, bf2) = val;
                                        eri(bf3, bf4, bf2, bf1) = val;
                                        eri(bf4, bf3, bf1, bf2) = val;
                                        eri(bf4, bf3, bf2, bf1) = val;
                                    }
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

