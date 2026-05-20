#include "Libint2IntegralProvider.h"

#if SCFCXX_ENABLE_LIBINT2_PROVIDER

#include <libint2.hpp>

#include <algorithm>
#include <stdexcept>

struct Libint2IntegralProvider::Impl {
    Molecule molecule;
    std::vector<std::string> basis_by_atom;
    IntegralBuildOptions options;

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

        // Build the libint2 atom list.
        std::vector<libint2::Atom> atoms;
        atoms.reserve(molecule.atoms.size());
        nuclear_charges.reserve(molecule.atoms.size());
        for (auto const& a : molecule.atoms) {
            const auto& libint_atom = a.atom;
            atoms.push_back(libint_atom);

            nuclear_charges.push_back({
                static_cast<double>(libint_atom.atomic_number),
                {libint_atom.x, libint_atom.y, libint_atom.z}});
        }

        // Determine a basis-set name. If per-atom basis is provided,
        // we require all entries are identical (for now).
        std::string basis_name = "";
        if (!this->basis_by_atom.empty()) {
            basis_name = this->basis_by_atom.front();
            if (!std::all_of(this->basis_by_atom.begin(), this->basis_by_atom.end(),
                             [&](const std::string& b) { return b == basis_name; })) {
                throw std::invalid_argument(
                    "Libint2IntegralProvider currently only supports a single basis set name for all atoms.");
            }
        }
        if (basis_name.empty()) {
            basis_name = "sto-3g";
        }

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

double Libint2IntegralProvider::ComputeERI(int, int, int, int) const { throw std::runtime_error("Not implemented"); }


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
                FillSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), kin_results[0]);
            }
            auto nuc_results = nuclear_engine.compute(shell1, shell2);
            if (nuc_results[0]) {
                FillSymmetricMatrix(hcore, i0, j0, shell1.size(), shell2.size(), nuc_results[0]);
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
                    const auto* buf = coulomb_engine.compute(basis[s1], basis[s2], basis[s3], basis[s4]);
                    if (buf && buf[0] != nullptr) {
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
                            const size_t bf1 = i0 + f1;
                            for (size_t f2 = 0; f2 < n2; ++f2) {
                                const size_t bf2 = j0 + f2;
                                for (size_t f3 = 0; f3 < n3; ++f3) {
                                    const size_t bf3 = k0 + f3;
                                    for (size_t f4 = 0; f4 < n4; ++f4) {
                                        const size_t bf4 = l0 + f4;
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

IntegralDerivatives Libint2IntegralProvider::ComputeFirstDerivatives() const {
    const auto& basis = impl_->basis;
    const auto nbf = static_cast<int>(basis.nbf());
    const int natoms = static_cast<int>(impl_->molecule.atoms.size());

    const int ncoords = 3 * natoms;
    IntegralDerivatives out;
    out.d_overlap.assign(ncoords, T2(nbf, nbf));
    out.d_hcore.assign(ncoords, T2(nbf, nbf));
    out.d_eri.assign(ncoords, T4(nbf, nbf, nbf, nbf));

    for (auto& m : out.d_overlap) m.setZero();
    for (auto& m : out.d_hcore) m.setZero();
    for (auto& m : out.d_eri) m.setZero();

    const auto max_nprim = basis.max_nprim();
    const auto max_l = basis.max_l();
    const int max_l_i = static_cast<int>(max_l);

    libint2::Engine overlap_engine(libint2::Operator::overlap, max_nprim,
                                   max_l_i, 1);
    libint2::Engine kinetic_engine(libint2::Operator::kinetic, max_nprim,
                                   max_l_i, 1);
    libint2::Engine nuclear_engine(
        libint2::Operator::nuclear, max_nprim, max_l_i, 1,
        std::numeric_limits<double>::epsilon(), impl_->nuclear_charges);
    libint2::Engine eri_engine(libint2::Operator::coulomb, max_nprim, max_l_i, 1);

    const auto& shell2bf = basis.shell2bf();
    const auto& shell2atom = impl_->shell2atom;
    const auto nshell = basis.size();

    auto process_derivs = [&](const double* buf, int n1, int n2, size_t i0,
                              size_t j0, int atom_index, int coord, auto& target) {
        if (atom_index < 0 || atom_index >= natoms) return;
        auto& mat = target[atom_index * 3 + coord];
        for (int p = 0; p < n1; ++p) {
            for (int q = 0; q < n2; ++q) {
                const double val = buf[p * n2 + q];
                mat(i0 + p, j0 + q) += val;
                if (i0 + p != j0 + q) {
                    mat(j0 + q, i0 + p) += val;
                }
            }
        }
    };

    // One-body integral derivatives.
    for (size_t s1 = 0; s1 < nshell; ++s1) {
        for (size_t s2 = 0; s2 <= s1; ++s2) {
            const auto& shell1 = basis[s1];
            const auto& shell2 = basis[s2];
            const size_t i0 = shell2bf[s1];
            const size_t j0 = shell2bf[s2];
            const int n1 = static_cast<int>(shell1.size());
            const int n2 = static_cast<int>(shell2.size());

            const int atom1 = static_cast<int>(shell2atom[s1]);
            const int atom2 = static_cast<int>(shell2atom[s2]);

            auto handle_one_body = [&](libint2::Engine& engine,
                                       std::vector<T2>& target) {
                auto results = engine.compute1(shell1, shell2);
                const auto nsets = results.size();
                for (size_t idx = 1; idx < nsets; ++idx) {
                    const auto deriv = idx - 1;  // 0-based derivative index
                    const int center = static_cast<int>(deriv / 3);
                    const int coord = static_cast<int>(deriv % 3);
                    const int atom = (center == 0) ? atom1 : atom2;
                    if (results[idx]) {
                        process_derivs(results[idx], n1, n2, i0, j0, atom, coord,
                                       target);
                    }
                }
            };

            handle_one_body(overlap_engine, out.d_overlap);
            handle_one_body(kinetic_engine, out.d_hcore);
            handle_one_body(nuclear_engine, out.d_hcore);
        }
    }

    // Two-electron integral derivatives.
    for (size_t a = 0; a < nshell; ++a) {
        for (size_t b = 0; b < nshell; ++b) {
            for (size_t c = 0; c < nshell; ++c) {
                for (size_t d = 0; d < nshell; ++d) {
                    const auto& shell_a = basis[a];
                    const auto& shell_b = basis[b];
                    const auto& shell_c = basis[c];
                    const auto& shell_d = basis[d];

                    const size_t a0 = shell2bf[a];
                    const size_t b0 = shell2bf[b];
                    const size_t c0 = shell2bf[c];
                    const size_t d0 = shell2bf[d];
                    const int na = static_cast<int>(shell_a.size());
                    const int nb = static_cast<int>(shell_b.size());
                    const int nc = static_cast<int>(shell_c.size());
                    const int nd = static_cast<int>(shell_d.size());

                    const int atom_a = static_cast<int>(shell2atom[a]);
                    const int atom_b = static_cast<int>(shell2atom[b]);
                    const int atom_c = static_cast<int>(shell2atom[c]);
                    const int atom_d = static_cast<int>(shell2atom[d]);

                    auto eri_results =
                        eri_engine.compute2<libint2::Operator::coulomb,
                                            libint2::BraKet::xx_xx,
                                            1>(shell_a, shell_b, shell_c,
                                              shell_d);
                    const auto nsets = eri_results.size();
                    if (nsets == 0) continue;

                    const size_t expected_derivs = 3 * 4;  // 4 centers (a,b,c,d), each has x/y/z
                    const size_t offset =
                        (nsets == expected_derivs + 1 && eri_results[0]) ? 1 : 0;

                    for (size_t idx = offset; idx < nsets; ++idx) {
                        const auto deriv = idx - offset;
                        const int center = static_cast<int>(deriv / 3);
                        const int coord = static_cast<int>(deriv % 3);
                        int atom = -1;
                        switch (center) {
                        case 0:
                            atom = atom_a;
                            break;
                        case 1:
                            atom = atom_b;
                            break;
                        case 2:
                            atom = atom_c;
                            break;
                        case 3:
                            atom = atom_d;
                            break;
                        default:
                            break;
                        }
                        if (atom < 0 || atom >= natoms) continue;

                        const double* buf = eri_results[idx];
                        if (!buf) continue;
                        auto& tensor = out.d_eri[atom * 3 + coord];

                        for (int ia = 0; ia < na; ++ia) {
                            for (int ib = 0; ib < nb; ++ib) {
                                for (int ic = 0; ic < nc; ++ic) {
                                    for (int id = 0; id < nd; ++id) {
                                        const size_t idxx = (((static_cast<size_t>(ia) * nb + ib) * nc + ic) * nd) + id;
                                        tensor(a0 + ia, b0 + ib, c0 + ic, d0 + id) +=
                                            buf[idxx];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

#endif
