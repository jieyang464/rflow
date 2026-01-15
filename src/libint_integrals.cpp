#include "integral_provider.h"

#include <libint2.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

// RAII guard to ensure libint2 is initialized before use and finalized at program end
struct LibintGuard {
  LibintGuard() { libint2::initialize(); }
  ~LibintGuard() { libint2::finalize(); }
};

static std::vector<libint2::Atom> build_libint_atoms(const IMolecule& mol) {
  std::vector<libint2::Atom> ai;
  ai.reserve(mol.atoms().size());
  for (const Atom& a : mol.atoms()) {
    libint2::Atom la;
    la.atomic_number = a.Z;
    la.x = a.x;
    la.y = a.y;
    la.z = a.z;
    ai.push_back(la);
  }
  return ai;
}

// Helper to choose basis: if basis_file_path empty, check ./basis/<basis_name>.g94 relative path first
static libint2::BasisSet choose_basis(const IMolecule& mol, const std::vector<libint2::Atom>& ai) {
  // Require either an explicit file path or a non-empty basis name
  if (mol.basis_file_path().empty() && mol.basis_name().empty()) {
    throw std::runtime_error("Basis set not specified: provide basis_name or basis_file_path");
  }
  // If user provided an explicit file path, set LIBINT_DATA_PATH to its directory and use basename
  if (!mol.basis_file_path().empty()) {
    std::string path = mol.basis_file_path();
    size_t pos = path.find_last_of('/');
    std::string dir = (pos == std::string::npos) ? std::string(".") : path.substr(0, pos);
    std::string file = (pos == std::string::npos) ? path : path.substr(pos + 1);
    std::string name = file;
    // strip .g94 if present
    if (name.size() > 4 && name.substr(name.size() - 4) == ".g94") {
      name = name.substr(0, name.size() - 4);
    }
    setenv("LIBINT_DATA_PATH", dir.c_str(), 1);
    return libint2::BasisSet(name, ai, /*throw_if_no_match=*/true);
  }

  // try relative path ./basis/<name>.g94
  std::string rel = std::string("basis/") + mol.basis_name() + ".g94";
  if (std::ifstream(rel)) {
    setenv("LIBINT_DATA_PATH", "./basis", 1);
    return libint2::BasisSet(mol.basis_name(), ai, /*throw_if_no_match=*/true);
  }

  // fallback to built-in (system) data path
  return libint2::BasisSet(mol.basis_name(), ai, /*throw_if_no_match=*/false);
}

}

class LibintIntegralProvider final : public IIntegralProvider {
public:
  size_t nbf(const IMolecule& mol) const override {
    static LibintGuard guard;
    std::vector<libint2::Atom> ai = build_libint_atoms(mol);
    libint2::BasisSet basis = choose_basis(mol, ai);
    return basis.nbf();
  }

  scfcpp::Matrix overlap(const IMolecule& mol) const override {
    static LibintGuard guard;
    return scfcpp::to_matrix(compute_one_electron(mol, libint2::Operator::overlap));
  }

  scfcpp::Matrix kinetic(const IMolecule& mol) const override {
    static LibintGuard guard;
    return scfcpp::to_matrix(compute_one_electron(mol, libint2::Operator::kinetic));
  }

  scfcpp::Matrix nuclear_attraction(const IMolecule& mol) const override {
    static LibintGuard guard;
    std::vector<libint2::Atom> ai = build_libint_atoms(mol);
    libint2::BasisSet basis = choose_basis(mol, ai);
    const size_t nbasis = basis.nbf();
    scfcpp::DenseMatrix V(nbasis, nbasis);
    V.setZero();

    libint2::Engine V_engine(libint2::Operator::nuclear, basis.max_nprim(), basis.max_l());
    V_engine.set_params(libint2::make_point_charges(ai));
    const libint2::Engine::target_ptr_vec& Vbufs = V_engine.results();

    size_t row_offset = 0;
    for (size_t s1 = 0; s1 != basis.size(); ++s1) {
      const libint2::Shell& sh1 = basis[s1];
      const size_t n1 = sh1.size();
      size_t col_offset = 0;
      for (size_t s2 = 0; s2 != basis.size(); ++s2) {
        const libint2::Shell& sh2 = basis[s2];
        const size_t n2 = sh2.size();
        V_engine.compute(sh1, sh2);
        const double* buf = Vbufs[0];
        if (buf == nullptr) {
          col_offset += n2;
          continue;
        }
        for (size_t f1 = 0; f1 != n1; ++f1) {
          for (size_t f2 = 0; f2 != n2; ++f2) {
            V(row_offset + f1, col_offset + f2) = buf[f1 * n2 + f2];
          }
        }
        col_offset += n2;
      }
      row_offset += n1;
    }
    return scfcpp::to_matrix(V);
  }

  scfcpp::Matrix hcore(const IMolecule& mol) const override {
    static LibintGuard guard;
    scfcpp::DenseMatrix T = scfcpp::to_dense(kinetic(mol));
    scfcpp::DenseMatrix V = scfcpp::to_dense(nuclear_attraction(mol));
    if (T.rows() != V.rows() || T.cols() != V.cols()) {
      throw std::runtime_error("T and V dimension mismatch");
    }
    return scfcpp::to_matrix(T + V);
  }

  double eri_abs_sum(const IMolecule& mol) const override {
    static LibintGuard guard;
    std::vector<libint2::Atom> ai = build_libint_atoms(mol);
    libint2::BasisSet basis = choose_basis(mol, ai);
    libint2::Engine eri_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l());
    const libint2::Engine::target_ptr_vec& ERIbufs = eri_engine.results();

    double eri_sum = 0.0;
    for (size_t s1 = 0; s1 != basis.size(); ++s1) {
      const libint2::Shell& sh1 = basis[s1];
      for (size_t s2 = 0; s2 != basis.size(); ++s2) {
        const libint2::Shell& sh2 = basis[s2];
        for (size_t s3 = 0; s3 != basis.size(); ++s3) {
          const libint2::Shell& sh3 = basis[s3];
          for (size_t s4 = 0; s4 != basis.size(); ++s4) {
            const libint2::Shell& sh4 = basis[s4];
            eri_engine.compute(sh1, sh2, sh3, sh4);
            const double* buf = ERIbufs[0];
            if (!buf) continue;
            const size_t n1 = sh1.size(), n2 = sh2.size();
            const size_t n3 = sh3.size(), n4 = sh4.size();
            const size_t n1234 = n1 * n2 * n3 * n4;
            for (size_t i = 0; i < n1234; ++i) eri_sum += std::abs(buf[i]);
          }
        }
      }
    }
    return eri_sum;
  }

  std::vector<double> eri_tensor(const IMolecule& mol) const override {
    static LibintGuard guard;
    std::vector<libint2::Atom> ai = build_libint_atoms(mol);
    libint2::BasisSet basis = choose_basis(mol, ai);
    const size_t nbf = basis.nbf();
    std::vector<double> eri(nbf * nbf * nbf * nbf, 0.0);

    struct EriIndex {
      size_t nbf;
      size_t operator()(size_t mu, size_t nu, size_t lam, size_t sig) const {
        return ((mu * nbf + nu) * nbf + lam) * nbf + sig;
      }
    } eri_index{nbf};

    std::vector<size_t> shell_offsets(basis.size() + 1, 0);
    for (size_t s = 0; s < basis.size(); ++s) {
      shell_offsets[s + 1] = shell_offsets[s] + basis[s].size();
    }

    libint2::Engine eri_engine(libint2::Operator::coulomb, basis.max_nprim(), basis.max_l());
    const libint2::Engine::target_ptr_vec& ERIbufs = eri_engine.results();
    for (size_t s1 = 0; s1 != basis.size(); ++s1) {
      const libint2::Shell& sh1 = basis[s1];
      const size_t n1 = sh1.size();
      for (size_t s2 = 0; s2 != basis.size(); ++s2) {
        const libint2::Shell& sh2 = basis[s2];
        const size_t n2 = sh2.size();
        for (size_t s3 = 0; s3 != basis.size(); ++s3) {
          const libint2::Shell& sh3 = basis[s3];
          const size_t n3 = sh3.size();
          for (size_t s4 = 0; s4 != basis.size(); ++s4) {
            const libint2::Shell& sh4 = basis[s4];
            const size_t n4 = sh4.size();
            eri_engine.compute(sh1, sh2, sh3, sh4);
            const double* buf = ERIbufs[0];
            if (!buf) continue;
            size_t off1 = shell_offsets[s1];
            size_t off2 = shell_offsets[s2];
            size_t off3 = shell_offsets[s3];
            size_t off4 = shell_offsets[s4];
            for (size_t f1 = 0; f1 != n1; ++f1)
              for (size_t f2 = 0; f2 != n2; ++f2)
                for (size_t f3 = 0; f3 != n3; ++f3)
                  for (size_t f4 = 0; f4 != n4; ++f4) {
                    size_t mu = off1 + f1;
                    size_t nu = off2 + f2;
                    size_t lam = off3 + f3;
                    size_t sig = off4 + f4;
                    size_t idx = eri_index(mu, nu, lam, sig);
                    size_t buf_idx = ((f1 * n2 + f2) * n3 + f3) * n4 + f4;
                    eri[idx] = buf[buf_idx];
                  }
          }
        }
      }
    }
    return eri;
  }

private:
  scfcpp::DenseMatrix compute_one_electron(const IMolecule& mol, libint2::Operator op) const {
    static LibintGuard guard;
    std::vector<libint2::Atom> ai = build_libint_atoms(mol);
    libint2::BasisSet basis = choose_basis(mol, ai);
    const size_t nbasis = basis.nbf();
    scfcpp::DenseMatrix M(nbasis, nbasis);
    M.setZero();

    libint2::Engine engine(op, basis.max_nprim(), basis.max_l());
    const libint2::Engine::target_ptr_vec& bufs = engine.results();

    size_t row_offset = 0;
    for (size_t s1 = 0; s1 != basis.size(); ++s1) {
      const libint2::Shell& sh1 = basis[s1];
      const size_t n1 = sh1.size();
      size_t col_offset = 0;
      for (size_t s2 = 0; s2 != basis.size(); ++s2) {
        const libint2::Shell& sh2 = basis[s2];
        const size_t n2 = sh2.size();
        engine.compute(sh1, sh2);
        const double* buf = bufs[0];
        if (buf == nullptr) {
          col_offset += n2;
          continue;
        }
        for (size_t f1 = 0; f1 != n1; ++f1) {
          for (size_t f2 = 0; f2 != n2; ++f2) {
            M(row_offset + f1, col_offset + f2) = buf[f1 * n2 + f2];
          }
        }
        col_offset += n2;
      }
      row_offset += n1;
    }
    return M;
  }
};

std::shared_ptr<IIntegralProvider> make_libint_integral_provider() {
  return std::make_shared<LibintIntegralProvider>();
}
