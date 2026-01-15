#pragma once

#include <string>
#include <utility>
#include <vector>

struct Atom {
  int Z;           // atomic number
  double x, y, z;  // Cartesian coordinates in Bohr (atomic units)
  // Optional: per-atom basis label (not used in BasisSet construction here; kept for future extension)
  std::string basis_label;
};

class IMolecule {
public:
  virtual ~IMolecule() = default;
  virtual const std::vector<Atom>& atoms() const = 0;
  virtual int net_charge() const = 0;
  virtual const std::string& basis_name() const = 0;
  virtual const std::string& basis_file_path() const = 0;
  virtual int multiplicity() const = 0;
  virtual int electron_count() const = 0;
  virtual std::pair<int,int> alpha_beta_electrons() const = 0;
};

class Molecule final : public IMolecule {
public:
  Molecule() = default;

  Molecule(std::vector<Atom> atoms, int net_charge,
           std::string basis_name,
           std::string basis_file_path = "")
      : atoms_(std::move(atoms)), net_charge_(net_charge),
        basis_name_(std::move(basis_name)), basis_file_path_(std::move(basis_file_path)) {}

  const std::vector<Atom>& atoms() const override { return atoms_; }
  int net_charge() const override { return net_charge_; }
  const std::string& basis_name() const override { return basis_name_; }
  const std::string& basis_file_path() const override { return basis_file_path_; }
  int multiplicity() const override { return multiplicity_; }

  // Spin multiplicity (2S+1). Default singlet (1). Must be positive integer.
  void set_multiplicity(int mult);

  // Total number of electrons N = sum Z - net_charge
  int electron_count() const override;

  // Alpha/Beta electron counts derived from N and multiplicity: Na=(N+M-1)/2, Nb=(N-M+1)/2
  // Throws if parity or sign constraints are violated.
  std::pair<int,int> alpha_beta_electrons() const override;

  // Set a uniform basis for all atoms (common practice). This sets the global
  // basis_name_/basis_file_path_ and mirrors the name into each atom's basis_label.
  void set_uniform_basis(const std::string& name, const std::string& file_path = "");

private:
  std::vector<Atom> atoms_;
  int net_charge_ = 0;
  std::string basis_name_ = ""; // must be provided by caller (e.g., in tests)
  std::string basis_file_path_;

  int multiplicity_ = 1; // 2S+1 (singlet default)
};
