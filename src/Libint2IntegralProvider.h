#pragma once

#include <memory>
#include <string>
#include <vector>

#include "IIntegralProvider.h"
#include "basis.h"
#include "types.h"

struct IntegralBuildOptions {
    // Cache the full ERI tensor on first use so that the element accessor
    // ComputeERI(i,j,k,l) is O(1).  Turning this off makes each element
    // accessor call rebuild the whole tensor and is only useful for testing.
    bool cache_eri = true;
};

// Stateful provider: geometry, basis and libint2 engine setup are expensive and
// are therefore built once and reused.
class Libint2IntegralProvider final : public IIntegralProvider {
public:
    // `basis_by_atom` may be empty (use each atom's own AtomWithBasis::basis_set),
    // hold a single name applied to every atom, or give one name per atom.
    explicit Libint2IntegralProvider(const Molecule& molecule,
                                     std::vector<std::string> basis_by_atom = {},
                                     IntegralBuildOptions options = {});
    ~Libint2IntegralProvider();

    Libint2IntegralProvider(const Libint2IntegralProvider&);
    Libint2IntegralProvider& operator=(const Libint2IntegralProvider&);
    Libint2IntegralProvider(Libint2IntegralProvider&&) noexcept;
    Libint2IntegralProvider& operator=(Libint2IntegralProvider&&) noexcept;

    void SetMolecule(const Molecule& molecule);
    void SetBasisByAtom(std::vector<std::string> basis_by_atom);
    void SetOptions(IntegralBuildOptions options);

    int NumBasisFunctions() const override;
    T2 ComputeHcore() const override;
    T2 ComputeOverlap() const override;
    T4 ComputeERI() const override;
    double ComputeERI(int i, int j, int k, int l) const override;
    double ComputeNuclearRepulsionEnergy() const override;
    std::vector<Atom> GetAtoms() const override;
    BasisShells GetBasisShells() const override;

    // Name of the basis set actually used for each atom.
    std::vector<std::string> BasisNames() const;

private:
    struct Impl;  // pimpl keeps libint2 headers out of most translation units
    std::unique_ptr<Impl> impl_;
};
