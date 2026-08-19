#pragma once

#include <memory>
#include <string>
#include <vector>

#include "IIntegralProvider.h"
#include "types.h"

// Options that typically control libint2 engine setup and what to compute.
struct IntegralBuildOptions {
    bool compute_overlap = true;
    bool compute_hcore = true;
    bool compute_eri = true;
    int derivative_order = 0;  // 0 = energies, 1 = gradients
    bool use_schwarz_screening = true;
};

// Stateful provider: realistic for libint2 because geometry, basis, and
// screening/cache setup are expensive and reused.
class Libint2IntegralProvider final : public IIntegralProvider {
public:
    Libint2IntegralProvider(const Molecule& molecule,
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

    T2 ComputeHcore() const override;
    T2 ComputeOverlap() const override;
    T4 ComputeERI() const override;
    double ComputeERI(int i, int j, int k, int l) const override;
    double ComputeNuclearRepulsionEnergy() const override;

private:
    struct Impl;  // pimpl keeps libint2 headers out of most translation units
    std::unique_ptr<Impl> impl_;
};
