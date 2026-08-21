#pragma once

#include <memory>

#include "jk_builder.h"
#include "types.h"
#include "dft/dft_helper.h"
#include "xc/xc_config.h"

struct FockBuildInput {
    const T2& Da;
    const T2& Db;
};

struct FockBuildResult {
    T2 Fa;
    T2 Fb;
    double electronic_energy{0.0};
};

class IFockBuilder {
public:
    virtual ~IFockBuilder() = default;
    virtual void build(const FockBuildInput& in, FockBuildResult& out) = 0;
};

// Unrestricted Hartree-Fock:
//   F_sigma = H + J[Da+Db] - K[D_sigma]
//   E_elec  = Tr[D H] + 1/2 Tr[D J] - 1/2 (Tr[Da Ka] + Tr[Db Kb])
class UHFBuilder : public IFockBuilder {
public:
    // `hcore` is borrowed, not copied: there is exactly one Hcore, owned by the
    // scope that owns the integral provider (main, or one iteration of a future
    // geometry-optimization loop).  It must outlive this builder.
    UHFBuilder(std::unique_ptr<IJKBuilder> jk_builder, const T2& hcore);

    void build(const FockBuildInput& in, FockBuildResult& out) override;

private:
    std::unique_ptr<IJKBuilder> jk_builder_;
    const T2& hcore_;

    // Backend-specific buffers sized by nbf, reused across SCF iterations.
    struct Scratch {
        T2 J;
        T2 Ka;
        T2 Kb;
    } scratch_;
};

// Unrestricted Kohn-Sham.  This class is the *only* place DFT enters the code:
// the integral providers, JK builds, SCF loop and density updaters are all
// functional-agnostic and are shared unchanged with the Hartree-Fock path.
//
//   F_sigma = H + J[Da+Db] - ax K[D_sigma] + Vxc_sigma
//   E_elec  = Tr[D H] + 1/2 Tr[D J] - ax/2 (Tr[Da Ka] + Tr[Db Kb]) + E_xc
class UKSBuilder : public IFockBuilder {
public:
    // jk_builder:  Coulomb and exact-exchange matrices.
    // hcore:       borrowed from the owning scope, as for UHFBuilder.
    // grid:        quadrature grid and cached AO values.  Built once, from the
    //              geometry and basis, and reused every iteration -- evaluating
    //              phi_mu(r) is the dominant DFT cost, so rebuilding it per cycle
    //              would dominate the run.  May be empty when the functional
    //              needs no grid (pure Hartree-Fock through this builder).
    // functional:  which XC treatment to use, and how much exact exchange to mix.
    UKSBuilder(std::unique_ptr<IJKBuilder> jk_builder,
               const T2& hcore,
               dft::DftGridContext grid,
               xc::FunctionalSpec functional);

    void build(const FockBuildInput& in, FockBuildResult& out) override;

    const xc::FunctionalSpec& functional() const { return functional_; }
    // Exchange-correlation energy from the most recent build().
    double last_exc() const { return scratch_.Exc; }
    // Tr(Da Vxc_a) + Tr(Db Vxc_b) from the most recent build().  Not used by the
    // energy expression above; kept as a diagnostic, because comparing it with a
    // finite difference of Exc is how you check that Vxc really is the
    // functional derivative of Exc.
    double last_tr_pvxc() const { return scratch_.tr_PVxc; }

private:
    std::unique_ptr<IJKBuilder> jk_builder_;
    const T2& hcore_;
    dft::DftGridContext grid_;
    xc::FunctionalSpec functional_;

    struct Scratch {
        T2 J;
        T2 Ka;
        T2 Kb;
        T2 Vxc_a;
        T2 Vxc_b;
        double Exc{0.0};
        double tr_PVxc{0.0};
    } scratch_;
};
