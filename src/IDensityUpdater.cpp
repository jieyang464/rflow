#include "IDensityUpdater.h"
#include "linalg.h"

#include <stdexcept>

CommutatorDensityUpdater::CommutatorDensityUpdater(int bch_order, double step)
    : bch_order_(bch_order), step_(step) {}

void CommutatorDensityUpdater::UpdateDensity(const T2& F, const T2& S, T2& D) const {
    // Calculate R = FDS - SDF
    T2 R = ComputeDIISError(F, D, S);
    
    // Apply step size
    R = R * step_;

    // Update D using S-metric BCH expansion: D_new = e^{-RS} D e^SR
    D = BCHTransform(R, D, S, bch_order_);

    // Purifier: McWeeny purification D = 3 D S D - 2 D S D S D
    // To ensure idempotency (D S D = D)
    T2 D_S_D = MatMul(MatMul(D, S), D);
    D = D_S_D * 3.0 - MatMul(MatMul(D_S_D, S), D) * 2.0;
}

FockDiagonalizationDensityUpdater::FockDiagonalizationDensityUpdater(double damping)
    : damping_(damping) {
    if (!(damping > 0.0) || damping > 1.0) {
        throw std::invalid_argument(
            "FockDiagonalizationDensityUpdater: damping must lie in (0, 1].");
    }
}

void FockDiagonalizationDensityUpdater::UpdateDensity(const T2& F, const T2& S, T2& D) const {
    // For this, we need the number of electrons. Since D is updated in place,
    // this updater usually expects n_elec or we compute it from Trace(DS)?
    double n_elec_d = Trace(MatMul(D, S));
    int n_elec = static_cast<int>(n_elec_d + 0.5); // Round to nearest integer

    DiagResult diag = DiagonalizeInSMetric(F, S);
    
    const Eigen::Index n = diag.C.dimension(0);
    T2 C_occ(n, n_elec);
    C_occ.setZero();
    for (Eigen::Index i = 0; i < n; ++i) {
        for (int m = 0; m < n_elec; ++m) {
            C_occ(i, m) = diag.C(i, m);
        }
    }

    const Eigen::array<Eigen::IndexPair<int>, 1> cols = {Eigen::IndexPair<int>(1, 1)};
    const T2 D_new = C_occ.contract(C_occ, cols);

    // Both D and D_new have trace n_elec in the S metric, so any convex mixture
    // does too -- damping cannot lose electrons.
    D = (damping_ == 1.0) ? D_new : (D_new * damping_ + D * (1.0 - damping_));
}
