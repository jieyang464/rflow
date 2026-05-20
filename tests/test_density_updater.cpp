#include <iostream>
#include <cmath>
#include "IDensityUpdater.h"
#include "linalg.h"
#include "types.h"

//unit tests for the density updators
int main() {
    T2 S(2, 2);
    S(0, 0) = 1.0; S(0, 1) = 0.4508;
    S(1, 0) = 0.4508; S(1, 1) = 1.0;

    T2 F(2, 2);
    F(0, 0) = -2.6527; F(0, 1) = -1.3472;
    F(1, 0) = -1.3472; F(1, 1) = -1.7318;

    // Get exact D from Fock diagonalization
    T2 D_exact(2, 2);
    // For FockDiagonalizationDensityUpdater, we need D to have Trace(DS) set up 
    // to give the correct number of electrons. We want 1 electron trace(DS)=1?
    // In Szabo Ostlund He-H, n_elec = 2, but D is a spatial orbital density matrix? No, unrestricted or restricted?
    // Let's assume restricted, 2 electrons = 1 doubly occupied spatial orbital.
    // The updater trace(DS) is n_elec. Let's initialize D_exact such that Trace(DS) = 1
    D_exact.setZero();
    D_exact(0, 0) = 1.0; // Initial guess for trace: D_exact*S has trace 1.0
    
    FockDiagonalizationDensityUpdater diag_updater;
    diag_updater.UpdateDensity(F, S, D_exact);
    
    std::cout << "Exact D from Diagonalization:\n";
    std::cout << D_exact(0,0) << " " << D_exact(0,1) << "\n";
    std::cout << D_exact(1,0) << " " << D_exact(1,1) << "\n";

    T2 R_exact = ComputeDIISError(F, D_exact, S);
    std::cout << "Error for exact D: " << MaxAbsElement(R_exact) << "\n\n";

    // Now test Commutator Density Updater
    T2 D_comm(2, 2);
    D_comm.setZero();
    D_comm(0, 0) = 1.0; // simple idempotent guess since S(0,0)=1

    CommutatorDensityUpdater comm_updater(4, 0.5); // use step=0.5 to be safe
    std::cout << "Starting Commutator Updates:\n";
    for(int i=0; i<50; ++i) {
        comm_updater.UpdateDensity(F, S, D_comm);
        T2 R = ComputeDIISError(F, D_comm, S);
        if (MaxAbsElement(R) < 1e-7) {
            std::cout << "Converged in " << i+1 << " iterations.\n";
            break;
        }
    }

    std::cout << "\nD from Commutator:\n";
    std::cout << D_comm(0,0) << " " << D_comm(0,1) << "\n";
    std::cout << D_comm(1,0) << " " << D_comm(1,1) << "\n";

    T2 diff = D_comm - D_exact;
    std::cout << "\nMax difference between D_comm and D_exact: " << MaxAbsElement(diff) << "\n";

    return 0;
}
