#pragma once

#include "types.h"

class IFockBuilder;
class IDensityUpdater;

struct SCFSettings {
    int Na{0};
    int Nb{0};
    int max_iter{100};
    double energy_tol{1e-8};
    // Convergence also requires the orbital-rotation gradient ||FDS - SDF||_max
    double gradient_tol{1e-6};
    // After convergence, diagonalize the converged Fock matrix once to produce
    // canonical orbitals if anything downstream needs them.
    bool finalize_orbitals{true};
    bool verbose{false};
};

struct DensityMatrices {
    T2 Da;
    T2 Db;
};

struct FockMatrices {
    T2 Fa;
    T2 Fb;
};

struct MOCoefficients {
    T2 Ca;
    T2 Cb;
};

struct FockEigenvalues {
    T2 eps_a;
    T2 eps_b;
};

struct SCFResults {
    DensityMatrices densityMatrices;
    FockMatrices    fockMatrices;
    MOCoefficients  moCoefficients;
    FockEigenvalues fockEigenvalues;
    int             iteration{0};
    bool            converged{false};
    double          energy{0.0};             // electronic_energy + nuclear_repulsion
    double          electronic_energy{0.0};
    double          nuclear_repulsion{0.0};
    double          orbital_gradient_norm{0.0};
    double          spin_squared{0.0};
};

// Core-Hamiltonian guess.  
// pass Na = Nb + 2 to break spin symmetry by hand, if not then it will be a 
// restricted guess (Ca = Cb) and the SCF will converge to a restricted solution.
void GenerateInitialGuess(const T2& Hcore, const T2& S, int Na, int Nb,
                          SCFResults& scfResults); //currently still done by diagonalizing Hcore 

void SCFLoop(const SCFSettings& settings,
             IFockBuilder& fock_builder,
             const IDensityUpdater& updater,
             const T2& S,
             SCFResults& scfResults);

double ComputeSpinSquared(const T2& Da, const T2& Db, const T2& S); //   <S^2> = Sz(Sz+1) + Nb - Tr(Da S Db S),  Sz = (Na - Nb)/2

void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao,
                    T2& h_mo, T4& eri_mo);
