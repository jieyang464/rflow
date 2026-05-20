#pragma once
#include "types.h"

class IFockBuilder;

struct SCFSettings {
    int Na{0};
    int Nb{0};
    int max_iter{100};
    double energy_tol{1e-8};
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
    double          energy{0.0};
    double          nuclear_repulsion{0.0};
    double          spin_squared{0.0};
};

//can be used to manually break symmetry by adding one alpha electron and removing one beta electron from the initial guess
void GenerateInitialGuess(const T2& Hcore, const T2& S, int Na, int Nb,
                          SCFResults& scfResults); 

void SCFLoop(const SCFSettings& settings,
             IFockBuilder& fock_builder,
             const class IDensityUpdater& updater,
             const T2& S, const T2& Hcore,
             SCFResults& scfResults);

double ComputeSpinSquared(const T2& Da, const T2& Db, const T2& S); //<S**2>

void TransformAO2MO(const T2& C, const T2& h_ao, const T4& eri_ao,
                    T2& h_mo, T4& eri_mo);
