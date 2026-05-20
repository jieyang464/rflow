#include "SzaboHeHIntegral.h"
#include <iostream>

T2 SzaboHeHIntegralProvider::ComputeHcore() const {
    T2 hcore;
    hcore.resize(2, 2);
    hcore.setZero();
    hcore(0, 0) = -2.652744703;
    hcore(0, 1) = -1.347205024;
    hcore(1, 0) = -1.347205024;
    hcore(1, 1) = -1.731828436;
    return hcore;
}

T2 SzaboHeHIntegralProvider::ComputeOverlap() const {
    T2 overlap;
    overlap.resize(2, 2);
    overlap.setZero();
    overlap(0, 0) = 1.0;
    overlap(0, 1) = 0.4508;
    overlap(1, 0) = 0.4508;
    overlap(1, 1) = 1.0;
    return overlap;
}

T4 SzaboHeHIntegralProvider::ComputeERI() const {
    T4 eri;
    eri.resize(2, 2, 2, 2);
    eri.setZero();
    eri(0, 0, 0, 0) = 1.307152;
    eri(0, 0, 0, 1) = 0.437279;
    eri(0, 0, 1, 0) = 0.437279;
    eri(0, 0, 1, 1) = 0.605703;
    eri(0, 1, 0, 0) = 0.437279;
    eri(0, 1, 0, 1) = 0.177267;
    eri(0, 1, 1, 0) = 0.177267;
    eri(0, 1, 1, 1) = 0.311795;
    eri(1, 0, 0, 0) = 0.437279;
    eri(1, 0, 0, 1) = 0.177267;
    eri(1, 0, 1, 0) = 0.177267;
    eri(1, 0, 1, 1) = 0.311795;
    eri(1, 1, 0, 0) = 0.605703;
    eri(1, 1, 0, 1) = 0.311795;
    eri(1, 1, 1, 0) = 0.311795;
    eri(1, 1, 1, 1) = 0.774608;
    return eri;
}

double SzaboHeHIntegralProvider::ComputeERI(int i, int j, int k, int l) const {
    // For SzaboHeH this feels inefficient, but fits the interface
    return ComputeERI()(i, j, k, l);
}

IntegralDerivatives SzaboHeHIntegralProvider::ComputeFirstDerivatives() const {
    IntegralDerivatives out;
    std::cout << "SzaboHeHIntegralDerivatives: returning zero derivatives for testing\n";   
    out.d_hcore.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& m : out.d_hcore) {
        m.resize(2, 2);
        m.setZero();
    }

    out.d_eri.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& t : out.d_eri) {
        t.resize(2, 2, 2, 2);
        t.setZero();
    }

    out.d_overlap.resize(6);  // 3 coords/atom * 2 atoms
    for (auto& m : out.d_overlap) {
        m.resize(2, 2);
        m.setZero();
    }

    return out;
}
