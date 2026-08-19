#pragma once

#include <cctype>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>

using T2 = Eigen::Tensor<double, 2>;  // 2D tensor for Hcore
using T3 = Eigen::Tensor<double, 3>;
using T4 = Eigen::Tensor<double, 4>;  // 4D tensor for ERI
using T5 = Eigen::Tensor<double, 5>;
using T6 = Eigen::Tensor<double, 6>;

struct Atom {
    std::string symbol;
    int atomic_number{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
};

struct AtomWithBasis {
    Atom atom;
    std::string basis_set{"STO-3G"};
};

struct Molecule {
    std::vector<AtomWithBasis> atoms;

};
