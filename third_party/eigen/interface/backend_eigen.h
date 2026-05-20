#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

struct EigenBackend {

    template <typename T, size_t Rank, bool IsSparse>
    struct Storage;

    template <typename T, size_t Rank>
    struct Storage<T, Rank, false> {
        using type = Eigen::Tensor<T, Rank>;
    };

    template <typename T>
    struct Storage<T, 2, false> {
        using type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    };

    template <typename T>
    struct Storage<T, 2, true> {
        using type = Eigen::SparseMatrix<T, Eigen::RowMajor>; 
    };

    template <typename T>
    struct Storage<T, 1, false> {
        using type = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    };

    template <typename T>
    static void diagonalize_dense(const typename Storage<T, 2, false>::type& mat, 
                                  typename Storage<T, 1, false>::type& evals, 
                                  typename Storage<T, 2, false>::type& evecs) 
    {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solver(mat);
        evals = solver.eigenvalues();
        evecs = solver.eigenvectors();
    }

    template <typename T>
    static void diagonalize_sparse(const typename Storage<T, 2, true>::type& mat, 
                                   typename Storage<T, 1, false>::type& evals, 
                                   typename Storage<T, 2, false>::type& evecs) 
    {
        std::cout << "[Backend] Warning: Sparse diagonalization requires Spectra.\n";
    }

    template <typename T, size_t RankA, size_t RankB>
    static void contract_tensors(const typename Storage<T, RankA, false>::type& A, 
                                 const typename Storage<T, RankB, false>::type& B, 
                                 typename Storage<T, RankA + RankB - 2, false>::type& C) 
    {
        std::cout << "Contracting dense tensors (generic dummy)...\n";
    }

    // UPDATED rule: Sparse T2 contracted with Dense T4 outputs a **Dense T2**
    // (Simulating a Coulomb/Exchange 2-index contraction over indices 2 and 3 of T4)
    template <typename T>
    static void contract_sparse_matrix_tensor4(const typename Storage<T, 2, true>::type& sparse,
                                               const typename Storage<T, 4, false>::type& t4,
                                               typename Storage<T, 2, false>::type& result) 
    {
        int dim0 = t4.dimension(0);
        int dim1 = t4.dimension(1);

        result.resize(dim0, dim1);
        result.setZero();

        // 2-index Contraction Loop: result(i,j) += sparse(r,c) * t4(i,j,r,c)
        for (int k = 0; k < sparse.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(sparse, k); it; ++it) {
                T val = it.value();
                int r = it.row(); // Contract with dimension 2
                int c = it.col(); // Contract with dimension 3

                for (int i = 0; i < dim0; ++i) {
                    for (int j = 0; j < dim1; ++j) {
                        result(i, j) += val * t4(i, j, r, c);
                    }
                }
            }
        }
    }
};
