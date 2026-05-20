#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>

// =====================================================================
// EIGEN BACKEND TYPES
// Maps our interface requirements to specific Eigen types.
// =====================================================================
namespace backend {
    template <typename T, size_t Rank, bool IsSparse>
    struct Storage;

    // 1. Generic Dense Tensor (Rank N)
    template <typename T, size_t Rank>
    struct Storage<T, Rank, false> {
        using type = Eigen::Tensor<T, Rank>;
    };

    // 2. Dense Matrix (Rank 2) -> Needs MatrixXd for diagonalization
    template <typename T>
    struct Storage<T, 2, false> {
        using type = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    };

    // 3. Dense Vector (Rank 1)
    template <typename T>
    struct Storage<T, 1, false> {
        using type = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    };

    // 4. Sparse Matrix (Rank 2)
    template <typename T>
    struct Storage<T, 2, true> {
        using type = Eigen::SparseMatrix<T, Eigen::RowMajor>; // RowMajor is often better for contraction loops
    };

    // 5. Sparse Vector (Rank 1)
    template <typename T>
    struct Storage<T, 1, true> {
        using type = Eigen::SparseVector<T>;
    };
}

// =====================================================================
// EIGEN BACKEND OPERATIONS (Hides the horrible loops!)
// =====================================================================
namespace backend {

    // --- Diagonalize Dense Matrix ---
    template <typename T>
    void diagonalize_dense(const typename Storage<T, 2, false>::type& mat, 
                           typename Storage<T, 1, false>::type& evals, 
                           typename Storage<T, 2, false>::type& evecs) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solver(mat);
        evals = solver.eigenvalues();
        evecs = solver.eigenvectors();
    }

    // --- Diagonalize Sparse Matrix (Dummy placeholder for Spectra/Arpack) ---
    template <typename T>
    void diagonalize_sparse(const typename Storage<T, 2, true>::type& mat, 
                            typename Storage<T, 1, false>::type& evals, 
                            typename Storage<T, 2, false>::type& evecs) {
        std::cout << "[Backend] Warning: Sparse diagonalization requires an external library like Spectra.\n";
        // E.g., Spectra::SymEigsSolver<...>
    }

    // --- Contract Sparse Matrix (Rank 2) with Dense Tensor (Rank 4) ---
    // Hides the horrible raw inner-iterator loop implementation out of the user's sight.
    template <typename T>
    void contract_sparse_matrix_tensor4(const typename Storage<T, 2, true>::type& sparse,
                                        const typename Storage<T, 4, false>::type& t4,
                                        typename Storage<T, 4, false>::type& result) 
    {
        // Assuming we contract sparse cols with t4's 2nd dimension (index 1)
        int dim0 = t4.dimension(0), dim1 = t4.dimension(1);
        int dim2 = t4.dimension(2), dim3 = t4.dimension(3);
        int new_dim1 = sparse.rows();

        result.resize(dim0, new_dim1, dim2, dim3);
        result.setZero();

        // The optimized iterator loop
        for (int k = 0; k < sparse.outerSize(); ++k) {
            for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(sparse, k); it; ++it) {
                T val = it.value();
                int r = it.row(); // target dim in result
                int c = it.col(); // source dim in t4

                for (int i = 0; i < dim0; ++i) {
                    for (int j = 0; j < dim2; ++j) {
                        for (int l = 0; l < dim3; ++l) {
                            result(i, r, j, l) += val * t4(i, c, j, l);
                        }
                    }
                }
            }
        }
    }

    // (More generic contract functions would go here...)
}
