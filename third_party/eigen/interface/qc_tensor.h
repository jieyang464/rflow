#pragma once
#include "qc_concepts.h"

// =====================================================================
// COMPILE-TIME BACKEND DISPATCHER (Controlled by Compiler Flags!)
// =====================================================================

#if defined(USE_EIGEN_BACKEND)
    #include "backend_eigen.h"
    using ActiveBackend = EigenBackend;
#elif defined(USE_CUSTOM_BACKEND)
    #include "backend_custom.h"
    using ActiveBackend = CustomBackend;
#else
    #error "FATAL: No tensor backend defined! Compile with -DUSE_EIGEN_BACKEND or -DUSE_CUSTOM_BACKEND"
#endif

// Validate the active backend against our Concept contract
static_assert(TensorBackend<ActiveBackend, double>, "Selected backend does not fulfill TensorBackend concept!");

template <typename T, size_t Rank, bool IsSparse = false>
class Tensor {
public:
    using StorageType = typename ActiveBackend::template Storage<T, Rank, IsSparse>::type;
    StorageType internal_data;

    Tensor() = default;

    template <size_t OtherRank>
    Tensor<T, Rank + OtherRank - 2, false> contract(const Tensor<T, OtherRank, false>& other) const {
        Tensor<T, Rank + OtherRank - 2, false> result;
        ActiveBackend::contract_tensors<T, Rank, OtherRank>(this->internal_data, other.internal_data, result.internal_data);
        return result;
    }

    Tensor<T, 2, false> contract(const Tensor<T, 4, false>& other) const 
    requires (Rank == 2 && IsSparse) 
    {
        Tensor<T, 2, false> result;
        ActiveBackend::contract_sparse_matrix_tensor4<T>(this->internal_data, other.internal_data, result.internal_data);
        return result;
    }

    void diagonalize(Tensor<T, 1, false>& eigenvalues, Tensor<T, 2, false>& eigenvectors) const 
    requires (Rank == 2 && !IsSparse) 
    {
        ActiveBackend::diagonalize_dense<T>(this->internal_data, eigenvalues.internal_data, eigenvectors.internal_data);
    }

    void diagonalize(Tensor<T, 1, false>& eigenvalues, Tensor<T, 2, false>& eigenvectors) const 
    requires (Rank == 2 && IsSparse) 
    {
        ActiveBackend::diagonalize_sparse<T>(this->internal_data, eigenvalues.internal_data, eigenvectors.internal_data);
    }
};

namespace qc {
    template <typename T> using T1       = Tensor<T, 1, false>;
    template <typename T> using T2       = Tensor<T, 2, false>;
    template <typename T> using T2sparse = Tensor<T, 2, true>;
    template <typename T> using T3       = Tensor<T, 3, false>;
    template <typename T> using T4       = Tensor<T, 4, false>;
}
