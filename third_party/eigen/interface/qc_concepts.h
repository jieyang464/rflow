#pragma once
#include <concepts>
#include <cstddef>

// =====================================================================
// C++20 CONCEPTS: THE BACKEND CONTRACT
// Defines exactly what types and static methods a Backend struct must provide.
// =====================================================================

template <typename Backend, typename T>
concept TensorBackend = requires {
    // 1. Must define a template struct/alias 'Storage' that takes <T, Rank, IsSparse>
    // Note: Checking templates in concepts can be tricky, so we check specific instantiations to ensure validity.
    typename Backend::template Storage<T, 4, false>::type; // Raw Tensor
    typename Backend::template Storage<T, 2, false>::type; // Dense Matrix
    typename Backend::template Storage<T, 2, true>::type;  // Sparse Matrix
    typename Backend::template Storage<T, 1, false>::type; // Dense Vector
};
