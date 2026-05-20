#pragma once
#include <string_view>
#include <array>
#include <stdexcept>
#include <vector>
#include <string>
#include <iostream>

// =====================================================================
// EINSUM PARSER (Compile-time or Runtime)
// Translates "xjyl", "pqyx" -> "jpql" into Eigen contraction pairs and shuffle maps.
// =====================================================================

struct EinsumPlan {
    std::vector<std::pair<int, int>> contract_pairs; // e.g. {{0,3}, {2,2}}
    std::vector<int> shuffle_map;                    // e.g. {0, 2, 3, 1}
    bool needs_shuffle;
};

// A simple runtime parser for the strings. 
// (In C++20 this can be made string_view + constexpr for zero overhead)
inline EinsumPlan parse_einsum(std::string_view strA, std::string_view strB, std::string_view strOut) {
    EinsumPlan plan;
    std::string uncontracted_A = "";
    std::string uncontracted_B = "";
    
    // 1. Find Contractions
    for (int i = 0; i < strA.size(); ++i) {
        bool contracted = false;
        for (int j = 0; j < strB.size(); ++j) {
            if (strA[i] == strB[j]) {
                plan.contract_pairs.push_back({i, j});
                contracted = true;
                break;
            }
        }
        if (!contracted) uncontracted_A += strA[i];
    }
    
    // Find uncontracted in B
    for (int j = 0; j < strB.size(); ++j) {
        bool contracted = false;
        for (int i = 0; i < strA.size(); ++i) {
            if (strA[i] == strB[j]) { contracted = true; break; }
        }
        if (!contracted) uncontracted_B += strB[j];
    }

    // 2. Eigen's Default Output Order
    std::string eigen_natural_output = uncontracted_A + uncontracted_B;

    // Error Checking
    if (eigen_natural_output.size() != strOut.size()) {
        throw std::runtime_error("Einsum error: output string length mismatch.");
    }

    // 3. Determine Shuffle
    plan.needs_shuffle = (eigen_natural_output != strOut);
    if (plan.needs_shuffle) {
        for (char c : strOut) {
            size_t pos = eigen_natural_output.find(c);
            if (pos == std::string::npos) {
                throw std::runtime_error(std::string("Einsum error: output char ") + c + " not found in inputs.");
            }
            plan.shuffle_map.push_back(pos);
        }
    }

    return plan;
}
