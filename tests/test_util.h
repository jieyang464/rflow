// test_util.h — Minimal assertion harness.
//
// Deliberately dependency-free: this project vendors Eigen and links libint2,
// and adding a unit-test framework on top of that is more build surface than a
// hundred lines of checks are worth.
#pragma once

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "types.h"

namespace testing {

class TestRunner {
public:
    explicit TestRunner(std::string suite) : suite_(std::move(suite)) {
        std::printf("=== %s ===\n", suite_.c_str());
    }

    void Section(const std::string& name) {
        std::printf("\n-- %s\n", name.c_str());
    }

    void Check(bool ok, const std::string& what, const std::string& detail = "") {
        ++checks_;
        if (ok) {
            std::printf("   ok    %s\n", what.c_str());
        } else {
            ++failures_;
            std::printf("   FAIL  %s%s%s\n", what.c_str(),
                        detail.empty() ? "" : "  |  ", detail.c_str());
        }
    }

    // Absolute-tolerance comparison; prints the actual deviation either way so a
    // passing run still shows how much headroom there is.
    void Near(double got, double want, double tol, const std::string& what) {
        const double diff = std::fabs(got - want);
        ++checks_;
        if (diff <= tol && std::isfinite(diff)) {
            std::printf("   ok    %-52s  |got - want| = %.3e  (tol %.1e)\n",
                        what.c_str(), diff, tol);
        } else {
            ++failures_;
            std::printf("   FAIL  %-52s  got %.12g  want %.12g  diff %.3e  (tol %.1e)\n",
                        what.c_str(), got, want, diff, tol);
        }
    }

    void Report(const std::string& label, double value) {
        std::printf("   info  %-52s  %.12g\n", label.c_str(), value);
    }

    int Summary() {
        std::printf("\n%s: %d checks, %d failed -> %s\n\n", suite_.c_str(), checks_,
                    failures_, failures_ == 0 ? "PASS" : "FAIL");
        return failures_ == 0 ? 0 : 1;
    }

    int failures() const { return failures_; }

private:
    std::string suite_;
    int checks_{0};
    int failures_{0};
};

inline double MaxAbsDiff(const T2& a, const T2& b) {
    if (a.dimension(0) != b.dimension(0) || a.dimension(1) != b.dimension(1)) {
        return std::numeric_limits<double>::infinity();
    }
    double worst = 0.0;
    for (Eigen::Index i = 0; i < a.dimension(0); ++i)
        for (Eigen::Index j = 0; j < a.dimension(1); ++j)
            worst = std::max(worst, std::fabs(a(i, j) - b(i, j)));
    return worst;
}

inline double MaxAbsDiff(const T4& a, const T4& b) {
    for (int d = 0; d < 4; ++d) {
        if (a.dimension(d) != b.dimension(d)) return std::numeric_limits<double>::infinity();
    }
    double worst = 0.0;
    for (Eigen::Index i = 0; i < a.dimension(0); ++i)
        for (Eigen::Index j = 0; j < a.dimension(1); ++j)
            for (Eigen::Index k = 0; k < a.dimension(2); ++k)
                for (Eigen::Index l = 0; l < a.dimension(3); ++l)
                    worst = std::max(worst, std::fabs(a(i, j, k, l) - b(i, j, k, l)));
    return worst;
}

inline void PrintMatrix(const std::string& label, const T2& m) {
    std::printf("   %s:\n", label.c_str());
    for (Eigen::Index i = 0; i < m.dimension(0); ++i) {
        std::printf("     ");
        for (Eigen::Index j = 0; j < m.dimension(1); ++j) std::printf("% .9f  ", m(i, j));
        std::printf("\n");
    }
}

}  // namespace testing
