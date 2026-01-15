// Optional LibXC wrapper for local XC evaluation (LDA/GGA)
#pragma once

struct XCInput {
    double rho;      // electron density at point
    double sigma;    // |grad rho|^2 (for GGA); ignored for LDA
};

struct XCOutput {
    double eps_xc;   // energy density per particle
    double v_rho;    // d(epsilon)/d rho (LDA term)
    double v_sigma;  // d(epsilon)/d sigma (GGA term); zero for LDA
};

enum class XCFunctionalKind { LDA, GGA };

// Evaluate XC locally using libxc if enabled at build time; otherwise throws.
// functional_id: integer from LibXC functional list (e.g., XC_LDA_X, etc.)
XCOutput evaluate_xc_point(int functional_id, XCFunctionalKind kind, const XCInput& in);
