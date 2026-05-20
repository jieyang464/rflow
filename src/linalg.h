// linalg.h — Linear algebra helpers (header-only, inline)
#pragma once
#include "types.h"

// ── Matrix multiply: C = A * B ──────────────────────────────────────────────
inline T2 MatMul(const T2& A, const T2& B) {
    const Eigen::array<Eigen::IndexPair<int>, 1> dims = {
        Eigen::IndexPair<int>(1, 0)};
    return A.contract(B, dims);
}

// ── Commutator: [A, B] = AB − BA ───────────────────────────────────────────
inline T2 Commutator(const T2& A, const T2& B) {
    return MatMul(A, B) - MatMul(B, A);
}

// ── S-metric commutator: [A, B]_S = ASB − BSA ──────────────────────────────
inline T2 SMetricCommutator(const T2& A, const T2& B, const T2& S) {
    return MatMul(MatMul(A, S), B) - MatMul(MatMul(B, S), A);
}

// ── DIIS error matrix: e = FDS − SDF ───────────────────────────────────────
inline T2 ComputeDIISError(const T2& F, const T2& D, const T2& S) {
    return MatMul(MatMul(F, D), S) - MatMul(MatMul(S, D), F);
}

// ── Max absolute element  (convergence metric) ─────────────────────────────
inline double MaxAbsElement(const T2& A) {
    Eigen::Tensor<double, 0> m = A.abs().maximum();
    return m(0);
}

// ── Trace: tr(A) ───────────────────────────────────────────────────────────
inline double Trace(const T2& A) {
    const Eigen::Index n = A.dimension(0);
    double sum = 0.0;
    for (Eigen::Index i = 0; i < n; ++i) sum += A(i, i);
    return sum;
}

// ── Identity matrix ────────────────────────────────────────────────────────
inline T2 Eye(Eigen::Index n) {
    T2 I(n, n);
    I.setZero();
    for (Eigen::Index i = 0; i < n; ++i) I(i, i) = 1.0;
    return I;
}

// ── Matrix exponential via Taylor series ────────────────────────────────────
//    e^X = I + X + X²/2! + X³/3! + … + X^order/order!
inline T2 MatExp(const T2& X, int order) {
    const Eigen::Index n = X.dimension(0);
    T2 result = Eye(n);      // accumulator  (starts at I)
    T2 term   = Eye(n);      // term_k = X^k / k!

    for (int k = 1; k <= order; ++k) {
        term   = MatMul(term, X) * (1.0 / k);   // X^k / k!
        result = result + term;
    }
    return result;
}

// ── BCH similarity transform: e^{−X} A e^{X} ──────────────────────────────
//
//    Using the adjoint‑action expansion:
//      e^{−X} A e^{X} = Σ_{k=0}^{∞}  (−1)^k / k!  ad_X^k(A)
//    where  ad_X^0(A) = A,  ad_X^k(A) = [X, ad_X^{k−1}(A)]
//
//    Truncated at the given order.
//
inline T2 BCHTransform(const T2& X, const T2& A, int order) {
    T2     result = A;
    T2     nested = A;          // ad_X^k(A)
    double sign   = -1.0;      // (−1)^k starts at k=1 → −1
    double fact   = 1.0;

    for (int k = 1; k <= order; ++k) {
        nested = Commutator(X, nested);       // ad_X^k(A)
        fact  *= k;
        result = result + nested * (sign / fact);
        sign   = -sign;
    }
    return result;
}

// ── S-metric BCH transform: e^{−XS} A e^{SX} ──────────────────────────────
//
//    Uses the S-metric adjoint action:
//      e^{−XS} A e^{SX} = Σ_{k=0}^{∞}  (−1)^k / k!  ad_{X,S}^k(A)
//    where  ad_{X,S}^0(A) = A,
//           ad_{X,S}^k(A) = [X, ad_{X,S}^{k−1}(A)]_S
//
//    Truncated at the given order.
//
inline T2 BCHTransform(const T2& X, const T2& A, const T2& S, int order) {
    T2     result = A;
    T2     nested = A;          // ad_{X,S}^k(A)
    double sign   = -1.0;      // (−1)^k starts at k=1 → −1
    double fact   = 1.0;

    for (int k = 1; k <= order; ++k) {
        nested = SMetricCommutator(X, nested, S);
        fact  *= k;
        result = result + nested * (sign / fact);
        sign   = -sign;
    }
    return result;
}

// ── Taylor similarity transform: e^{−X} A e^{X} ───────────────────────────
//    Compute e^{−X} and e^{X} explicitly, then multiply.
//    More expensive than BCH (3 MatExp + 2 MatMul) but can be more
//    accurate at low orders for larger ‖X‖.
//
inline T2 TaylorTransform(const T2& X, const T2& A, int order) {
    T2 expNegX = MatExp(X * (-1.0), order);
    T2 expX    = MatExp(X, order);
    return MatMul(MatMul(expNegX, A), expX);
}

// ── Result of generalized eigenvalue solve ─────────────────────────────────
struct DiagResult {
    T2 C;      // Coefficient matrix (eigenvectors)
    T2 eps;    // Eigenvalues as diagonal matrix
};

// ── Generalized eigenvalue problem: FC = SCε ────────────────────────────────
//    Diagonalize F in the metric of S.
inline DiagResult DiagonalizeInSMetric(const T2& F, const T2& S) {
    const Eigen::Index n = F.dimension(0);
    using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                 Eigen::ColMajor>;
    const Eigen::Map<const Matrix> f_view(F.data(), n, n);
    const Eigen::Map<const Matrix> s_view(S.data(), n, n);

    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver;
    solver.compute(Matrix(f_view), Matrix(s_view), Eigen::ComputeEigenvectors);

    DiagResult out;
    out.C = T2(n, n);
    out.eps = T2(n, n);
    out.eps.setZero();
    for (int i = 0; i < n; ++i)
        out.eps(i, i) = solver.eigenvalues()(i);
    std::copy(solver.eigenvectors().data(),
              solver.eigenvectors().data() + (n * n), out.C.data());
    return out;
}

// ── Build density matrix: D = C_occ · C_occᵀ ───────────────────────────────
//    Fills the first n_elec columns of C into C_occ and contracts.
inline T2 BuildDensityMatrix(const T2& C, int n_elec) {
    const Eigen::Index n = C.dimension(0);

    T2 C_occ(n, n);
    C_occ.setZero();
    for (Eigen::Index i = 0; i < n; ++i)
        for (int m = 0; m < n_elec; ++m)
            C_occ(i, m) = C(i, m);

    // D = C_occ * C_occ^T
    const Eigen::array<Eigen::IndexPair<int>, 1> cols = {
        Eigen::IndexPair<int>(1, 1)};
    return C_occ.contract(C_occ, cols);
}
