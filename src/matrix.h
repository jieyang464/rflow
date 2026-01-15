#pragma once

#include <Eigen/Dense>
#ifdef SCFCPP_USE_EIGEN_SPARSE
#include <Eigen/Sparse>
#endif

namespace scfcpp {

#ifdef SCFCPP_USE_EIGEN_SPARSE
using Matrix = Eigen::SparseMatrix<double>;
using DenseMatrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Vector3 = Eigen::Vector3d;

inline DenseMatrix to_dense(const Matrix& m) { return DenseMatrix(m); }
inline Matrix to_matrix(const DenseMatrix& m) { return m.sparseView(); }
#else
using Matrix = Eigen::MatrixXd;
using DenseMatrix = Eigen::MatrixXd;
using Vector = Eigen::VectorXd;
using Vector3 = Eigen::Vector3d;

inline DenseMatrix to_dense(const Matrix& m) { return m; }
inline Matrix to_matrix(const DenseMatrix& m) { return m; }
#endif

}
