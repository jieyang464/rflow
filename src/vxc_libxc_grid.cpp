#include "vxc_libxc_grid.h"
#include "vxc_evaluator.h"
#include "integral_provider.h"
#include "libxc_wrapper.h"
#include <stdexcept>

VxcFunctor make_libxc_vxc_functor() {
  return [](const VxcInputs& in) -> VxcResult {
    VxcResult out;
    if (!in.mol || !in.integrals || !in.Pa || !in.Pb) return out;

    try {
      const size_t nbf = in.integrals->nbf(*in.mol);
      scfcpp::DenseMatrix Va_dense(nbf, nbf);
      scfcpp::DenseMatrix Vb_dense(nbf, nbf);
      Va_dense.setZero();
      Vb_dense.setZero();
      out.Exc = 0.0;
      out.tr_PVxc = 0.0;

      bool want_x = (in.xc.exchange == "LDA_X");
      bool want_c = (in.xc.correlation == "LDA_C_VWN" || in.xc.correlation == "LDA_C_VWN_RPA");
      // Overlap-based diagonal approximation (illustrative, no grid):
      // q_mu = (P * S)_{mu,mu}
      scfcpp::DenseMatrix S = scfcpp::to_dense(in.integrals->overlap(*in.mol));
      scfcpp::DenseMatrix P = scfcpp::to_dense(*in.Pa) + scfcpp::to_dense(*in.Pb);
      scfcpp::DenseMatrix PS = P * S;
      for (size_t mu = 0; mu < nbf; ++mu) {
        double rho = std::max(0.0, PS(mu, mu));
        double eps_xc = 0.0;
        double v_rho = 0.0;
        if (rho > 1e-16) {
#ifdef USE_LIBXC
          if (want_x) {
            XCOutput xo = evaluate_xc_point(XC_LDA_X, XCFunctionalKind::LDA, XCInput{rho, 0.0});
            eps_xc += xo.eps_xc; v_rho += xo.v_rho;
          }
          if (want_c) {
            XCOutput co = evaluate_xc_point(XC_LDA_C_VWN, XCFunctionalKind::LDA, XCInput{rho, 0.0});
            eps_xc += co.eps_xc; v_rho += co.v_rho;
          }
#else
          if (want_x) {
            const double Cx = -0.75 * std::pow(3.0/M_PI, 1.0/3.0);
            double epsx = Cx * std::cbrt(rho);
            double vx   = (4.0/3.0) * epsx;
            eps_xc += epsx; v_rho += vx;
          }
          // VWN correlation not implemented without libxc (contributes zero)
#endif
        }
        out.Exc += rho * eps_xc;
        Va_dense(mu, mu) = v_rho;
        Vb_dense(mu, mu) = v_rho;
        out.tr_PVxc += v_rho * P(mu, mu);
      }
      out.Va = scfcpp::to_matrix(Va_dense);
      out.Vb = scfcpp::to_matrix(Vb_dense);
    } catch (const std::exception&) {
      // On any failure (e.g., no libxc), return zeros
      if (in.mol && in.integrals) {
        const size_t nbf = in.integrals->nbf(*in.mol);
        out.Va = scfcpp::Matrix(nbf, nbf);
        out.Vb = scfcpp::Matrix(nbf, nbf);
        out.Va.setZero();
        out.Vb.setZero();
      } else {
        out.Va = scfcpp::Matrix();
        out.Vb = scfcpp::Matrix();
      }
      out.Exc = 0.0;
      out.tr_PVxc = 0.0;
    }

    return out;
  };
}
