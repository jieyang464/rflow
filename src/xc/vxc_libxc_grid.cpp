#include "xc/vxc_libxc_grid.h"
#include "scf.h"
#include "xc/libxc_wrapper.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {

constexpr double kPi = 3.14159265358979323846;

double OverlapWeightedPopulation(const T2& D, const T2& S, Eigen::Index mu) {
  const Eigen::Index nbf = D.dimension(0);
  double population = 0.0;
  for (Eigen::Index nu = 0; nu < nbf; ++nu) {
    population += D(mu, nu) * S(nu, mu);
  }
  return std::max(0.0, population);
}

void AccumulateBuiltinLsdaExchange(double rho_a, double rho_b,
                                   double& exc, double& vxa, double& vxb) {
  const double kCx = -1.5 * std::cbrt(3.0 / (4.0 * kPi));

  auto rho43 = [](double rho) {
    return rho > 0.0 ? std::pow(rho, 4.0 / 3.0) : 0.0;
  };
  auto rho13 = [](double rho) {
    return rho > 0.0 ? std::cbrt(rho) : 0.0;
  };

  exc += kCx * (rho43(rho_a) + rho43(rho_b));
  vxa += (4.0 / 3.0) * kCx * rho13(rho_a);
  vxb += (4.0 / 3.0) * kCx * rho13(rho_b);
}

void AccumulateLibxcLda(int functional_id, double rho_a, double rho_b,
                        double& exc, double& vxa, double& vxb) {
  const SpinXCOutput xc = evaluate_spin_lda_point(functional_id, SpinXCInput{rho_a, rho_b});
  exc += (rho_a + rho_b) * xc.eps_xc;
  vxa += xc.v_rho_a;
  vxb += xc.v_rho_b;
}

}  // namespace

VxcFunctor make_libxc_vxc_functor() {
  return [](const UksDensityInput& in, UksVxcOutput& out) {
    const T2& S = in.overlap;
    const Eigen::Index nbf = S.dimension(0);
    if (S.dimension(1) != nbf || in.Da.dimension(0) != nbf || in.Da.dimension(1) != nbf ||
        in.Db.dimension(0) != nbf || in.Db.dimension(1) != nbf) {
      throw std::invalid_argument("UKS Vxc functor requires square Da/Db/S matrices of identical shape.");
    }

    for (Eigen::Index mu = 0; mu < nbf; ++mu) {
      const double rho_a = OverlapWeightedPopulation(in.Da, S, mu);
      const double rho_b = OverlapWeightedPopulation(in.Db, S, mu);

      double exc_mu = 0.0;
      double vxa_mu = 0.0;
      double vxb_mu = 0.0;

      try {
        if (in.xc_functional_id == xc::BuiltinFunctionalId::LsdaExchange) {
          AccumulateBuiltinLsdaExchange(rho_a, rho_b, exc_mu, vxa_mu, vxb_mu);
        } else if (in.xc_functional_id != xc::BuiltinFunctionalId::None) {
          AccumulateLibxcLda(in.xc_functional_id, rho_a, rho_b,
                             exc_mu, vxa_mu, vxb_mu);
        }
      } catch (const std::exception&) {
        // Preserve zero XC contributions if libxc is unavailable or rejects the id.
      }

      out.Exc += exc_mu;
      out.Vxca(mu, mu) += vxa_mu;
      out.Vxcb(mu, mu) += vxb_mu;
      // Use the overlap-weighted populations `rho_a`/`rho_b` (Mulliken)
      // to accumulate Tr(P Vxc, approximated on the grid).
      out.tr_PVxc += rho_a * vxa_mu + rho_b * vxb_mu;
    }
  };
}
