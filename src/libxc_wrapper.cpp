#include "libxc_wrapper.h"
#include <stdexcept>

#ifdef USE_LIBXC
#include <xc.h>
#endif

XCOutput evaluate_xc_point(int functional_id, XCFunctionalKind kind, const XCInput& in) {
#ifndef USE_LIBXC
    (void)functional_id; (void)kind; (void)in;
    throw std::runtime_error("LibXC not enabled. Define USE_LIBXC and link libxc to use LDA/GGA.");
#else
    XCOutput out{0.0, 0.0, 0.0};

    xc_func_type func;
    if (xc_func_init(&func, functional_id, XC_UNPOLARIZED) != 0) {
        throw std::runtime_error("xc_func_init failed for given functional id");
    }

    const int n = 1;
    if (kind == XCFunctionalKind::LDA) {
        double rho = std::max(in.rho, 0.0);
        double eps[n];
        double vrho[n];
        xc_lda_exc_vxc(&func, n, &rho, eps, vrho);
        out.eps_xc = eps[0];
        out.v_rho  = vrho[0];
        out.v_sigma = 0.0;
    } else if (kind == XCFunctionalKind::GGA) {
        double rho = std::max(in.rho, 0.0);
        double sigma = std::max(in.sigma, 0.0);
        double eps[n];
        double vrho[n];
        double vsigma[n];
        xc_gga_exc_vxc(&func, n, &rho, &sigma, eps, vrho, vsigma);
        out.eps_xc = eps[0];
        out.v_rho  = vrho[0];
        out.v_sigma = vsigma[0];
    } else {
        xc_func_end(&func);
        throw std::runtime_error("Unsupported functional kind in evaluate_xc_point");
    }

    xc_func_end(&func);
    return out;
#endif
}
