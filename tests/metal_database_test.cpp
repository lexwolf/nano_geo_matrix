#define CUP_BACKEND_QUASI_STATIC

#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/single.hpp>
#include "cup.hpp"

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool approx_equal(const std::complex<double>& a,
                  const std::complex<double>& b,
                  double tol = 1.0e-12)
{
    return std::abs(a - b) < tol;
}

void print_eps_row(const std::string& label,
                   nanosphere& ns,
                   double omega_eV)
{
    const std::complex<double> eps = ns.metal(omega_eV);

    std::cout << std::setw(12) << label
              << "  omega = " << std::fixed << std::setprecision(6) << omega_eV
              << " eV"
              << "   eps = ("
              << std::setw(14) << std::setprecision(8) << std::real(eps)
              << ", "
              << std::setw(14) << std::setprecision(8) << std::imag(eps)
              << ")\n";
}

} // namespace

int main()
{
    nanosphere ns_default;
    nanosphere ns_jc;
    nanosphere ns_unical;
    nanosphere ns_drude;

    ns_default.init();
    ns_jc.init();
    ns_unical.init();
    ns_drude.init();

    // 1) default spline silver should be equivalent to explicit JC
    ns_default.set_metal("silver", "spline");
    ns_jc.set_metal("silver", "spline", 0, "jc");

    const std::vector<double> probe_omegas = {
        1.50, 2.00, 2.50, 3.00, 3.50
    };

    bool default_matches_jc = true;

    std::cout << "\n=== Checking default silver spline vs explicit JC ===\n";
    for (double omega : probe_omegas) {
        const std::complex<double> eps_default = ns_default.metal(omega);
        const std::complex<double> eps_jc      = ns_jc.metal(omega);

        const bool ok = approx_equal(eps_default, eps_jc);
        default_matches_jc = default_matches_jc && ok;

        std::cout << "omega = " << std::fixed << std::setprecision(6) << omega
                  << " eV   "
                  << (ok ? "[OK]   " : "[FAIL] ")
                  << "default = (" << std::setprecision(8)
                  << std::real(eps_default) << ", " << std::imag(eps_default) << ")"
                  << "   jc = ("
                  << std::real(eps_jc) << ", " << std::imag(eps_jc) << ")\n";
    }

    if (!default_matches_jc) {
        std::cerr << "\nERROR: default spline silver is not matching explicit jc.\n";
        return EXIT_FAILURE;
    }

    // 2) compare JC and UNICAL side by side
    ns_unical.set_metal("silver", "spline", 0, "unical");

    std::cout << "\n=== Silver spline datasets: JC vs UNICAL ===\n";
    for (double omega : probe_omegas) {
        const std::complex<double> eps_jc      = ns_jc.metal(omega);
        const std::complex<double> eps_unical  = ns_unical.metal(omega);
        const std::complex<double> delta       = eps_unical - eps_jc;

        std::cout << "omega = " << std::fixed << std::setprecision(6) << omega << " eV\n";
        std::cout << "   jc      = (" << std::setprecision(8)
                  << std::real(eps_jc) << ", " << std::imag(eps_jc) << ")\n";
        std::cout << "   unical  = (" << std::real(eps_unical) << ", "
                  << std::imag(eps_unical) << ")\n";
        std::cout << "   delta   = (" << std::real(delta) << ", "
                  << std::imag(delta) << ")\n";
    }

    // 3) drude still works
    ns_drude.set_metal("silver", "drude");

    std::cout << "\n=== Silver drude reference ===\n";
    for (double omega : probe_omegas) {
        print_eps_row("drude", ns_drude, omega);
    }

    // 4) warning path: unical + sel != 0
    std::cout << "\n=== Triggering warning path: UNICAL with sel = +1 ===\n";
    nanosphere ns_warn;
    ns_warn.init();
    ns_warn.set_metal("silver", "spline", +1, "unical");

    for (double omega : probe_omegas) {
        print_eps_row("unical+1", ns_warn, omega);
    }

    std::cout << "\nAll smoke tests completed.\n";
    return EXIT_SUCCESS;
}