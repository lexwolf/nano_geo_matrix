#define CUP_BACKEND_QUASI_STATIC

#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/single.hpp>
#include "cup.hpp"

#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

namespace {

void print_eps_row(const char* label,
                   nanosphere& ns,
                   double omega_eV,
                   bool expect_zero_imag,
                   bool& ok)
{
    const std::complex<double> eps = ns.dielectric(omega_eV);
    const double imag_eps = std::imag(eps);
    const bool imag_ok = expect_zero_imag ? (std::abs(imag_eps) < 1.0e-6)
                                          : (imag_eps > 0.0);
    ok = ok && imag_ok;

    std::cout << std::setw(8) << label
              << "  omega = " << std::fixed << std::setprecision(6) << omega_eV
              << " eV"
              << "   eps = ("
              << std::setw(14) << std::setprecision(8) << std::real(eps)
              << ", "
              << std::setw(14) << std::imag(eps)
              << ")   "
              << (imag_ok ? "[OK]" : "[FAIL]") << '\n';
}

} // namespace

int main()
{
    // Smoke test only:
    // - verifies dielectric tables load
    // - verifies interpolation works
    // - verifies glass is nearly lossless and ITO is lossy
    // This is not a full numerical validation against reference data.
    nanosphere glass;
    nanosphere ito;

    glass.init();
    ito.init();

    glass.set_dielectric("glass");
    ito.set_dielectric("ito");

    const std::vector<double> probe_omegas = {1.50, 2.00, 3.00};

    bool ok = true;

    std::cout << "\n=== Glass dielectric spline ===\n";
    for (double omega : probe_omegas) {
        print_eps_row("glass", glass, omega, true, ok);
    }

    std::cout << "\n=== ITO dielectric spline ===\n";
    for (double omega : probe_omegas) {
        print_eps_row("ito", ito, omega, false, ok);
    }

    if (!ok) {
        std::cerr << "\nERROR: dielectric database test failed.\n";
        return EXIT_FAILURE;
    }

    std::cout << "\nDielectric database smoke test completed.\n";
    return EXIT_SUCCESS;
}
