#include "stopping_arithmetic.hpp"

#include <cfenv>
#include <iostream>
#include <limits>

namespace {
int failures = 0;
void check(bool result, const char* message)
{
    if (!result) {
        std::cerr << message << '\n';
        ++failures;
    }
}
}

int main()
{
    using namespace cup_stopping_detail;
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::feclearexcept(FE_ALL_EXCEPT);
    for (double invalid : {0., -1., inf, -inf, nan}) {
        check(!relative_change_below(invalid, 1., 1.e-6), "invalid current enabled Fourier collection");
        check(!relative_change_below(invalid, invalid, 1.e-6), "invalid pair enabled Fourier collection");
        check(!amplitude_above(1.e30, invalid, 5.e26), "invalid E0 triggered cutoff");
        check(!amplitude_above(inf, invalid, 5.e26), "invalid E0 allowed infinite numerator division");
    }
    for (double invalid : {inf, -inf, nan}) {
        check(!relative_change_below(1., invalid, 1.e-6), "nonfinite previous enabled Fourier collection");
    }
    check(!relative_change_below(0., 0., 1.e-10), "zero Fourier intensity accepted");
    check(!(std::fetestexcept(FE_DIVBYZERO | FE_INVALID)), "invalid division in rejected inputs");

    // This legacy guarded comparison remains solely for Fourier collection.
    // Rabi/block-average guards were removed with the obsolete detector.
    for (double tolerance : {1.e-30, 1.e-20, 1.e-10, 1.e-7, 1.e-6}) {
        for (double now : {1., 2., 1.e-5}) {
            for (double before : {0., 1., 2., now, now*(1.-1.e-8)}) {
                check(relative_change_below(now, before, tolerance) ==
                      (std::fabs((now-before)/now) < tolerance), "finite comparison changed");
            }
        }
    }
    for (double field : {0.5, 1., 2.}) {
        for (double intensity : {1., 5.e26, 1.e30}) {
            check(amplitude_above(intensity, field, 5.e26) ==
                  (intensity/field > 5.e26), "finite amplitude cutoff changed");
        }
    }
    check(!(std::fetestexcept(FE_DIVBYZERO | FE_INVALID)), "invalid division in finite inputs");
    std::cout << "Stopping arithmetic failures: " << failures << '\n';
    return failures ? 1 : 0;
}
