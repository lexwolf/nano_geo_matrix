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
        check(!relative_change_below(invalid, 1., 1.e-6), "invalid current converged");
        check(!relative_change_below(invalid, invalid, 1.e-6), "invalid pair converged");
        check(!averaging_window_complete(invalid, 1000.), "invalid omeR completed a window");
        check(!averaging_window_complete(1., invalid), "invalid elapsed completed a window");
        check(!amplitude_above(1.e30, invalid, 5.e26), "invalid E0 triggered cutoff");
        check(!amplitude_above(inf, invalid, 5.e26), "invalid E0 allowed infinite numerator division");

        double current = 2., previous = 2.;
        invalidate_average_if_needed(invalid, 1000., 2000., current, previous);
        check(!relative_change_below(current, previous, 1.e-6), "invalid omeR left stale averages");
        current = previous = 2.;
        invalidate_average_if_needed(1., invalid, 2000., current, previous);
        check(!relative_change_below(current, previous, 1.e-6), "invalid elapsed left stale averages");
        current = previous = 2.;
        store_average(2000., invalid, current, previous);
        check(!relative_change_below(current, previous, 1.e-6), "invalid elapsed produced an average");
    }
    for (double invalid : {inf, -inf, nan}) {
        check(!relative_change_below(1., invalid, 1.e-6), "nonfinite previous converged");
        double current = 2., previous = 2.;
        invalidate_average_if_needed(1., 1000., invalid, current, previous);
        check(current == 0. && previous == 0., "nonfinite total left stale averages");
        current = previous = 2.;
        store_average(invalid, 1000., current, previous);
        check(!relative_change_below(current, previous, 1.e-6), "nonfinite total produced convergence");
    }
    check(!relative_change_below(0., 0., 1.e-6), "zero averages converged");
    double current = 2., previous = 2.;
    store_average(0., 1000., current, previous);
    check(!relative_change_below(current, previous, 1.e-6), "zero block average converged");
    check(!(std::fetestexcept(FE_DIVBYZERO | FE_INVALID)), "invalid division in rejected inputs");

    // All tolerances used by either backend retain their original formula.
    for (double tolerance : {1.e-30, 1.e-20, 1.e-10, 1.e-7, 1.e-6}) {
        for (double now : {1., 2., 1.e-5}) {
            for (double before : {0., 1., 2., now, now*(1.-1.e-8)}) {
                check(relative_change_below(now, before, tolerance) ==
                      (std::fabs((now-before)/now) < tolerance), "finite comparison changed");
            }
        }
    }
    for (double frequency : {0.5, 1., 2.}) {
        for (double elapsed : {1., 500., 1000., 2000.}) {
            check(averaging_window_complete(frequency, elapsed) ==
                  (elapsed >= 1000/frequency), "finite window changed");
        }
    }
    current = 2.; previous = 1.;
    invalidate_average_if_needed(1., 1000., 3000., current, previous);
    check(current == 2. && previous == 1., "valid block invalidated averages");
    store_average(3000., 1000., current, previous);
    check(current == 3. && previous == 2., "finite average changed");
    for (double field : {0.5, 1., 2.}) {
        for (double intensity : {1., 5.e26, 1.e30}) {
            check(amplitude_above(intensity, field, 5.e26) ==
                  (intensity/field > 5.e26), "finite amplitude cutoff changed");
        }
    }
    // Overflowing a valid quotient must also discard both old averages.
    store_average(std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::min(), current, previous);
    check(current == 0. && previous == 0., "overflowed average left stale averages");
    check(!(std::fetestexcept(FE_DIVBYZERO | FE_INVALID)), "invalid division in finite inputs");
    std::cout << "Stopping arithmetic failures: " << failures << '\n';
    return failures ? 1 : 0;
}
