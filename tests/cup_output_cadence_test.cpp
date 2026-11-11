#include "automatic_output.hpp"

#include <cfenv>
#include <iostream>
#include <limits>

int main()
{
    using cup_output_detail::automatic_cadence;
    int failures = 0;
    const auto check = [&](bool ok) { if (!ok) ++failures; };
    const double frequency = 1.60217733000103e-19*9.6*2.*std::acos(-1.)/6.626068e-34;
    const double step_ps = 1.e12/frequency;
    const auto default_cadence = automatic_cadence(0.1, frequency, 1.);
    check(default_cadence.steps == 1459);
    check(default_cadence.actual_ps >= 0.1 && default_cadence.actual_ps < 0.1 + step_ps);
    const auto tiny = automatic_cadence(step_ps/10., frequency, 1.);
    check(tiny.steps == 1 && std::abs(tiny.actual_ps-step_ps) < 1.e-18);
    // Exact integral-step inputs do not gain an extra step.
    check(automatic_cadence(0.125, 8.e12, 1.).steps == 1);
    check(automatic_cadence(0.25, 8.e12, 1.).steps == 2);
    for (double invalid : {0., -1., std::numeric_limits<double>::infinity(),
                           std::numeric_limits<double>::quiet_NaN()}) {
        for (int position = 0; position < 3; ++position) {
            bool threw = false;
            try {
                automatic_cadence(position == 0 ? invalid : 0.1,
                                  position == 1 ? invalid : frequency,
                                  position == 2 ? invalid : 1.);
            } catch (const std::invalid_argument&) { threw = true; }
            check(threw);
        }
    }
    std::feclearexcept(FE_ALL_EXCEPT);
    for (double huge : {1.e100, std::numeric_limits<double>::max()}) {
        bool threw = false;
        try { automatic_cadence(huge, frequency, 1.); }
        catch (const std::overflow_error&) { threw = true; }
        check(threw);
    }
    check(!std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW));
    check(!cup_output_detail::scheduled(0, 3));
    check(!cup_output_detail::scheduled(2, 3));
    check(cup_output_detail::scheduled(3, 3));
    check(cup_output_detail::scheduled(6, 3));
    std::cout << "Cadence failures=" << failures << ", default actual interval="
              << std::setprecision(17) << default_cadence.actual_ps << " ps\n";
    return failures ? 1 : 0;
}
