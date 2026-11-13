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
    // Use the same row prefix writer as all four solvers, then append value
    // columns at the stream's original precision and width.
    const auto row = [&](double time, bool automatic, int precision) {
        std::ostringstream output;
        output.precision(precision);
        cup_output_detail::write_timestamp(output, time, automatic);
        check(output.precision() == precision);
        output << '\t' << std::setw(13) << 1.23456789012345
               << '\t' << std::setw(13) << -9.87654321098765;
        return output.str();
    };
    const auto legacy_row = [](double time, int precision) {
        std::ostringstream output;
        output.precision(precision);
        output << "  " << std::setw(8) << std::setiosflags(std::ios::left) << time
               << '\t' << std::setw(13) << 1.23456789012345
               << '\t' << std::setw(13) << -9.87654321098765;
        return output.str();
    };
    const auto parsed_time = [](const std::string& text) {
        return std::stod(text.substr(0, text.find('\t')));
    };
    const auto check_pair = [&](double scheduled_time, double terminal_time) {
        const auto first = row(scheduled_time, true, 6);
        const auto last = row(terminal_time, true, 6);
        check(first != last);
        check(parsed_time(first) == scheduled_time);
        check(parsed_time(last) == terminal_time);
        check(parsed_time(first) < parsed_time(last));
    };
    const int scheduled_step = 101*default_cadence.steps;
    check(cup_output_detail::scheduled(scheduled_step, default_cadence.steps));
    check(!cup_output_detail::scheduled(scheduled_step+1, default_cadence.steps));
    // Exactly the timestamp calculation used in the integration loops.
    const double scheduled_time = scheduled_step*1.e12/frequency;
    const double terminal_time = (scheduled_step+1)*1.e12/frequency;
    check(std::abs(scheduled_time-10.103477308065345) < 1.e-12);
    check(std::abs(terminal_time-10.10354587175883) < 1.e-12);
    check(parsed_time(legacy_row(scheduled_time, 6)) ==
          parsed_time(legacy_row(terminal_time, 6))); // Previously both 10.1035.
    check_pair(scheduled_time, terminal_time);
    const int large_step = 999000*default_cadence.steps;
    check_pair(large_step*1.e12/frequency, (large_step+1)*1.e12/frequency);
    // A sub-timestep request saves every integration step. Also check adjacent
    // representable doubles at the ceiling, where six digits lose even more.
    check(tiny.steps == 1);
    check_pair(large_step*1.e12/frequency, (large_step+tiny.steps)*1.e12/frequency);
    check_pair(100000., std::nextafter(100000., std::numeric_limits<double>::infinity()));
    for (double time : {0., scheduled_time, terminal_time, 100000.}) {
        for (int precision : {6, 9}) {
            const auto automatic = row(time, true, precision);
            const auto legacy = legacy_row(time, precision);
            check(automatic.substr(automatic.find('\t')) == legacy.substr(legacy.find('\t')));
            check(row(time, false, precision) == legacy); // Fixed-duration bytes.
        }
    }
    std::cout << "Cadence failures=" << failures << ", default actual interval="
              << std::setprecision(17) << default_cadence.actual_ps << " ps\n";
    return failures ? 1 : 0;
}
