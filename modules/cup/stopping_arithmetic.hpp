#pragma once

#include <cmath>

namespace cup_stopping_detail {

inline bool finite_positive(double value)
{
    return std::isfinite(value) && value > 0.;
}

inline bool relative_change_below(double current, double previous, double tolerance)
{
    return finite_positive(current) && std::isfinite(previous) &&
           std::fabs((current - previous)/current) < tolerance;
}

inline bool amplitude_above(double intensity, double field, double threshold)
{
    return finite_positive(field) && intensity/field > threshold;
}

inline void invalidate_average_if_needed(double frequency, double elapsed, double total,
                                         double& current, double& previous)
{
    if (!finite_positive(frequency) || !finite_positive(elapsed) ||
        !std::isfinite(total)) {
        current = previous = 0.;
    }
}

inline bool averaging_window_complete(double frequency, double elapsed)
{
    return finite_positive(frequency) && finite_positive(elapsed) &&
           elapsed >= 1000/frequency;
}

inline void store_average(double total, double elapsed, double& current, double& previous)
{
    if (finite_positive(elapsed) && std::isfinite(total)) {
        const double average = total/elapsed;
        if (std::isfinite(average)) {
            previous = current;
            current = average;
            return;
        }
    }
    // Zero is unavailable to relative_change_below; discard stale averages too.
    current = previous = 0.;
}

} // namespace cup_stopping_detail
