#include "convergence.hpp"
#include <functional>
#include <iostream>
#include <string>

using cup_convergence::Monitor;
using cup_convergence::Options;
int failures = 0;
void check(bool ok, const char* message) { if (!ok) { ++failures; std::cerr << message << '\n'; } }
struct TraceResult { bool accepted; cup_convergence::Result result; std::size_t peak_storage; };
TraceResult trace(const std::function<double(double)>& intensity, double end = 120.,
                  Options options = {}, bool periodic = true, double dt = .02) {
    Monitor m(options, periodic); std::size_t storage = 0;
    for (int i = 0; i*dt <= end; ++i) {
        const bool done = m.update(i*dt, intensity(i*dt));
        storage = std::max(storage, m.stored_records());
        if (done) return {true, m.result(), storage};
    }
    return {false, {}, storage};
}
int main() {
    const double pi = std::acos(-1.);
    for (double level : {0., 1., 1.e-30, 1.e30}) {
        auto r = trace([=](double) { return level; });
        check(r.accepted && std::string(r.result.classification) == "sustained_flat", "constant not flat");
        check(r.result.qualifying_end_ps == 30. && r.result.confirmation_end_ps == 40., "flat confirmation interval incorrect");
        check(r.result.mean_intensity == level && r.result.peak_to_peak == 0., "flat statistics incorrect");
    }
    auto relaxed = trace([](double t) { return 2.+std::exp(-t/2.); });
    check(relaxed.accepted && relaxed.result.start_ps >= 20., "relaxation accepted prematurely");
    auto periodic = trace([=](double t) { return 10.+2.*std::sin(2.*pi*t/4.); });
    check(periodic.accepted && std::string(periodic.result.classification) == "sustained_periodic", "stable periodic not accepted");
    check(periodic.result.qualifying_end_ps-periodic.result.start_ps >= 30. &&
          periodic.result.confirmation_end_ps-periodic.result.qualifying_end_ps >= 7.99,
          "periodic evidence or confirmation too short");
    check(std::abs(periodic.result.mean_period_ps-4.) < .001 &&
          std::abs(periodic.result.mean_intensity-10.) < .001 &&
          std::abs(periodic.result.mean_amplitude-4.) < .001, "periodic statistics incorrect");
    auto background = trace([=](double t) { return 1.e6+80.*std::sin(2.*pi*t/4.); });
    check(background.accepted && std::string(background.result.classification) == "sustained_periodic", "small resolved oscillation on large background lost");
    for (double sign : {-1., 1.}) {
        auto evolving = trace([=](double t) { return 10.+std::exp(sign*t/100.)*std::sin(2.*pi*t/4.); });
        check(!evolving.accepted, "evolving amplitude qualified prematurely");
    }
    check(!trace([=](double t) { return 10.+.01*t+std::sin(2.*pi*t/4.); }).accepted, "slow mean drift accepted");
    check(!trace([=](double t) { return 10.+std::sin(2.*pi*(t/4.+.0005*t*t)); }).accepted, "changing period accepted");
    check(!trace([](double t) { return t < 35. ? 1. : 1.+(t-35.); }, 80.).accepted, "confirmation failed to reject temporary plateau");
    auto tiny = trace([=](double t) { return 10.+1.e-5*std::sin(31.*t)+1.e-6*std::sin(57.*t); });
    check(tiny.accepted && std::string(tiny.result.classification) == "sustained_flat", "sub-resolution extrema treated as cycles");
    check(!trace([=](double t) { return 10.+std::sin(2.*pi*t/4.); }, 29.).accepted, "insufficient periodic duration accepted");
    Options many; many.periodic_cycles = 50;
    check(!trace([=](double t) { return 10.+std::sin(2.*pi*t/4.); }, 80., many).accepted, "insufficient cycles accepted");
    auto plateaus = trace([](double t) {
        const double phase = std::fmod(t, 4.);
        return phase < 1. ? 3. : phase < 2. ? 3.-2.*(phase-1.) : phase < 3. ? 1. : 1.+2.*(phase-3.);
    });
    check(plateaus.accepted && std::string(plateaus.result.classification) == "sustained_periodic", "plateau extrema not deterministic");
    Monitor invalid;
    for (int t = 0; t < 39; ++t) invalid.update(t, 1.);
    check(!invalid.update(39., std::numeric_limits<double>::quiet_NaN()), "nonfinite sample accepted");
    for (int t = 40; t < 80; ++t) check(!invalid.update(t, 1.), "invalid evidence reused");
    check(invalid.update(80., 1.), "fresh evidence after NaN not accepted");
    check(!invalid.update(81., std::numeric_limits<double>::infinity()) && invalid.stored_records() == 0, "infinity did not invalidate accepted evidence");
    check(!invalid.update(82., -1.), "negative intensity accepted");
    Monitor irregular;
    for (double t : {0., .7, 7.3, 10.2, 21., 39.9}) check(!irregular.update(t, 2.), "irregular timestamps shortened windows");
    check(irregular.update(40.2, 2.) && irregular.result().confirmation_end_ps == 40., "window boundary interpolation failed");
    Options absolute; absolute.relative_flatness = 0.; absolute.absolute_flatness = .01;
    check(trace([=](double t) { return .001+.0001*std::sin(t); }, 50., absolute).accepted, "absolute intensity tolerance not used");
    Options collective; collective.flat_window_ps = 1.; collective.absolute_flatness = 1.1; collective.relative_flatness = 0.;
    check(!trace([](double t) { return t; }, 20., collective, false).accepted, "collective range not enforced");
    Options weighted = collective; weighted.absolute_flatness = 10.;
    Monitor linear(weighted, false);
    for (double t : {0., .3, 1.4, 2.1, 3.2, 4.2}) linear.update(t, t);
    check(linear.result().mean_intensity == 2. && linear.result().confirmation_end_ps == 4., "means not time weighted at exact boundaries");
    Options bounded; bounded.max_periodic_cycles = 16;
    auto rapid = trace([=](double t) { return 10.+std::sin(2.*pi*t/.1); }, 80., bounded, true, .005);
    check(!rapid.accepted && rapid.peak_storage <= 16+4, "storage bound or insufficient cycle span violated");
    // Nanoshell-like weak seed: no whole-run maximum hides exponential growth.
    auto weak = trace([](double t) { const double x = 1.e-12*std::exp(t/2.); return x/(1.+x); }, 160.);
    check(weak.accepted && weak.result.start_ps > 50., "weak delayed seed accepted prematurely");
    auto delayed_cycles = trace([=](double t) {
        const double envelope = t < 20. ? std::exp((t-60.)/3.) : std::min(1., std::exp((t-60.)/3.));
        return envelope*(2.+std::sin(2.*pi*t/5.));
    }, 160.);
    check(delayed_cycles.accepted && std::string(delayed_cycles.result.classification) == "sustained_periodic" &&
          delayed_cycles.result.start_ps >= 60., "delayed periodic envelope accepted prematurely");
    for (double scale : {1.e-12, 1.e12}) {
        const auto scaled = trace([=](double t) { return scale*(10.+std::sin(2.*pi*t/4.)); });
        check(scaled.accepted && std::string(scaled.result.classification) == "sustained_periodic", "periodic scale dependence");
    }
    check(!trace([=](double t) { return 10.+.001*t+std::sin(2.*pi*t/4.); }).accepted,
          "maxima/minima spread ignored below cycle-mean tolerance");
    check(!trace([](double t) {
        const double duty = .2+.005*t;
        return std::fmod(t, 4.) < duty*4. ? 2. : 1.;
    }, 100.).accepted, "changing duty-cycle means accepted");
    const auto failed_periodic_confirmation = trace([=](double t) {
        return 10.+(t < 35. ? 1. : 2.)*std::sin(2.*pi*t/4.);
    });
    check(failed_periodic_confirmation.accepted && failed_periodic_confirmation.result.start_ps >= 35.,
          "failed periodic confirmation retained old evidence");
    Monitor broken_cycle;
    bool completed = false;
    for (int i = 0; i <= 5000 && !completed; ++i) {
        const double t = i*.02;
        completed = broken_cycle.update(t, i == 1500 ? std::numeric_limits<double>::quiet_NaN() : 10.+std::sin(2.*pi*t/4.));
    }
    check(completed && broken_cycle.result().start_ps > 30., "periodic evidence crossed nonfinite gap");
    check(!trace([=](double t) { return 10.+std::sin(2.*pi*t/4.); }, 120., {}, false).accepted,
          "analytical flat-only monitor accepted periodicity");
    const auto changed_scale = trace([=](double t) {
        return t < 20. ? 1.e200 : 1.e-200*(10.+std::sin(2.*pi*t/4.));
    });
    check(changed_scale.accepted && std::abs(changed_scale.result.mean_intensity/1.e-199-1.) < .001,
          "cycle integrals contaminated by large prehistory");
    for (int field = 0; field < 16; ++field) {
        Options bad;
        switch (field) {
            case 0: bad.flat_window_ps = 0.; break;
            case 1: bad.relative_flatness = 0.; break;
            case 2: bad.absolute_flatness = -1.; break;
            case 3: bad.relative_period = std::numeric_limits<double>::infinity(); break;
            case 4: bad.periodic_cycles = 1; break;
            case 5: bad.flat_confirmation_windows = -1; break;
            case 6: bad.max_periodic_cycles = 5; break;
            case 7: bad.periodic_min_span_ps = std::numeric_limits<double>::quiet_NaN(); break;
            case 8: bad.flat_windows = 0; break;
            case 9: bad.relative_flatness = std::numeric_limits<double>::quiet_NaN(); break;
            case 10: bad.relative_amplitude = 0.; break;
            case 11: bad.relative_cycle_mean = -1.; break;
            case 12: bad.absolute_flatness = std::numeric_limits<double>::infinity(); break;
            case 13: bad.max_periodic_cycles = 65537; break;
            case 14: bad.flat_confirmation_windows = 4096; break;
            case 15: bad.periodic_confirmation_cycles = std::numeric_limits<int>::max(); break;
        }
        bool threw = false; try { Monitor m(bad); } catch (const std::invalid_argument&) { threw = true; }
        check(threw, "invalid options accepted");
    }
    std::cout << "Detector failures=" << failures << "; flat=[0,30]+confirmation to 40 ps; periodic=["
              << periodic.result.start_ps << ',' << periodic.result.qualifying_end_ps << "]+confirmation to "
              << periodic.result.confirmation_end_ps << "; weak-seed start=" << weak.result.start_ps << " ps\n";
    return failures ? 1 : 0;
}
