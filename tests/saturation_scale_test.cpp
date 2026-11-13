#include <nano_geo_matrix/core/saturation_scale.hpp>

#include <iostream>
#include <limits>

using namespace nano_geo_matrix::saturation;

int main() {
    int failures = 0;
    const auto check = [&](bool ok, const char* message) {
        if (!ok) { ++failures; std::cerr << message << '\n'; }
    };
    const auto near = [](double a, double b) {
        return std::abs(a-b) <= 2.e-12 * std::max(std::abs(a), std::abs(b));
    };
    const auto rejects = [&](Inputs input) {
        try { calculate(input); } catch (const std::invalid_argument&) { return true; }
        return false;
    };
    const Inputs silver{9.6, .15, .13625607};
    const auto scale = calculate(silver);
    check(near(scale.tau2_normalized, 128.), "silver tau2");
    check(near(scale.tau1_normalized, 640.), "silver tau1");
    check(near(scale.Isat_normalized, .011467379031260773), "audited silver saturation intensity");
    check(near(scale.Esat_normalized*scale.Esat_normalized, scale.Isat_normalized), "field/intensity relationship");
    check(scale.inputs.tau1_policy == Tau1Policy::require_match && !scale.mismatch_overridden,
          "default authority/provenance");

    // Independent physical construction, rather than repeating the eliminated formula.
    constexpr double hbar=1.054571817e-34, ev=1.602176634e-19, eps0=8.8541878128e-12;
    for (double plasma : {9.6, 3.2}) {
        for (double gain : {.03, .13625607, -.2}) {
            for (double ratio : {2., 5., 8.}) {
                Inputs in{plasma, .15, gain, ratio};
                const auto s=calculate(in);
                for (double mu : {10.*3.33564e-30, 20.*3.33564e-30}) {
                    const double tau2_si=2.*hbar/(in.linewidth_eV*ev), tau1_si=ratio*tau2_si;
                    const double density=3.*hbar*eps0*std::abs(gain)/(mu*mu*tau2_si);
                    const double unit=std::sqrt(density*hbar*(plasma*ev/hbar)/eps0);
                    const double esat_si=hbar/mu*std::sqrt(3./(tau1_si*tau2_si));
                    check(near(s.Esat_normalized, esat_si/unit), "independent SI field round trip");
                }
                check(near(s.effective_tau1_over_tau2, ratio), "custom relaxation ratio");
            }
        }
    }
    Inputs supplied=silver;
    supplied.supplied_tau1_normalized=640.;
    check(calculate(supplied).supplied_tau1_matches, "matching tau1 accepted");
    supplied.supplied_tau1_normalized=640.*(1.+.5e-12);
    check(calculate(supplied).supplied_tau1_matches, "roundoff tolerated");
    supplied.supplied_tau1_normalized=640.*(1.+4.e-12);
    check(rejects(supplied), "meaningful mismatch rejected");
    // The historical missing-plasma conversion must fail without an explicit override.
    supplied.supplied_tau1_normalized=10./silver.linewidth_eV;
    check(rejects(supplied), "legacy tau1 rejected by default");
    supplied.tau1_policy=Tau1Policy::energy_inputs_authoritative;
    const auto energies=calculate(supplied);
    check(near(energies.Isat_normalized, scale.Isat_normalized) && energies.mismatch_overridden &&
          !energies.supplied_tau1_matches && energies.inputs.supplied_tau1_normalized==supplied.supplied_tau1_normalized,
          "energy authority retains conflicting supplied value");
    supplied.tau1_policy=Tau1Policy::supplied_tau1_authoritative;
    const auto override=calculate(supplied);
    check(near(override.Isat_normalized/scale.Isat_normalized, 9.6), "explicit legacy authority reproduces factor");
    check(near(override.energy_derived_tau1_normalized, 640.) && override.mismatch_overridden &&
          near(override.effective_tau1_over_tau2, 5./9.6), "override records changed effective relaxation ratio");
    check(std::string(policy_name(override.inputs.tau1_policy))=="supplied_tau1_authoritative", "serializable policy name");

    const double inf=std::numeric_limits<double>::infinity(), nan=std::numeric_limits<double>::quiet_NaN();
    for (double bad : {0., -1., inf, nan}) {
        Inputs in=silver; in.plasma_energy_eV=bad; check(rejects(in), "invalid plasma");
        in=silver; in.linewidth_eV=bad; check(rejects(in), "invalid linewidth");
        in=silver; in.tau1_over_tau2=bad; check(rejects(in), "invalid ratio");
        in=silver; in.supplied_tau1_normalized=bad; check(rejects(in), "invalid supplied tau1");
    }
    for (double bad : {0., inf, -inf, nan}) {
        Inputs in=silver; in.gain=bad; check(rejects(in), "invalid gain");
    }
    for (auto policy : {Tau1Policy::energy_inputs_authoritative, Tau1Policy::supplied_tau1_authoritative}) {
        Inputs in=silver; in.tau1_policy=policy; check(rejects(in), "override requires supplied value");
    }
    Inputs unknown=silver; unknown.tau1_policy=static_cast<Tau1Policy>(99);
    check(rejects(unknown), "unknown policy rejected");
    const auto range_rejected = [](Inputs in) {
        try { calculate(in); } catch (const std::overflow_error&) { return true; }
        return false;
    };
    Inputs extreme=silver;
    extreme.plasma_energy_eV=std::numeric_limits<double>::max();
    extreme.linewidth_eV=std::numeric_limits<double>::min();
    check(range_rejected(extreme), "unrepresentable relaxation rejected");
    extreme=silver; extreme.gain=std::numeric_limits<double>::denorm_min();
    check(range_rejected(extreme), "unrepresentable saturation rejected");
    // Multiplication would overflow in double even though the final scale is representable.
    extreme={1., 2., 1.e200, 1.e120};
    const auto tiny=calculate(extreme);
    check(tiny.Isat_normalized>0. && std::isfinite(tiny.Esat_normalized), "intermediate overflow avoided");
    extreme.tau1_over_tau2=1.e200;
    check(range_rejected(extreme), "underflow to zero rejected");
    if (!failures) std::cout << "Saturation scales, SI cross-checks, authority policies and invalid inputs: PASS\n";
    return failures ? 1 : 0;
}
