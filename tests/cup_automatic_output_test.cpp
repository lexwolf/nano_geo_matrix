#include <chrono>
#include <filesystem>
#include <map>
#include <stdexcept>
#ifdef CUP_TEST_MIE
#define CUP_BACKEND_MIE
#include <nano_geo_matrix/bessel/myBessel_bl.hpp>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/mie/geometry/gimme_p.hpp>
#include <nano_geo_matrix/mie/geometry/single.hpp>
#else
#define CUP_BACKEND_QUASI_STATIC
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/single.hpp>
#endif
#include "cup.hpp"

namespace {
namespace fs = std::filesystem;
void require(bool ok, const std::string& message)
{
    if (!ok) throw std::runtime_error(message);
}
struct File {
    std::string text;
    std::map<std::string, std::string> metadata;
    std::vector<std::vector<double>> rows;
};
File read(const fs::path& path, std::size_t columns)
{
    File result;
    std::ifstream input(path);
    require(bool(input), "Missing output: " + path.string());
    std::string line;
    while (std::getline(input, line)) {
        result.text += line + '\n';
        if (line.empty()) continue;
        if (line[0] == '#') {
            const auto split = line.find(" = ");
            if (split != std::string::npos)
                result.metadata[line.substr(2, split-2)] = line.substr(split+3);
            continue;
        }
        std::istringstream values(line);
        std::vector<double> row;
        std::string token;
        while (values >> token) {
            const double value = std::stod(token);
            require(std::isfinite(value), "Nonfinite output: " + path.string());
            row.push_back(value);
        }
        require(row.size() == columns, "Changed columns: " + path.string());
        if (!result.rows.empty())
            require(row[0] > result.rows.back()[0], "Duplicate or decreasing timestamp");
        result.rows.push_back(row);
    }
    require(!result.rows.empty(), "No data rows: " + path.string());
    return result;
}
struct Run {
    std::complex<double> state;
    std::vector<File> files;
};

Run run(const fs::path& directory, bool analytical, double duration,
        double interval = 0.1, double gain = .01, int order = 1,
        double energy_scale = 1., bool use_default = false, double pump_time = 0.)
{
    fs::create_directories(directory / "num");
    fs::create_directories(directory / "anl");
    nanosphere s{};
    s.init();
    s.set_metal("silver", "drude");
    s.set_active("lorentz");
    s.G = gain; s.Dome = .1*energy_scale; s.ome_g = 3.5*energy_scale;
    s.Ome_p *= energy_scale;
    s.ome_p *= energy_scale;
    s.ome_p2 = s.ome_p*s.ome_p;
    s.Gam_d *= energy_scale;
    // Scaling energy and radius inversely keeps normalized Mie geometry intact.
    s.a = (order == 1 ? 10. : 50.)/energy_scale;
    s.nOme = 0.;
    char model[] = "drude", metal[] = "silver", host[] = "vacuum";
    const double frequency = eV2j*s.Ome_p*2.*M_PI/h;
    const double step_ps = 1.e12/frequency;
    Run result;
#ifdef CUP_TEST_MIE
#define ORDER_ARG order,
#else
#define ORDER_ARG
    (void)order;
#endif
    if (analytical) {
        if (use_default)
            result.state = s.analytical(ORDER_ARG model, metal, host, 1., 3.5*energy_scale,
                                       duration, pump_time, nullptr, 0., directory.c_str());
        else
            result.state = s.analytical(ORDER_ARG model, metal, host, 1., 3.5*energy_scale,
                                       duration, pump_time, nullptr, 0., directory.c_str(), interval);
    } else {
        if (use_default)
            result.state = s.numerical(ORDER_ARG model, metal, host, 1., 3.5*energy_scale,
                                      duration, pump_time, nullptr, 0., "step", false, directory.c_str());
        else
            result.state = s.numerical(ORDER_ARG model, metal, host, 1., 3.5*energy_scale,
                                      duration, pump_time, nullptr, 0., "step", false, directory.c_str(), interval);
    }
    require(std::isfinite(result.state.real()) && std::isfinite(result.state.imag()), "Nonfinite returned state");
    const std::string prefix = analytical ? "anl" : "num";
#ifdef CUP_TEST_MIE
    if (order == 1) result.files.push_back(read(directory / (prefix+"time.dat"), 3));
    const std::string suffix = order == 1 ? "01.dat" : "02.dat";
    result.files.push_back(read(directory / prefix / ("a"+suffix), 3));
    result.files.push_back(read(directory / prefix / ("func"+suffix), 11));
#else
    result.files.push_back(read(directory / (prefix+"time.dat"), analytical ? 5 : 9));
    result.files.push_back(read(directory / (prefix+"func.dat"), 8));
#endif
    for (const auto& file : result.files) {
        require(file.rows.size() == result.files[0].rows.size(), "Files have different row counts");
        for (std::size_t j = 0; j < file.rows.size(); ++j)
            require(file.rows[j][0] == result.files[0].rows[j][0], "Files have different timestamps");
        if (duration == 0.) {
            const auto cadence = cup_output_detail::automatic_cadence(interval, frequency, 1.);
            require(std::stod(file.metadata.at("requested_output_interval_ps")) == interval, "Wrong requested interval");
            require(std::stod(file.metadata.at("actual_output_interval_ps")) == cadence.actual_ps, "Wrong rounded interval");
            require(std::stoi(file.metadata.at("output_interval_steps")) == cadence.steps, "Wrong step count");
            const auto terminal_step = std::llround(file.rows.back()[0]/step_ps);
            const auto scheduled_rows = terminal_step/cadence.steps;
            const auto expected_rows = scheduled_rows + (terminal_step % cadence.steps != 0 || terminal_step == 0 ? 1 : 0);
            require(file.rows.size() == static_cast<std::size_t>(expected_rows), "Missing or duplicate terminal row");
            for (long long j = 0; j < scheduled_rows; ++j) {
                const double expected_time = (j+1)*cadence.actual_ps;
                require(std::abs(file.rows[j][0]-expected_time) < 5.e-6*std::max(1., expected_time), "Wrong regular cadence");
            }
            require(!file.metadata.at("termination_reason").empty(), "Missing termination reason");
        } else {
            require(file.metadata.count("actual_output_interval_ps") == 0 &&
                    file.metadata.count("termination_reason") == 0, "Fixed-run metadata changed");
            const int steps = static_cast<int>(duration*frequency*1.e-12);
            const int cap = analytical ? 500000 :
#ifdef CUP_TEST_MIE
                500000;
#else
                50000;
#endif
            const int stride = std::max(1, steps/cap);
            require(file.rows.size() == static_cast<std::size_t>((steps+1)/stride), "Fixed row count changed");
            for (std::size_t j = 0; j < file.rows.size(); ++j) {
                const double expected_time = ((j+1)*stride-1)*step_ps;
                require(std::abs(file.rows[j][0]-expected_time) < 5.e-6*std::max(1., expected_time), "Fixed row scheduling changed");
            }
        }
    }
    return result;
}
void same_stop(const Run& a, const Run& b)
{
    require(a.state == b.state, "Output spacing changed returned state");
    for (std::size_t j = 0; j < a.files.size(); ++j)
        require(a.files[j].rows.back() == b.files[j].rows.back(), "Output spacing changed terminal state/time");
    require(a.files[0].metadata.at("termination_reason") == b.files[0].metadata.at("termination_reason"), "Output spacing changed termination reason");
}
}

int main()
{
    const auto directory = fs::temp_directory_path() /
        ("ngm-auto-output-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    struct Cleanup { fs::path path; ~Cleanup() { fs::remove_all(path); } } cleanup{directory};
    try {
        const double step_ps = 1.e12/(eV2j*9.6*2.*M_PI/h);
        for (bool analytical : {false, true}) {
            const std::string name = analytical ? "analytical" : "numerical";
            bool rejected = false;
            try { run(directory/(name+"-invalid"), analytical, 0., 0.); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Solver accepted zero automatic interval");
            rejected = false;
            try { run(directory/(name+"-overflow"), analytical, 0., 1.e100); }
            catch (const std::overflow_error&) { rejected = true; }
            require(rejected, "Solver accepted overflowing automatic interval");
            const auto normal = run(directory/(name+"-default"), analytical, 0., .1, .01, 1, 1., true);
            const auto alternate = run(directory/(name+"-alternate"), analytical, 0., .07);
            same_stop(normal, alternate);
            const auto every_step = run(directory/(name+"-every-step"), analytical, 0., step_ps/2.);
            same_stop(normal, every_step);
            require(normal.files[0].metadata.at("termination_reason") ==
                    (analytical ? "existing_instantaneous_intensity_convergence" : "existing_block_average_convergence"), "Unexpected convergence criterion");
            if (!analytical) {
                const auto delayed = run(directory/(name+"-spacing"), false, 0., .1, .01, 1, 1., false, .2);
                require(delayed.files[0].rows.size() >= 4, "Too few default-cadence rows to verify spacing");
            }
            const auto fixed = run(directory/(name+"-fixed"), analytical, .003);
            const auto fixed_other = run(directory/(name+"-fixed-other"), analytical, .003, -1.);
            require(fixed.state == fixed_other.state, "Fixed state changed");
            for (std::size_t j = 0; j < fixed.files.size(); ++j)
                require(fixed.files[j].text == fixed_other.files[j].text, "Fixed output changed with unused interval");
            // Scale all energies down to exhaust the unchanged 100000 ps ceiling
            // in 15 iterations, preserving normalized equations and Mie geometry.
            const auto exhausted = run(directory/(name+"-maximum"), analytical, 0., 15000., .01, 1, 1.e-8);
            require(exhausted.files[0].metadata.at("termination_reason") == "maximum_duration_exhaustion", "Maximum duration reason missing");
            const double slow_step = step_ps/1.e-8;
            require(std::abs(exhausted.files[0].rows.back()[0] - std::floor(100000./slow_step)*slow_step) < 1., "Maximum terminal time incorrect");
#ifdef CUP_TEST_MIE
            const auto higher_order = run(directory/(name+"-order2"), analytical, 0., .1, .01, 2);
            require(higher_order.files.size() == 2, "Missing order-specific output");
#endif
            std::cout << name << ": terminal=" << normal.files[0].rows.back()[0]
                      << " ps, reason=" << normal.files[0].metadata.at("termination_reason")
                      << ", rows=" << normal.files[0].rows.size() << '\n';
        }
        const auto early = run(directory/"early-amplitude", true, 0., .1, 4.);
        require(early.files[0].metadata.at("termination_reason") == "amplitude_cutoff", "Expected amplitude cutoff");
        require(early.files[0].rows.size() == 1 && early.files[0].rows[0][0] < .1, "Early terminal row not saved");
        std::cout << "early analytical amplitude cutoff=" << early.files[0].rows[0][0] << " ps\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
