#define CUP_BACKEND_QUASI_STATIC
#include <chrono>
#include <filesystem>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/single.hpp>
#include "cup.hpp"

// Exercise a growing solution beyond the legacy amplitude cutoff without overflow.
int main()
{
    namespace fs = std::filesystem;
    const auto dir = fs::temp_directory_path() /
        ("ngm-qs-duration-" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    fs::create_directories(dir);
    struct Cleanup {
        fs::path path;
        ~Cleanup() { fs::remove_all(path); }
    } cleanup{dir};

    nanosphere s{};
    s.init();
    s.set_metal("silver", "drude");
    s.set_active("lorentz");
    s.G = 1.;
    s.Dome = 0.1;
    s.ome_g = 3.5;
    char model[] = "drude", metal[] = "silver", host[] = "vacuum";
    constexpr double field = 1., duration = 0.2;
    const double frequency = eV2j*s.Ome_p*2.*M_PI/h;
    const double cadence = 1.e12/frequency;
    const auto run = [&](double requested) {
        return s.analytical(model, metal, host, field, 3.5, requested,
                            0., nullptr, 0., dir.c_str());
    };
    const auto fixed = run(duration);
    bool ok = std::isfinite(fixed.real()) && std::isfinite(fixed.imag()) &&
              std::isfinite(std::norm(fixed)) && std::norm(fixed)/field > 5.e26;
    double last = -1., first_crossing = -1.;
    std::size_t rows = 0;
    for (const char* file : {"anltime.dat", "anlfunc.dat"}) {
        std::ifstream input(dir / file);
        std::string line;
        double endpoint = -1.;
        std::size_t count = 0;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            std::vector<double> values;
            std::string token;
            while (row >> token) {
                const double value = std::stod(token);
                ok = ok && std::isfinite(value);
                values.push_back(value);
            }
            const bool dynamics = std::string(file) == "anltime.dat";
            if (values.size() != (dynamics ? 5u : 8u)) return 1;
            endpoint = values[0];
            ++count;
            if (dynamics && first_crossing < 0. &&
                (values[3]*values[3] + values[4]*values[4])/field > 5.e26)
                first_crossing = endpoint;
        }
        // Output uses six significant digits; allow its rounding as well as dt.
        ok = ok && count > 0 && std::abs(endpoint-duration) <= cadence + 1.e-6;
        last = endpoint;
        rows = count;
    }
    ok = ok && first_crossing > 0. && first_crossing < duration - cadence;
    const auto automatic = run(0.);
    const double cutoff_value = std::norm(automatic)/field;
    // The returned state must remain just beyond the existing amplitude gate,
    // independently of the automatic output cadence.
    ok = ok && std::isfinite(cutoff_value) && cutoff_value > 5.e26 &&
         cutoff_value < 6.e26;
    std::cout << "fixed endpoint=" << last << " ps, rows=" << rows
              << ", first amplitude crossing=" << first_crossing << " ps"
              << ", final norm=" << std::norm(fixed)
              << ", automatic cutoff norm=" << cutoff_value << '\n';
    return ok ? 0 : 1;
}
