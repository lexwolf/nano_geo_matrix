#define CUP_BACKEND_QUASI_STATIC
#include <chrono>
#include <filesystem>
#include <nano_geo_matrix/core/mathNN.hpp>
#include <nano_geo_matrix/quasi_static/geometry/nanoshell.hpp>
#include "cup.hpp"

int main() {
    namespace fs = std::filesystem;
    const auto base = fs::temp_directory_path() / ("ngm-shell-convergence-"+
        std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    struct Cleanup { fs::path path; ~Cleanup() { fs::remove_all(path); } } cleanup{base};
    for (bool analytical : {false, true}) {
        const auto out = base/(analytical ? "analytical" : "numerical"); fs::create_directories(out);
        nanosphere s{}; s.init(); s.set_metal("silver", "drude"); s.set_active("lorentz");
        s.G = .01; s.Dome = .1; s.ome_g = 3.5; s.a = 10.; s.nOme = 0.;
        char model[] = "drude", metal[] = "silver", host[] = "vacuum", solvent[] = "vacuum";
        const auto state = analytical
            ? s.analytical(model, metal, host, 1.e-12, 3.5, 0., .5, solvent, .5, out.c_str())
            : s.numerical(model, metal, host, 1.e-12, 3.5, 0., .5, solvent, .5, "step", true, out.c_str());
        if (!std::isfinite(state.real()) || !std::isfinite(state.imag())) return 1;
        std::ifstream file(out/(analytical ? "anltime.dat" : "numtime.dat"));
        std::string line; double start = -1., end = -1.; bool converged = false;
        while (std::getline(file, line)) {
            if (line.find("# convergence_start_ps = ") == 0) start = std::stod(line.substr(line.find('=')+1));
            if (line.find("# convergence_confirmation_end_ps = ") == 0) end = std::stod(line.substr(line.find('=')+1));
            if (line == "# termination_reason = sustained_flat") converged = true;
        }
        if (!converged || start < .5 || end-start < 39.999) return 1;
        std::cout << (analytical ? "analytical" : "numerical") << " nanoshell weak seed E0=1e-12: ["
                  << start << ',' << end << "] ps, norm=" << std::norm(state) << '\n';
    }
    return 0;
}
