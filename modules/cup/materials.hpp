/*
 * cup/materials.H
 * Shared material/permittivity routines for cup.H and cupN.H.
 *
 * IMPORTANT: this header must be included AFTER class nanosphere is declared.
 */
#pragma once
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

inline void nanosphere::set_metal(const char* mtl_cstr, const char* mdl_cstr, int sel)
{
    namespace fs = std::filesystem;
    using std::string_view;

    const string_view mtl = (mtl_cstr != nullptr) ? string_view{mtl_cstr} : string_view{};
    const string_view mdl = (mdl_cstr != nullptr) ? string_view{mdl_cstr} : string_view{};

    // store selection:
    // -1: LOW (Re-Dr, Im-Di), 0: BASE, +1: HIGH (Re+Dr, Im+Di)
    const int s = (sel > 0) ? +1 : (sel < 0 ? -1 : 0);

    auto fail = [](const std::string& msg) -> void {
        std::cerr << "Error: " << msg << '\n';
        std::exit(EXIT_FAILURE);
    };

    auto warn = [](const std::string& msg) -> void {
        std::cerr << "Warning: " << msg << '\n';
    };

    struct metal_spec {
        string_view name;
        double Ome_p_ev;
        double Gam_d_ev;
        double eps_inf_val;
        const char* jc_filename;
        const char* allowed_models;
    };

    static constexpr metal_spec gold{
        "gold",
        8.91,
        0.0759,
        9.0685,
        "goldJCeV.dat",
        "drude, drude-lorentz, spline"
    };

    static constexpr metal_spec silver{
        "silver",
        9.6,
        0.0228,
        5.3,
        "silverJCeV.dat",
        "drude, spline"
    };

    const metal_spec* spec = nullptr;

    if (mtl == "gold") {
        spec = &gold;
    } else if (mtl == "silver") {
        spec = &silver;
    } else {
        fail(std::string{mtl} + " not found, options are: gold, silver");
    }

    // If the class already has a cleanup routine for previous spline/GSL state,
    // call it here before reinitializing. For now we preserve the existing ecology.
    rows = 0;
    spln = 0;

    Ome_p   = spec->Ome_p_ev;
    Gam_d   = spec->Gam_d_ev;
    eps_inf = spec->eps_inf_val;

    if (mdl == "drude") {
        if (sel != 0) {
            warn(
                "'drude' model does not support JC uncertainty "
                "(sel=" + std::to_string(sel) + "). Ignoring uncertainty selection."
            );
        }
    }
    else if (mdl == "spline") {
        spln = 1;

        const fs::path jc_rel =
            fs::path{"data"} / "materials" / "metals" / spec->jc_filename;

        std::ifstream inp;
        fs::path jcfile;
        for (const fs::path& prefix : {fs::path{"."}, fs::path{".."}, fs::path{"../.."}}) {
            jcfile = prefix / jc_rel;
            inp.open(jcfile);
            if (inp) {
                break;
            }
            inp.clear();
        }
        if (!inp) {
            fail("could not open " + jc_rel.string() + " (tried ./, ../, ../../)");
        }

        std::vector<double> omem_vec;
        std::vector<double> reps_vec;
        std::vector<double> ieps_vec;

        double lam = 0.0;
        double om  = 0.0;
        double re  = 0.0;
        double im  = 0.0;
        double Dr  = 0.0;
        double Di  = 0.0;

        std::string line;
        while (std::getline(inp, line)) {
            std::istringstream iss(line);
            if (iss >> lam >> om >> re >> im >> Dr >> Di) {
                omem_vec.push_back(om);
                reps_vec.push_back(re + s * Dr);
                ieps_vec.push_back(im + s * Di);
            }
        }

        rows = static_cast<int>(omem_vec.size());
        if (rows == 0) {
            fail("empty JC file: " + jcfile.string());
        }

        omem = new double[rows];
        double* reps = new double[rows];
        double* ieps = new double[rows];

        for (std::size_t i = 0; i < rows; ++i) {
            omem[i] = omem_vec[i];
            reps[i] = reps_vec[i];
            ieps[i] = ieps_vec[i];
        }

        acc   = gsl_interp_accel_alloc();
        reeps = gsl_spline_alloc(gsl_interp_cspline, rows);
        imeps = gsl_spline_alloc(gsl_interp_cspline, rows);

        gsl_spline_init(reeps, omem, reps, rows);
        gsl_spline_init(imeps, omem, ieps, rows);
    }
    else {
        fail(
            "you chose '" + std::string{mdl} +
            "' while options for " + std::string{spec->name} +
            " are: " + spec->allowed_models
        );
    }

    ome_p  = Ome_p;
    ome_p2 = ome_p * ome_p;
}


inline double nanosphere::set_host(const char* hst){
  double eps_0;
  if(strcmp(hst,"silica")==0) eps_0 = 2.1316; // Si02
  else if(strcmp(hst,"solgel")==0) eps_0 = 2.1331;//3.8
  else if(strcmp(hst,"vacuum")==0) eps_0 = 1.0;
  else if(strcmp(hst,"water")==0) eps_0 = 1.7689;
  else if(strcmp(hst,"glass")==0) eps_0 = 1.2247;
  else if(strcmp(hst,"PMMA")==0) eps_0 = 2.2201;
  else if(strcmp(hst,"ethanol")==0) eps_0 = 1.8496;//1.1662;
  else{
std::cout<<"option are: silica, solgel, vacuum, water, glass, PMMA, ethanol"<<std::endl;
exit(1);
}
  return eps_0;
  }


inline void nanosphere::set_active(const char* mod){
  if(strcmp(mod,"lorentz")==0) act=1.;
    else if(strcmp(mod,"flat")==0) act=0.;
    else{
      std::cout<<"option are: flat, lorentz"<<std::endl;
      exit(1);
      }
  }


inline std::complex<double> nanosphere::metal(double ome) {
    std::complex<double> eps;

    if (spln == 1) {
        double x_min = omem[0];
        double x_max = omem[rows - 1];

        if (std::real(ome) < x_min || std::real(ome) > x_max) {
            std::cerr << "⚠️ GSL interpolation error: omega = " << std::real(ome)
                      << " is outside the interpolation range [" << x_min << ", " << x_max << "]\n";
            std::cerr << "Aborting safely to avoid crash.\n";
            std::exit(EXIT_FAILURE);
        }

        eps.real(gsl_spline_eval(reeps, ome, acc));
        eps.imag(gsl_spline_eval(imeps, ome, acc));
    } else {
        eps = eps_inf - ome_p2 / (ome * (ome + img * Gam_d));
    } 
    ceps_inf = eps + ome_p2 / (ome * (ome + img * Gam_d));
    return eps;
}


inline std::complex<double> nanosphere::confinement(double ome){
  std::complex<double> dlt;
  dlt=ome_p2/(ome*(ome+img*Gam_d))-ome_p2/(ome*(ome+img*(Gam_d+Gam)));
  return dlt;
  }


inline std::complex<double> nanosphere::active(double ome, double epsh){
  std::complex<double> eps;
      double nG = fabs(G);
      eps = std::complex<double>(epsh, 0.) + act*nG*Dome/(2.*(ome - ome_g) + img*Dome) - (1 - act)*nG*img;
      return eps;
    }


inline std::complex<double> nanosphere::set_GamG(double G, double tau2){
std::complex<double> GamG;
GamG   = -img*fabs(G)/tau2;
return GamG;
}


inline double * nanosphere::normalized_variables (){
double *nv, tau1, tau2, gamd, ome_g_norm;
nv = new double[4];
		
tau2 = 2.*Ome_p/Dome;
tau1 = 5.*tau2;

gamd = .5*Gam_d/Ome_p;
ome_g_norm = ome_g/Ome_p;

nv[0] = tau1;
nv[1] = tau2;
nv[2] = gamd;
nv[3] = ome_g_norm;

return nv;
}


inline std::complex<double> * nanosphere::set_ome_dep_vrbls(double ome, double ome_g, double tau2, double gamd){
std::complex<double> OmeG, GamM, OmeM;
std::complex<double> *odv;
odv = new std::complex<double> [3];

OmeG = img*(ome-ome_g)-1./tau2;
GamM = 1./(2.*(gamd-img*ome));
OmeM = (ome*(ome+2.*img*gamd))*GamM; 

odv[0] = OmeG;
odv[1] = OmeM;
odv[2] = GamM;

return odv;
}
