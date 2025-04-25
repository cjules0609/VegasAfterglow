//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "utilities.h"
struct LightCurveData {
    double nu{0};
    Array t;
    Array Fv_obs;
    Array Fv_err;
    Array Fv_model;

    double calcChiSquare() const;
};

struct SpectrumData {
    double t{0};
    Array nu;
    Array Fv_obs;
    Array Fv_err;
    Array Fv_model;

    double calcChiSquare() const;
};

struct MultiBandData {
    using List = std::vector<double>;

    std::vector<LightCurveData> light_curve;
    std::vector<SpectrumData> spectrum;

    double calcChiSquare() const;
    void addObsLightCurve(double nu, List const& t, List const& Fv_obs, List const& Fv_err);
    void addObsSpectrum(double t, List const& nu, List const& Fv_obs, List const& Fv_err);
};

struct Params {
    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double theta_v{0};
    double theta_w{con::pi / 2};
    double p{2.3};
    double eps_e{0.1};
    double eps_B{0.01};
    double n_ism{1};
    double A_star{0.01};
    double xi{1};
    double k_jet{2};
};

struct ConfigParams {
    double lumi_dist{1e26};
    double z{0};
    std::string medium{"ism"};
    std::string jet{"tophat"};
    size_t t_grid{24};
    size_t phi_grid{24};
    size_t theta_grid{24};
    double rtol{1e-4};
};

struct MultiBandModel {
    using List = std::vector<double>;
    using Grid = std::vector<std::vector<double>>;
    MultiBandModel() = delete;
    MultiBandModel(MultiBandData const& data);

    void configure(ConfigParams const& param);
    double chiSquare(Params const& param);
    Grid lightCurve(Params const& param, List const& t, List const& nu);
    Grid spectrum(Params const& param, List const& nu, List const& t);

   private:
    void buildSystem(Params const& param, Array const& t_eval, Observer& obs, SynElectronGrid& electrons,
                     SynPhotonGrid& photons);
    MultiBandData obs_data;
    ConfigParams config;
    // SynElectronGrid electrons;
    // SynPhotonGrid photons;
    // Observer obs;
    Array t_eval;
};