//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <vector>

#include "macros.h"
#include "mesh.h"
#include "utilities.h"

struct LightCurveData {
    using List = std::vector<double>;
    double nu{0};
    List t;
    List Fv_obs;
    List Fv_err;
    List Fv_model;

    double calcChiSquare() const;
};

struct SpectrumData {
    using List = std::vector<double>;
    double t{0};
    List nu;
    List Fv_obs;
    List Fv_err;
    List Fv_model;

    double calcChiSquare() const;
};

struct MultiBandData {
    using List = std::vector<double>;
    std::vector<LightCurveData> light_curve;
    List t_grid;
    List lc_band;

    std::vector<SpectrumData> spectrum;
    List nu_grid;
    List spectrum_t;

    double calcChiSquare() const;
    void genEstimatePoints();
    void addObsLightCurve(double nu, List const& t, List const& Fv_obs, List const& Fv_err);
    void addObsSpectrum(double t, List const& nu, List const& Fv_obs, List const& Fv_err);
    void addModelLightCurve(MeshGrid const& light_curves);
    void addModelSpectrum(MeshGrid const& spectra);
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
    size_t t_grid{64};
    size_t phi_grid{64};
    size_t theta_grid{64};
    double rtol{1e-6};
};

struct MultiBandModel {
    MultiBandModel() = delete;
    MultiBandModel(MultiBandData const& data) : obs_data(data) { obs_data.genEstimatePoints(); }

    void configure(ConfigParams const& param);
    double chiSquare(Params const& param);
    std::vector<std::vector<double>> lightCurve(Params const& param, std::vector<double> const& t,
                                                std::vector<double> const& nu);
    std::vector<std::vector<double>> spectrum(Params const& param, std::vector<double> const& nu,
                                              std::vector<double> const& t);
    std::vector<double> chiSquareBatch(std::vector<Params> const& param_batch);

   private:
    MultiBandData obs_data;
    ConfigParams config;
};
