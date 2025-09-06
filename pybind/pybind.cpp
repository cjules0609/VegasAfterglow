//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#define FORCE_IMPORT_ARRAY  // numpy C api loading, must before any xtensor-python headers
#include "pybind.h"

#include "mcmc.h"
#include "pymodel.h"

PYBIND11_MODULE(VegasAfterglowC, m) {
    xt::import_numpy();
    // Jet bindings
    py::object zero2d_fn = py::cpp_function(func::zero_2d);
    py::object zero3d_fn = py::cpp_function(func::zero_3d);

    //========================================================================================================
    //                                 Model bindings
    //========================================================================================================
    py::class_<PyMagnetar>(m, "Magnetar").def(py::init<Real, Real, Real>(), py::arg("L0"), py::arg("t0"), py::arg("q"));

    m.def("TophatJet", &PyTophatJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("GaussianJet", &PyGaussianJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("PowerLawJet", &PyPowerLawJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"), py::arg("k_e"),
          py::arg("k_g"), py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("TwoComponentJet", &PyTwoComponentJet, py::arg("theta_c"), py::arg("E_iso_c"), py::arg("Gamma0_c"),
          py::arg("theta_w"), py::arg("E_iso_w"), py::arg("Gamma0_w"), py::arg("spreading") = false,
          py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("StepPowerLawJet", &PyStepPowerLawJet, py::arg("theta_c"), py::arg("E_iso_c"), py::arg("Gamma0_c"),
          py::arg("E_iso_w"), py::arg("Gamma0_w"), py::arg("k_e"), py::arg("k_g"), py::arg("spreading") = false,
          py::arg("duration") = 1, py::arg("magnetar") = py::none());

    py::class_<Ejecta>(m, "Ejecta")
        .def(py::init<BinaryFunc, BinaryFunc, BinaryFunc, TernaryFunc, TernaryFunc, bool, Real>(), py::arg("E_iso"),
             py::arg("Gamma0"), py::arg("sigma0") = zero2d_fn, py::arg("E_dot") = zero3d_fn,
             py::arg("M_dot") = zero3d_fn, py::arg("spreading") = false, py::arg("duration") = 1);

    // Medium bindings
    m.def("ISM", &PyISM, py::arg("n_ism"));

    m.def("Wind", &PyWind, py::arg("A_star"), py::arg("n_ism") = 0, py::arg("n_0") = con::inf);

    py::class_<Medium>(m, "Medium").def(py::init<TernaryFunc>(), py::arg("rho"));

    // Observer bindings
    py::class_<PyObserver>(m, "Observer")
        .def(py::init<Real, Real, Real, Real>(), py::arg("lumi_dist"), py::arg("z"), py::arg("theta_obs"),
             py::arg("phi_obs") = 0);

    // Radiation bindings
    py::class_<PyRadiation>(m, "Radiation")
        .def(py::init<Real, Real, Real, Real, bool, bool, bool>(), py::arg("eps_e"), py::arg("eps_B"), py::arg("p"),
             py::arg("xi_e") = 1, py::arg("IC_cooling") = false, py::arg("SSC") = false, py::arg("KN") = false);

    // Model bindings
    py::class_<PyModel>(m, "Model")
        .def(py::init<Ejecta, Medium, PyObserver, PyRadiation, std::optional<PyRadiation>, std::tuple<Real, Real, Real>,
                      Real, bool>(),
             py::arg("jet"), py::arg("medium"), py::arg("observer"), py::arg("forward_rad"),
             py::arg("reverse_rad") = py::none(), py::arg("resolutions") = std::make_tuple(0.3, 1, 10),
             py::arg("rtol") = 1e-5, py::arg("axisymmetric") = true)
        .def("specific_flux", &PyModel::specific_flux, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())
        .def("specific_flux_series", &PyModel::specific_flux_series, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())
        .def("specific_flux_series_with_expo", &PyModel::specific_flux_series_with_expo, py::arg("t"), py::arg("nu"),
             py::arg("expo_time"), py::arg("num_points") = 10, py::call_guard<py::gil_scoped_release>())
        .def("details", &PyModel::details, py::arg("t_min"), py::arg("t_max"),
             py::call_guard<py::gil_scoped_release>());

    //========================================================================================================
    //                                 MCMC bindings
    //========================================================================================================
    // Parameters for MCMC modeling
    py::class_<Params>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("theta_v", &Params::theta_v)

        .def_readwrite("n_ism", &Params::n_ism)
        .def_readwrite("n0", &Params::n0)
        .def_readwrite("A_star", &Params::A_star)

        .def_readwrite("E_iso", &Params::E_iso)
        .def_readwrite("Gamma0", &Params::Gamma0)
        .def_readwrite("theta_c", &Params::theta_c)
        .def_readwrite("k_e", &Params::k_e)
        .def_readwrite("k_g", &Params::k_g)
        .def_readwrite("tau", &Params::duration)

        .def_readwrite("E_iso_w", &Params::E_iso_w)
        .def_readwrite("Gamma0_w", &Params::Gamma0_w)
        .def_readwrite("theta_w", &Params::theta_w)

        .def_readwrite("L0", &Params::L0)
        .def_readwrite("t0", &Params::t0)
        .def_readwrite("q", &Params::q)

        .def_readwrite("p", &Params::p)
        .def_readwrite("eps_e", &Params::eps_e)
        .def_readwrite("eps_B", &Params::eps_B)
        .def_readwrite("xi_e", &Params::xi_e)

        .def_readwrite("p_r", &Params::p_r)
        .def_readwrite("eps_e_r", &Params::eps_e_r)
        .def_readwrite("eps_B_r", &Params::eps_B_r)
        .def_readwrite("xi_e_r", &Params::xi_e_r);

    // Parameters for modeling that are not used in the MCMC
    py::class_<ConfigParams>(m, "Setups")
        .def(py::init<>())
        .def_readwrite("lumi_dist", &ConfigParams::lumi_dist)
        .def_readwrite("z", &ConfigParams::z)
        .def_readwrite("medium", &ConfigParams::medium)
        .def_readwrite("jet", &ConfigParams::jet)
        .def_readwrite("t_resol", &ConfigParams::t_resol)
        .def_readwrite("phi_resol", &ConfigParams::phi_resol)
        .def_readwrite("theta_resol", &ConfigParams::theta_resol)
        .def_readwrite("rtol", &ConfigParams::rtol)
        .def_readwrite("rvs_shock", &ConfigParams::rvs_shock)
        .def_readwrite("fwd_SSC", &ConfigParams::fwd_SSC)
        .def_readwrite("rvs_SSC", &ConfigParams::rvs_SSC)
        .def_readwrite("IC_cooling", &ConfigParams::IC_cooling)
        .def_readwrite("KN", &ConfigParams::KN)
        .def_readwrite("magnetar", &ConfigParams::magnetar);

    // MultiBandData bindings
    py::class_<MultiBandData>(m, "ObsData")
        .def(py::init<>())
        .def("add_light_curve", &MultiBandData::add_light_curve, py::arg("nu_cgs"), py::arg("t_cgs"),
             py::arg("Fnu_cgs"), py::arg("Fnu_err"), py::arg("weights") = py::none())
        .def("add_spectrum", &MultiBandData::add_spectrum, py::arg("t_cgs"), py::arg("nu_cgs"), py::arg("Fnu_cgs"),
             py::arg("Fnu_err"), py::arg("weights") = py::none());

    // MultiBandModel bindings
    py::class_<MultiBandModel>(m, "VegasMC")
        .def(py::init<MultiBandData const &>(), py::arg("obs_data"))
        .def("set", &MultiBandModel::configure, py::arg("param"))
        .def("estimate_chi2", &MultiBandModel::estimate_chi2, py::arg("param"),
             py::call_guard<py::gil_scoped_release>())
        .def("specific_flux", &MultiBandModel::specific_flux, py::arg("param"), py::arg("t_cgs"), py::arg("nu_cgs"),
             py::call_guard<py::gil_scoped_release>());
}
