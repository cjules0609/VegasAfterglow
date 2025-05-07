//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mcmc.h"
#include "emission.h"

namespace py = pybind11;
PYBIND11_MODULE(VegasAfterglowC, m) {
    // Parameters for MCMC modeling
    py::class_<Params>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("E_iso", &Params::E_iso)
        .def_readwrite("Gamma0", &Params::Gamma0)
        .def_readwrite("theta_c", &Params::theta_c)
        .def_readwrite("theta_v", &Params::theta_v)
        .def_readwrite("theta_w", &Params::theta_w)
        .def_readwrite("p", &Params::p)
        .def_readwrite("eps_e", &Params::eps_e)
        .def_readwrite("eps_B", &Params::eps_B)
        .def_readwrite("n_ism", &Params::n_ism)
        .def_readwrite("A_star", &Params::A_star)
        .def_readwrite("xi", &Params::xi)
        .def_readwrite("k_jet", &Params::k_jet)
        .def("__repr__", [](const Params &p) {
            return "<Params E_iso=" + std::to_string(p.E_iso) + ", Gamma0=" + std::to_string(p.Gamma0) +
                   ", theta_c=" + std::to_string(p.theta_c) + ", theta_v=" + std::to_string(p.theta_v) +
                   ", theta_w=" + std::to_string(p.theta_w) + ", p=" + std::to_string(p.p) +
                   ", eps_e=" + std::to_string(p.eps_e) + ", eps_B=" + std::to_string(p.eps_B) +
                   ", n_ism=" + std::to_string(p.n_ism) + ", A_star=" + std::to_string(p.A_star) +
                   ", xi=" + std::to_string(p.xi) + ", k_jet=" + std::to_string(p.k_jet) + ">";
        });
    // Parameters for modeling that are not used in the MCMC
    py::class_<ConfigParams>(m, "Setups")
        .def(py::init<>())
        .def_readwrite("lumi_dist", &ConfigParams::lumi_dist)
        .def_readwrite("z", &ConfigParams::z)
        .def_readwrite("medium", &ConfigParams::medium)
        .def_readwrite("jet", &ConfigParams::jet)
        .def_readwrite("t_grid", &ConfigParams::t_grid)
        .def_readwrite("phi_grid", &ConfigParams::phi_grid)
        .def_readwrite("theta_grid", &ConfigParams::theta_grid)
        .def_readwrite("rtol", &ConfigParams::rtol)
        .def("__repr__", [](const ConfigParams &c) {
            return "<ConfigParams lumi_dist=" + std::to_string(c.lumi_dist) + ", z=" + std::to_string(c.z) +
                   ", medium='" + c.medium + "', jet='" + c.jet + "', t_grid_num=" + std::to_string(c.t_grid) +
                   ", phi_grid_num=" + std::to_string(c.phi_grid) + ", theta_grid_num=" + std::to_string(c.theta_grid) +
                   ", rtol=" + std::to_string(c.rtol) + ">";
        });

    // MultiBandData bindings
    py::class_<MultiBandData>(m, "ObsData")
        .def(py::init<>())
        .def("add_light_curve", &MultiBandData::add_light_curve, py::arg("nu_cgs"), py::arg("t_cgs"),
             py::arg("Fnu_cgs"), py::arg("Fnu_err"))
        .def("add_spectrum", &MultiBandData::add_spectrum, py::arg("t_cgs"), py::arg("nu_cgs"), py::arg("Fnu_cgs"),
             py::arg("Fnu_err"));

    // MultiBandModel bindings
    py::class_<MultiBandModel>(m, "VegasMC")
        .def(py::init<MultiBandData const &>(), py::arg("obs_data"))
        .def("set", &MultiBandModel::configure, py::arg("param"))
        .def("estimate_chi2", &MultiBandModel::estimate_chi2, py::arg("param"),
             py::call_guard<py::gil_scoped_release>())
        .def("light_curves", &MultiBandModel::light_curves, py::arg("param"), py::arg("t_cgs"), py::arg("nu_cgs"),
             py::call_guard<py::gil_scoped_release>())
        .def("spectra", &MultiBandModel::spectra, py::arg("param"), py::arg("nu_cgs"), py::arg("t_cgs"),
             py::call_guard<py::gil_scoped_release>());

    py::class_<Emission>(m, "Emission")
        .def(py::init<const Params&, const ConfigParams&>())
        .def("generate_lc", &Emission::generate_lc, py::arg("nu"), py::arg("t"))
        .def("generate_spec", &Emission::generate_spec, py::arg("nu"), py::arg("t"));

}