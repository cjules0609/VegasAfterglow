//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "inverse-compton.h"

#include <cmath>
#include <iostream>
#include <thread>

#include "macros.h"
#include "utilities.h"

/********************************************************************************************************************
 * FUNCTION: createICPhotonGrid
 * DESCRIPTION: Creates and returns an ICPhotonGrid with dimensions (phi_size x theta_size x r_size)
 *              using boost::multi_array with the specified extents.
 ********************************************************************************************************************/
ICPhotonGrid createICPhotonGrid(size_t phi_size, size_t theta_size, size_t r_size) {
    return ICPhotonGrid(boost::extents[phi_size][theta_size][r_size]);
}

/********************************************************************************************************************
 * INLINE FUNCTION: order
 * DESCRIPTION: Returns true if the three arguments are in strictly increasing order.
 ********************************************************************************************************************/
inline bool order(double a, double b, double c) { return a < b && b < c; }

/********************************************************************************************************************
 * METHOD: ICPhoton::I_nu
 * DESCRIPTION: Computes the photon intensity at frequency nu using logarithmic-logarithmic interpolation
 *              (specifically, loglogInterpEqSpaced) on the IC photon spectrum data.
 ********************************************************************************************************************/
double ICPhoton::I_nu(double nu) const { return loglogInterpEqSpaced(nu, this->nu_IC_, this->j_nu_, true, true); }

/********************************************************************************************************************
 * INLINE FUNCTION: eta_rad
 * DESCRIPTION: Computes the radiative efficiency parameter (ηₑ) given minimum electron Lorentz factors.
 *              If gamma_c is less than gamma_m, it returns 1; otherwise, it returns (gamma_c/gamma_m)^(2-p).
 ********************************************************************************************************************/
inline double eta_rad(double gamma_m, double gamma_c, double p) {
    return gamma_c < gamma_m ? 1 : std::pow(gamma_c / gamma_m, (2 - p));
}

/********************************************************************************************************************
 * FUNCTION: effectiveYThomson
 * DESCRIPTION: Computes the effective Compton Y parameter in the Thomson regime.
 *              It iteratively solves for Y until convergence using the relation:
 *                  Y0 = (sqrt(1+4b) - 1)/2,
 *              where b = (ηₑ * eps_e / eps_B). The electron cooling parameters are updated during each iteration.
 ********************************************************************************************************************/
double effectiveYThomson(double B, double t_com, double eps_e, double eps_B, SynElectrons const& e) {
    double eta_e = eta_rad(e.gamma_m, e.gamma_c, e.p);
    double b = eta_e * eps_e / eps_B;
    double Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    double Y1 = 2 * Y0;
    for (; std::fabs((Y1 - Y0) / Y0) > 1e-5;) {
        Y1 = Y0;
        double gamma_c = syn_gamma_c(t_com, B, e.Ys, e.p);
        eta_e = eta_rad(e.gamma_m, gamma_c, e.p);
        b = eta_e * eps_e / eps_B;
        Y0 = (std::sqrt(1 + 4 * b) - 1) / 2;
    }
    return Y0;
}

/********************************************************************************************************************
 * FUNCTION: genICPhotons
 * DESCRIPTION: Generates a grid of ICPhoton objects (ICPhotonGrid) from the provided SynElectronGrid and
 *              SynPhotonGrid. For each grid cell, it calls the gen() method of the ICPhoton to compute the IC
 *              photon spectrum based on the local electron and synchrotron photon properties.
 ********************************************************************************************************************/
ICPhotonGrid genICPhotons(SynElectronGrid const& e, SynPhotonGrid const& ph) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];
    ICPhotonGrid IC_ph = createICPhotonGrid(phi_size, theta_size, r_size);

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                // Generate the IC photon spectrum for each grid cell.
                IC_ph[i][j][k].gen(e[i][j][k], ph[i][j][k]);
            }
        }
    }
    return IC_ph;
}

/********************************************************************************************************************
 * FUNCTION: eCoolingThomson
 * DESCRIPTION: Applies electron cooling in the Thomson regime.
 *              For each cell in the SynElectronGrid, it computes the effective Y parameter using effectiveYThomson,
 *              clears the current inverse Compton Y parameters (Ys), and stores the computed Y_T.
 *              Finally, it updates the electrons based on the new Y parameter.
 ********************************************************************************************************************/
void eCoolingThomson(SynElectronGrid& e, SynPhotonGrid const& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];

    for (size_t i = 0; i < phi_size; i++) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                double Y_T =
                    effectiveYThomson(shock.B[i][j][k], shock.t_com[i][j][k], shock.eps_e, shock.eps_B, e[i][j][k]);

                // Clear previous inverse Compton Y parameters and store the new value.
                e[i][j][k].Ys.clear();
                e[i][j][k].Ys.emplace_back(Y_T);
            }
        }
    }
    updateElectrons4Y(e, shock);
}

/********************************************************************************************************************
 * FUNCTION: eCoolingKleinNishina
 * DESCRIPTION: Applies electron cooling in the Klein-Nishina regime.
 *              Similar to eCoolingThomson, but for each cell, it creates an InverseComptonY object with additional
 *              parameters from the synchrotron photon grid.
 ********************************************************************************************************************/
void eCoolingKleinNishina(SynElectronGrid& e, SynPhotonGrid const& ph, Shock const& shock) {
    size_t phi_size = e.shape()[0];
    size_t theta_size = e.shape()[1];
    size_t r_size = e.shape()[2];
    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < r_size; ++k) {
                double Y_T =
                    effectiveYThomson(shock.B[i][j][k], shock.t_com[i][j][k], shock.eps_e, shock.eps_B, e[i][j][k]);
                // Clear existing Ys and emplace a new InverseComptonY with additional synchrotron frequency parameters.
                e[i][j][k].Ys.clear();
                e[i][j][k].Ys.emplace_back(ph[i][j][k].nu_m, ph[i][j][k].nu_c, shock.B[i][j][k], Y_T);
            }
        }
    }
    updateElectrons4Y(e, shock);
}