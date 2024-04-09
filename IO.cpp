#include "IO.h"

#include <fstream>
#include <iostream>
void print_array(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void write2file(Coord const& coord, std::string const& filename) {
    std::ofstream file_r(filename + "_r.txt");
    std::ofstream file_theta(filename + "_theta.txt");
    std::ofstream file_phi(filename + "_phi.txt");

    for (size_t i = 0; i < coord.r.size(); ++i) {
        file_r << coord.r[i] << " ";
    }

    file_r << std::endl;

    for (size_t i = 0; i < coord.theta.size(); ++i) {
        file_theta << coord.theta[i] << " ";
    }

    file_theta << std::endl;

    for (size_t i = 0; i < coord.phi.size(); ++i) {
        file_phi << coord.phi[i] << " ";
    }

    file_phi << std::endl;
}

void write2file(Shock const& shock, std::string const& filename) {
    std::ofstream file(filename + "_Gamma.txt");
    std::ofstream file_B(filename + "_B.txt");
    std::ofstream file_D_com(filename + "_D_com.txt");
    std::ofstream file_t_com(filename + "_t_com.txt");

    for (size_t j = 0; j < shock.Gamma.size(); ++j) {
        for (size_t k = 0; k < shock.Gamma[j].size(); ++k) {
            file << shock.Gamma[j][k] << " ";
            file_B << shock.B[j][k] << " ";
            file_D_com << shock.D_com[j][k] << " ";
            file_t_com << shock.t_com[j][k] << " ";
        }
        file << std::endl;
        file_B << std::endl;
        file_D_com << std::endl;
        file_t_com << std::endl;
    }
}

void write2file(SynPhotonsMesh const& syn_rad, std::string const& filename) {
    std::ofstream file_I_peak(filename + "_I_nu_peak.txt");
    std::ofstream file_nu_peak(filename + "_nu_E_peak.txt");
    std::ofstream file_nu_m(filename + "_nu_m.txt");
    std::ofstream file_nu_c(filename + "_nu_c.txt");
    std::ofstream file_nu_a(filename + "_nu_a.txt");
    std::ofstream file_nu_M(filename + "_nu_Max.txt");
    for (size_t i = 0; i < syn_rad.size(); ++i) {
        for (size_t j = 0; j < syn_rad[i].size(); ++j) {
            file_I_peak << syn_rad[i][j].j_nu_peak << " ";
            file_nu_peak << syn_rad[i][j].nu_E_peak << " ";
            file_nu_m << syn_rad[i][j].nu_m << " ";
            file_nu_c << syn_rad[i][j].nu_c << " ";
            file_nu_a << syn_rad[i][j].nu_a << " ";
            file_nu_M << syn_rad[i][j].nu_M << " ";
        }
        file_I_peak << std::endl;
        file_nu_peak << std::endl;
        file_nu_m << std::endl;
        file_nu_c << std::endl;
        file_nu_a << std::endl;
        file_nu_M << std::endl;
    }
}

void write2file(SynElectronsMesh const& syn_rad, std::string const& filename) {
    std::ofstream file_n_tot(filename + "_n_tot.txt");
    std::ofstream file_gamma_peak(filename + "_gamma_N_peak.txt");
    std::ofstream file_gamma_m(filename + "_gamma_m.txt");
    std::ofstream file_gamma_c(filename + "_gamma_c.txt");
    std::ofstream file_gamma_a(filename + "_gamma_a.txt");
    std::ofstream file_gamma_M(filename + "_gamma_Max.txt");
    for (size_t i = 0; i < syn_rad.size(); ++i) {
        for (size_t j = 0; j < syn_rad[i].size(); ++j) {
            file_n_tot << syn_rad[i][j].n_tot << " ";
            file_gamma_peak << syn_rad[i][j].gamma_N_peak << " ";
            file_gamma_m << syn_rad[i][j].gamma_m << " ";
            file_gamma_c << syn_rad[i][j].gamma_c << " ";
            file_gamma_a << syn_rad[i][j].gamma_a << " ";
            file_gamma_M << syn_rad[i][j].gamma_M << " ";
        }
        file_n_tot << std::endl;
        file_gamma_peak << std::endl;
        file_gamma_m << std::endl;
        file_gamma_c << std::endl;
        file_gamma_a << std::endl;
        file_gamma_M << std::endl;
    }
}

void write2file(MeshGrid3d const& array, std::string const& filename) {
    std::ofstream file(filename + ".txt");
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] << " ";
            }
            file << std::endl;
        }
    }
}
void write2file(MeshGrid const& grid, std::string const& filename) {
    std::ofstream file(filename + ".txt");
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] << " ";
        }
        file << std::endl;
    }
}

void write2file(Array const& array, std::string const& filename) {
    std::ofstream file(filename + ".txt");
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] << " ";
    }
}