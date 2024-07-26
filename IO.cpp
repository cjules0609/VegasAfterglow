#include "IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

#include "macros.h"
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

    if (!file_r || !file_theta || !file_phi) {
        std::cerr << "Error opening files " << filename << "_*.txt" << std::endl;
        return;
    }

    file_r.precision(16);
    file_theta.precision(16);
    file_phi.precision(16);

    for (size_t i = 0; i < coord.r.size(); ++i) {
        file_r << coord.r[i] / con::cm << " ";
    }

    for (size_t i = 0; i < coord.theta.size(); ++i) {
        file_theta << coord.theta[i] << " ";
    }

    for (size_t i = 0; i < coord.phi.size(); ++i) {
        file_phi << coord.phi[i] << " ";
    }

    file_r << std::endl;
    file_theta << std::endl;
    file_phi << std::endl;
}

void write2file(Shock const& shock, std::string const& filename) {
    std::ofstream file(filename + "_Gamma.txt");
    std::ofstream file_B(filename + "_B.txt");
    std::ofstream file_D_com(filename + "_D_com.txt");
    std::ofstream file_t_com(filename + "_t_com.txt");
    std::ofstream file_np(filename + "_n_p.txt");
    std::ofstream file_e_th(filename + "_e_th.txt");

    if (!file || !file_B || !file_D_com || !file_t_com || !file_np) {
        std::cerr << "Error opening files " << filename << "_*.txt" << std::endl;
        return;
    }

    file.precision(16);
    file_B.precision(16);
    file_D_com.precision(16);
    file_t_com.precision(16);
    file_np.precision(16);
    file_e_th.precision(16);

    for (size_t j = 0; j < shock.Gamma.size(); ++j) {
        for (size_t k = 0; k < shock.Gamma[j].size(); ++k) {
            file << shock.Gamma[j][k] << " ";
            file_B << shock.B[j][k] << " ";  // TODO:add B unit
            file_D_com << shock.width_eff[j][k] / con::cm << " ";
            file_t_com << shock.t_com[j][k] / con::sec << " ";
            file_np << shock.n_p[j][k] * (con::cm * con::cm * con::cm) << " ";
            file_e_th << shock.e_th[j][k] / con::erg << " ";
        }
        file << '\n';
        file_B << '\n';
        file_D_com << '\n';
        file_t_com << '\n';
        file_np << '\n';
        file_e_th << '\n';
    }
}

void write2file(SynPhotonsMesh const& syn_rad, std::string const& filename) {
    std::ofstream file_I_peak(filename + "_I_nu_peak.txt");
    std::ofstream file_nu_peak(filename + "_nu_E_peak.txt");
    std::ofstream file_nu_m(filename + "_nu_m.txt");
    std::ofstream file_nu_c(filename + "_nu_c.txt");
    std::ofstream file_nu_a(filename + "_nu_a.txt");
    std::ofstream file_nu_M(filename + "_nu_Max.txt");

    if (!file_I_peak || !file_nu_peak || !file_nu_m || !file_nu_c || !file_nu_a || !file_nu_M) {
        std::cerr << "Error opening files " << filename << "_*.txt" << std::endl;
        return;
    }

    file_I_peak.precision(16);
    file_nu_peak.precision(16);
    file_nu_m.precision(16);
    file_nu_c.precision(16);
    file_nu_a.precision(16);
    file_nu_M.precision(16);

    for (size_t i = 0; i < syn_rad.size(); ++i) {
        for (size_t j = 0; j < syn_rad[i].size(); ++j) {
            file_I_peak << syn_rad[i][j].L_nu_peak / (con::erg / con::sec) << " ";
            file_nu_peak << syn_rad[i][j].nu_E_peak / con::Hz << " ";
            file_nu_m << syn_rad[i][j].nu_m / con::Hz << " ";
            file_nu_c << syn_rad[i][j].nu_c / con::Hz << " ";
            file_nu_a << syn_rad[i][j].nu_a / con::Hz << " ";
            file_nu_M << syn_rad[i][j].nu_M / con::Hz << " ";
        }
        file_I_peak << '\n';
        file_nu_peak << '\n';
        file_nu_m << '\n';
        file_nu_c << '\n';
        file_nu_a << '\n';
        file_nu_M << '\n';
    }
}

void write2file(SynElectronsMesh const& syn_rad, std::string const& filename) {
    std::ofstream file_n_tot(filename + "_n_tot.txt");
    std::ofstream file_gamma_peak(filename + "_gamma_N_peak.txt");
    std::ofstream file_gamma_m(filename + "_gamma_m.txt");
    std::ofstream file_gamma_c(filename + "_gamma_c.txt");
    std::ofstream file_gamma_a(filename + "_gamma_a.txt");
    std::ofstream file_gamma_M(filename + "_gamma_Max.txt");
    std::ofstream file_N_tot(filename + "_N.txt");

    if (!file_n_tot || !file_gamma_peak || !file_gamma_m || !file_gamma_c || !file_gamma_a || !file_gamma_M ||
        !file_N_tot) {
        std::cerr << "Error opening files " << std::endl;
        return;
    }

    file_n_tot.precision(16);
    file_gamma_peak.precision(16);
    file_gamma_m.precision(16);
    file_gamma_c.precision(16);
    file_gamma_a.precision(16);
    file_gamma_M.precision(16);
    file_N_tot.precision(16);

    for (size_t i = 0; i < syn_rad.size(); ++i) {
        for (size_t j = 0; j < syn_rad[i].size(); ++j) {
            file_n_tot << syn_rad[i][j].n_tot * (con::cm * con::cm * con::cm) << " ";
            file_gamma_peak << syn_rad[i][j].gamma_N_peak << " ";
            file_gamma_m << syn_rad[i][j].gamma_m << " ";
            file_gamma_c << syn_rad[i][j].gamma_c << " ";
            file_gamma_a << syn_rad[i][j].gamma_a << " ";
            file_gamma_M << syn_rad[i][j].gamma_M << " ";
            file_N_tot << syn_rad[i][j].N_tot << " ";
        }
        file_n_tot << '\n';
        file_gamma_peak << '\n';
        file_gamma_m << '\n';
        file_gamma_c << '\n';
        file_gamma_a << '\n';
        file_gamma_M << '\n';
        file_N_tot << '\n';
    }

    for (size_t i = 0; i < syn_rad[0][0].Ys.size(); ++i) {
        std::ofstream file_Y(filename + "_Yt-" + std::to_string(i) + ".txt");
        std::ofstream file_gamma_hat_m(filename + "_gamma_hat_m-" + std::to_string(i) + ".txt");
        std::ofstream file_gamma_hat_c(filename + "_gamma_hat_c-" + std::to_string(i) + ".txt");

        if (!file_Y || !file_gamma_hat_m || !file_gamma_hat_c) {
            std::cerr << "Error opening files " << std::endl;
            return;
        }

        file_Y.precision(16);
        file_gamma_hat_m.precision(16);
        file_gamma_hat_c.precision(16);

        for (size_t i = 0; i < syn_rad.size(); ++i) {
            for (size_t j = 0; j < syn_rad[i].size(); ++j) {
                file_N_tot << syn_rad[i][j].N_tot << " ";
                file_Y << syn_rad[i][j].Ys[i].Y_T << " ";
                file_gamma_hat_m << syn_rad[i][j].Ys[i].gamma_hat_m << " ";
                file_gamma_hat_c << syn_rad[i][j].Ys[i].gamma_hat_c << " ";
            }
            file_Y << '\n';
            file_gamma_hat_m << '\n';
            file_gamma_hat_c << '\n';
        }
    }
}

void write2file(MeshGrid3d const& array, std::string const& filename) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }

    file.precision(16);
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] << " ";
            }
            file << '\n';
        }
    }
}
void write2file(MeshGrid const& grid, std::string const& filename) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }

    file.precision(16);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] << " ";
        }
        file << '\n';
    }
}

void write2file(Array const& array, std::string const& filename) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }

    file.precision(16);
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] << " ";
    }
}

void write2file(MeshGrid3d const& array, std::string const& filename, double unit) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }
    file.precision(16);
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] / unit << " ";
            }
            file << '\n';
        }
    }
}
void write2file(MeshGrid const& grid, std::string const& filename, double unit) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }
    file.precision(16);
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] / unit << " ";
        }
        file << '\n';
    }
}

void write2file(Array const& array, std::string const& filename, double unit) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt" << std::endl;
        return;
    }
    file.precision(16);
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] / unit << " ";
    }
}