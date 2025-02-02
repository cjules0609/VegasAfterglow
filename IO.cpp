#include "IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

#include "macros.h"
void printArray(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void output(Coord const& coord, std::string const& filename) {
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

void writeGrid(std::string filename, auto const& names, auto const& data, auto const& units) {
    for (size_t l = 0; l < names.size(); ++l) {
        std::ofstream file(filename + "_" + names[l] + ".txt");
        if (!file) {
            std::cerr << "Error opening files " << filename + "_" + names[l] + ".txt" << std::endl;
            return;
        }
        file.precision(16);

        for (size_t i = 0; i < (*data[l]).size(); ++i) {
            for (size_t j = 0; j < (*data[l])[i].size(); ++j) {
                for (size_t k = 0; k < (*data[l])[i][j].size(); ++k) {
                    file << (*data[l])[i][j][k] / units[l] << " ";
                }
                file << '\n';
            }
            file << '\n';
        }
    }
}
void output(Shock const& shock, std::string const& filename) {
    std::array<std::string, 5> strs = {"Gamma", "B", "t_com", "t_eng", "Sigma"};
    std::array<MeshGrid3d const*, 5> data = {&(shock.Gamma_rel), &(shock.B), &(shock.t_com), &(shock.t_eng),
                                             &(shock.column_num_den)};
    std::array<double, 5> units = {1, 1, con::sec, con::sec, 1 / con::cm2};

    writeGrid(filename, strs, data, units);
}

void output(PromptPhotonsGrid const& prompt_pj, std::string const& filename) {}

void output(SynPhotonGrid const& ph, std::string const& filename) {}

void output(SynElectronGrid const& syn_rad, std::string const& filename) {}

void output(MeshGrid3d const& array, std::string const& filename) {
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

void output(MeshGrid const& grid, std::string const& filename) {
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

void output(Array const& array, std::string const& filename) {
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

void output(MeshGrid3d const& array, std::string const& filename, double unit) {
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
void output(MeshGrid const& grid, std::string const& filename, double unit) {
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

void output(Array const& array, std::string const& filename, double unit) {
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