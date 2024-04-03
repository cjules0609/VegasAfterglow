#include "IO.h"

#include <fstream>
#include <iostream>
void print_array(Array const& arr) {
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    std::cout << std::endl;
}

void write2file(SynElectronMesh const& syn_rad, std::string const& filename) {
    std::ofstream file_I_peak(filename + "_I_nu_peak.txt");
    std::ofstream file_nu_peak(filename + "_nu_E_peak.txt");
    std::ofstream file_nu_m(filename + "_nu_m.txt");
    std::ofstream file_nu_c(filename + "_nu_c.txt");
    std::ofstream file_nu_a(filename + "_nu_a.txt");
    std::ofstream file_nu_M(filename + "_nu_Max.txt");
    for (size_t i = 0; i < syn_rad.size(); ++i) {
        for (size_t j = 0; j < syn_rad[i].size(); ++j) {
            file_I_peak << syn_rad[i][j].I_nu_peak << " ";
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