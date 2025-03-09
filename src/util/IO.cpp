//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "IO.h"

#include <fstream>
#include <iomanip>
#include <iostream>

#include "macros.h"
/********************************************************************************************************************
 * FUNCTION: printArray
 * DESCRIPTION: Prints all elements of a 1D Array to the standard output.
 ********************************************************************************************************************/
void printArray(Array const& arr) {
    // Iterate through each element of the array and output it.
    for (auto const& a : arr) {
        std::cout << a << " ";
    }
    // End the line after printing all elements.
    std::cout << std::endl;
}

/********************************************************************************************************************
 * FUNCTION: output (Coord version)
 * DESCRIPTION: Outputs the coordinate arrays (r, theta, phi) from a Coord object to separate text files.
 *              The radius values are normalized by con::cm.
 ********************************************************************************************************************/
void output(Coord const& coord, std::string const& filename, size_t precision) {
    // Open files for r, theta, and phi data.
    std::ofstream file_r(filename + "_t.txt");
    std::ofstream file_theta(filename + "_theta.txt");
    std::ofstream file_phi(filename + "_phi.txt");

    if (!file_r || !file_theta || !file_phi) {
        std::cerr << "Error opening files " << filename << "_*.txt\n";
        return;
    }

    file_r.precision(precision);
    file_theta.precision(precision);
    file_phi.precision(precision);

    // Write radius values normalized by con::cm.
    for (size_t i = 0; i < coord.t.size(); ++i) {
        file_r << coord.t[i] / con::sec << " ";
    }

    // Write theta values.
    for (size_t i = 0; i < coord.theta.size(); ++i) {
        file_theta << coord.theta[i] << " ";
    }

    // Write phi values.
    for (size_t i = 0; i < coord.phi.size(); ++i) {
        file_phi << coord.phi[i] << " ";
    }

    // End each file with a newline.
    file_r << "\n";
    file_theta << "\n";
    file_phi << "\n";
}

/********************************************************************************************************************
 * FUNCTION: writeGrid
 * DESCRIPTION: Writes multiple 3D grid data arrays to files. For each name provided in 'names', a file is created
 *              with that name appended to the base filename. The data are normalized by the corresponding unit.
 ********************************************************************************************************************/

template <typename T, typename U, typename V>
void writeGrid(std::string filename, T const& names, U const& data, V const& units, size_t precision) {
    for (size_t l = 0; l < names.size(); ++l) {
        std::ofstream file(filename + "_" + names[l] + ".txt");
        if (!file) {
            std::cerr << "Error opening file " << filename + "_" + names[l] + ".txt\n";
            return;
        }

        file.precision(precision);

        // Loop over the 3D grid dimensions and output the normalized values.
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

/********************************************************************************************************************
 * FUNCTION: output (Shock version)
 * DESCRIPTION: Outputs various fields of a Shock object to files by using writeGrid.
 *              The output includes Gamma_rel, B, t_com, t_eng, and column_num_den (labeled as "Sigma").
 ********************************************************************************************************************/
void output(Shock const& shock, std::string const& filename, size_t precision) {
    std::array<std::string, 7> strs = {"Gamma", "Gamma_rel", "B", "t_com", "r", "theta", "Sigma"};
    std::array<MeshGrid3d const*, 7> data = {&(shock.Gamma),         &(shock.Gamma_rel), &(shock.B),
                                             &(shock.t_com),         &(shock.r),         &(shock.theta),
                                             &(shock.column_num_den)};
    std::array<Real, 7> units = {1, 1, con::Gauss, con::sec, con::cm, 1, 1 / con::cm2};

    writeGrid(filename, strs, data, units, precision);
}

/********************************************************************************************************************
 * FUNCTION: output (PromptPhotonsGrid version)
 * DESCRIPTION: (Empty implementation) Intended to output a PromptPhotonsGrid to a file.
 ********************************************************************************************************************/
void output(PromptPhotonsGrid const& ph, std::string const& filename, size_t precision) {
    std::ofstream file(filename + ".txt");
    if (!file) {
        std::cerr << "Error opening file " << filename + ".txt\n";
        return;
    }
    file.precision(precision);
    for (size_t i = 0; i < ph.size(); ++i) {
        for (size_t j = 0; j < ph[i].size(); ++j) {
            for (size_t k = 0; k < ph[i][j].size(); ++k) {
                file << ph[i][j][k].E_nu_peak << " ";
            }
            file << '\n';
        }
        file << '\n';
    }
}

/********************************************************************************************************************
 * FUNCTION: output (SynPhotonGrid version)
 * DESCRIPTION: (Empty implementation) Intended to output a SynPhotonGrid to a file.
 ********************************************************************************************************************/
void output(SynPhotonGrid const& ph, std::string const& filename, size_t precision) {
    std::array<std::string, 3> strs = {"nu_a", "nu_m", "nu_c"};
    std::array<std::ofstream, 3> files;
    for (size_t i = 0; i < strs.size(); ++i) {
        files[i].open(filename + "_" + strs[i] + ".txt");
        if (!files[i]) {
            std::cerr << "Error opening file " << filename + "_" + strs[i] << ".txt\n";
            return;
        }
        files[i].precision(precision);
    }

    for (size_t i = 0; i < ph.size(); ++i) {
        for (size_t j = 0; j < ph[i].size(); ++j) {
            for (size_t k = 0; k < ph[i][j].size(); ++k) {
                files[0] << ph[i][j][k].nu_a / con::Hz << " ";
                files[1] << ph[i][j][k].nu_m / con::Hz << " ";
                files[2] << ph[i][j][k].nu_c / con::Hz << " ";
            }
            files[0] << '\n';
            files[1] << '\n';
            files[2] << '\n';
        }
        files[0] << '\n';
        files[1] << '\n';
        files[2] << '\n';
    }
}

/********************************************************************************************************************
 * FUNCTION: output (SynElectronGrid version)
 * DESCRIPTION: (Empty implementation) Intended to output a SynElectronGrid to a file.
 ********************************************************************************************************************/
void output(SynElectronGrid const& e, std::string const& filename, size_t precision) {
    std::array<std::string, 4> strs = {"gamma_a", "gamma_m", "gamma_c", "I_peak"};
    std::array<std::ofstream, 4> files;
    for (size_t i = 0; i < strs.size(); ++i) {
        files[i].open(filename + "_" + strs[i] + ".txt");
        if (!files[i]) {
            std::cerr << "Error opening file " << filename + "_" + strs[i] << ".txt\n";
            return;
        }
        files[i].precision(precision);
    }

    for (size_t i = 0; i < e.size(); ++i) {
        for (size_t j = 0; j < e[i].size(); ++j) {
            for (size_t k = 0; k < e[i][j].size(); ++k) {
                files[0] << e[i][j][k].gamma_a << " ";
                files[1] << e[i][j][k].gamma_m << " ";
                files[2] << e[i][j][k].gamma_c << " ";
                files[3] << e[i][j][k].I_nu_peak << " ";
            }
            files[0] << '\n';
            files[1] << '\n';
            files[2] << '\n';
            files[3] << '\n';
        }
        files[0] << '\n';
        files[1] << '\n';
        files[2] << '\n';
        files[3] << '\n';
    }
}

/********************************************************************************************************************
 * FUNCTION: output (MeshGrid3d version without unit)
 * DESCRIPTION: Outputs a 3D MeshGrid to a file. The filename is appended with ".txt".
 ********************************************************************************************************************/
void output(MeshGrid3d const& array, std::string const& filename, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Iterate over each element in the 3D grid and write to file.
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] << " ";
            }
            file << '\n';
        }
    }
}

/********************************************************************************************************************
 * FUNCTION: output (MeshGrid version without unit)
 * DESCRIPTION: Outputs a 2D MeshGrid to a file.
 ********************************************************************************************************************/
void output(MeshGrid const& grid, std::string const& filename, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Loop through the grid and output each element.
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] << " ";
        }
        file << '\n';
    }
}

/********************************************************************************************************************
 * FUNCTION: output (Array version without unit)
 * DESCRIPTION: Outputs a 1D Array to a file.
 ********************************************************************************************************************/
void output(Array const& array, std::string const& filename, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Output each element of the array.
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] << " ";
    }
}

/********************************************************************************************************************
 * FUNCTION: output (MeshGrid3d version with unit)
 * DESCRIPTION: Outputs a 3D MeshGrid to a file with each value normalized by the given unit.
 ********************************************************************************************************************/
void output(MeshGrid3d const& array, std::string const& filename, Real unit, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Normalize each element by 'unit' before output.
    for (size_t i = 0; i < array.size(); ++i) {
        for (size_t j = 0; j < array[i].size(); ++j) {
            for (size_t k = 0; k < array[i][j].size(); ++k) {
                file << array[i][j][k] / unit << " ";
            }
            file << '\n';
        }
    }
}

/********************************************************************************************************************
 * FUNCTION: output (MeshGrid version with unit)
 * DESCRIPTION: Outputs a 2D MeshGrid to a file with each value normalized by the given unit.
 ********************************************************************************************************************/
void output(MeshGrid const& grid, std::string const& filename, Real unit, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Output each normalized element.
    for (size_t i = 0; i < grid.size(); ++i) {
        for (size_t j = 0; j < grid[i].size(); ++j) {
            file << grid[i][j] / unit << " ";
        }
        file << '\n';
    }
}

/********************************************************************************************************************
 * FUNCTION: output (Array version with unit)
 * DESCRIPTION: Outputs a 1D Array to a file with each value normalized by the given unit.
 ********************************************************************************************************************/
void output(Array const& array, std::string const& filename, Real unit, size_t precision) {
    std::ofstream file(filename + ".txt");

    if (!file) {
        std::cerr << "Error opening file " << filename << ".txt\n";
        return;
    }
    file.precision(precision);
    // Output each element normalized by 'unit'.
    for (size_t i = 0; i < array.size(); ++i) {
        file << array[i] / unit << " ";
    }
}
