//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <string>

#ifndef _WIN32
#include "xtensor-io/xnpz.hpp"
#endif

#include "mesh.h"
#include "prompt.h"
#include "shock.h"
#include "synchrotron.h"
/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Output and Printing Functions
 * DESCRIPTION: These function prototypes declare various utility routines for printing and outputting data.
 *              They include functions to print a 1D Array to the console and to output different types of grids
 *              (SynPhotonGrid, SynElectronGrid, PromptPhotonsGrid, Shock, Coord, MeshGrid3d, MeshGrid, and Array)
 *              to a file, optionally with a specified unit.
 ********************************************************************************************************************/
// Write Array to a CSV file with an optional unit for value scaling
void write_csv(std::string const& filename, Array const& array, Real unit = 1.0);
// Write MeshGrid to a CSV file with an optional unit for value scaling
void write_csv(std::string const& filename, MeshGrid const& grid, Real unit = 1.0);
// Write MeshGrid3d to a CSV file with an optional unit for value scaling
void write_csv(std::string const& filename, MeshGrid3d const& grid3d, Real unit = 1.0);

#ifndef _WIN32
// Write SynPhotonGrid to an NPZ file
void write_npz(std::string const& filename, SynPhotonGrid const& syn_ph);
// Write SynElectronGrid to an NPZ file
void write_npz(std::string const& filename, SynElectronGrid const& syn_e);
// Write Shock to an NPZ file
void write_npz(std::string const& filename, Shock const& shock);
// Write Coord to an NPZ file
void write_npz(std::string const& filename, Coord const& coord);

// Write an array to an NPZ file with an optional unit for value scaling
template <typename T>
void write_npz(std::string const& filename, const T& array, Real unit = 1.0) {
    xt::dump_npz(filename + ".npz", "array", xt::eval(array / unit), false, false);
}

// Helper recursive function for writing multiple arrays to a single NPZ file
template <typename T, typename... Rest>
void write_npz_recursive(std::string const& filename, bool first, std::string const& name, const T& array,
                         const Rest&... rest) {
    auto arr = xt::eval(array);                                 // ensure evaluated
    xt::dump_npz(filename + ".npz", name, arr, false, !first);  // append after first write

    if constexpr (sizeof...(rest) > 0) {
        write_npz_recursive(filename, false, rest...);  // continue with rest
    }
}

// Write multiple named arrays to a single NPZ file
// Usage: write_npz("filename", "name1", array1, "name2", array2, ...)
template <typename... Args>
void write_npz(std::string const& filename, Args const&... args) {
    static_assert(sizeof...(args) % 2 == 0, "Arguments must be pairs: name1, array1, name2, array2, ...");

    write_npz_recursive(filename, true, args...);
}
#endif

// Write an array to an NPY file with an optional unit for value scaling
template <typename T>
void write_npy(std::string const& filename, const T& array, Real unit = 1.0) {
    xt::dump_npy(filename + ".npy", xt::eval(array / unit));
}