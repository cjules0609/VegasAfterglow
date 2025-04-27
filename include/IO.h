//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

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

#include <string>

void write_csv(std::string const& filename, Array const& array, Real unit = 1.0);
void write_csv(std::string const& filename, MeshGrid const& grid, Real unit = 1.0);
void write_csv(std::string const& filename, MeshGrid3d const& grid3d, Real unit = 1.0);
void write_npz(std::string const& filename, SynPhotonGrid const& syn_ph);
void write_npz(std::string const& filename, SynElectronGrid const& syn_e);
void write_npz(std::string const& filename, Shock const& shock);
void write_npz(std::string const& filename, Coord const& coord);

template <typename T>
void write_npy(std::string const& filename, const T& array, Real unit = 1.0) {
    xt::dump_npy(filename + ".npy", xt::eval(array / unit));
}

template <typename T>
void write_npz(std::string const& filename, const T& array, Real unit = 1.0) {
    xt::dump_npz(filename + ".npz", "array", xt::eval(array / unit), false, false);
}

template <typename T, typename... Rest>
void write_npz_recursive(std::string const& filename, bool first, std::string const& name, const T& array,
                         const Rest&... rest) {
    auto arr = xt::eval(array);                                 // ensure evaluated
    xt::dump_npz(filename + ".npz", name, arr, false, !first);  // append after first write

    if constexpr (sizeof...(rest) > 0) {
        write_npz_recursive(filename, false, rest...);  // continue with rest
    }
}

template <typename... Args>
void write_npz(std::string const& filename, Args const&... args) {
    static_assert(sizeof...(args) % 2 == 0, "Arguments must be pairs: name1, array1, name2, array2, ...");

    write_npz_recursive(filename, true, args...);
}
