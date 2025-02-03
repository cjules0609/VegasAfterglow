//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#ifndef _IO_H_
#define _IO_H_

#include "mesh.h"
#include "prompt.h"
#include "shock.h"
#include "synchrotron.h"
void printArray(Array const& arr);
void output(SynPhotonGrid const& syn_ph, std::string const& filename);
void output(SynElectronGrid const& syn_e, std::string const& filename);
void output(PromptPhotonsGrid const& prompt_pj, std::string const& filename);
void output(Shock const& shock, std::string const& filename);
void output(Coord const& coord, std::string const& filename);
void output(MeshGrid3d const& array, std::string const& filename);
void output(MeshGrid const& grid, std::string const& filename);
void output(Array const& array, std::string const& filename);
void output(MeshGrid3d const& array, std::string const& filename, double unit);
void output(MeshGrid const& grid, std::string const& filename, double unit);
void output(Array const& array, std::string const& filename, double unit);

inline void print() {  // Base case
    std::cout << std::endl;
}

template <typename First, typename... Rest>
void print(First first, Rest... rest) {
    std::cout << first << " ";
    print(rest...);  // Recursive call
}
#endif