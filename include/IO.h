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
/********************************************************************************************************************
 * FUNCTION PROTOTYPES: Output and Printing Functions
 * DESCRIPTION: These function prototypes declare various utility routines for printing and outputting data.
 *              They include functions to print a 1D Array to the console and to output different types of grids
 *              (SynPhotonGrid, SynElectronGrid, PromptPhotonsGrid, Shock, Coord, MeshGrid3d, MeshGrid, and Array)
 *              to a file, optionally with a specified unit.
 ********************************************************************************************************************/
void printArray(Array const& arr);  // Prints the elements of a 1D Array.

void output(SynPhotonGrid const& syn_ph, std::string const& filename);         // Outputs a SynPhotonGrid to a file.
void output(SynElectronGrid const& syn_e, std::string const& filename);        // Outputs a SynElectronGrid to a file.
void output(PromptPhotonsGrid const& prompt_pj, std::string const& filename);  // Outputs a PromptPhotonsGrid to a file.
void output(Shock const& shock, std::string const& filename);                  // Outputs a Shock object to a file.
void output(Coord const& coord, std::string const& filename);                  // Outputs a Coord object to a file.
void output(MeshGrid3d const& array, std::string const& filename);             // Outputs a 3D MeshGrid to a file.
void output(MeshGrid const& grid, std::string const& filename);                // Outputs a 2D MeshGrid to a file.
void output(Array const& array, std::string const& filename);                  // Outputs a 1D Array to a file.
void output(MeshGrid3d const& array, std::string const& filename, Real unit);  // Outputs a 3D MeshGrid with a unit.
void output(MeshGrid const& grid, std::string const& filename, Real unit);     // Outputs a 2D MeshGrid with a unit.
void output(Array const& array, std::string const& filename, Real unit);       // Outputs a 1D Array with a unit.

/********************************************************************************************************************
 * FUNCTION: print (Base Case)
 * DESCRIPTION: A base-case overload for the variadic print function. It simply outputs a newline.
 ********************************************************************************************************************/
inline void print() {  // Base case: terminates the recursion.
    std::cout << std::endl;
}

/********************************************************************************************************************
 * TEMPLATE FUNCTION: print (Variadic)
 * DESCRIPTION: Prints multiple arguments to standard output, separated by a space.
 *              This template function recursively prints all provided arguments.
 ********************************************************************************************************************/
template <typename First, typename... Rest>
void print(First first, Rest... rest) {
    std::cout << first << " ";
    print(rest...);  // Recursive call to print the remaining arguments.
}
#endif