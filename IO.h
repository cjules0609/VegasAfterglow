#ifndef _IO_H_
#define _IO_H_

#include "mesh.h"
#include "synchrotron.h"
void printArray(Array const& arr);
void write2file(SynRadMesh const& syn_rad, std::string const& filename);
void write2file(MeshGrid3d const& array, std::string const& filename);
void write2file(MeshGrid const& grid, std::string const& filename);
void write2file(Array const& array, std::string const& filename);

#endif