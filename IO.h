#ifndef _IO_H_
#define _IO_H_

#include "mesh.h"
#include "shock.h"
#include "synchrotron.h"
void print_array(Array const& arr);
void write2file(SynPhotonsMesh const& syn_ph, std::string const& filename);
void write2file(SynElectronsMesh const& syn_e, std::string const& filename);

void write2file(Shock const& shock, std::string const& filename);
void write2file(Coord const& coord, std::string const& filename);

void write2file(MeshGrid3d const& array, std::string const& filename);
void write2file(MeshGrid const& grid, std::string const& filename);
void write2file(Array const& array, std::string const& filename);


void write2file(MeshGrid3d const& array, std::string const& filename, double unit);
void write2file(MeshGrid const& grid, std::string const& filename, double unit);
void write2file(Array const& array, std::string const& filename, double unit);
#endif