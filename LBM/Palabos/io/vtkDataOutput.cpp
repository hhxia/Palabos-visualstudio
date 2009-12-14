/* This file is part of the Palabos library.
 * Copyright (C) 2009 Jonas Latt
 * E-mail contact: jonas@lbmethod.org
 * The most recent release of Palabos can be downloaded at 
 * <http://www.lbmethod.org/palabos/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "parallelism/mpiManager.h"
#include "core/serializer.h"
#include "core/serializer.hh"
#include "io/vtkDataOutput.h"
#include "io/vtkDataOutput.hh"
#include "io/serializerIO.h"
#include "io/serializerIO.hh"
#include "io/base64.h"
#include "io/base64.hh"

namespace plb {
    
////////// class VtkDataWriter3D ////////////////////////////////////////

VtkDataWriter3D::VtkDataWriter3D(std::string const& fileName_)
    : fileName(fileName_),
      ostr(0)
{
    if (global::mpi().isMainProcessor()) {
        ostr = new std::ofstream(fileName.c_str());
        if (!(*ostr)) {
            std::cerr << "could not open file " <<  fileName << "\n";
            return;
        }
    }
}

VtkDataWriter3D::~VtkDataWriter3D() {
    delete ostr;
}

void VtkDataWriter3D::writeHeader(Box3D domain, Array<double,3> origin, double deltaX)
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "<?xml version=\"1.0\"?>\n";
        (*ostr) << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        (*ostr) << "<ImageData WholeExtent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\" "
                << "Origin=\""
                << origin[0] << " " << origin[1] << " " << origin[2] << "\" "
                << "Spacing=\""
                << deltaX << " " << deltaX << " " << deltaX << "\">\n";
    }
}

void VtkDataWriter3D::startPiece(Box3D domain) {
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "<Piece Extent=\""
                << domain.x0 << " " << domain.x1 << " "
                << domain.y0 << " " << domain.y1 << " "
                << domain.z0 << " " << domain.z1 << "\">\n";
        (*ostr) << "<PointData>\n";
    }
}

void VtkDataWriter3D::endPiece() {
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "</PointData>\n";
        (*ostr) << "</Piece>\n";
    }
}

void VtkDataWriter3D::writeFooter() {
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "</ImageData>\n";
        (*ostr) << "</VTKFile>\n";
    }
}

template<>
std::string VtkTypeNames<bool>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<char>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned char>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<short int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned short int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<long int>::getBaseName() {
    return "Int";
}

template<>
std::string VtkTypeNames<unsigned long int>::getBaseName() {
    return "UInt";
}

template<>
std::string VtkTypeNames<float>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<double>::getBaseName() {
    return "Float";
}

template<>
std::string VtkTypeNames<long double>::getBaseName() {
    return "Float";
}

}
