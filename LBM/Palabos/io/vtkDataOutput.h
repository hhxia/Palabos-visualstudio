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


#ifndef VTK_DATA_OUTPUT_H
#define VTK_DATA_OUTPUT_H

#include "core/globalDefs.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "core/serializer.h"
#include "core/dataFieldBase3D.h"
#include "core/array.h"

namespace plb {

class VtkDataWriter3D {
public:
    VtkDataWriter3D(std::string const& fileName_);
    ~VtkDataWriter3D();
    void writeHeader(Box3D domain, Array<double,3> origin, double deltaX);
    void startPiece(Box3D domain);
    void endPiece();
    void writeFooter();
    template <typename T>
    void writeDataField( DataSerializer<T> const* serializer,
                         std::string const& name, T scalingFactor, plint nDim );
private:
    VtkDataWriter3D(VtkDataWriter3D const& rhs);
    VtkDataWriter3D operator=(VtkDataWriter3D const& rhs);
private:
    std::string fileName;
    std::ofstream *ostr;
};

template<typename T>
class VtkImageOutput2D {
public:
    VtkImageOutput2D(std::string fName, T deltaX_=(T)1);
    VtkImageOutput2D(std::string fName, T deltaX_, Array<T,2> offset);
    ~VtkImageOutput2D();
    template<typename TConv>
    void writeData(ScalarFieldBase2D<T> const& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1);
    template<plint n, typename TConv>
    void writeData(TensorFieldBase2D<T,n> const& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
private:
    void writeHeader(plint nx_, plint ny_);
    void writeFooter();
private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    T deltaX;
    Array<T,2> offset;
    bool headerWritten;
    plint nx, ny;
};

template<typename T>
class VtkImageOutput3D {
public:
    VtkImageOutput3D(std::string fName, T deltaX_=(T)1);
    VtkImageOutput3D(std::string fName, T deltaX_, Array<T,3> offset);
    ~VtkImageOutput3D();
    template<typename TConv>
    void writeData(ScalarFieldBase3D<T> const& scalarField,
                   std::string scalarFieldName, TConv scalingFactor=(T)1);
    template<plint n, typename TConv>
    void writeData(TensorFieldBase3D<T,n> const& tensorField,
                   std::string tensorFieldName, TConv scalingFactor=(T)1);
private:
    void writeHeader(plint nx_, plint ny_, plint nz_);
    void writeFooter();
private:
    std::string fullName;
    VtkDataWriter3D vtkOut;
    T deltaX;
    Array<T,3> offset;
    bool headerWritten;
    plint nx, ny, nz;
};


template<typename T> void writeVTKData3D (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarFieldBase3D<T> const& scalarField,
        std::string const& vectorFieldName,
        TensorFieldBase3D<T,3> const& vectorField,
        T deltaX, T deltaT );

} // namespace plb

#endif
