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

#ifndef VTK_DATA_OUTPUT_HH
#define VTK_DATA_OUTPUT_HH

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/vtkDataOutput.h"
#include "io/serializerIO.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace plb {


////////// struct VtkTypeNames ///////////////////////////////////

template<typename T>
class VtkTypeNames {
public:
    static std::string getName();
private:
    static std::string getBaseName();
};

template<typename T>
std::string VtkTypeNames<T>::getName() {
    std::stringstream sstream;
    sstream << getBaseName();
    sstream << 8 * sizeof(T);

    std::string tName;
    sstream >> tName;
    return tName;
}


////////// class VtkDataWriter3D ////////////////////////////////////////

template<typename T>
void VtkDataWriter3D::writeDataField(DataSerializer<T> const* serializer,
                                     std::string const& name, T scalingFactor, plint nDim)
{
    if (global::mpi().isMainProcessor()) {
        (*ostr) << "<DataArray type=\"" << VtkTypeNames<T>::getName()
                << "\" Name=\"" << name
                << "\" format=\"binary\" encoding=\"base64";
        if (nDim>1) {
            (*ostr) << "\" NumberOfComponents=\"" << nDim;
        }
        (*ostr) << "\">\n";
    }

    // Undocumented requirement of the vtk xml file format:
    // in front of every binary blob, base64 or raw-binary, appended or not, 
    // there is an UInt32 length indicator, giving the size of the binary blob in bytes;
    // when using base64 encoding, that length header must be encoded separately;
    // there must be no newline between the encoded length indicator and the encoded data block.
    //
    // those properties are properly handled by the serializer2ostr function, if pluint plint is
    // equal to UInt32. If not, you are on your own.

    bool enforceUint=true; // VTK uses "unsigned" to indicate the size of data, even on a 64-bit machine.
    serializerToBase64Stream<T>(new ScalingSerializer<T>(serializer, scalingFactor), *ostr, enforceUint);

    if (global::mpi().isMainProcessor()) {
        (*ostr) << "\n</DataArray>\n";
    }
}


////////// class VtkImageOutput2D ////////////////////////////////////

template<typename T>
VtkImageOutput2D<T>::VtkImageOutput2D(std::string fName, T deltaX_)
    : fullName ( global::directories().getVtkOutDir() + fName+".vti" ),
      vtkOut( fullName ),
      deltaX(deltaX_),
      offset(T(),T()),
      headerWritten( false )
{ }

template<typename T>
VtkImageOutput2D<T>::VtkImageOutput2D(std::string fName, T deltaX_, Array<T,2> offset_)
    : fullName ( global::directories().getVtkOutDir() + fName+".vti" ),
      vtkOut( fullName ),
      deltaX(deltaX_),
      offset(offset_),
      headerWritten( false )
{ }

template<typename T>
VtkImageOutput2D<T>::~VtkImageOutput2D() {
    writeFooter();
}

template<typename T>
void VtkImageOutput2D<T>::writeHeader(plint nx_, plint ny_) {
    if (headerWritten) {
        PLB_PRECONDITION(nx == nx_);
        PLB_PRECONDITION(ny == ny_);
    }
    else {
        nx = nx_;
        ny = ny_;
        vtkOut.writeHeader(Box3D(0,nx-1,0,ny-1,0,0), Array<T,3>(offset[0],offset[1],T()), deltaX);
        vtkOut.startPiece(Box3D(0,nx-1,0,ny-1,0,0));
        headerWritten = true;
    }
}

template<typename T>
void VtkImageOutput2D<T>::writeFooter() {
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template<typename T>
template<typename TConv>
void VtkImageOutput2D<T>::writeData( ScalarFieldBase2D<T> const& scalarField,
                                     std::string scalarFieldName, TConv scalingFactor )
{
    writeHeader(scalarField.getNx(), scalarField.getNy());
    vtkOut.writeDataField (
            new TypeConversionSerializer<T,TConv> (
                scalarField.getBlockSerializer(scalarField.getBoundingBox(), IndexOrdering::backward) ),
            scalarFieldName, scalingFactor, 1);
}

template<typename T>
template<plint n, typename TConv>
void VtkImageOutput2D<T>::writeData( TensorFieldBase2D<T,n> const& tensorField,
                                     std::string tensorFieldName, TConv scalingFactor )
{
    writeHeader(tensorField.getNx(), tensorField.getNy());
    vtkOut.writeDataField (
            new TypeConversionSerializer<T,TConv> (
                tensorField.getBlockSerializer(tensorField.getBoundingBox(), IndexOrdering::backward) ),
            tensorFieldName, scalingFactor, n);
}

////////// class VtkImageOutput3D ////////////////////////////////////

template<typename T>
VtkImageOutput3D<T>::VtkImageOutput3D(std::string fName, T deltaX_)
    : fullName ( global::directories().getVtkOutDir() + fName+".vti" ),
      vtkOut( fullName ),
      deltaX(deltaX_),
      offset(T(),T(),T()),
      headerWritten( false )
{ }

template<typename T>
VtkImageOutput3D<T>::VtkImageOutput3D(std::string fName, T deltaX_, Array<T,3> offset_)
    : fullName ( global::directories().getVtkOutDir() + fName+".vti" ),
      vtkOut( fullName ),
      deltaX(deltaX_),
      offset(offset_),
      headerWritten( false )
{ }

template<typename T>
VtkImageOutput3D<T>::~VtkImageOutput3D() {
    writeFooter();
}

template<typename T>
void VtkImageOutput3D<T>::writeHeader(plint nx_, plint ny_, plint nz_) {
    if (headerWritten) {
        PLB_PRECONDITION(nx == nx_);
        PLB_PRECONDITION(ny == ny_);
        PLB_PRECONDITION(nz == nz_);
    }
    else {
        nx = nx_;
        ny = ny_;
        nz = nz_;
        vtkOut.writeHeader(Box3D(0,nx-1,0,ny-1,0,nz-1), offset, deltaX);
        vtkOut.startPiece(Box3D(0,nx-1,0,ny-1,0,nz-1));
        headerWritten = true;
    }
}

template<typename T>
void VtkImageOutput3D<T>::writeFooter() {
    if (headerWritten) {
        vtkOut.endPiece();
        vtkOut.writeFooter();
        headerWritten = false;
    }
}

template<typename T>
template<typename TConv>
void VtkImageOutput3D<T>::writeData( ScalarFieldBase3D<T> const& scalarField,
                                     std::string scalarFieldName, TConv scalingFactor )
{
    writeHeader(scalarField.getNx(), scalarField.getNy(), scalarField.getNz());
    vtkOut.writeDataField (
            new TypeConversionSerializer<T,TConv> (
                scalarField.getBlockSerializer(scalarField.getBoundingBox(), IndexOrdering::backward) ),
            scalarFieldName, scalingFactor, 1 );
}

template<typename T>
template<plint n, typename TConv>
void VtkImageOutput3D<T>::writeData( TensorFieldBase3D<T,n> const& tensorField,
                                     std::string tensorFieldName, TConv scalingFactor )
{
    writeHeader(tensorField.getNx(), tensorField.getNy(), tensorField.getNz());
    vtkOut.writeDataField (
            new TypeConversionSerializer<T,TConv> (
                tensorField.getBlockSerializer(tensorField.getBoundingBox(), IndexOrdering::backward) ),
            tensorFieldName, scalingFactor, n );
}


////////// Free Functions //////////////////////////////////////////////

template<typename T> void writeVTKData3D (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarFieldBase3D<T> const& scalarField,
        std::string const& vectorFieldName,
        TensorFieldBase3D<T,3> const& vectorField,
        T deltaX, T deltaT )
{
    std::string fullName = global::directories().getVtkOutDir() + fName+".vti";
    VtkDataWriter3D vtkOut(fullName);
    plint nx = scalarField.getNx();
    plint ny = scalarField.getNy();
    plint nz = scalarField.getNz();
    vtkOut.writeHeader(Box3D(0,nx-1,0,ny-1,0,nz-1), Array<T,3>(T(),T(),T()), deltaX);

    vtkOut.startPiece(Box3D(0,nx-1,0,ny-1,0,nz-1));
    vtkOut.writeDataField (
            scalarField.getBlockSerializer(scalarField.getBoundingBox(),IndexOrdering::backward),
            scalarFieldName, 1./deltaT, 1 );

    vtkOut.writeDataField (
            vectorField.getBlockSerializer(scalarField.getBoundingBox(),IndexOrdering::backward),
            vectorFieldName, deltaX/deltaT, 3 );
    vtkOut.endPiece();

    vtkOut.writeFooter();
}


}  // namespace plb

#endif


