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
    
template
void VtkDataWriter3D::writeDataField (
        DataSerializer<char> const* serializer, std::string const& name, 
        char scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<unsigned char> const* serializer, std::string const& name, 
        unsigned char scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<short int> const* serializer, std::string const& name, 
        short int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<unsigned short int> const* serializer, std::string const& name, 
        unsigned short int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<int> const* serializer, std::string const& name, 
        int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<unsigned int> const* serializer, std::string const& name, 
        unsigned int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<long int> const* serializer, std::string const& name, 
        long int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<unsigned long int> const* serializer, std::string const& name, 
        unsigned long int scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<float> const* serializer, std::string const& name, 
        float scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<double> const* serializer, std::string const& name, 
        double scalingFactor, plint nDim );

template
void VtkDataWriter3D::writeDataField (
        DataSerializer<long double> const* serializer, std::string const& name, 
        long double scalingFactor, plint nDim );
        
template void writeVTKData3D<double> (
        std::string const& fName,
        std::string const& scalarFieldName,
        ScalarFieldBase3D<double> const& scalarField,
        std::string const& vectorFieldName,
        TensorFieldBase3D<double,3> const& vectorField,
        double deltaX, double deltaT );

template class VtkImageOutput2D<double>;

template void VtkImageOutput2D<double>::writeData<double> (
        ScalarFieldBase2D<double> const& scalarField,
        std::string scalarFieldName, double scalingFactor );

template void VtkImageOutput2D<double>::writeData<float> (
        ScalarFieldBase2D<double> const& scalarField,
        std::string scalarFieldName, float scalingFactor );

template void VtkImageOutput2D<double>::writeData<2, double> (
        TensorFieldBase2D<double,2> const& tensorField,
        std::string tensorFieldName, double scalingFactor );

template void VtkImageOutput2D<double>::writeData<2, float> (
        TensorFieldBase2D<double,2> const& tensorField,
        std::string tensorFieldName, float scalingFactor );

template void VtkImageOutput2D<double>::writeData<9, double> (
        TensorFieldBase2D<double,9> const& tensorField,
        std::string tensorFieldName, double scalingFactor );

template void VtkImageOutput2D<double>::writeData<9, float> (
        TensorFieldBase2D<double,9> const& tensorField,
        std::string tensorFieldName, float scalingFactor );


template class VtkImageOutput3D<double>;

template void VtkImageOutput3D<double>::writeData<double> (
        ScalarFieldBase3D<double> const& scalarField,
        std::string scalarFieldName, double scalingFactor );

template void VtkImageOutput3D<double>::writeData<float> (
        ScalarFieldBase3D<double> const& scalarField,
        std::string scalarFieldName, float scalingFactor );

template void VtkImageOutput3D<double>::writeData<3, double> (
        TensorFieldBase3D<double,3> const& tensorField,
        std::string tensorFieldName, double scalingFactor );

template void VtkImageOutput3D<double>::writeData<3, float> (
        TensorFieldBase3D<double,3> const& tensorField,
        std::string tensorFieldName, float scalingFactor );

template void VtkImageOutput3D<double>::writeData<19, double> (
        TensorFieldBase3D<double,19> const& tensorField,
        std::string tensorFieldName, double scalingFactor );

template void VtkImageOutput3D<double>::writeData<19, float> (
        TensorFieldBase3D<double,19> const& tensorField,
        std::string tensorFieldName, float scalingFactor );

}
