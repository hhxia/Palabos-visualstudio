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

/** \file
 * Helper functions for data field initialization -- header file.
 */
#ifndef DATA_FIELD_INITIALIZER_3D_H
#define DATA_FIELD_INITIALIZER_3D_H

#include "core/globalDefs.h"
#include "core/dataFieldBase3D.h"
#include "atomicBlock/dataField3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"

namespace plb {

/// Initialize scalar-field with the same constant value on each cell.
template<typename T>
void setToConstant(ScalarFieldBase3D<T>& field, Box3D domain, T value);

/// Initialize tensor-field with the same constant tensor/vector on each cell.
template<typename T, int nDim>
void setToConstant( TensorFieldBase3D<T,nDim>& field, Box3D domain, 
                    Array<T,nDim> const& value );

/// Initialize scalar-field with the a value from a function.
template<typename T, class Function>
void setToFunction(ScalarFieldBase3D<T>& field, Box3D domain, Function f);

/// Initialize tensor-field with a vector/tensor value from a function.
template<typename T, int nDim, class Function>
void setToFunction( TensorFieldBase3D<T,nDim>& field, Box3D domain, 
                    Function f );

/// Assign the component "index" of its space coordinate to each cell.
template<typename T>
void setToCoordinate(ScalarFieldBase3D<T>& field, Box3D domain, plint index);

/// Assign its space coordinate to each cell.
template<typename T>
void setToCoordinates(TensorFieldBase3D<T,3>& field,  Box3D domain);

}  // namespace plb

#endif  // DATA_FIELD_INITIALIZER_3D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "simulationSetup/dataFieldInitializerGenerics2D.h"
