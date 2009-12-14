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
 * Base class for scalar, vector and tensor fields for 3D data analysis -- header file.
 */

#ifndef DATA_FIELD_BASE_3D_H
#define DATA_FIELD_BASE_3D_H

#include "core/globalDefs.h"
#include "core/dataFieldBase2D.h"
#include "core/block3D.h"
#include "core/geometry3D.h"
#include "core/array.h"

namespace plb {

/// Interface for the variants of 3D scalar, vector and tensor fields.
template<typename T>
class ScalarFieldBase3D : virtual public Block3D<T> {
public:
    virtual ~ScalarFieldBase3D() { }
public:
    virtual void reset() =0;
    virtual T& get(plint iX, plint iY, plint iZ) =0;
    virtual T const& get(plint iX, plint iY, plint iZ) const =0;
};

template<typename T, int nDim>
class TensorFieldBase3D : virtual public Block3D<T> {
public:
    virtual ~TensorFieldBase3D() { }
public:
    virtual void reset() =0;
    virtual Array<T,nDim>& get(plint iX, plint iY, plint iZ) =0;
    virtual Array<T,nDim> const& get(plint iX, plint iY, plint iZ) const =0;
};

}


#endif
