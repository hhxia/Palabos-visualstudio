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
 * Helper functions for domain initialization -- header file.
 */
#ifndef DATA_ANALYSIS_3D_H
#define DATA_ANALYSIS_3D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase3D.h"
#include "core/dataFieldBase3D.h"


namespace plb {

/* *************** Data-analysis functions for BlockLattice ********** */

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase3D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase3D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase3D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, BoolMask boolMask);

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase3D<T,Descriptor>& lattice, BoolMask boolMask);


/* *************** Data-analysis functions for ScalarField *********** */

template<typename T>
T computeAverage(ScalarFieldBase3D<T>& scalarField, Box3D domain);

template<typename T>
T computeAverage(ScalarFieldBase3D<T>& scalarField);


template<typename T>
T computeMin(ScalarFieldBase3D<T>& scalarField, Box3D domain);

template<typename T>
T computeMin(ScalarFieldBase3D<T>& scalarField);


template<typename T>
T computeMax(ScalarFieldBase3D<T>& scalarField, Box3D domain);

template<typename T>
T computeMax(ScalarFieldBase3D<T>& scalarField);


template<typename T>
T computeBoundedAverage(ScalarFieldBase3D<T>& scalarField, Box3D domain);

template<typename T>
T computeBoundedAverage(ScalarFieldBase3D<T>& scalarField);


template<typename T, class BoolMask> 
plint count(ScalarFieldBase3D<T>& field, Box3D domain, BoolMask boolMask);

template<typename T, class BoolMask> 
plint count(ScalarFieldBase3D<T>& field, BoolMask boolMask);


/* *************** Data-analysis functions for TensorField *********** */

template<typename T, int nDim, class BoolMask> 
plint count(TensorFieldBase3D<T,nDim>& field, Box3D domain, BoolMask boolMask);

template<typename T, int nDim, class BoolMask> 
plint count(TensorFieldBase3D<T,nDim>& field, BoolMask boolMask);


}  // namespace plb

#endif  // DATA_ANALYSIS_3D_H
