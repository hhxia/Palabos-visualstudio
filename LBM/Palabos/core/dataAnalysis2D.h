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
#ifndef DATA_ANALYSIS_2D_H
#define DATA_ANALYSIS_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "core/dataFieldBase2D.h"


namespace plb {

/* *************** Data-analysis functions for BlockLattice ********** */

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageDensity(BlockLatticeBase2D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageRhoBar(BlockLatticeBase2D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain);

template<typename T, template<typename U> class Descriptor> 
T computeAverageEnergy(BlockLatticeBase2D<T,Descriptor>& lattice);


template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, BoolMask boolMask);

template<typename T, template<typename U> class Descriptor, class BoolMask> 
plint count(BlockLatticeBase2D<T,Descriptor>& lattice, BoolMask boolMask);


/* *************** Data-analysis functions for ScalarField *********** */

template<typename T>
T computeAverage(ScalarFieldBase2D<T>& scalarField, Box2D domain);

template<typename T>
T computeAverage(ScalarFieldBase2D<T>& scalarField);


template<typename T>
T computeMin(ScalarFieldBase2D<T>& scalarField, Box2D domain);

template<typename T>
T computeMin(ScalarFieldBase2D<T>& scalarField);


template<typename T>
T computeMax(ScalarFieldBase2D<T>& scalarField, Box2D domain);

template<typename T>
T computeMax(ScalarFieldBase2D<T>& scalarField);


template<typename T>
T computeBoundedAverage(ScalarFieldBase2D<T>& scalarField, Box2D domain);

template<typename T>
T computeBoundedAverage(ScalarFieldBase2D<T>& scalarField);


template<typename T, class BoolMask> 
plint count(ScalarFieldBase2D<T>& field, Box2D domain, BoolMask boolMask);

template<typename T, class BoolMask> 
plint count(ScalarFieldBase2D<T>& field, BoolMask boolMask);


/* *************** Data-analysis functions for TensorField *********** */

template<typename T, int nDim, class BoolMask> 
plint count(TensorFieldBase2D<T,nDim>& field, Box2D domain, BoolMask boolMask);

template<typename T, int nDim, class BoolMask> 
plint count(TensorFieldBase2D<T,nDim>& field, BoolMask boolMask);

}  // namespace plb

#endif  // DATA_ANALYSIS_2D_H
