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
 * Helper functions for domain initialization -- generic implementation.
 */
#include "core/globalDefs.h"
#include "multiBlock/multiLatticeInitializer3D.h"
#include "multiBlock/multiDataCouplingWrapper3D.h"
#include "simulationSetup/latticeInitializerFunctionals3D.h"

#ifndef MULTI_LATTICE_INITIALIZER_3D_HH
#define MULTI_LATTICE_INITIALIZER_3D_HH

namespace plb {

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& boolMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromMaskFunctional3D<T,Descriptor>(dynamics,whichFlag),
            domain, lattice, boolMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& boolMask,
                     Dynamics<T,Descriptor>* dynamics, bool whichFlag )
{
    defineDynamics(lattice, boolMask, lattice.getBoundingBox(), dynamics, whichFlag);
}


template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& intMask,
                     Box3D domain, Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    applyProcessingFunctional (
            new DynamicsFromIntMaskFunctional3D<T,Descriptor>(dynamics,whichFlag),
            domain, lattice, intMask );
}

template<typename T, template<typename U> class Descriptor>
void defineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& intMask,
                     Dynamics<T,Descriptor>* dynamics, int whichFlag )
{
    defineDynamics(lattice, intMask, lattice.getBoundingBox(), dynamics, whichFlag);
}

}  // namespace plb

#endif  // MULTI_LATTICE_INITIALIZER_3D_HH
