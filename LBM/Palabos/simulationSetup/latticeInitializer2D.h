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
#ifndef LATTICE_INITIALIZER_2D_H
#define LATTICE_INITIALIZER_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "core/dynamics.h"
#include "simulationSetup/latticeInitializerFunctionals2D.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
void apply(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, OneCellFunctional2D<T,Descriptor>* f);

template<typename T, template<class U> class Descriptor>
void applyIndexed(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, OneCellIndexedFunctional2D<T,Descriptor>* f);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D boundingBox,
                    DomainFunctional2D* domain, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, plint iX, plint iY, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, DotList2D const& dotList, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Array<T,2> velocity);

template<typename T, template<class U> class Descriptor, class Function>
void setBoundaryVelocity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Function f);

template<typename T, template<class U> class Descriptor>
void setBoundaryDensity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T rho);

template<typename T, template<class U> class Descriptor, class Function>
void setBoundaryDensity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Function f);

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T density, Array<T,2> velocity);

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice,
                             Box2D boundingBox, DomainFunctional2D* domain,
                             T density, Array<T,2> velocity);

template<typename T, template<class U> class Descriptor, class Function>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Function f);

template<typename T, template<class U> class Descriptor, class Function>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice,
                             Box2D boundingBox, DomainFunctional2D* domain, Function f);

template<typename T, template<class U> class Descriptor>
void stripeOffDensityOffset(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T deltaRho);


template<typename T, template<class U> class Descriptor>
void setCompositeDynamics( BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics );

}  // namespace plb

#endif  // LATTICE_INITIALIZER_2D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "simulationSetup/latticeInitializerGenerics2D.h"
