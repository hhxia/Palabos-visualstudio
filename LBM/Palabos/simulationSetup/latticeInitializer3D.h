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
#ifndef LATTICE_INITIALIZER_3D_H
#define LATTICE_INITIALIZER_3D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "core/dynamics.h"
#include "simulationSetup/latticeInitializerFunctionals3D.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
void apply(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, OneCellFunctional3D<T,Descriptor>* f);

template<typename T, template<class U> class Descriptor>
void applyIndexed(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, OneCellIndexedFunctional3D<T,Descriptor>* f);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D boundingBox,
                    DomainFunctional3D* domain, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase3D<T,Descriptor>& lattice, DotList3D const& dotList, Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase3D<T,Descriptor>& lattice, plint iX, plint iY, plint iZ,
                    Dynamics<T,Descriptor>* dynamics);

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, Array<T,3> velocity);

template<typename T, template<class U> class Descriptor, class Function>
void setBoundaryVelocity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, Function f);

template<typename T, template<class U> class Descriptor>
void setBoundaryDensity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, T rho);

template<typename T, template<class U> class Descriptor, class Function>
void setBoundaryDensity(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, Function f);

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, T rho, Array<T,3> velocity);

template<typename T, template<class U> class Descriptor, class Function>
void initializeAtEquilibrium(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, Function f);

template<typename T, template<class U> class Descriptor>
void stripeOffDensityOffset(BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain, T deltaRho);


template<typename T, template<class U> class Descriptor>
void setCompositeDynamics( BlockLatticeBase3D<T,Descriptor>& lattice, Box3D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics );

}  // namespace plb

#endif  // LATTICE_INITIALIZER_3D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "simulationSetup/latticeInitializerGenerics3D.h"
