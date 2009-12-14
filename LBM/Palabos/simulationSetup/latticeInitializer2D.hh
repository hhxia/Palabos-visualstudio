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
#ifndef LATTICE_INITIALIZER_2D_HH
#define LATTICE_INITIALIZER_2D_HH

#include "simulationSetup/latticeInitializer2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
void apply(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, OneCellFunctional2D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericLatticeFunctional2D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void applyIndexed(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, OneCellIndexedFunctional2D<T,Descriptor>* f) {
    applyProcessingFunctional(new GenericIndexedLatticeFunctional2D<T,Descriptor>(f), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateDynamicsFunctional2D<T,Descriptor>(dynamics), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D boundingBox,
                    DomainFunctional2D* domain, Dynamics<T,Descriptor>* dynamics) {
    applyProcessingFunctional (
        new InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>(dynamics, domain),
        boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice,
                    DotList2D const& dotList, Dynamics<T,Descriptor>* dynamics)
{
    applyProcessingFunctional (
        new InstantiateDotDynamicsFunctional2D<T,Descriptor>(dynamics), dotList, lattice );
}

template<typename T, template<class U> class Descriptor>
void defineDynamics(BlockLatticeBase2D<T,Descriptor>& lattice,
                    plint iX, plint iY, Dynamics<T,Descriptor>* dynamics)
{
    DotList2D pos; pos.addDot(Dot2D(iX,iY));
    defineDynamics(lattice, pos, dynamics);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, Array<T,2> velocity) {
    applyProcessingFunctional(new SetConstBoundaryVelocityFunctional2D<T,Descriptor>(velocity), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setBoundaryDensity(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T rho) {
    applyProcessingFunctional(new SetConstBoundaryDensityFunctional2D<T,Descriptor>(rho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T rho, Array<T,2> velocity) {
    applyProcessingFunctional(new IniConstEquilibriumFunctional2D<T,Descriptor>(rho, velocity), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void stripeOffDensityOffset(BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain, T deltaRho) {
    applyProcessingFunctional(new StripeOffDensityOffsetFunctional2D<T,Descriptor>(deltaRho), domain, lattice);
}

template<typename T, template<class U> class Descriptor>
void setCompositeDynamics( BlockLatticeBase2D<T,Descriptor>& lattice, Box2D domain,
                           CompositeDynamics<T,Descriptor>* compositeDynamics )
{
    applyProcessingFunctional (
            new InstantiateCompositeDynamicsFunctional2D<T,Descriptor>(compositeDynamics),
            domain, lattice );
}

}  // namespace plb

#endif  // BLOCK_OPERATIONS_2D_HH
