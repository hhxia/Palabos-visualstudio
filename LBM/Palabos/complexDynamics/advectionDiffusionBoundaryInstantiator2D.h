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

/* Main author: Orestis Malaspinas
 */

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_2D_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/geometry2D.h"
#include "simulationSetup/latticeInitializer2D.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator2D : 
        public OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Descriptor> 
{
public:
    AdvectionDiffusionBoundaryConditionInstantiator2D();

    void addTemperatureBoundary0N(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureBoundary0P(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureBoundary1N(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureBoundary1P(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);

    void addTemperatureCornerNN(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureCornerNP(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureCornerPN(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice);
    void addTemperatureCornerPP(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice);

private:
    template<int direction, int orientation>
        void addTemperatureBoundary(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice);
    template<int normalX, int normalY>
        void addTemperatureCorner(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice);
};

///////// class AdvectionDiffusionBoundaryConditionInstantiator2D ////////////////////////

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
        AdvectionDiffusionBoundaryConditionInstantiator2D()
{ }

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    PLB_PRECONDITION(domain.x0==domain.x1 || domain.y0==domain.y1);

    setCompositeDynamics (
            lattice, domain,
            BoundaryManager::template getAdvectionDiffusionBoundaryDynamics<direction,orientation>(new NoDynamics<T,Descriptor>) );

    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getAdvectionDiffusionBoundaryProcessor<direction,orientation>(domain);
    if (generator) {
        lattice.addInternalProcessor(*generator);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int xNormal, int yNormal>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureCorner(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    setCompositeDynamics (
            lattice, Box2D(x,x,y,y),
            BoundaryManager::template getAdvectionDiffusionCornerDynamics<xNormal,yNormal>(new NoDynamics<T,Descriptor>) );

    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getAdvectionDiffusionCornerProcessor<xNormal,yNormal>(x, y);
    if (generator) {
        lattice.addInternalProcessor(*generator);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary0N(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<0,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary0P(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<0,1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary1N(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary1P(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<1,1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNN(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1,-1>(x,y, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNP(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1,1>(x,y, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPN(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureCorner<1,-1>(x,y, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPP(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice)
{
    addTemperatureCorner<1,1>(x,y, lattice);
}

}

#endif
