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

/** \file A helper for initialising 3D boundaries -- header file.  */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H
#define ADVECTION_DIFFUSION_BOUNDARY_INSTANTIATOR_3D_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/geometry3D.h"
#include "simulationSetup/latticeInitializer3D.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
class AdvectionDiffusionBoundaryConditionInstantiator3D : 
		public OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor> {
public:
    AdvectionDiffusionBoundaryConditionInstantiator3D();

    void addTemperatureBoundary0N (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureBoundary0P (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureBoundary1N (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureBoundary1P (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureBoundary2N (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureBoundary2P (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );

    void addTemperatureEdge0NN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge0NP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge0PN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge0PP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge1NN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge1NP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge1PN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge1PP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge2NN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge2NP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge2PN (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureEdge2PP (
            Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );

    void addTemperatureCornerNNN (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerNNP (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerNPN (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerNPP (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerPNN (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerPNP (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerPPN (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
    void addTemperatureCornerPPP (
            plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );

private:
    template<int direction, int orientation>
        void addTemperatureBoundary (
                Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    template<int plane, int normal1, int normal2>
        void addTemperatureEdge (
                Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice );
    template<int normalX, int normalY, int normalZ>
        void addTemperatureCorner (
                plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice );
};


///////// class AdvectionDiffusionBoundaryConditionInstantiator3D ////////////////////////

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::AdvectionDiffusionBoundaryConditionInstantiator3D()
{ }

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int direction, int orientation>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    PLB_PRECONDITION( domain.x0==domain.x1 || domain.y0==domain.y1 || domain.z0==domain.z1 );

    setCompositeDynamics (
            lattice, domain,
            BoundaryManager::template getTemperatureBoundaryDynamics<direction,orientation>(new NoDynamics<T,Descriptor>) );

    DataProcessorGenerator3D<T>* generator
        = BoundaryManager::template getTemperatureBoundaryProcessor<direction,orientation>(domain);
    if (generator) {
        lattice.addInternalProcessor(*generator);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int plane, int normal1, int normal2>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    PLB_PRECONDITION(
            ( domain.x0==domain.x1 && domain.y0==domain.y1 ) ||
            ( domain.x0==domain.x1 && domain.z0==domain.z1 ) ||
            ( domain.y0==domain.y1 && domain.z0==domain.z1 ) );

    setCompositeDynamics (
            lattice, domain,
            BoundaryManager::template getTemperatureEdgeDynamics<plane,normal1,normal2>(new NoDynamics<T,Descriptor>) );

    DataProcessorGenerator3D<T>* generator
        = BoundaryManager::template getTemperatureEdgeProcessor<plane,normal1,normal2>(domain);
    if (generator) {
        lattice.addInternalProcessor(*generator);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCorner(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    setCompositeDynamics (
            lattice, Box3D(x,x, y,y, z,z),
            BoundaryManager::template  getTemperatureCornerDynamics<xNormal,yNormal,zNormal>(new NoDynamics<T,Descriptor>) );

    DataProcessorGenerator3D<T>* generator
        = BoundaryManager::template getTemperatureCornerProcessor<xNormal,yNormal,zNormal>(x, y, z);
    if (generator) {
        lattice.addInternalProcessor(*generator);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary0N(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<0,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary0P(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<0,1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary1N(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary1P(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary2N(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<2,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureBoundary2P(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureBoundary<2, 1>(domain, lattice);
}


template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge0NN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<0,-1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge0NP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<0,-1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge0PN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<0, 1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge0PP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<0, 1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge1NN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<1,-1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge1NP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<1,-1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge1PN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<1, 1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge1PP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<1, 1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge2NN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<2,-1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge2NP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<2,-1, 1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge2PN(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<2, 1,-1>(domain, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureEdge2PP(Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureEdge<2, 1, 1>(domain, lattice);
}



template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNNN(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1,-1,-1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNNP(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1,-1, 1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNPN(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1, 1,-1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerNPP(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner<-1, 1, 1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPNN(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner< 1,-1,-1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPNP(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner< 1,-1, 1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPPN(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner< 1, 1,-1>(x,y,z, lattice);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void AdvectionDiffusionBoundaryConditionInstantiator3D<T,Descriptor,BoundaryManager>::
    addTemperatureCornerPPP(plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice)
{
    addTemperatureCorner< 1, 1, 1>(x,y,z, lattice);
}

}  // namespace plb

#endif
