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

/** \file
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH

#include "complexDynamics/advectionDiffusionBoundaries.h"
#include "complexDynamics/advectionDiffusionBoundaryCondition3D.h"
#include "complexDynamics/advectionDiffusionBoundaryInstantiator3D.h"

namespace plb {
template<typename T, template<typename U> class Descriptor>
void OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>::setTemperatureConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice )
{
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();

    addTemperatureBoundary0N(Box3D(   0,   0,   1,ny-2,   1,nz-2), lattice);
    addTemperatureBoundary0P(Box3D(nx-1,nx-1,   1,ny-2,   1,nz-2), lattice);
    addTemperatureBoundary1N(Box3D(   1,nx-2,   0,   0,   1,nz-2), lattice);
    addTemperatureBoundary1P(Box3D(   1,nx-2,ny-1,ny-1,   1,nz-2), lattice);
    addTemperatureBoundary2N(Box3D(   1,nx-2,   1,ny-2,   0,   0), lattice);
    addTemperatureBoundary2P(Box3D(   1,nx-2,   1,ny-2,nz-1,nz-1), lattice);

    addTemperatureEdge0NN(Box3D(   1,nx-2,   0,   0,   0,   0), lattice);
    addTemperatureEdge0NP(Box3D(   1,nx-2,   0,   0,nz-1,nz-1), lattice);
    addTemperatureEdge0PN(Box3D(   1,nx-2,ny-1,ny-1,   0,   0), lattice);
    addTemperatureEdge0PP(Box3D(   1,nx-2,ny-1,ny-1,nz-1,nz-1), lattice);

    addTemperatureEdge1NN(Box3D(   0,   0,   1,ny-2,   0,   0), lattice);
    addTemperatureEdge1NP(Box3D(nx-1,nx-1,   1,ny-2,   0,   0), lattice);
    addTemperatureEdge1PN(Box3D(   0,   0,   1,ny-2,nz-1,nz-1), lattice);
    addTemperatureEdge1PP(Box3D(nx-1,nx-1,   1,ny-2,nz-1,nz-1), lattice);

    addTemperatureEdge2NN(Box3D(   0,   0,   0,   0,   1,nz-2), lattice);
    addTemperatureEdge2NP(Box3D(   0,   0,ny-1,ny-1,   1,nz-2), lattice);
    addTemperatureEdge2PN(Box3D(nx-1,nx-1,   0,   0,   1,nz-2), lattice);
    addTemperatureEdge2PP(Box3D(nx-1,nx-1,ny-1,ny-1,   1,nz-2), lattice);

    addTemperatureCornerNNN(   0,   0,   0, lattice);
    addTemperatureCornerNNP(   0,   0,nz-1, lattice);
    addTemperatureCornerNPN(   0,ny-1,   0, lattice);
    addTemperatureCornerNPP(   0,ny-1,nz-1, lattice);
    addTemperatureCornerPNN(nx-1,   0,   0, lattice);
    addTemperatureCornerPNP(nx-1,   0,nz-1, lattice);
    addTemperatureCornerPPN(nx-1,ny-1,   0, lattice);
    addTemperatureCornerPPP(nx-1,ny-1,nz-1, lattice);
}

template<typename T, template<typename U> class Descriptor>
class AdvectionDiffusionBoundaryManager3D {
public:

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getTemperatureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getTemperatureBoundaryProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getTemperatureEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getTemperatureEdgeProcessor(Box3D domain);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getTemperatureCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getTemperatureCornerProcessor(plint x, plint y, plint z);

};


////////// AdvectionDiffusionBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
    getTemperatureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new AdvectionDiffusionBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
        getTemperatureBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
        getTemperatureEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new AdvectionDiffusionEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
        getTemperatureEdgeProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
        getTemperatureCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new AdvectionDiffusionCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    AdvectionDiffusionBoundaryManager3D<T,Descriptor>::
        getTemperatureCornerProcessor(plint x, plint y, plint z)
{
    return 0;
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Descriptor>* createLocalAdvectionDiffusionBoundaryCondition3D() {
    return new AdvectionDiffusionBoundaryConditionInstantiator3D <
                   T, Descriptor,
                   AdvectionDiffusionBoundaryManager3D<T,Descriptor> > ();
}

}  // namespace plb

#endif
