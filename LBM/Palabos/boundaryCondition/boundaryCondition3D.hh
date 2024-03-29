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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef BOUNDARY_CONDITION_3D_HH
#define BOUNDARY_CONDITION_3D_HH

#include "boundaryCondition/boundaryCondition3D.h"
#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/boundaryInstantiator3D.h"
#include "core/blockSurface3D.h"
#include "core/plbDebug.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        boundary::BcType bcType )
{
    setVelocityConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), bcType );
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        Box3D applicationDomain,
        boundary::BcType bcType )
{
    setVelocityConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), applicationDomain, bcType );
}


template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        Box3D block, Box3D applicationDomain,
        boundary::BcType bcType )
{
    plint boundaryWidth = 1;
    BlockSurface3D surf(block, boundaryWidth);
    Box3D intersection;
    if (intersect(surf.surface0N(), applicationDomain, intersection)) {
        addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface0P(), applicationDomain, intersection)) {
        addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1N(), applicationDomain, intersection)) {
        addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1P(), applicationDomain, intersection)) {
        addVelocityBoundary1P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2N(), applicationDomain, intersection)) {
        addVelocityBoundary2N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2P(), applicationDomain, intersection)) {
        addVelocityBoundary2P(intersection, lattice, bcType);
    }

    if (intersect(surf.edge0NN(), applicationDomain, intersection)) {
        addExternalVelocityEdge0NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0NP(), applicationDomain, intersection)) {
        addExternalVelocityEdge0NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PN(), applicationDomain, intersection)) {
        addExternalVelocityEdge0PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0PP(), applicationDomain, intersection)) {
        addExternalVelocityEdge0PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge1NN(), applicationDomain, intersection)) {
        addExternalVelocityEdge1NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1NP(), applicationDomain, intersection)) {
        addExternalVelocityEdge1NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PN(), applicationDomain, intersection)) {
        addExternalVelocityEdge1PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1PP(), applicationDomain, intersection)) {
        addExternalVelocityEdge1PP(intersection, lattice, bcType);
    }

    if (intersect(surf.edge2NN(), applicationDomain, intersection)) {
        addExternalVelocityEdge2NN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2NP(), applicationDomain, intersection)) {
        addExternalVelocityEdge2NP(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PN(), applicationDomain, intersection)) {
        addExternalVelocityEdge2PN(intersection, lattice, bcType);
    }
    if (intersect(surf.edge2PP(), applicationDomain, intersection)) {
        addExternalVelocityEdge2PP(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNNN(), applicationDomain, intersection)) {
        addExternalVelocityCornerNNN (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNNP(), applicationDomain, intersection)) {
        addExternalVelocityCornerNNP (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPN(), applicationDomain, intersection)) {
        addExternalVelocityCornerNPN (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerNPP(), applicationDomain, intersection)) {
        addExternalVelocityCornerNPP (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNN(), applicationDomain, intersection)) {
        addExternalVelocityCornerPNN (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPNP(), applicationDomain, intersection)) {
        addExternalVelocityCornerPNP (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPN(), applicationDomain, intersection)) {
        addExternalVelocityCornerPPN (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
    if (intersect(surf.cornerPPP(), applicationDomain, intersection)) {
        addExternalVelocityCornerPPP (
                intersection.x0, intersection.y0, intersection.z0, lattice, bcType);
    }
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        boundary::BcType bcType )
{
    setPressureConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), bcType );
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        Box3D applicationDomain,
        boundary::BcType bcType )
{
    setPressureConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), applicationDomain, bcType );
}


template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition3D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase3D<T,Descriptor>& lattice,
        Box3D block, Box3D applicationDomain,
        boundary::BcType bcType )
{
    plint boundaryWidth = 1;
    BlockSurface3D surf(block, boundaryWidth);
    Box3D intersection;
    if (intersect(surf.surface0N(), applicationDomain, intersection)) {
        addPressureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface0P(), applicationDomain, intersection)) {
        addPressureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1N(), applicationDomain, intersection)) {
        addPressureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface1P(), applicationDomain, intersection)) {
        addPressureBoundary1P(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2N(), applicationDomain, intersection)) {
        addPressureBoundary2N(intersection, lattice, bcType);
    }
    if (intersect(surf.surface2P(), applicationDomain, intersection)) {
        addPressureBoundary2P(intersection, lattice, bcType);
    }

    if (intersect(surf.edge0NN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge0NP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge0PN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge0PP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }

    if (intersect(surf.edge1NN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge1NP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge1PN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge1PP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }

    if (intersect(surf.edge2NN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge2NP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge2PN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.edge2PP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }

    if (intersect(surf.cornerNNN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerNNP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerNPN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerNPP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPNN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPNP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPPN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPPP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
}

template<typename T, template<typename U> class Descriptor>
class RegularizedBoundaryManager3D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getVelocityBoundaryProcessor(Box3D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getPressureBoundaryProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityEdgeProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityEdgeProcessor(Box3D domain);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y, plint z);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y, plint z);
};

template<typename T, template<typename U> class Descriptor>
class EquilibriumBoundaryManager3D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getVelocityBoundaryProcessor(Box3D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getPressureBoundaryProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);

    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityEdgeProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityEdgeProcessor(Box3D domain);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y, plint z);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y, plint z);
};

template<typename T, template<typename U> class Descriptor>
class InterpolationBoundaryManager3D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getVelocityBoundaryProcessor(Box3D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator3D<T>*
            getPressureBoundaryProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityEdgeProcessor(Box3D domain);

    template<int plane, int normal1, int normal2>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityEdgeProcessor(Box3D domain);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y, plint z);

    template<int xNormal, int yNormal, int zNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal>
        static DataProcessorGenerator3D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y, plint z);
};


////////// RegularizedBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* RegularizedBoundaryManager3D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* RegularizedBoundaryManager3D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedDensityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getPressureBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeProcessor(Box3D domain)
{
    return new OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return new OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    RegularizedBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return 0;
}


////////// EquilibriumBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* EquilibriumBoundaryManager3D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* EquilibriumBoundaryManager3D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedDensityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getPressureBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeProcessor(Box3D domain)
{
    return new OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return new OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    EquilibriumBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return 0;
}


////////// InterpolationBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager3D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new VelocityDirichletBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box3D domain)
{
     return new PlaneFdBoundaryProcessorGenerator3D
                   <T,Descriptor, direction,orientation>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager3D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new DensityDirichletBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getPressureBoundaryProcessor(Box3D domain)
{
    return new PlaneFdBoundaryProcessorGenerator3D
                   <T,Descriptor, direction,orientation>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeProcessor(Box3D domain)
{
    return new OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return new OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    InterpolationBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return 0;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createLocalBoundaryCondition3D() {
    return new BoundaryConditionInstantiator3D <
                   T, Descriptor,
                   RegularizedBoundaryManager3D<T,Descriptor> > ();
}

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createEquilibriumBoundaryCondition3D() {
    return new BoundaryConditionInstantiator3D <
                   T, Descriptor,
                   EquilibriumBoundaryManager3D<T,Descriptor> > ();
}

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createInterpBoundaryCondition3D() {
    return new BoundaryConditionInstantiator3D <
                   T, Descriptor,
                   InterpolationBoundaryManager3D<T,Descriptor> > ();
}

}  // namespace plb

#endif
