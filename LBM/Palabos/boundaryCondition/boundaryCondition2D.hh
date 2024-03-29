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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef BOUNDARY_CONDITION_2D_HH
#define BOUNDARY_CONDITION_2D_HH

#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/regularizedBoundaryDynamics2D.h"
#include "boundaryCondition/equilibriumBoundaryDynamics.h"
#include "boundaryCondition/boundaryInstantiator2D.h"
#include "core/blockSurface2D.h"
#include "core/plbDebug.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        boundary::BcType bcType )
{
    setVelocityConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), bcType );
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        Box2D applicationDomain,
        boundary::BcType bcType )
{
    setVelocityConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), applicationDomain, bcType );
}


template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setVelocityConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        Box2D block, Box2D applicationDomain,
        boundary::BcType bcType )
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addVelocityBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addVelocityBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addVelocityBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addVelocityBoundary1P(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        addExternalVelocityCornerNN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        addExternalVelocityCornerNP(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        addExternalVelocityCornerPN(intersection.x0, intersection.y0, lattice, bcType);
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        addExternalVelocityCornerPP(intersection.x0, intersection.y0, lattice, bcType);
    }
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        boundary::BcType bcType )
{
    setPressureConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), bcType );
}

template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        Box2D applicationDomain,
        boundary::BcType bcType )
{
    setPressureConditionOnBlockBoundaries (
            lattice, lattice.getBoundingBox(), applicationDomain, bcType );
}


template<typename T, template<typename U> class Descriptor>
void OnLatticeBoundaryCondition2D<T,Descriptor>::setPressureConditionOnBlockBoundaries (
        BlockLatticeBase2D<T,Descriptor>& lattice,
        Box2D block, Box2D applicationDomain,
        boundary::BcType bcType )
{
    plint boundaryWidth = 1;
    BlockSurface2D surf(block, boundaryWidth);
    Box2D intersection;
    if (intersect(surf.edge0N(), applicationDomain, intersection)) {
        addPressureBoundary0N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge0P(), applicationDomain, intersection)) {
        addPressureBoundary0P(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1N(), applicationDomain, intersection)) {
        addPressureBoundary1N(intersection, lattice, bcType);
    }
    if (intersect(surf.edge1P(), applicationDomain, intersection)) {
        addPressureBoundary1P(intersection, lattice, bcType);
    }

    if (intersect(surf.cornerNN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerNP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPN(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
    if (intersect(surf.cornerPP(), applicationDomain, intersection)) {
        PLB_ASSERT( false );
    }
}

template<typename T, template<typename U> class Descriptor>
class RegularizedBoundaryManager2D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getVelocityBoundaryProcessor(Box2D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getPressureBoundaryProcessor(Box2D domain);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y);
};

template<typename T, template<typename U> class Descriptor>
class EquilibriumBoundaryManager2D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getVelocityBoundaryProcessor(Box2D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getPressureBoundaryProcessor(Box2D domain);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y);
};

template<typename T, template<typename U> class Descriptor>
class InterpolationBoundaryManager2D {
public:
    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getVelocityBoundaryProcessor(Box2D domain);

    template<int direction, int orientation>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation>
        static DataProcessorGenerator2D<T>*
            getPressureBoundaryProcessor(Box2D domain);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getExternalVelocityCornerProcessor(plint x, plint y);

    template<int xNormal, int yNormal>
        static BoundaryCompositeDynamics<T,Descriptor>*
            getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal>
        static DataProcessorGenerator2D<T>*
            getInternalVelocityCornerProcessor(plint x, plint y);
};


////////// RegularizedBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box2D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedDensityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryProcessor(Box2D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* RegularizedBoundaryManager2D<T,Descriptor>::
    getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y)
{
    return new OuterVelocityCornerProcessorGenerator2D<T,Descriptor, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    RegularizedBoundaryManager2D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T,Descriptor, xNormal, yNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    RegularizedBoundaryManager2D<T,Descriptor>::getInternalVelocityCornerProcessor
        (plint x, plint y)
{
    return 0;
}


////////// EquilibriumBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new EquilibriumVelocityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box2D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new EquilibriumDensityBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryProcessor(Box2D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* EquilibriumBoundaryManager2D<T,Descriptor>::
    getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y)
{
    return new OuterVelocityCornerProcessorGenerator2D<T,Descriptor, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    EquilibriumBoundaryManager2D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T,Descriptor, xNormal, yNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    EquilibriumBoundaryManager2D<T,Descriptor>::getInternalVelocityCornerProcessor
        (plint x, plint y)
{
    return 0;
}


////////// InterpolationBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager2D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new VelocityDirichletBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    InterpolationBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box2D domain)
{
    return new StraightFdBoundaryProcessorGenerator2D
                   < T,Descriptor,direction,orientation  >  (domain);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager2D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new DensityDirichletBoundaryDynamics
                   <T,Descriptor,direction,orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    InterpolationBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryProcessor(Box2D domain)
{
    return new StraightFdBoundaryProcessorGenerator2D
                   < T,Descriptor,direction,orientation >  (domain);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager2D<T,Descriptor>::
    getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    InterpolationBoundaryManager2D<T,Descriptor>::getExternalVelocityCornerProcessor(plint x, plint y)
{
    return new OuterVelocityCornerProcessorGenerator2D<T,Descriptor, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* InterpolationBoundaryManager2D<T,Descriptor>::
    getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T,Descriptor, xNormal, yNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    InterpolationBoundaryManager2D<T,Descriptor>::getInternalVelocityCornerProcessor (plint x, plint y)
{
    return 0;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createLocalBoundaryCondition2D() {
    return new BoundaryConditionInstantiator2D <
                   T, Descriptor,
                   RegularizedBoundaryManager2D<T,Descriptor> >;
}

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createEquilibriumBoundaryCondition2D() {
    return new BoundaryConditionInstantiator2D <
                   T, Descriptor,
                   EquilibriumBoundaryManager2D<T,Descriptor> >;
}

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createInterpBoundaryCondition2D() {
    return new BoundaryConditionInstantiator2D <
                   T, Descriptor,
                   InterpolationBoundaryManager2D<T,Descriptor > >;
}

}  // namespace plb

#endif
