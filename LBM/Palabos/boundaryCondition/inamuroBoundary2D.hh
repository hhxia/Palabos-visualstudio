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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef INAMURO_BOUNDARY_2D_HH
#define INAMURO_BOUNDARY_2D_HH

#include "boundaryCondition/inamuroBoundary2D.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.h"
#include "boundaryCondition/inamuroAnalyticalDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics2D.h"
#include "boundaryCondition/boundaryInstantiator2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class InamuroBoundaryManager2D {
public:
    template<int direction, int orientation> static BoundaryCompositeDynamics<T,Descriptor>*
        getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation> static DataProcessorGenerator2D<T>*
        getVelocityBoundaryProcessor(Box2D domain);

    template<int direction, int orientation> static BoundaryCompositeDynamics<T,Descriptor>*
        getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation> static DataProcessorGenerator2D<T>*
        getPressureBoundaryProcessor(Box2D domain);

    template<int xNormal, int yNormal> static BoundaryCompositeDynamics<T,Descriptor>*
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal> static DataProcessorGenerator2D<T>*
        getExternalVelocityCornerProcessor(plint x, plint y);

    template<int xNormal, int yNormal> static BoundaryCompositeDynamics<T,Descriptor>*
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal> static DataProcessorGenerator2D<T>*
        getInternalVelocityCornerProcessor(plint x, plint y);
};

////////// InamuroBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InamuroBoundaryManager2D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new InamuroAnalyticalVelocityDynamics<T,Descriptor, direction, orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    InamuroBoundaryManager2D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box2D domain)
{
    return 0;
}

 
template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* InamuroBoundaryManager2D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new InamuroAnalyticalPressureDynamics<T,Descriptor, direction, orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator2D<T>*
    InamuroBoundaryManager2D<T,Descriptor>::
        getPressureBoundaryProcessor(Box2D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* InamuroBoundaryManager2D<T,Descriptor>::
    getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    InamuroBoundaryManager2D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y)
{
    return new OuterVelocityCornerProcessorGenerator2D<T,Descriptor, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
BoundaryCompositeDynamics<T,Descriptor>* InamuroBoundaryManager2D<T,Descriptor>::
    getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics2D<T,Descriptor, xNormal, yNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal>
DataProcessorGenerator2D<T>*
    InamuroBoundaryManager2D<T,Descriptor>::getInternalVelocityCornerProcessor
        (plint x, plint y)
{
    return 0;
}

////////// Factory function //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createInamuroBoundaryCondition2D() {
    return new BoundaryConditionInstantiator2D <
                   T, Descriptor,
                   InamuroBoundaryManager2D<T,Descriptor> > ();
}

}  // namespace plb

#endif
