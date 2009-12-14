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

#ifndef ZOU_HE_BOUNDARY_3D_HH
#define ZOU_HE_BOUNDARY_3D_HH

#include "boundaryCondition/zouHeBoundary3D.h"
#include "boundaryCondition/zouHeDynamics.h"
#include "boundaryCondition/zouHeDynamics.hh"
#include "boundaryCondition/regularizedBoundaryDynamics3D.h"
#include "boundaryCondition/boundaryInstantiator3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ZouHeBoundaryManager3D {
public:
    template<int direction, int orientation> static BoundaryCompositeDynamics<T,Descriptor>*
        getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation> static DataProcessorGenerator3D<T>*
        getVelocityBoundaryProcessor(Box3D domain);

    template<int direction, int orientation> static BoundaryCompositeDynamics<T,Descriptor>*
        getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int direction, int orientation> static DataProcessorGenerator3D<T>*
        getPressureBoundaryProcessor(Box3D domain);

    template<int plane, int normal1, int normal2> static BoundaryCompositeDynamics<T,Descriptor>*
        getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2> static DataProcessorGenerator3D<T>*
        getExternalVelocityEdgeProcessor(Box3D domain);

    template<int plane, int normal1, int normal2> static BoundaryCompositeDynamics<T,Descriptor>*
        getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int plane, int normal1, int normal2> static DataProcessorGenerator3D<T>*
        getInternalVelocityEdgeProcessor(Box3D domain);

    template<int xNormal, int yNormal, int zNormal> static BoundaryCompositeDynamics<T,Descriptor>*
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal> static DataProcessorGenerator3D<T>*
        getExternalVelocityCornerProcessor(plint x, plint y, plint z);

    template<int xNormal, int yNormal, int zNormal> static BoundaryCompositeDynamics<T,Descriptor>*
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics);
    template<int xNormal, int yNormal, int zNormal> static DataProcessorGenerator3D<T>*
        getInternalVelocityCornerProcessor(plint x, plint y, plint z);
};

////////// ZouHeBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* ZouHeBoundaryManager3D<T,Descriptor>::
    getVelocityBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new ZouHeVelocityDynamics<T,Descriptor, direction, orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getVelocityBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
BoundaryCompositeDynamics<T,Descriptor>* ZouHeBoundaryManager3D<T,Descriptor>::
    getPressureBoundaryDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new ZouHePressureDynamics<T,Descriptor, direction, orientation>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int direction, int orientation>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getPressureBoundaryProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getExternalVelocityEdgeProcessor(Box3D domain)
{
    return new OuterVelocityEdgeProcessorGenerator3D<T,Descriptor, plane,normal1,normal2>(domain);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
BoundaryCompositeDynamics<T,Descriptor>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerEdgeDynamics3D<T,Descriptor,plane,normal1,normal2>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int plane, int normal1, int normal2>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getInternalVelocityEdgeProcessor(Box3D domain)
{
    return 0;
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new StoreVelocityDynamics<T,Descriptor>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getExternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return new OuterVelocityCornerProcessorGenerator3D<T,Descriptor, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
BoundaryCompositeDynamics<T,Descriptor>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerDynamics(Dynamics<T,Descriptor>* baseDynamics)
{
    return new RegularizedVelocityInnerCornerDynamics3D<T,Descriptor,xNormal,yNormal,zNormal>(baseDynamics);
}

template<typename T, template<typename U> class Descriptor>
template<int xNormal, int yNormal, int zNormal>
DataProcessorGenerator3D<T>*
    ZouHeBoundaryManager3D<T,Descriptor>::
        getInternalVelocityCornerProcessor(plint x, plint y, plint z)
{
    return 0;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createZouHeBoundaryCondition3D() {
    return new BoundaryConditionInstantiator3D <
                   T, Descriptor,
                   ZouHeBoundaryManager3D<T,Descriptor> > ();
}

}  // namespace plb

#endif
