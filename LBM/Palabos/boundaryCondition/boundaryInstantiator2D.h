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

/** \file A helper for initialising 2D boundaries -- header file.  */
#ifndef BOUNDARY_INSTANTIATOR_2D_H
#define BOUNDARY_INSTANTIATOR_2D_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/geometry2D.h"
#include "boundaryCondition/boundaryCondition2D.h"
#include "boundaryCondition/neumannCondition2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "simulationSetup/latticeInitializer2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
class BoundaryConditionInstantiator2D : public OnLatticeBoundaryCondition2D<T,Descriptor> {
public:
    BoundaryConditionInstantiator2D();

    void addVelocityBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addVelocityBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addVelocityBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addVelocityBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                               boundary::BcType bcType );

    void addPressureBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addPressureBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addPressureBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );
    void addPressureBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                boundary::BcType bcType );

    void addExternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addExternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addExternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addExternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );

    void addInternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addInternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addInternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
    void addInternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                      boundary::BcType bcType );
private:
    template<int direction, int orientation>
        void addVelocityBoundary( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                  boundary::BcType bcType );
    template<int direction, int orientation>
        void addPressureBoundary( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                  boundary::BcType bcType );
    template<int normalX, int normalY>
        void addExternalVelocityCorner( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType );
    template<int normalX, int normalY>
        void addInternalVelocityCorner( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType );
};

///////// class BoundaryConditionInstantiator2D ////////////////////////

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::BoundaryConditionInstantiator2D()
{ }

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addVelocityBoundary(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0==domain.x1 || domain.y0==domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction==0) ? orientation : 0,
        normalY = (direction==1) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics (
            lattice, domain,
            BoundaryManager::template
                getVelocityBoundaryDynamics<direction,orientation>(new NoDynamics<T,Descriptor>) );

    // In case an outflow condition is used, start by instantiating a data processor which copies
    //   all velocity value from the previous lattice cell.
    if (bcType==boundary::outflow || bcType==boundary::neumann) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyVelocityFunctional2D<T,Descriptor, normalX, normalY>,
                domain, lattice, processorLevel );
    }
    // In case a normal outflow condition is used, start by instantiating a data processor which copies
    //   the normal velocity value from the previous lattice cell, and sets the other components to zero.
    if (bcType==boundary::normalOutflow) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyNormalVelocityFunctional2D<T,Descriptor, normalX, normalY>,
                domain, lattice, processorLevel );
    }
    else
    // In case a freeslip condition is used, start by instantiating a data processor which copies
    //   the tangential velocity values from the previous lattice cell.
    if (bcType==boundary::freeslip) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyTangentialVelocityFunctional2D<T,Descriptor, normalX, normalY>,
                domain, lattice, processorLevel );
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data processor.
    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getVelocityBoundaryProcessor<direction,orientation>(domain);
    if (generator) {
        // In case the boundary is of type Neumann, the processor must be instantiated at level 1,
        //   because it must be executed after the copy-density data processor.
        plint processorLevel = bcType==boundary::dirichlet ? 0 : 1;
        lattice.addInternalProcessor(*generator, processorLevel);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int direction, int orientation>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addPressureBoundary(Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice, boundary::BcType bcType)
{
    PLB_PRECONDITION(domain.x0==domain.x1 || domain.y0==domain.y1);

    // Convert (direction,orientation) description of the normal vector into a (normalX,normalY)
    //   description, as it is requried by the data processor for Neumann boundaries.
    enum {
        normalX = (direction==0) ? orientation : 0,
        normalY = (direction==1) ? orientation : 0
    };

    // Instantiate the dynamics of the boundary as a composite dynamics, based on the one currently
    //   residing on the lattice.
    setCompositeDynamics (
            lattice, domain,
            BoundaryManager::template
                getPressureBoundaryDynamics<direction,orientation>(new NoDynamics<T,Descriptor>) );

    // In case a Neumann condition is used, start by instantiating a data processor which copies
    //   the density value from the previous lattice cell.
    if (bcType==boundary::neumann) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyDensityFunctional2D<T,Descriptor, normalX, normalY>,
                domain, lattice, processorLevel );
    }

    // If the boundary condition has a non-local component, instantiate a corresponding data processor.
    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getPressureBoundaryProcessor<direction,orientation>(domain);
    if (generator) {
        // In case the boundary is of type Neumann, the processor must be instantiated at level 1,
        //   because it must be executed after the copy-density data processor.
        plint processorLevel = bcType==boundary::dirichlet ? 0 : 1;
        lattice.addInternalProcessor(*generator, processorLevel);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addExternalVelocityCorner(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice, boundary::BcType bcType)
{
    setCompositeDynamics (
            lattice, Box2D(x,x,y,y),
            BoundaryManager::template
                getExternalVelocityCornerDynamics<xNormal,yNormal>(new NoDynamics<T,Descriptor>) );

    if (bcType==boundary::neumann || bcType==boundary::outflow) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }
    else
    if (bcType==boundary::freeslip) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyTangentialVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }
    else
    if (bcType==boundary::normalOutflow) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyNormalVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }


    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getExternalVelocityCornerProcessor<xNormal,yNormal>(x, y);
    if (generator) {
        // In case the boundary is of type Neumann, the processor must be instantiated at level 1,
        //   because it must be executed after the copy-density data processor.
        plint processorLevel = bcType==boundary::dirichlet ? 0 : 1;
        lattice.addInternalProcessor(*generator, processorLevel);
        delete generator;
    }
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
template<int xNormal, int yNormal>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addInternalVelocityCorner(plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice, boundary::BcType bcType)
{
    setCompositeDynamics (
            lattice, Box2D(x,x,y,y),
            BoundaryManager::template
                getInternalVelocityCornerDynamics<xNormal,yNormal>(new NoDynamics<T,Descriptor>) );

    if (bcType==boundary::neumann || bcType==boundary::outflow) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }
    else
    if (bcType==boundary::freeslip) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyTangentialVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }
    else
    if (bcType==boundary::normalOutflow) {
        // This data processor must be at level 1, because it also copies values on the envelope,
        //   and must therefore be executed only after a communication step.
        plint processorLevel = 1;
        integrateProcessingFunctional (
                new CopyNormalVelocityFunctional2D<T,Descriptor, xNormal,yNormal>,
                Box2D(x,x,y,y), lattice, processorLevel );
    }

    DataProcessorGenerator2D<T>* generator
        = BoundaryManager::template getInternalVelocityCornerProcessor<xNormal,yNormal>(x, y);
    if (generator) {
        // In case the boundary is of type Neumann, the processor must be instantiated at level 1,
        //   because it must be executed after the copy-density data processor.
        plint processorLevel = bcType==boundary::dirichlet ? 0 : 1;
        lattice.addInternalProcessor(*generator, processorLevel);
        delete generator;
    }
}


template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addVelocityBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addVelocityBoundary<0,-1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addVelocityBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addVelocityBoundary<0,1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addVelocityBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addVelocityBoundary<1,-1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addVelocityBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addVelocityBoundary<1,1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addPressureBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addPressureBoundary<0,-1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addPressureBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addPressureBoundary<0,1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addPressureBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addPressureBoundary<1,-1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addPressureBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                           boundary::BcType bcType )
{
    addPressureBoundary<1,1>(domain, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addExternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addExternalVelocityCorner<-1,-1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addExternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addExternalVelocityCorner<-1,1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addExternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addExternalVelocityCorner<1,-1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addExternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addExternalVelocityCorner<1,1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addInternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addInternalVelocityCorner<-1,-1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addInternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addInternalVelocityCorner<-1,1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addInternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addInternalVelocityCorner<1,-1>(x,y, lattice, bcType);
}

template<typename T, template<typename U> class Descriptor, class BoundaryManager>
void BoundaryConditionInstantiator2D<T,Descriptor,BoundaryManager>::
    addInternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                 boundary::BcType bcType )
{
    addInternalVelocityCorner<1,1>(x,y, lattice, bcType);
}

}  // namespace plb


#endif
