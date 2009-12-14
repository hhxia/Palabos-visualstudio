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

#ifndef BOUNDARY_CONDITION_2D_H
#define BOUNDARY_CONDITION_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "boundaryCondition/boundaryCondition.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor2D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class OnLatticeBoundaryCondition2D {
public:
    virtual ~OnLatticeBoundaryCondition2D() { }

    virtual void addVelocityBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addPressureBoundary0N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary0P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary1N( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary1P( Box2D domain, BlockLatticeBase2D<T,Descriptor>& lattice,
                                        boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addExternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addInternalVelocityCornerNN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerNP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPN( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPP( plint x, plint y, BlockLatticeBase2D<T,Descriptor>& lattice,
                                              boundary::BcType bcType=boundary::dirichlet ) =0;

    /// Set velocity/Neumann condition on outer boundaries of the lattice.
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set velocity/Neumann condition on a sub-domain, on the outer boundaries of
    ///   the lattice.
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block resp. multi-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     **/
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                Box2D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set velocity/Neumann condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain.
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                Box2D block, Box2D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set Pressure condition on outer boundaries of the lattice.
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set Pressure condition on a sub-domain, on the outer boundaries of
    ///   the lattice.
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block resp. multi-block. For boundaries inside the domain, use
     *  the method which takes two Box2D arguments.
     *  Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                Box2D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set Pressure condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain.
    /** Attention: pressure conditions are implemented for edges only. On corners,
     *  this function has no effect.
     **/
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase2D<T,Descriptor>& lattice,
                                                Box2D block, Box2D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

};


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createLocalBoundaryCondition2D();

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createEquilibriumBoundaryCondition2D();

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition2D<T,Descriptor>* createInterpBoundaryCondition2D();

}  // namespace plb

#endif
