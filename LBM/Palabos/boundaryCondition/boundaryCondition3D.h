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

/** \file A helper for initialising 3D boundaries -- header file.  */

#ifndef BOUNDARY_CONDITION_3D_H
#define BOUNDARY_CONDITION_3D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase3D.h"
#include "boundaryCondition/boundaryCondition.h"
#include "boundaryCondition/finiteDifferenceBoundaryProcessor3D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class OnLatticeBoundaryCondition3D {
public:
    virtual ~OnLatticeBoundaryCondition3D() { }

    virtual void addVelocityBoundary0N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary0P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary1N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary1P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary2N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addVelocityBoundary2P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addPressureBoundary0N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary0P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary1N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary1P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary2N (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addPressureBoundary2P (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addExternalVelocityEdge0NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge0NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge0PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge0PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge1NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge1NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge1PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge1PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge2NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge2NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge2PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityEdge2PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addInternalVelocityEdge0NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge0NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge0PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge0PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge1NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge1NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge1PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge1PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge2NN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge2NP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge2PN (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityEdge2PP (
              Box3D domain, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addExternalVelocityCornerNNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerNNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerNPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerNPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addExternalVelocityCornerPPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    virtual void addInternalVelocityCornerNNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerNNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerNPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerNPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPNN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPNP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPPN (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;
    virtual void addInternalVelocityCornerPPP (
              plint x, plint y, plint z, BlockLatticeBase3D<T,Descriptor>& lattice,
              boundary::BcType bcType=boundary::dirichlet ) =0;

    /// Set velocity/Neumann condition on outer boundaries of the lattice.
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set velocity/Neumann condition on a sub-domain, on the outer boundaries of
    ///   the lattice.
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block resp. multi-block. For boundaries inside the domain, use
     *  the method which takes two Box3D arguments.
     **/
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                Box3D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set velocity/Neumann condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain.
    void setVelocityConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                Box3D block, Box3D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );


    /// Set Pressure condition on outer boundaries of the lattice.
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set Pressure condition on a sub-domain, on the outer boundaries of
    ///   the lattice.
    /** Attention: this function only has an effect when it is used on the outer surface
     *  of the atomic-block resp. multi-block. For boundaries inside the domain, use
     *  the method which takes two Box3D arguments.
     **/
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                Box3D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

    /// Set Pressure condition on the block boundaries, but only on places which
    ///    intersect with the area of applicationDomain.
    void setPressureConditionOnBlockBoundaries( BlockLatticeBase3D<T,Descriptor>& lattice,
                                                Box3D block, Box3D applicationDomain,
                                                boundary::BcType bcType=boundary::dirichlet );

};


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createLocalBoundaryCondition3D();

template<typename T, template<typename U> class Descriptor>
OnLatticeBoundaryCondition3D<T,Descriptor>* createEquilibriumBoundaryCondition3D();

template<typename T, template<typename U> class Descriptor>
    OnLatticeBoundaryCondition3D<T,Descriptor>* createInterpBoundaryCondition3D();

}  // namespace plb

#endif
