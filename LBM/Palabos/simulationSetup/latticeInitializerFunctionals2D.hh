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
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef LATTICE_INITIALIZER_FUNCTIONALS_2D_HH
#define LATTICE_INITIALIZER_FUNCTIONALS_2D_HH

#include "simulationSetup/latticeInitializer2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
OneCellFunctional2D<T,Descriptor>::~OneCellFunctional2D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}


template<typename T, template<class U> class Descriptor>
void OneCellFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale)
{ }


template<typename T, template<class U> class Descriptor>
OneCellIndexedFunctional2D<T,Descriptor>::~OneCellIndexedFunctional2D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional2D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale)
{ }


template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional2D<T,Descriptor>::GenericLatticeFunctional2D (
        OneCellFunctional2D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional2D<T,Descriptor>::GenericLatticeFunctional2D (
        GenericLatticeFunctional2D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional2D<T,Descriptor>::~GenericLatticeFunctional2D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional2D<T,Descriptor>& GenericLatticeFunctional2D<T,Descriptor>::operator= (
        GenericLatticeFunctional2D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            f->execute(lattice.get(iX,iY));
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional2D<T,Descriptor>*
    GenericLatticeFunctional2D<T,Descriptor>::clone() const
{
    return new GenericLatticeFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    f->getModificationPattern(isWritten);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericLatticeFunctional2D<T,Descriptor>::appliesTo() const {
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale) {
    f->rescale(dxScale, dtScale);
}



template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T,Descriptor>::GenericIndexedLatticeFunctional2D (
        OneCellIndexedFunctional2D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T,Descriptor>::GenericIndexedLatticeFunctional2D (
        GenericIndexedLatticeFunctional2D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T,Descriptor>::~GenericIndexedLatticeFunctional2D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T,Descriptor>&
    GenericIndexedLatticeFunctional2D<T,Descriptor>::operator= (
        GenericIndexedLatticeFunctional2D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    Dot2D relativeOffset = lattice.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            f->execute ( iX+relativeOffset.x,
                         iY+relativeOffset.y,
                         lattice.get(iX,iY) );
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional2D<T,Descriptor>*
    GenericIndexedLatticeFunctional2D<T,Descriptor>::clone() const
{
    return new GenericIndexedLatticeFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    f->getModificationPattern(isWritten);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericIndexedLatticeFunctional2D<T,Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    f->rescale(dxScale, dtScale);
}


/* *************** Class InstantiateDynamicsFunctional2D ************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T,Descriptor>::InstantiateDynamicsFunctional2D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T,Descriptor>::InstantiateDynamicsFunctional2D (
        InstantiateDynamicsFunctional2D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T,Descriptor>&
    InstantiateDynamicsFunctional2D<T,Descriptor>::operator= (
        InstantiateDynamicsFunctional2D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T,Descriptor>::~InstantiateDynamicsFunctional2D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.attributeDynamics(iX,iY, dynamics->clone());
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsFunctional2D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional2D<T,Descriptor>*
    InstantiateDynamicsFunctional2D<T,Descriptor>::clone() const 
{
    return new InstantiateDynamicsFunctional2D<T,Descriptor>(*this);
}


/* ************* Class InstantiateComplexDomainDynamicsFunctional2D ** */

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional2D (
        Dynamics<T,Descriptor>* dynamics_, DomainFunctional2D* domain_ )
    : dynamics(dynamics_),
      domain(domain_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional2D (
        InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      domain(rhs.domain->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>&
    InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::operator= (
        InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    delete domain; domain = rhs.domain->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::~InstantiateComplexDomainDynamicsFunctional2D()
{
    delete dynamics;
    delete domain;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::process (
        Box2D boundingBox, BlockLattice2D<T,Descriptor>& lattice )
{
    Dot2D relativeOffset = lattice.getLocation();
    for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
        for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
            if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y)) {
                lattice.attributeDynamics(iX,iY, dynamics->clone());
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>*
    InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>::clone() const 
{
    return new InstantiateComplexDomainDynamicsFunctional2D<T,Descriptor>(*this);
}


/* ************* Class InstantiateDotDynamicsFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T,Descriptor>::InstantiateDotDynamicsFunctional2D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T,Descriptor>::InstantiateDotDynamicsFunctional2D (
        InstantiateDotDynamicsFunctional2D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T,Descriptor>&
    InstantiateDotDynamicsFunctional2D<T,Descriptor>::operator= (
        InstantiateDotDynamicsFunctional2D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T,Descriptor>::~InstantiateDotDynamicsFunctional2D()
{
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional2D<T,Descriptor>::process (
        DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
        Dot2D const& dot = dotList.getDot(iDot);
        lattice.attributeDynamics(dot.x, dot.y, dynamics->clone());
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDotDynamicsFunctional2D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional2D<T,Descriptor>*
    InstantiateDotDynamicsFunctional2D<T,Descriptor>::clone() const 
{
    return new InstantiateDotDynamicsFunctional2D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromMaskFunctional2D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T,Descriptor>::DynamicsFromMaskFunctional2D (
        Dynamics<T,Descriptor>* dynamics_, bool whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T,Descriptor>::DynamicsFromMaskFunctional2D (
        DynamicsFromMaskFunctional2D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T,Descriptor>&
    DynamicsFromMaskFunctional2D<T,Descriptor>::operator= (
        DynamicsFromMaskFunctional2D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T,Descriptor>::~DynamicsFromMaskFunctional2D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                      ScalarField2D<T>& mask )
{
    Dot2D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            bool flag = (bool) util::roundToInt (
                           mask.get(iX+offset.x, iY+offset.y) );
            if ( util::boolIsEqual(flag, whichFlag) ) {
                lattice.attributeDynamics(iX,iY, dynamics->clone());
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromMaskFunctional2D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional2D<T,Descriptor>*
    DynamicsFromMaskFunctional2D<T,Descriptor>::clone() const 
{
    return new DynamicsFromMaskFunctional2D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromIntMaskFunctional2D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T,Descriptor>::DynamicsFromIntMaskFunctional2D (
        Dynamics<T,Descriptor>* dynamics_, int whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T,Descriptor>::DynamicsFromIntMaskFunctional2D (
        DynamicsFromIntMaskFunctional2D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T,Descriptor>&
    DynamicsFromIntMaskFunctional2D<T,Descriptor>::operator= (
        DynamicsFromIntMaskFunctional2D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T,Descriptor>::~DynamicsFromIntMaskFunctional2D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                      ScalarField2D<T>& mask )
{
    Dot2D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            bool flag = (bool) util::roundToInt (
                           mask.get(iX+offset.x, iY+offset.y) );
            if ( flag == whichFlag ) {
                lattice.attributeDynamics(iX,iY, dynamics->clone());
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromIntMaskFunctional2D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
    isWritten[1] = false;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional2D<T,Descriptor>*
    DynamicsFromIntMaskFunctional2D<T,Descriptor>::clone() const 
{
    return new DynamicsFromIntMaskFunctional2D<T,Descriptor>(*this);
}


/* ************* Class SetConstBoundaryVelocityFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional2D<T,Descriptor>::SetConstBoundaryVelocityFunctional2D (
        Array<T,Descriptor<T>::d> velocity )
    : u(velocity)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).defineVelocity(u);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional2D<T,Descriptor>*
    SetConstBoundaryVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityFunctional2D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    T latticeVelocityScale = dtScale / dxScale;
    u[0] *= latticeVelocityScale;
    u[1] *= latticeVelocityScale;
}


/* ************* Class SetConstBoundaryDensityFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional2D<T,Descriptor>::SetConstBoundaryDensityFunctional2D(T rho_)
    : rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).defineDensity(rho);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional2D<T,Descriptor>*
    SetConstBoundaryDensityFunctional2D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryDensityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryDensityFunctional2D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale )
{
    // Nothing to do: density is scale invariant.
}


/* ************* Class IniConstEquilibriumFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional2D<T,Descriptor>::IniConstEquilibriumFunctional2D (
        T density_, Array<T,Descriptor<T>::d> velocity )
    : rhoBar(Descriptor<T>::rhoBar(density_)),
      j     (density_*velocity[0], density_*velocity[1]),
      jSqr  (VectorTemplate<T,Descriptor>::normSqr(j))
{ }

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                lattice.get(iX,iY)[iPop] =
                    lattice.get(iX,iY).computeEquilibrium(iPop, rhoBar, j, jSqr);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional2D<T,Descriptor>*
    IniConstEquilibriumFunctional2D<T,Descriptor>::clone() const
{
    return new IniConstEquilibriumFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumFunctional2D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale)
{
    T latticeVelocityScale = dtScale / dxScale;
    j[0] *= latticeVelocityScale;
    j[1] *= latticeVelocityScale;
    jSqr *= latticeVelocityScale*latticeVelocityScale;
}


/* ************* Class StripeOffDensityOffsetFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional2D<T,Descriptor>::StripeOffDensityOffsetFunctional2D(T deltaRho_)
    : deltaRho(deltaRho_)
{ }

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor>& cell         = lattice.get(iX,iY);
            Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
            plint orderOfDecomposition = 0;
            std::vector<T> rawData;
            dynamics.decompose(cell, rawData, orderOfDecomposition);
            T& rhoBar = rawData[0];
            rhoBar -= deltaRho;
            dynamics.recompose(cell, rawData, orderOfDecomposition);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional2D<T,Descriptor>*
    StripeOffDensityOffsetFunctional2D<T,Descriptor>::clone() const
{
    return new StripeOffDensityOffsetFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT StripeOffDensityOffsetFunctional2D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten) const
{
    isWritten[0] = true;
}

// Density is scale invariant.
template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional2D<T,Descriptor>::rescale(T dxScale, T dtScale) { }


/* ************* Class InstantiateCompositeDynamicsFunctional2D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::InstantiateCompositeDynamicsFunctional2D (
        CompositeDynamics<T,Descriptor>* compositeDynamics_ )
    : compositeDynamics(compositeDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::InstantiateCompositeDynamicsFunctional2D (
        InstantiateCompositeDynamicsFunctional2D<T,Descriptor> const& rhs )
    : compositeDynamics(rhs.compositeDynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T,Descriptor>&
    InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::operator= (
        InstantiateCompositeDynamicsFunctional2D<T,Descriptor> const& rhs )
{
    delete compositeDynamics; compositeDynamics = rhs.compositeDynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::~InstantiateCompositeDynamicsFunctional2D()
{
    delete compositeDynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Dynamics<T,Descriptor>* baseDynamics = lattice.get(iX,iY).getDynamics().cloneComposeable();
            Dynamics<T,Descriptor>* dynamics = compositeDynamics->cloneWithNewBase(baseDynamics);
            lattice.attributeDynamics(iX,iY, dynamics);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::appliesTo() const
{
    // Composite dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::getModificationPattern (
        std::vector<bool>& isWritten ) const
{
    isWritten[0] = true;
}

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional2D<T,Descriptor>*
    InstantiateCompositeDynamicsFunctional2D<T,Descriptor>::clone() const 
{
    return new InstantiateCompositeDynamicsFunctional2D<T,Descriptor>(*this);
}

}  // namespace plb

#endif  // LATTICE_INITIALIZER_FUNCTIONALS_2D_HH
